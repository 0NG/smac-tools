#include <iostream>
#include <fstream>
#include <iomanip>
#include <atomic>
#include <stdexcept>
#include <algorithm>

#include <cxxopts.hpp>

#include "ortools_kit.hpp"
#include "sbox_espresso.hpp"
#include "aes.hpp"

using namespace std;

constexpr int NROUND = 6;
constexpr int NTHREAD = 10;
constexpr int THRESHOLD = 18;
constexpr bool ENABLELOG = true;
constexpr int TIMEOUT = -1;

static vector<BoolVec> dp;
static vector<BoolVec> state;
static vector<BoolVec> intermediate;
static vector<BoolVar> actives;
static vector<BoolVec> SA12;
static vector<BoolVec> LSA23;

static auto ddt = compute_ddt();
static auto Lbm = compute_lbm();

void print_solution(const sat::CpSolverResponse &response, const BoolVec &x)
{
    for (int i = 0; i < 16; ++i) {
        int tmp = 0;
        for (int j = 0; j < 8; ++j)
            tmp = (tmp << 1) | SolutionIntegerValue(response, x[8 * i + j]);
        cout << (tmp > 0) << ", ";
    }
    return;
}

void print_solution_bytes(const sat::CpSolverResponse &response, const BoolVec &x)
{
    for (int i = 0; i < 16; ++i) {
        int tmp = 0;
        for (int j = 0; j < 8; ++j)
            tmp = (tmp << 1) | SolutionIntegerValue(response, x[8 * i + j]);
        cout << std::hex << std::setw(2) << std::setfill('0') << tmp << std::dec << ",";
    }
    return;
}

std::vector<std::array<int, 16>> gen_byte_trail(const sat::CpSolverResponse &response)
{
    std::vector<std::array<int, 16>> dpActive;
    for (auto &dpi : dp) {
        std::array<int, 16> tmp;
        for (int i = 0; i < 16; ++i) {
            int byte = 0;
            for (int j = 0; j < 8; ++j)
                byte = (byte << 1) | SolutionIntegerValue(response, dpi[8 * i + j]);
            tmp[i] = (byte > 0);
        }
        dpActive.push_back(tmp);
    }
    return dpActive;
}

std::array<int, 16> gen_pattern(const sat::CpSolverResponse &response, const std::array<uint8_t, 16> &sigma)
{
    // Y = sigma(X) should be treated as: Y[i] = X[sigma[i]]
    std::array<int, 16> pattern;
    std::fill(pattern.begin(), pattern.end(), -1);

    for (int i = 0; i < state.size(); i += 3) {
        for (int j = 0; j < 16; ++j) {
            int byte = 0;
            for (int k = 0; k < 8; ++k) {
                byte = (byte << 1) | SolutionIntegerValue(response, state[i][8 * j + k]);
            }

            if (byte) {
                if (pattern[j] != -1 && pattern[j] != sigma[j]) throw std::logic_error("There is a bug in generating patterns");

                pattern[j] = sigma[j];
            }
        }
    }

    return pattern;
}

BoolVec bmXbv(sat::CpModelBuilder &model, const bm<128> &m, BoolVec &bv)
{
    BoolVec r = NewBoolVec(model, 128);

    for (int row = 0; row < 128; ++row) {
        BoolVec tmp{{model.TrueVar(), r[row]}};
        for (int col = 0; col < 128; ++col)
            if (m[row][col])
                tmp.push_back(bv[col]);
        model.AddBoolXor(tmp);
    }

    return r;
}

void bmXbv(sat::CpModelBuilder &model, const bm<128> &m, BoolVec &bv, BoolVec &r)
{
    for (int row = 0; row < 128; ++row) {
        BoolVec tmp{{model.TrueVar(), r[row]}};
        for (int col = 0; col < 128; ++col)
            if (m[row][col])
                tmp.push_back(bv[col]);
        model.AddBoolXor(tmp);
    }

    return;
}

static BoolVec mui_bv(sat::CpModelBuilder &model, BoolVec &bv)
{
    BoolVec r = NewBoolVec(model, 16);

    for (int i = 0; i < 16; ++i) {
        BoolVec tmp1;
        BoolVec tmp2;
        for (int j = 0; j < 8; ++j) {
            tmp1.push_back(bv[i * 8 + j]);
            tmp2.push_back(bv[i * 8 + j].Not());

            // MILP
            model.AddGreaterOrEqual(r[i], bv[i * 8 + j]);
        }

        model.AddBoolOr(tmp1).OnlyEnforceIf(r[i]);
        model.AddBoolAnd(tmp2).OnlyEnforceIf(r[i].Not());

        // MILP
        model.AddLessOrEqual(r[i], LinearExpr::Sum(tmp1));
    }

    return r;
}

inline auto bytes2bits(const std::vector<unsigned short> bytes)
{
    std::vector<int64_t> bits;
    const auto cnt = bytes.size();

    for (int i = 0; i < cnt; ++i)
        for (int j = 0; j < 8; ++j)
            bits.push_back((bytes[i] >> (7 - j)) & 1);

    return bits;
}

void add_espresso_ddt(sat::CpModelBuilder &model, BoolVec &xy)
{
    BoolVec iterms;
    for (auto &eIterm : espressoBf) {
        BoolVec tmp0;
        BoolVec tmp1;
        auto literal = model.NewBoolVar();

        for (int i = 0; i < 16; ++i)
            if (eIterm[i] == 0) {
                tmp0.push_back(xy[i].Not());
                tmp1.push_back(xy[i]);
            } else if (eIterm[i] == 1) {
                tmp0.push_back(xy[i]);
                tmp1.push_back(xy[i].Not());
            }

        model.AddBoolAnd(tmp0).OnlyEnforceIf(literal);
        model.AddBoolOr(tmp1).OnlyEnforceIf(literal.Not());
        iterms.push_back(literal);
    }

    model.AddBoolOr(iterms);
    return;
}

void add_ddt_active(sat::CpModelBuilder &model, BoolVec &x, BoolVec &y)
{
    for (int i = 0; i < 16; ++i) {
        BoolVec tmpIn;
        BoolVec tmpOut;
        for (int j = 0; j < 8; ++j) tmpIn.push_back(x[8 * i + j]);
        for (int j = 0; j < 8; ++j) tmpOut.push_back(y[8 * i + j]);

        BoolVec tmpInNot;
        BoolVec tmpOutNot;
        for (int j = 0; j < 8; ++j) tmpInNot.push_back(x[8 * i + j].Not());
        for (int j = 0; j < 8; ++j) tmpOutNot.push_back(y[8 * i + j].Not());

        auto isActive =  model.NewBoolVar();
        model.AddBoolOr(tmpIn).OnlyEnforceIf(isActive);
        model.AddBoolAnd(tmpInNot).OnlyEnforceIf(isActive.Not());
        model.AddBoolOr(tmpOut).OnlyEnforceIf(isActive);
        model.AddBoolAnd(tmpOutNot).OnlyEnforceIf(isActive.Not());

        actives.push_back(isActive);
    }
    return;
}

void add_ddt(sat::CpModelBuilder &model, BoolVec &x, BoolVec &y)
{
    /*
    std::vector<std::vector<int64_t>> values;
    for (int i = 0; i < 256; ++i)
        for (int j = 0; j < 256; ++j)
            if (ddt[i][j]) {
                values.push_back(bytes2bits({static_cast<unsigned short>(i), static_cast<unsigned short>(j)}));
            }
    */

    for (int i = 0; i < 16; ++i) {
        BoolVec tmp;
        for (int j = 0; j < 8; ++j) tmp.push_back(x[8 * i + j]);
        for (int j = 0; j < 8; ++j) tmp.push_back(y[8 * i + j]);

        //BVAssign(model, tmp, values);
        add_espresso_ddt(model, tmp);
    }
    return;
}

static void add_simple_l(sat::CpModelBuilder &model, BoolVec &bv, BoolVec &r)
{
    auto muiBv = mui_bv(model, bv);
    auto muiR = mui_bv(model, r);
    for (int i = 0; i < 4; ++i) {
        BoolVec tmp;
        auto d = model.NewBoolVar();
        for (int j = 0; j < 4; ++j) {
            model.AddGreaterOrEqual(d, muiBv[SRp[4 * i + j]]);
            model.AddGreaterOrEqual(d, muiR[4 * i + j]);
            tmp.push_back(muiBv[SRp[4 * i + j]]);
            tmp.push_back(muiR[4 * i + j]);
        }
        model.AddGreaterOrEqual(LinearExpr::Sum(tmp), 5 * d);
    }
    return;
}

std::vector<std::array<int, 16>> search_active(const int nRound, const std::array<uint8_t, 16> &sigma, const int nThread, const int threshold, const int logLevel, const std::vector<int> &dpi0, const int timeLimit)
{
    dp.clear();
    state.clear();
    intermediate.clear();
    actives.clear();
    SA12.clear();
    LSA23.clear();

    SatParameters parameters;
    parameters.set_num_search_workers(nThread);
    parameters.set_log_search_progress(logLevel == 1);
    if (timeLimit != -1) parameters.set_max_time_in_seconds(timeLimit);

    CpModelBuilder cp_model;
    auto pbm = perm2bm(sigma);

    {
        auto tmp1 = NewBoolVec(cp_model, 128);
        auto tmp2 = NewBoolVec(cp_model, 128);
        LSA23.push_back(tmp1);
        LSA23.push_back(tmp2);
        for (int i = 0; i < 128; ++i) {
            cp_model.AddEquality(tmp1[i], 0);
            cp_model.AddEquality(tmp2[i], 0);
        }
    }

    // first round
    {
        auto dp0 = NewBoolVec(cp_model, 128);
        dp.push_back(dp0);
        cp_model.AddBoolOr(dp0);

        auto sdp0 = bmXbv(cp_model, pbm, dp0);

        state.push_back(sdp0);
        state.push_back(dp0);
        state.push_back(dp0);
    }

    for (int i = 1; i < nRound - 1; ++i) {
        auto dpi = NewBoolVec(cp_model, 128);
        dp.push_back(dpi);

        auto d0 = NewBoolVec(cp_model, 128);
        auto d1 = NewBoolVec(cp_model, 128);
        add_ddt_active(cp_model, state[3 * (i - 1) + 0], d0);
        add_ddt_active(cp_model, state[3 * (i - 1) + 1], d1);

        if (i == (nRound - 1) - 1) {
            SA12.push_back(d0);
            SA12.push_back(d1);
        }

        auto d2 = NewBoolVec(cp_model, 128);
        auto d3 = NewBoolVec(cp_model, 128);
        bmXbv(cp_model, Lbm, d0, d2);
        add_simple_l(cp_model, d0, d2);
        bmXbv(cp_model, Lbm, d1, d3);
        add_simple_l(cp_model, d1, d3);

        LSA23.push_back(d2);
        LSA23.push_back(d3);

        auto s1 = NewBoolVec(cp_model, 128);
        auto s2 = NewBoolVec(cp_model, 128);
        auto s3 = NewBoolVec(cp_model, 128);

        BVXor(cp_model, d2, dpi, s2);
        BVXor(cp_model, d3, dpi, s3);

        auto d4 = NewBoolVec(cp_model, 128);
        BVXor3(cp_model, LSA23[2 * (i - 1) + 0], LSA23[2 * (i - 1) + 1], dpi, d4);
        bmXbv(cp_model, pbm, d4, s1);

        state.push_back(s1);
        state.push_back(s2);
        state.push_back(s3);
    }

    // last round (i == nRound - 1), dpi is L^-1 times the real dpi
    if constexpr (1) {
        const int i = nRound - 1;
        auto dpi = NewBoolVec(cp_model, 128);
        dp.push_back(dpi);

        add_ddt_active(cp_model, state[3 * (i - 1) + 0], dpi);
        add_ddt_active(cp_model, state[3 * (i - 1) + 1], dpi);

        BVXor(cp_model, SA12[0], SA12[1], dpi);

        cp_model.AddBoolOr(dpi);
    }

    for (auto i : dpi0) {
        if (i == -1) break;
        for (int b = 0; b < 128; ++b)
            cp_model.AddEquality(dp[i][b], 0);
    }

    auto obj = cp_model.NewIntVar(Domain(0, actives.size()));
    cp_model.AddEquality(obj, LinearExpr::Sum(actives));

    cp_model.Minimize(obj);

    /*===========================================================================*/
    Model model;
    std::atomic<bool> stopped(false);
    model.GetOrCreate<TimeLimit>()->RegisterExternalBooleanAsLimit(&stopped);
    model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
        auto newObj = SolutionIntegerValue(r, obj);
        if (threshold != -1 && newObj < threshold) stopped = true;
    }));
    model.Add(NewSatParameters(parameters));
    /*===========================================================================*/

    auto model_built = cp_model.Build();

    const auto response = SolveCpModel(model_built, &model);
    const auto status = response.status();

    // CpSolverStatus::INFEASIBLE = 3
    std::vector<std::array<int, 16>> dpActive;
    if (status == CpSolverStatus::OPTIMAL || status == CpSolverStatus::FEASIBLE) {
        auto result = SolutionIntegerValue(response, obj);

        dpActive = gen_byte_trail(response);
        std::array<int, 16> tmp;
        tmp[0] = result;
        dpActive.push_back(tmp);

        if (logLevel == 2) {
            if (result < threshold) {
                const auto pattern = gen_pattern(response, sigma);
                cout << "pattern: "; for (auto &s : pattern) cout << s << ' '; cout << endl;
                return dpActive;
            }
            cout << "sigma: "; for (auto &s : sigma) cout << static_cast<unsigned int>(s) << ' '; cout << endl;
            return dpActive;
        }

        cout << "====================" << endl;
        cout << "status: " << status << endl;
        cout << "sigma: "; for (auto &s : sigma) cout << static_cast<unsigned int>(s) << ' '; cout << endl;
        cout << "obj  : " << result << endl;

        for (int i = 0; i < dp.size(); ++i) {
            cout << "mui(dp" << i << ")                                 : "; print_solution(response, dp[i]); cout << endl; 
        }

        cout << "wall time: " << response.wall_time() << endl;
    } else {
        cout << "====================" << endl;
        cout << "infeasible: "; for (auto &s : sigma) cout << static_cast<unsigned int>(s) << ' '; cout << endl;
    }

    return dpActive;
}

template<int nRound>
bool _verify_trail(const std::array<uint8_t, 16> &sigma, const int nThread, std::vector<std::array<int, 16>> &dpActive, const int objVal) { return false; }

template<>
bool _verify_trail<3>(const std::array<uint8_t, 16> &sigma, const int nThread, std::vector<std::array<int, 16>> &dpActive, const int objVal)
{
    auto pbm = perm2bm(sigma);

    SatParameters parameters;
    parameters.set_num_search_workers(nThread);
    parameters.set_log_search_progress(false);

    CpModelBuilder cp_model;

    auto iLdp2  = NewBoolVec(cp_model, 128);
    auto dp0    = NewBoolVec(cp_model, 128);
    auto sA10   = NewBoolVec(cp_model, 128);
    auto sA20   = NewBoolVec(cp_model, 128);
    auto sdp0   = bmXbv(cp_model, pbm, dp0);

    add_ddt(cp_model, sdp0, sA10);
    add_ddt(cp_model, dp0, sA20);
    BVXor(cp_model, sA10, sA20, iLdp2);

    auto dp1  = NewBoolVec(cp_model, 128);
    auto sdp1 = bmXbv(cp_model, pbm, dp1);

    add_ddt(cp_model, sdp1, iLdp2);

    auto d0 = bmXbv(cp_model, Lbm, sA10);
    add_simple_l(cp_model, sA10, d0);
    auto d1 = BVXor(cp_model, d0, dp1);

    //
    auto dt = bmXbv(cp_model, Lbm, sA20);
    add_simple_l(cp_model, sA20, dt);
    auto dtt = BVXor(cp_model, dt, dp1);
    //

    add_ddt(cp_model, d1, iLdp2);

    auto realdp = bmXbv(cp_model, Lbm, iLdp2);

    auto muia10 = mui_bv(cp_model, sA10);
    auto muia20 = mui_bv(cp_model, sA20);
    auto muisdp1 = mui_bv(cp_model, sdp1);
    auto muid1 = mui_bv(cp_model, d1);
    auto muiiLdp2 = mui_bv(cp_model, iLdp2);

    for (int i = 0; i < 16; ++i) cp_model.AddEquality(muisdp1[i], muiiLdp2[i]);
    for (int i = 0; i < 16; ++i) cp_model.AddEquality(muid1[i], muiiLdp2[i]);

    BoolVec tmpSum;
    for (int i = 0; i < 16; ++i) {
        tmpSum.push_back(muia10[i]);
        tmpSum.push_back(muia20[i]);
        tmpSum.push_back(muiiLdp2[i]);
        tmpSum.push_back(muiiLdp2[i]);
    }

    auto obj = cp_model.NewIntVar(Domain(0, tmpSum.size()));
    cp_model.AddEquality(obj, LinearExpr::Sum(tmpSum));

    auto muidp0 = mui_bv(cp_model, dp0);
    auto muidp1 = mui_bv(cp_model, dp1);
    auto muidp2 = mui_bv(cp_model, iLdp2);
    for (int i = 0; i < 16; ++i) {
        cp_model.AddEquality(muidp0[i], dpActive[0][i]);
        cp_model.AddEquality(muidp1[i], dpActive[1][i]);
        cp_model.AddEquality(muidp2[i], dpActive[2][i]);
    }

    cp_model.AddEquality(obj, objVal);

    auto model_built = cp_model.Build();

    const auto response = SolveWithParameters(model_built, parameters);
    const auto status = response.status();

    // CpSolverStatus::INFEASIBLE = 3
    if (status == CpSolverStatus::OPTIMAL || status == CpSolverStatus::FEASIBLE) {
        auto result = SolutionIntegerValue(response, obj);
        cout << "dp0: "; print_solution_bytes(response, dp0); cout << endl;
        cout << "dp1: "; print_solution_bytes(response, dp1); cout << endl;
        cout << "dp2: "; print_solution_bytes(response, realdp); cout << endl;
        cout << "Round 1:" << endl;
        cout << "sdp0: "; print_solution_bytes(response, sdp0); cout << endl;
        cout << "Round 2:" << endl;
        cout << "S(A1): "; print_solution_bytes(response, sA10); cout << endl;
        cout << "S(A2): "; print_solution_bytes(response, sA20); cout << endl;
        cout << "state1: "; print_solution_bytes(response, sdp1); cout << endl;
        cout << "state2: "; print_solution_bytes(response, d1); cout << endl;
        cout << "state3: "; print_solution_bytes(response, dtt); cout << endl;
        cout << "Round 3:" << endl;
        cout << "S(A1): "; print_solution_bytes(response, iLdp2); cout << endl;
        cout << "S(A2): "; print_solution_bytes(response, iLdp2); cout << endl;
        if (result == objVal) return true;
    }

    return false;
}

template<>
bool _verify_trail<4>(const std::array<uint8_t, 16> &sigma, const int nThread, std::vector<std::array<int, 16>> &dpActive, const int objVal)
{
    auto pbm = perm2bm(sigma);

    SatParameters parameters;
    parameters.set_num_search_workers(nThread);
    parameters.set_log_search_progress(false);

    CpModelBuilder cp_model;

    auto dp0 = NewBoolVec(cp_model, 128);
    auto dp1 = NewBoolVec(cp_model, 128);
    auto dp2 = NewBoolVec(cp_model, 128);
    auto dp3 = NewBoolVec(cp_model, 128);

    BoolVec nzExpr;
    for (int i = 0; i < 128; ++i) {
        nzExpr.push_back(dp0[i]);
        nzExpr.push_back(dp1[i]);
        nzExpr.push_back(dp2[i]);
        nzExpr.push_back(dp3[i]);
    }
    cp_model.AddGreaterOrEqual(LinearExpr::Sum(nzExpr), 1);

    auto d0 = NewBoolVec(cp_model, 128); // S(x + s(Pt1))
    auto d1 = NewBoolVec(cp_model, 128); // S(y + s(Pt0))
    //auto d2 = NewBoolVec(cp_model, 128); // LS(y + s(Pt0))
    auto d3 = NewBoolVec(cp_model, 128); // LS(y + s(Pt0)) + Pt1
    auto d4 = NewBoolVec(cp_model, 128); // S(LS(y + s(Pt0)) + Pt1)

    auto sdp0 = bmXbv(cp_model, pbm, dp0);
    auto sdp1 = bmXbv(cp_model, pbm, dp1);
    add_ddt(cp_model, sdp1, d0);
    add_ddt(cp_model, sdp0, d1);
    auto d2 = bmXbv(cp_model, Lbm, d1);
    add_simple_l(cp_model, d1, d2);
    BVXor(cp_model, d2, dp1, d3);
    add_ddt(cp_model, d3, d4);
    BVXor(cp_model, d0, d4, dp3);

    auto d5 = NewBoolVec(cp_model, 128); // S(z + Pt0)
    //auto d6 = NewBoolVec(cp_model, 128); // LS(y + s(Pt0))
    //auto d7 = NewBoolVec(cp_model, 128); // LS(z + Pt0)
    auto d8 = NewBoolVec(cp_model, 128); // LS(y + s(Pt0)) + LS(z + Pt0) + Pt2
    auto d9 = NewBoolVec(cp_model, 128); // S(LS(y + s(Pt0) + LS(z + Pt0) + Pt2)

    add_ddt(cp_model, dp0, d5);
    //auto d6 = bmXbv(cp_model, Lbm, d1);
    auto d7 = bmXbv(cp_model, Lbm, d5);
    add_simple_l(cp_model, d5, d7);
    //BVXor3(cp_model, d6, d7, dp2, d8);
    BVXor3(cp_model, d2, d7, dp2, d8);
    add_ddt(cp_model, d8, d9);
    bmXbv(cp_model, pbm, d9, dp3);

    //auto d10 = NewBoolVec(cp_model, 128); // LS(x + s(Pt1))
    auto d11 = NewBoolVec(cp_model, 128); // LS(x + s(Pt1)) + Pt2

    auto d10 = bmXbv(cp_model, Lbm, d0);
    add_simple_l(cp_model, d0, d10);
    BVXor(cp_model, d10, dp2, d11);
    add_ddt(cp_model, d11, dp3);

    auto muidp0 = mui_bv(cp_model, dp0);
    auto muidp1 = mui_bv(cp_model, dp1);
    auto muidp2 = mui_bv(cp_model, dp2);
    auto muidp3 = mui_bv(cp_model, dp3);
    auto muid3 = mui_bv(cp_model, d3);
    auto muid8 = mui_bv(cp_model, d8);
    auto muid11 = mui_bv(cp_model, d11);

    //for (int i = 0; i < 16; ++i) cp_model.AddEquality(muisdp1[i], muiiLdp2[i]);
    //for (int i = 0; i < 16; ++i) cp_model.AddEquality(muid1[i], muiiLdp2[i]);

    BoolVec tmpSum;
    for (int i = 0; i < 16; ++i) {
        tmpSum.push_back(muidp0[i]);
        tmpSum.push_back(muidp0[i]);
        tmpSum.push_back(muidp1[i]);
        tmpSum.push_back(muid3[i]);
        tmpSum.push_back(muid8[i]);
        tmpSum.push_back(muid11[i]);
    }

    auto obj = cp_model.NewIntVar(Domain(0, tmpSum.size()));
    cp_model.AddEquality(obj, LinearExpr::Sum(tmpSum));

    cp_model.AddBoolOr(dp0);
    cp_model.AddBoolOr(dp3);

    for (int i = 0; i < 16; ++i) {
        cp_model.AddEquality(muidp0[i], dpActive[0][i]);
        cp_model.AddEquality(muidp1[i], dpActive[1][i]);
        cp_model.AddEquality(muidp2[i], dpActive[2][i]);
        cp_model.AddEquality(muidp3[i], dpActive[3][i]);
    }
    cp_model.AddEquality(obj, objVal);

    auto model_built = cp_model.Build();

    const auto response = SolveWithParameters(model_built, parameters);
    const auto status = response.status();

    // CpSolverStatus::INFEASIBLE = 3
    if (status == CpSolverStatus::OPTIMAL || status == CpSolverStatus::FEASIBLE) {
        auto result = SolutionIntegerValue(response, obj);
        if (result == objVal) return true;
    }

    return false;
}

template<>
bool _verify_trail<5>(const std::array<uint8_t, 16> &sigma, const int nThread, std::vector<std::array<int, 16>> &dpActive, const int objVal)
{
    auto pbm = perm2bm(sigma);

    SatParameters parameters;
    parameters.set_num_search_workers(nThread);
    parameters.set_log_search_progress(false);

    CpModelBuilder cp_model;

    auto dp0 = NewBoolVec(cp_model, 128);
    auto dp1 = NewBoolVec(cp_model, 128);
    auto dp2 = NewBoolVec(cp_model, 128);
    auto dp3 = NewBoolVec(cp_model, 128);
    auto dp4 = NewBoolVec(cp_model, 128); // L^{-1}dp4

    BoolVec nzExpr;
    for (int i = 0; i < 128; ++i) {
        nzExpr.push_back(dp0[i]);
        nzExpr.push_back(dp1[i]);
        nzExpr.push_back(dp2[i]);
        nzExpr.push_back(dp3[i]);
        nzExpr.push_back(dp4[i]);
    }
    cp_model.AddGreaterOrEqual(LinearExpr::Sum(nzExpr), 1);

    auto d0  = NewBoolVec(cp_model, 128); // S(y + s(dp0))
    auto d1  = NewBoolVec(cp_model, 128); // S(z + dp0)
    auto d2  = NewBoolVec(cp_model, 128); // LS(y + s(dp0))
    auto d3  = NewBoolVec(cp_model, 128); // LS(z + dp0)
    auto d4  = NewBoolVec(cp_model, 128); // LS(y + s(dp0)) + LS(z + dp0) + dp2
    auto d5  = NewBoolVec(cp_model, 128); // S(LS(y + s(dp0)) + LS(z + dp0) + dp2)
    auto d6  = NewBoolVec(cp_model, 128); // s(S(LS(y + s(dp0)) + LS(z + dp0) + dp2))

    auto d7  = NewBoolVec(cp_model, 128); // S(x + s(dp1))
    auto d8  = NewBoolVec(cp_model, 128); // LS(x + s(dp1))
    auto d9  = NewBoolVec(cp_model, 128); // LS(x + s(dp1)) + dp2
    auto d10 = NewBoolVec(cp_model, 128); // S(LS(x + s(dp1)) + dp2)

    auto sdp0 = bmXbv(cp_model, pbm, dp0);
    auto sdp1 = bmXbv(cp_model, pbm, dp1);

    add_ddt(cp_model, sdp0, d0);
    add_ddt(cp_model, dp0, d1);
    bmXbv(cp_model, Lbm, d0, d2);
    bmXbv(cp_model, Lbm, d1, d3);
    BVXor3(cp_model, d2, d3, dp2, d4);
    add_ddt(cp_model, d4, d5);
    bmXbv(cp_model, pbm, d5, d6);

    add_ddt(cp_model, sdp1, d7);
    bmXbv(cp_model, Lbm, d7, d8);
    BVXor(cp_model, d8, dp2, d9);
    add_ddt(cp_model, d9, d10);

    BVXor(cp_model, d6, d10, dp4);

    auto d11  = NewBoolVec(cp_model, 128); // LS(y + s(dp0)) + dp1
    auto d12  = NewBoolVec(cp_model, 128); // S(LS(y + s(dp0)) + dp1)
    auto d13  = NewBoolVec(cp_model, 128); // LS(LS(y + s(dp0)) + dp1)
    auto d14  = NewBoolVec(cp_model, 128); // LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3
    auto d15  = NewBoolVec(cp_model, 128); // S(LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3)

    BVXor(cp_model, d2, dp1, d11);
    add_ddt(cp_model, d11, d12);
    bmXbv(cp_model, Lbm, d12, d13);
    BVXor3(cp_model, d8, d13, dp3, d14);
    add_ddt(cp_model, d14, d15);
    bmXbv(cp_model, pbm, d15, dp4);

    auto d16  = NewBoolVec(cp_model, 128); // Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2))
    auto d17  = NewBoolVec(cp_model, 128); // Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2)) + dp3

    bmXbv(cp_model, Lbm, d6, d16);
    BVXor(cp_model, d16, dp3, d17);
    add_ddt(cp_model, d17, dp4);

    auto muidp0 = mui_bv(cp_model, dp0);
    auto muidp1 = mui_bv(cp_model, dp1);
    auto muidp2 = mui_bv(cp_model, dp2);
    auto muidp3 = mui_bv(cp_model, dp3);
    auto muidp4 = mui_bv(cp_model, dp4);
    auto muid4  = mui_bv(cp_model, d4);
    auto muid9  = mui_bv(cp_model, d9);
    auto muid11 = mui_bv(cp_model, d11);
    auto muid14 = mui_bv(cp_model, d14);
    auto muid17 = mui_bv(cp_model, d17);

    BoolVec tmpSum;
    for (int i = 0; i < 16; ++i) {
        tmpSum.push_back(muidp0[i]);
        tmpSum.push_back(muidp0[i]);
        tmpSum.push_back(muidp1[i]);
        tmpSum.push_back(muid4[i]);
        tmpSum.push_back(muid9[i]);
        tmpSum.push_back(muid11[i]);
        tmpSum.push_back(muid14[i]);
        tmpSum.push_back(muid17[i]);
    }

    auto obj = cp_model.NewIntVar(Domain(0, tmpSum.size()));
    cp_model.AddEquality(obj, LinearExpr::Sum(tmpSum));

    cp_model.AddBoolOr(dp0);
    cp_model.AddBoolOr(dp4);

    for (int i = 0; i < 16; ++i) {
        cp_model.AddEquality(muidp0[i], dpActive[0][i]);
        cp_model.AddEquality(muidp1[i], dpActive[1][i]);
        cp_model.AddEquality(muidp2[i], dpActive[2][i]);
        cp_model.AddEquality(muidp3[i], dpActive[3][i]);
        cp_model.AddEquality(muidp4[i], dpActive[4][i]);
    }

    cp_model.AddEquality(obj, objVal);

    auto model_built = cp_model.Build();

    const auto response = SolveWithParameters(model_built, parameters);
    const auto status = response.status();

    // CpSolverStatus::INFEASIBLE = 3
    if (status == CpSolverStatus::OPTIMAL || status == CpSolverStatus::FEASIBLE) {
        auto result = SolutionIntegerValue(response, obj);
        if (result == objVal) return true;
    }

    return false;
}

template<>
bool _verify_trail<6>(const std::array<uint8_t, 16> &sigma, const int nThread, std::vector<std::array<int, 16>> &dpActive, const int objVal)
{
    auto pbm = perm2bm(sigma);

    SatParameters parameters;
    parameters.set_num_search_workers(nThread);
    parameters.set_log_search_progress(false);

    CpModelBuilder cp_model;

    auto dp0 = NewBoolVec(cp_model, 128);
    auto dp1 = NewBoolVec(cp_model, 128);
    auto dp2 = NewBoolVec(cp_model, 128);
    auto dp3 = NewBoolVec(cp_model, 128);
    auto dp4 = NewBoolVec(cp_model, 128);
    auto dp5 = NewBoolVec(cp_model, 128); // L^{-1}dp5

    /*
    BoolVec nzExpr;
    for (int i = 0; i < 128; ++i) {
        nzExpr.push_back(dp0[i]);
        nzExpr.push_back(dp1[i]);
        nzExpr.push_back(dp2[i]);
        nzExpr.push_back(dp3[i]);
        nzExpr.push_back(dp4[i]);
    }
    cp_model.AddGreaterOrEqual(LinearExpr::Sum(nzExpr), 1);
    */

    auto d0  = NewBoolVec(cp_model, 128); // S(y + s(dp0))
    auto d1  = NewBoolVec(cp_model, 128); // S(z + dp0)
    auto d2  = NewBoolVec(cp_model, 128); // LS(y + s(dp0))
    auto d3  = NewBoolVec(cp_model, 128); // LS(z + dp0)
    auto d4  = NewBoolVec(cp_model, 128); // LS(y + s(dp0)) + LS(z + dp0) + dp2
    auto d5  = NewBoolVec(cp_model, 128); // S(LS(y + s(dp0)) + LS(z + dp0) + dp2)
    auto d6  = NewBoolVec(cp_model, 128); // s(S(LS(y + s(dp0)) + LS(z + dp0) + dp2))

    auto d7  = NewBoolVec(cp_model, 128); // S(x + s(dp1))
    auto d8  = NewBoolVec(cp_model, 128); // LS(x + s(dp1))
    auto d9  = NewBoolVec(cp_model, 128); // LS(x + s(dp1)) + dp2
    auto d10 = NewBoolVec(cp_model, 128); // S(LS(x + s(dp1)) + dp2)

    auto n0 = NewBoolVec(cp_model, 128); // Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2))
    auto n1 = NewBoolVec(cp_model, 128); // LS(LS(x + s(dp1)) + dp2)
    auto n2 = NewBoolVec(cp_model, 128); // Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2)) + LS(LS(x + s(dp1)) + dp2) + dp4
    auto n3 = NewBoolVec(cp_model, 128); // s(Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2)) + LS(LS(x + s(dp1)) + dp2) + dp4)

    auto sdp0 = bmXbv(cp_model, pbm, dp0);
    auto sdp1 = bmXbv(cp_model, pbm, dp1);

    add_ddt(cp_model, sdp0, d0);
    add_ddt(cp_model, dp0, d1);
    bmXbv(cp_model, Lbm, d0, d2);
    add_simple_l(cp_model, d0, d2);
    bmXbv(cp_model, Lbm, d1, d3);
    add_simple_l(cp_model, d1, d3);
    BVXor3(cp_model, d2, d3, dp2, d4);
    add_ddt(cp_model, d4, d5);
    bmXbv(cp_model, pbm, d5, d6);

    add_ddt(cp_model, sdp1, d7);
    bmXbv(cp_model, Lbm, d7, d8);
    add_simple_l(cp_model, d7, d8);
    BVXor(cp_model, d8, dp2, d9);
    add_ddt(cp_model, d9, d10);

    bmXbv(cp_model, Lbm, d6, n0);
    add_simple_l(cp_model, d6, n0);
    bmXbv(cp_model, Lbm, d10, n1);
    add_simple_l(cp_model, d10, n1);
    BVXor3(cp_model, n0, n1, dp4, n2);
    bmXbv(cp_model, pbm, n2, n3);

    add_ddt(cp_model, n3, dp5);

    auto d11  = NewBoolVec(cp_model, 128); // LS(y + s(dp0)) + dp1
    auto d12  = NewBoolVec(cp_model, 128); // S(LS(y + s(dp0)) + dp1)
    auto d13  = NewBoolVec(cp_model, 128); // LS(LS(y + s(dp0)) + dp1)
    auto d14  = NewBoolVec(cp_model, 128); // LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3
    auto d15  = NewBoolVec(cp_model, 128); // S(LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3)

    auto n4  = NewBoolVec(cp_model, 128); // s(S(LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3))
    auto n5  = NewBoolVec(cp_model, 128); // Ls(S(LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3))
    auto n6  = NewBoolVec(cp_model, 128); // Ls(S(LS(x + s(dp1)) + LS(LS(y + s(dp0)) + dp1) + dp3)) + dp4

    BVXor(cp_model, d2, dp1, d11);
    add_ddt(cp_model, d11, d12);
    bmXbv(cp_model, Lbm, d12, d13);
    add_simple_l(cp_model, d12, d13);
    BVXor3(cp_model, d8, d13, dp3, d14);
    add_ddt(cp_model, d14, d15);

    bmXbv(cp_model, pbm, d15, n4);
    bmXbv(cp_model, Lbm, n4, n5);
    add_simple_l(cp_model, n4, n5);
    BVXor(cp_model, n5, dp4, n6);
    add_ddt(cp_model, n6, dp5);

    auto d16  = NewBoolVec(cp_model, 128); // Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2))
    auto d17  = NewBoolVec(cp_model, 128); // Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2)) + dp3

    auto n7  = NewBoolVec(cp_model, 128); // S(Ls(S(LS(y + s(dp0)) + LS(z + dp0) + dp2)) + dp3)

    bmXbv(cp_model, Lbm, d6, d16);
    add_simple_l(cp_model, d6, d16);
    BVXor(cp_model, d16, dp3, d17);
    add_ddt(cp_model, d17, n7);

    BVXor(cp_model, n4, n7, dp5);

    auto muidp0 = mui_bv(cp_model, dp0);
    auto muidp1 = mui_bv(cp_model, dp1);
    auto muidp2 = mui_bv(cp_model, dp2);
    auto muidp3 = mui_bv(cp_model, dp3);
    auto muidp4 = mui_bv(cp_model, dp4);
    auto muidp5 = mui_bv(cp_model, dp5);
    auto muid4  = mui_bv(cp_model, d4);
    auto muid9  = mui_bv(cp_model, d9);
    auto muid11 = mui_bv(cp_model, d11);
    auto muid14 = mui_bv(cp_model, d14);
    auto muid17 = mui_bv(cp_model, d17);

    BoolVec tmpSum;
    for (int i = 0; i < 16; ++i) {
        tmpSum.push_back(muidp0[i]);
        tmpSum.push_back(muidp0[i]);
        tmpSum.push_back(muidp1[i]);
        tmpSum.push_back(muid4[i]);
        tmpSum.push_back(muid9[i]);
        tmpSum.push_back(muid11[i]);
        tmpSum.push_back(muid14[i]);
        tmpSum.push_back(muid17[i]);
        tmpSum.push_back(muidp5[i]);
        tmpSum.push_back(muidp5[i]);
    }

    auto obj = cp_model.NewIntVar(Domain(0, tmpSum.size()));
    cp_model.AddEquality(obj, LinearExpr::Sum(tmpSum));

    cp_model.AddBoolOr(dp0);
    cp_model.AddBoolOr(dp5);

    for (int i = 0; i < 16; ++i) {
        cp_model.AddEquality(muidp0[i], dpActive[0][i]);
        cp_model.AddEquality(muidp1[i], dpActive[1][i]);
        cp_model.AddEquality(muidp2[i], dpActive[2][i]);
        cp_model.AddEquality(muidp3[i], dpActive[3][i]);
        cp_model.AddEquality(muidp4[i], dpActive[4][i]);
        cp_model.AddEquality(muidp5[i], dpActive[5][i]);
    }

    cp_model.AddEquality(obj, objVal);

    auto model_built = cp_model.Build();

    const auto response = SolveWithParameters(model_built, parameters);
    const auto status = response.status();

    // CpSolverStatus::INFEASIBLE = 3
    if (status == CpSolverStatus::OPTIMAL || status == CpSolverStatus::FEASIBLE) {
        auto result = SolutionIntegerValue(response, obj);
        if (result == objVal) return true;
    }

    return false;
}

bool verify_trail(const int nRound, const std::array<uint8_t, 16> &sigma, const int nThread, std::vector<std::array<int, 16>> &dpActive, const int objVal)
{
    switch (nRound) {
        case 3:
            return _verify_trail<3>(sigma, nThread, dpActive, objVal);
            break;
        case 4:
            return _verify_trail<4>(sigma, nThread, dpActive, objVal);
            break;
        case 5:
            return _verify_trail<5>(sigma, nThread, dpActive, objVal);
            break;
        case 6:
            return _verify_trail<6>(sigma, nThread, dpActive, objVal);
            break;
        default:
            throw std::invalid_argument("no model for validation");
            break;
    }

    return false;
}

std::array<uint8_t, 16> parse_sigma(const char *s, const int start = 0)
{
    std::array<uint8_t, 16> sigma{{}};

    int si = 0;
    int p = 0;
    for (int i = start + 0; i < start + (10 + 6 * 2 + 15); ++i) {
        if (s[i] == ',') {
            sigma[si] = p;
            p = 0;
            ++si;
            continue;
        }

        p = 10 * p + (s[i] - 0x30);
    }
    sigma[si] = p;
    return sigma;
}

std::vector<std::array<uint8_t, 16>> parse_sigma_file(const char *s)
{
    ifstream iFile;
    iFile.open(s);

    std::vector<std::array<uint8_t, 16>> sigmaList;

    std::string line;
    while (getline(iFile, line)) {
        const auto len = line.length();
        int start = 0;
        int state = 0;
        while (state != 4) {
            switch (line[start])
            {
                case 's':
                    state += (state == 0);
                    state = (state == 1 ? state : 0);
                    break;
        
                case 'g':
                    state += (state == 1);
                    state = (state == 2 ? state : 2);
                    break;
        
                case '=':
                    state += (state == 2);
                    state = (state == 3 ? state : 3);
                    break;
        
                case '{':
                    state += (state == 3);
                    state = (state == 4 ? state : 4);
                    break;
                
                default:
                    state = 0;
                    break;
            }
    
            ++start;
            if (start == len) { throw std::invalid_argument("not a valid sigma file or no sigma in this file (sg={xxx})"); }
        }
    
        sigmaList.push_back(parse_sigma(line.c_str(), start));
    }

    iFile.close();
    return sigmaList;
}

int main(int argc, char *argv[])
{
    if (argc > 1) {

        cxxopts::Options options("./executable", "SMAC t+n for n > 3. Universal version for different rates.");
        options.add_options()
            ("s,sigma", "Single sigma (e.g. 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)", cxxopts::value<std::string>()->default_value(""))
            ("f,file", "Sigma file [/path/to/file]", cxxopts::value<std::string>()->default_value(""))
            ("n,round", "Number of rounds (clocks) [4, 5, 6, 7, ...]", cxxopts::value<int>())
            ("j,thread", "Number of threads", cxxopts::value<int>()->default_value("10"))
            ("b,threshold", "Stop when a solution less than this value is found. (default: no threshold)", cxxopts::value<int>()->default_value("-1"))
            ("l,log", "Log level: 0=only trails, 1=detailed search progress, 2=only sigmas that survive (default: 0) [0, 1, 2]", cxxopts::value<int>()->default_value("0"))
            ("t,time", "Timeout, time limit (in seconds, default: unlimited)", cxxopts::value<int>()->default_value("-1"))
            ("p,dpi", "The index(s) of dp's that is/are set to 0 (separated by commas if more than one are set, default: none of them is/are set to 0)", cxxopts::value<std::vector<int>>()->default_value("-1"))
            ("verify", "Verify the byte trail (default: don't verify)", cxxopts::value<bool>()->default_value("false"))
            ("h,help", "Help menu", cxxopts::value<bool>()->default_value("false"));

        auto results = options.parse(argc, argv);
        if (results["help"].as<bool>()) {
            cout << options.help() << endl;
            return 0;
        }

        const auto sigmaStr = results["sigma"].as<std::string>();
        const auto fileName = results["file"].as<std::string>();

        std::vector<std::array<uint8_t, 16>> sigmaList;
        if (sigmaStr.length() > 0) {
            const auto sigma = parse_sigma(sigmaStr.c_str());
            sigmaList.push_back(sigma);
        } else {
            sigmaList = parse_sigma_file(fileName.c_str());
        }

        const auto nRound = results["round"].as<int>();
        const auto nThread = results["thread"].as<int>();
        const auto threshold = results["threshold"].as<int>();
        const auto logLevel = results["log"].as<int>();
        const auto timeLimit = results["time"].as<int>();
        const auto dpi = results["dpi"].as<std::vector<int>>();
        const auto verify = results["verify"].as<bool>();

        for (auto i : dpi)
            if (i == 0 || i == nRound - 1) {
                cout << "dp0 and dp" << nRound - 1 << " cannot be set to 0." << endl;
                return 0;
            }

        if (logLevel != 2) cout << sigmaList.size() << " sigmas loaded." << endl;
        if (sigmaList.size() == 0) return 0;

        if (logLevel != 2) {
            cout << "Sigma: "; for (auto &s : sigmaList[0]) cout << static_cast<unsigned int>(s) << ' '; cout << endl;
            cout << "Sigma file: '" << fileName << "'." << endl;
            cout << "Check t + " << nRound << "." << endl;
            cout << "Using " << nThread << " threads." << endl;
            cout << "Stop when cost < " << threshold << " is found." << endl;
            cout << "Log level: " << logLevel << "." << endl;
            cout << "Timeout in " << timeLimit << " seconds." << endl;
            cout << "dp "; for(auto i : dpi) cout << i << " "; cout << "is/are set to 0." << endl;
        }

        for (const auto &sigma : sigmaList) {
            auto dpActive = search_active(nRound, sigma, nThread, threshold, logLevel, dpi, timeLimit);
            if (logLevel != 2 && verify) {
                const auto isConsistent = dpActive.size() > 1 && verify_trail(nRound, sigma, nThread, dpActive, dpActive[dpActive.size() - 1][0]);
                if (!isConsistent) {
                    cout << "====================" << endl;
                    cout << "inconsistent: "; for (auto &s : sigma) cout << static_cast<unsigned int>(s) << ' '; cout << endl;
                    return 0;
                } else {
                    cout << "verified" << endl;
                }
            }
        }
        return 0;
    }

    std::vector<std::array<uint8_t, 16>> sigmaList = {{
        {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}}, // s0
        {{0,7,14,11,4,13,10,1,8,15,6,3,12,5,2,9}}, // s1
        {{7,14,15,10,12,13,3,0,4,6,1,5,8,11,2,9}}  // s42
    }};

    return 0;
}
