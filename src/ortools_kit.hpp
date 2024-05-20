#ifndef ORTOOLS_KIT_HPP
#define ORTOOLS_KIT_HPP

#include <vector>
#include <string>
#include <cassert>

#include <ortools/sat/cp_model.h>
#include <ortools/util/time_limit.h>

using namespace operations_research;
using namespace operations_research::sat;

using BoolVec = std::vector<sat::BoolVar>;
using IntVec = std::vector<sat::IntVar>;

BoolVec NewBoolVec(sat::CpModelBuilder &model, const int size)
{
    // bv[0] is the lsb
    BoolVec bv;
    for (int i = 0; i < size; ++i) bv.push_back(model.NewBoolVar());
    return bv;
}


IntVec BV2IV(BoolVec &bv)
{
    const int len = bv.size();

    IntVec iv;
    for (int i = 0; i < len; ++i) iv.push_back(sat::IntVar(bv[i]));
    return iv;
}

void BVXor(sat::CpModelBuilder &model, BoolVec &bv0, BoolVec &bv1, BoolVec &bv2)
{
    /*
     * bv0 ^ bv1 = bv2
     */
    const int len = bv0.size();

    for (int i = 0; i < len; ++i) {
        model.AddBoolXor({bv0[i], bv1[i], bv2[i], model.TrueVar()});
    }

    return;
}

void BVAssignIf(sat::CpModelBuilder &model, BoolVec &bv, const std::vector<std::vector<int64_t>> &values, sat::IntVar b)
{
    auto iv = BV2IV(bv);
    iv.push_back(b);
    auto table = model.AddAllowedAssignments(iv);

    for (auto &value : values) {
        table.AddTuple(value);
    }

    return;
} 

void BVAssignIf(sat::CpModelBuilder &model, BoolVec &bv, const std::vector<std::vector<int64_t>> &values, sat::BoolVar b)
{
    BVAssignIf(model, bv, values, sat::IntVar(b));
    return;
} 

void BVAssign(sat::CpModelBuilder &model, BoolVec &bv, const std::vector<std::vector<int64_t>> &values)
{
    auto iv = BV2IV(bv);
    auto table = model.AddAllowedAssignments(iv);

    for (auto &value : values) {
        table.AddTuple(value);
    }

    return;
} 

BoolVec BVRor(const BoolVec &bv, const int rotation)
{
    const int len = bv.size();
    const int rn = rotation % len;

    BoolVec output;
    for (int i = rn; i < len; ++i) {
        output.push_back(bv[i]);
    }
    for (int i = 0; i < rn; ++i) {
        output.push_back(bv[i]);
    }
    return output;
} 

BoolVec BVRol(BoolVec &bv, const int rotation)
{
    const int len = bv.size();
    const int rn = rotation % len;

    return BVRor(bv, len - rn);
} 

BoolVec BVXor(sat::CpModelBuilder &model, BoolVec &bv0, BoolVec &bv1)
{
    const int len = bv0.size();
    BoolVec bv2 = NewBoolVec(model, len);

    for (int i = 0; i < len; ++i) {
        model.AddBoolXor({bv0[i], bv1[i], bv2[i], model.TrueVar()});
    }

    return bv2;
}

void BVXor3(sat::CpModelBuilder &model, BoolVec &bv0, BoolVec &bv1, BoolVec &bv2, BoolVec &bv3)
{
    const int len = bv0.size();

    for (int i = 0; i < len; ++i) {
        model.AddBoolXor({bv0[i], bv1[i], bv2[i], bv3[i], model.TrueVar()});
    }

    return;
}

#endif