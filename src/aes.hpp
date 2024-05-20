#ifndef AES_HPP
#define AES_HPP

#include <array>

using gf8v = std::array<bool, 8>;
using gf8m = std::array<gf8v, 8>;

template<int bl>
using bv = std::array<bool, bl>;

template<int bl>
using bm = std::array<bv<bl>, bl>;

static constexpr gf8m C1{{
    {{1, 0, 0, 0, 0, 0, 0, 0}},
    {{0, 1, 0, 0, 0, 0, 0, 0}},
    {{0, 0, 1, 0, 0, 0, 0, 0}},
    {{0, 0, 0, 1, 0, 0, 0, 0}},
    {{0, 0, 0, 0, 1, 0, 0, 0}},
    {{0, 0, 0, 0, 0, 1, 0, 0}},
    {{0, 0, 0, 0, 0, 0, 1, 0}},
    {{0, 0, 0, 0, 0, 0, 0, 1}}
}};

static constexpr gf8m C2{{
    {{0, 1, 0, 0, 0, 0, 0, 0}},
    {{0, 0, 1, 0, 0, 0, 0, 0}},
    {{0, 0, 0, 1, 0, 0, 0, 0}},
    {{1, 0, 0, 0, 1, 0, 0, 0}},
    {{1, 0, 0, 0, 0, 1, 0, 0}},
    {{0, 0, 0, 0, 0, 0, 1, 0}},
    {{1, 0, 0, 0, 0, 0, 0, 1}},
    {{1, 0, 0, 0, 0, 0, 0, 0}}
}};

static constexpr gf8m C3{{
    {{1, 1, 0, 0, 0, 0, 0, 0}},
    {{0, 1, 1, 0, 0, 0, 0, 0}},
    {{0, 0, 1, 1, 0, 0, 0, 0}},
    {{1, 0, 0, 1, 1, 0, 0, 0}},
    {{1, 0, 0, 0, 1, 1, 0, 0}},
    {{0, 0, 0, 0, 0, 1, 1, 0}},
    {{1, 0, 0, 0, 0, 0, 1, 1}},
    {{1, 0, 0, 0, 0, 0, 0, 1}}
}};

static constexpr std::array<gf8m, 3> Ci{{ C1, C2, C3 }};

static constexpr std::array<std::array<uint8_t, 4>, 4> MCm{{
    {{2, 3, 1, 1}},
    {{1, 2, 3, 1}},
    {{1, 1, 2, 3}},
    {{3, 1, 1, 2}}
}};

static constexpr std::array<uint8_t, 16> SRp{{ 0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11 }};

static constexpr std::array<uint8_t, 256> Sbox{{
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
    0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
    0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
    0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
    0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
    0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
    0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
    0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
    0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
    0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
    0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
}};

[[maybe_unused]]
static bm<128> compute_srbm()
{
    bm<128> SRbm{{}};
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 8; ++j) {
            auto onep = SRp[i] * 8 + j;
            SRbm[i * 8 + j][onep] = 1;
        }
    }
    return SRbm;
}

[[maybe_unused]]
static bm<128> compute_mcbm()
{
    bm<128> MCbm{{}};
    for (int i = 0; i < 4; ++i) {
        for (int mi = 0; mi < 4; ++mi) {
            for (int k = 0; k < 8; ++k) {
                for (int mj = 0; mj < 4; ++mj) {
                    for (int bi = 0; bi < 8; ++bi) {
                        MCbm[i * 32 + mi * 8 + k][i * 32 + mj * 8 + bi] = Ci[MCm[mi][mj] - 1][k][bi];
                    }
                }
            }
        }
    }
    return MCbm;
}

[[maybe_unused]]
static bm<128> compute_lbm()
{
    auto SRbm = compute_srbm();
    auto MCbm = compute_mcbm();
    bm<128> Lbm{{}};

    for (int i = 0; i < 128; ++i)
        for (int j = 0; j < 128; ++j)
            for (int k = 0; k < 128; ++k)
                Lbm[i][j] ^= (MCbm[i][k] && SRbm[k][j]);
    return Lbm;
}

[[maybe_unused]]
static bm<128> perm2bm(std::array<uint8_t, 16> perm)
{
    bm<128> Pbm{{}};
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 8; ++j) {
            Pbm[i * 8 + j][perm[i] * 8 + j] = 1;
        }
    }
    return Pbm;
}

[[maybe_unused]]
static bm<128> inverse_bm(bm<128> m)
{
    bm<128> ibm{{}};
    for (int i = 0; i < 128; ++i) ibm[i][i] = 1;

    for (int col = 0; col < 128; ++col) {
        for (int pivot = col; pivot < 128; ++pivot)
            if (m[pivot][col]) {
                bv<128> tmp = m[pivot];
                m[pivot] = m[col];
                m[col] = tmp;

                tmp = ibm[pivot];
                ibm[pivot] = ibm[col];
                ibm[col] = tmp;
            }

        for (int row = 0; row < 128; ++row) {
            if (row == col) continue;
            if (m[row][col]) {
                for (int i = 0; i < 128; ++i) {
                    m[row][i] ^= m[col][i];
                    ibm[row][i] ^= ibm[col][i];
                }
            }
        }
    }
    return ibm;
}

[[maybe_unused]]
static bm<128> bmXbm(const bm<128> &a, const bm<128> &b)
{
    bm<128> r{{}};
    for (int col = 0; col < 128; ++col)
        for (int row = 0; row < 128; ++row) {
            bool tmp = 0;
            for (int i = 0; i < 128; ++i)
                tmp ^= (a[row][i] && b[i][col]);
            r[row][col] = tmp;
        }
    return r;
}

bool inner_product(const std::array<bool, 8> &a, const std::array<bool, 8> &b)
{
    bool r = 0;
    for (int i = 0; i < 8; ++i)
        r ^= (a[i] && b[i]);
    return r;
}

template<int bl>
bool inner_product(const bv<bl> &a, const bv<bl> &b)
{
    bv<bl> r{{}};
    for (int i = 0; i < bl; ++i)
        r ^= (a[i] && b[i]);
    return r;
}

std::array<bool, 8> n2v(uint8_t n)
{
    std::array<bool, 8> v{};
    for (int i = 0; i < 8; ++i)
        v[i] = (n >> (7 - i)) & 1;
    return v;
}

uint8_t v2n(std::array<bool, 8> v)
{
    uint8_t n = 0;
    for (int i = 0; i < 8; ++i)
        n = (n << 1) | v[i];
    return n;
}

template<int mi>
uint8_t mul(uint8_t x)
{
    static_assert(mi >= 1 && mi <= 3);
    auto v = n2v(x);
    uint8_t r = 0;
    for (int i = 0; i < 8; ++i) {
        auto tmp = inner_product(Ci[mi - 1][i], v);
        r = (r << 1) | tmp;
    }
    return r;
}

[[maybe_unused]]
static auto compute_ddt()
{
    std::array<std::array<unsigned short, 256>, 256> ddt{{}};
    for (int x = 0; x < 256; ++x)
        for (int dx = 0; dx < 256; ++dx) {
            const auto dy = Sbox[x] ^ Sbox[x ^ dx];
            ++ddt[dx][dy];
        }
    return ddt;
}

#endif // AES_HPP
