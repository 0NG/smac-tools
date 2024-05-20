#include <iostream>
#include <immintrin.h>
#include <cstdint>
#include <cstring>

using std::cout;
using std::endl;

#include "smac_refcode.h"
#include "smac_testvectors.h"

int main()
{
    cout << "This is a demo of SMAC-1." << endl;

    int aad_sz[4] = {0, 1, 16, 19};
    int ct_sz[4]  = {0, 0, 16, 13};
    for (int i = 0; i < 4; ++i) {
        uint8_t tag[32];
        int tag_sz = 32;
        SMAC(keys[i], iv[i], aad[i], aad_sz[i], cipher[i], ct_sz[i], tag, tag_sz);

        for (int j = 0; j < tag_sz; ++j) {
            if (j < 16) {
                if (tag[j] != a2[i][j])
                    cout << "[wrong demo]" << endl;
            } else {
                if (tag[j] != a3[i][j - 16])
                    cout << "[wrong demo]" << endl;
            }
        }
    }

    cout << "[passed]" << endl;
    return 0;
}
