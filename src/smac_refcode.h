#define IS_SMAC1 1  /* Select the instance: SMAC-1 or SMAC-3/4 */
#define SIGMA (IS_SMAC1\
	? _mm_setr_epi8(0,7,14,11,4,13,10,1,8,15,6,3,12,5,2,9)\
	: _mm_setr_epi8(7,14,15,10,12,13,3,0,4,6,1,5,8,11,2,9))

#define load(ptr)     _mm_loadu_si128((__m128i*)(ptr))
#define store(ptr, x) _mm_storeu_si128((__m128i*)(ptr), x)
#define aes(a, k)     _mm_aesenc_si128(a, k)
#define sigma(x)      _mm_shuffle_epi8(x, SIGMA)
#define xor2(x, y)    _mm_xor_si128(x, y)
#define xor3(x, y, z)  xor2(xor2(x,y),z)

void SMAC_Compress(__m128i& A1, __m128i& A2, __m128i& A3, uint8_t* msg)
{	__m128i M = msg ? load(msg) : _mm_cvtsi32_si128(1);
	__m128i T = sigma(xor3(A2, A3, M));
	A3 = aes(A2, M);
	A2 = aes(A1, M);
	A1 = T;
}

void SMAC_InitFinal(__m128i& A1, __m128i& A2, __m128i& A3)
{	__m128i T1 = A1, T2 = A2, T3 = A3;
	for (int i = 0; i < 9; i++)
		SMAC_Compress(A1, A2, A3, NULL);
	A1 = xor2(A1, T1);
	A2 = xor2(A2, T2);
	A3 = xor2(A3, T3);
}

// (!) In this implementation, aad and ct must reserve 16 extra bytes
void SMAC(uint8_t key[32], uint8_t iv[16], uint8_t* aad, int aad_sz, uint8_t * ct, int ct_sz, uint8_t * tag, int tag_sz)
{	// initialise with the key and iv
	__m128i A1 = load(key + 16), A2 = load(key), A3 = load(iv);	
	SMAC_InitFinal(A1, A2, A3);
 
	// zeroise ending unaligned bytes, and add LEN-block to ct
	memset(aad + aad_sz, 0, 16);
	memset(ct + ct_sz, 0, 16);
	int aad_blocks = (aad_sz + 15) >> 4;
	int ct_blocks = (ct_sz + 15) >> 4;
	*(uint64_t*)(ct + (ct_blocks * 16) + 0) = aad_sz * 8;
	*(uint64_t*)(ct + (ct_blocks * 16) + 8) = ct_sz * 8;

	// compress full blocks, including the ending LEN-block
	for (int i = 0; i <= (aad_blocks + ct_blocks); i++)
	{	uint8_t* msg = i < aad_blocks ? (aad + i * 16) 
                        : (ct + (i - aad_blocks) * 16);
		SMAC_Compress(A1, A2, A3, msg);
		if (!IS_SMAC1 && (i % 3) == 2)
			SMAC_Compress(A1, A2, A3, NULL);
	}

	// finalise and derive the MAC value
	SMAC_InitFinal(A1, A2, A3);
	memcpy(tag, (uint8_t*)&A2, (tag_sz <= 16 ? tag_sz : 16));
	if (tag_sz > 16)
		memcpy(tag + 16, (uint8_t*)&A3, tag_sz - 16);
}