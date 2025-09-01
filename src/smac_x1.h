#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#define SMAC_VER 1  /* Select the instance: {1, 34, 12} for SMAC-{1, 3/4, 1/2} */

#define sigma \
		  (SMAC_VER == 1  ? _mm_setr_epi8(0, 7, 14, 11, 4, 13, 10, 1, 8, 15, 6, 3, 12, 5, 2, 9)\
		: (SMAC_VER == 34 ? _mm_setr_epi8(7, 14, 15, 10, 12, 13, 3, 0, 4, 6, 1, 5, 8, 11, 2, 9)\
		:                   _mm_setr_epi8(0, 11, 7, 14, 6, 4, 1, 15, 9, 3, 8, 5, 13, 2, 10, 12)))

#define const1 _mm_cvtsi32_si128(1)

#define COMPRESS1(MSG, A1, A2, A3, B1, B2, B3) /* A3 and B3 must be distinct registers */\
			A3 = _mm_ternarylogic_epi64(A3, A2, MSG, 0x96);\
			B3 = _mm_aesenc_si128(A2, MSG);\
			B2 = _mm_aesenc_si128(A1, MSG);\
			B1 = _mm_shuffle_epi8(A3, sigma)

static inline void SMAC_compress_max3_blocks(__m128i& A1, __m128i& A2, __m128i& A3, uint8_t* msg, int blocks)
{	__m128i M0, T3;
	for (int i = 0; i < blocks; i++)
	{	M0 = _mm_loadu_si128((__m128i*)msg + i);
		COMPRESS1(M0, A1, A2, A3, A1, A2, T3);
		A3 = T3;
		if (SMAC_VER == 12 || (SMAC_VER == 34 && i == 2))
		{	COMPRESS1(const1, A1, A2, A3, A1, A2, T3);
			A3 = T3;
		}
	}
}

static inline void SMAC_compress_full_blocks(__m128i &A1, __m128i& A2, __m128i& A3, uint8_t* msg, int blocks)
{
	__m128i M0, M1, M2, M3, T1, T2, T3;

#define step (SMAC_VER == 1 ? 4 : (SMAC_VER == 34 ? 3 : 2))

	for (; blocks >= step; blocks -= step, msg += step * 16)
	{	
		M0 = _mm_loadu_si128((__m128i*)(msg));
		COMPRESS1(M0, A1, A2, A3, T1, T2, T3);

		if (SMAC_VER == 12)
		{
			COMPRESS1(const1, T1, T2, T3, A1, A2, A3);
		}
		else
		{
			M1 = _mm_loadu_si128((__m128i*)(msg + 16));
			COMPRESS1(M1, T1, T2, T3, A1, A2, A3);
		}

		if (SMAC_VER == 12)
			M2 = _mm_loadu_si128((__m128i*)(msg + 16));
		else
			M2 = _mm_loadu_si128((__m128i*)(msg + 32));
		COMPRESS1(M2, A1, A2, A3, T1, T2, T3);


		if (SMAC_VER == 1)
		{
			M3 = _mm_loadu_si128((__m128i*)(msg + 48));
			COMPRESS1(M3, T1, T2, T3, A1, A2, A3);
		}
		else
		{
			COMPRESS1(const1, T1, T2, T3, A1, A2, A3);
		}
	}

	// remained unaligned blocks
	SMAC_compress_max3_blocks(A1, A2, A3, msg, blocks);
}

static inline void SMAC_InitFinal(__m128i& A1, __m128i& A2, __m128i& A3)
{
	__m128i B1, B2, B3;
	__m128i T1, T2, T3;

	B1 = A1;
	B2 = _mm_xor_si128(A2, const1);
	B3 = _mm_xor_si128(A3, const1);

	for (int i = 0; i < 8; i += 2)
	{	COMPRESS1(const1, A1, A2, A3, T1, T2, T3);
		COMPRESS1(const1, T1, T2, T3, A1, A2, A3);
	}

	// optimise the last round
	A3 = _mm_ternarylogic_epi64(A3, A2, const1, 0x96);
	T1 = _mm_aesenc_si128(A2, B3);
	A2 = _mm_aesenc_si128(A1, B2);
	A1 = _mm_shuffle_epi8(A3, sigma);
	A1 = _mm_xor_si128(A1, B1);
	A3 = T1;
}

void SMAC(uint8_t key[32], uint8_t iv[16], uint8_t* aad, int aad_sz, uint8_t* ct, int ct_sz, uint8_t* tag, int tag_sz)
{
	// Initialise with the key and iv
	__m128i A1 = _mm_loadu_si128((__m128i*)(key + 16));
	__m128i A2 = _mm_loadu_si128((__m128i*)(key));
	__m128i A3 = _mm_loadu_si128((__m128i*)(iv));
	SMAC_InitFinal(A1, A2, A3);

	// zeroise ending unaligned bytes, and add LEN-block
	memset(aad + aad_sz, 0, 16);
	memset(ct + ct_sz, 0, 16);
	int aad_blocks = (aad_sz + 15) >> 4;
	int ct_blocks = (ct_sz + 15) >> 4;
	*(uint64_t*)(ct + (ct_blocks * 16) + 0) = aad_sz * 8;
	*(uint64_t*)(ct + (ct_blocks * 16) + 8) = ct_sz * 8;
	ct_blocks++;

	// compress full blocks, including the ending LEN-block to ct
	if (SMAC_VER != 34)
	{	
		SMAC_compress_full_blocks(A1, A2, A3, aad, aad_blocks);
		SMAC_compress_full_blocks(A1, A2, A3, ct, ct_blocks);
	}
	else
	{
		int aad_rem = aad_blocks % 3;
		int aad_full = aad_blocks - aad_rem;
		uint8_t tmp[48];
		SMAC_compress_full_blocks(A1, A2, A3, aad, aad_full);
		if (aad_rem)
		{
			memcpy(tmp, aad + aad_full * 16, aad_rem * 16);
			int k = 3 - aad_rem;
			if (k > ct_blocks) k = ct_blocks;
			memcpy(tmp + aad_rem * 16, ct, k * 16);
			aad_rem += k;
			SMAC_compress_max3_blocks(A1, A2, A3, tmp, aad_rem);
			ct += k * 16;
			ct_blocks -= k;
			if(ct_blocks)
				SMAC_compress_full_blocks(A1, A2, A3, ct, ct_blocks);
		}
		else
			SMAC_compress_full_blocks(A1, A2, A3, ct, ct_blocks);
	}

	SMAC_InitFinal(A1, A2, A3);
	memcpy(tag, (uint8_t*)&A2, (tag_sz <= 16 ? tag_sz : 16));
	if (tag_sz > 16)
		memcpy(tag + 16, (uint8_t*)&A3, tag_sz - 16);

}

