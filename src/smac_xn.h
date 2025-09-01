
#define ROUND512(MSG, A1, A2, A3, B1, B2, B3)\
			A3 = _mm512_ternarylogic_epi64(A3, A2, MSG, 0x96);\
			B3 = _mm512_aesenc_epi128(A2, MSG);\
			B2 = _mm512_aesenc_epi128(A1, MSG);\
			B1 = _mm512_shuffle_epi8(A3, sigma512)

#define ROUND128(MSG, A1, A2, A3, B1, B2, B3)\
			A3 = _mm_ternarylogic_epi64(A3, A2, MSG, 0x96);\
			B3 = _mm_aesenc_si128(A2, MSG);\
			B2 = _mm_aesenc_si128(A1, MSG);\
			B1 = _mm_shuffle_epi8(A3, _mm512_castsi512_si128(sigma512))

// ------------------------------------------------------------
// Align data to m 16n-byte blocks (n must be divisible by 4). 
// Change the sizes to m == the number of full blocks
// ------------------------------------------------------------
template<int n> inline void Smac1xn_align_data(uint8_t* aad, int& aad_sz, uint8_t*& ct, int& ct_sz)
{
	const __m512i zero = _mm512_setzero_si512();
	_mm_storeu_si128((__m128i*)(aad + aad_sz), _mm512_castsi512_si128(zero));
	_mm_storeu_si128((__m128i*)(ct + ct_sz), _mm512_castsi512_si128(zero));

	int k = (ct_sz + 15) >> 4;
	((uint64_t*)(ct + k * 16))[0] = aad_sz << 3;
	((uint64_t*)(ct + k * 16))[1] = ct_sz << 3;
	aad_sz = (aad_sz + 15) >> 4;
	ct_sz = k + 1;

	k = aad_sz % n;
	if (k)
	{
		int t = n - k, flag = t > ct_sz;
		if (flag) t = ct_sz;

		if (n == 4)
			_mm512_storeu_si512(aad + aad_sz * 16 + 0, _mm512_loadu_si512(ct + 0));
		else if (n == 8)
		{
			_mm512_storeu_si512(aad + aad_sz * 16 + 0, _mm512_loadu_si512(ct + 0));
			_mm512_storeu_si512(aad + aad_sz * 16 + 64, _mm512_loadu_si512(ct + 64));
		}
		else
			memcpy(aad + aad_sz * 16, ct, t * 16);

		ct += t * 16;
		aad_sz += t;
		ct_sz -= t;

		if (flag)
			if (n == 4)
				_mm512_storeu_si512(aad + aad_sz * 16 + 0, zero);
			else if (n == 8)
			{
				_mm512_storeu_si512(aad + aad_sz * 16 + 0, zero);
				_mm512_storeu_si512(aad + aad_sz * 16 + 64, zero);
			}
			else
				memset(aad + aad_sz * 16, 0, n * 16);
	}

	if (n == 4)
		_mm512_storeu_si512(ct + ct_sz * 16 + 0, zero);
	else if (n == 8)
	{
		_mm512_storeu_si512(ct + ct_sz * 16 + 0, zero);
		_mm512_storeu_si512(ct + ct_sz * 16 + 64, zero);
	}
	else
		memset(ct + ct_sz * 16, 0, n * 16);

	aad_sz = (aad_sz + n - 1) / n;
	ct_sz = (ct_sz + n - 1) / n;
}


// =======================================================
// SMAC-1x4
// =======================================================
inline void Smac1x4_compress(__m512i& _A1, __m512i& _A2, __m512i& _A3, __m512i* msg, int blocks)
{
	__m512i A1 = _A1, A2 = _A2, A3 = _A3;
	const __m512i sigma512 = _mm512_broadcast_i32x4(_mm_setr_epi8(0, 7, 14, 11, 4, 13, 10, 1, 8, 15, 6, 3, 12, 5, 2, 9));
	int k;

	for (k = 0; k <= (blocks - 4); k += 4)
	{
		__m512i T1, T2, T3, M0, M1, M2, M3;
		M0 = _mm512_loadu_si512(msg + k + 0); ROUND512(M0, A1, A2, A3, T1, T2, T3);
		M1 = _mm512_loadu_si512(msg + k + 1); ROUND512(M1, T1, T2, T3, A1, A2, A3);
		M2 = _mm512_loadu_si512(msg + k + 2); ROUND512(M2, A1, A2, A3, T1, T2, T3);
		M3 = _mm512_loadu_si512(msg + k + 3); ROUND512(M3, T1, T2, T3, A1, A2, A3);
	}

	if (blocks & 2)
	{
		__m512i T1, T2, T3, M0, M1;
		M0 = _mm512_loadu_si512(msg + k + 0); ROUND512(M0, A1, A2, A3, T1, T2, T3);
		M1 = _mm512_loadu_si512(msg + k + 1); ROUND512(M1, T1, T2, T3, A1, A2, A3);
		k += 2;
	}

	if (blocks & 1)
	{
		__m512i M0 = _mm512_loadu_si512(msg + k + 0); ROUND512(M0, A1, A2, A3, _A1, _A2, _A3);
	}
	else
	{
		_A1 = A1, _A2 = A2, _A3 = A3;
	}
}


// ct[] and aad[] must reserve 128 bytes
void Smac1x4(uint8_t key[32], uint8_t iv[16], uint8_t* aad, int aad_sz, uint8_t* ct, int ct_sz, uint8_t* tag, int tag_sz)
{
	// Key-IV loading
	__m512i A1, A2, A3, A1r, A2r, A3r;
	iv[15] = (4 - 1) * 16;
	A3r = _mm512_broadcast_i32x4(_mm_loadu_si128((__m128i*)(iv)));
	A1r = A1 = _mm512_broadcast_i32x4(_mm_loadu_si128((__m128i*)(key + 16)));
	A2r = A2 = _mm512_broadcast_i32x4(_mm_loadu_si128((__m128i*)(key)));
	A3 = _mm512_xor_si512(A3r, _mm512_setr_epi64(0, 0ULL << 56, 0, 1ULL << 56, 0, 2ULL << 56, 0, 3ULL << 56));

	// Initialisation, 8 dummy clocks
	const __m512i sigma512 = _mm512_broadcast_i32x4(_mm_setr_epi8(0, 7, 14, 11, 4, 13, 10, 1, 8, 15, 6, 3, 12, 5, 2, 9));
	const __m512i const512 = _mm512_broadcast_i32x4(_mm_cvtsi32_si128(1));

	for (int i = 0; i < 8; i += 2)
	{
		__m512i T1, T2, T3;
		ROUND512(const512, A1, A2, A3, T1, T2, T3);
		ROUND512(const512, T1, T2, T3, A1, A2, A3);
	}
	__m512i T3;
	ROUND512(const512, A1, A2, A3, A1, A2, T3);
	A3 = T3;

	// Prepare AAD and CT blocks, add Lengths, padding and align
	Smac1xn_align_data<4>(aad, aad_sz, ct, ct_sz);

	// PRP-PRF switch
	A1 = _mm512_xor_si512(A1, A1r);
	A2 = _mm512_xor_si512(A2, A2r);
	A3 = _mm512_xor_si512(A3, A3r);

	// compress AAD & CT in 4x16 blocks
	Smac1x4_compress(A1, A2, A3, (__m512i*)aad, aad_sz);
	Smac1x4_compress(A1, A2, A3, (__m512i*)ct, ct_sz);

	// Finalisation 1
	for (int i = 0; i < 6; i += 2)
	{
		__m512i T1, T2, T3;
		ROUND512(const512, A1, A2, A3, T1, T2, T3);
		ROUND512(const512, T1, T2, T3, A1, A2, A3);
	}

	// XOR
	__m128i A1d = _mm_xor_si128(_mm512_castsi512_si128(A1), _mm512_extracti64x2_epi64(A1, 1));
	__m128i A2d = _mm_xor_si128(_mm512_castsi512_si128(A2), _mm512_extracti64x2_epi64(A2, 1));
	__m128i A3d = _mm_xor_si128(_mm512_castsi512_si128(A3), _mm512_extracti64x2_epi64(A3, 1));
	A1d = _mm_ternarylogic_epi64(A1d, _mm512_extracti64x2_epi64(A1, 2), _mm512_extracti64x2_epi64(A1, 3), 0x96);
	A2d = _mm_ternarylogic_epi64(A2d, _mm512_extracti64x2_epi64(A2, 2), _mm512_extracti64x2_epi64(A2, 3), 0x96);
	A3d = _mm_ternarylogic_epi64(A3d, _mm512_extracti64x2_epi64(A3, 2), _mm512_extracti64x2_epi64(A3, 3), 0x96);

	// Finalisation 2
	__m128i A2e = A2d, T1d, T2d, T3d;
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, A1d, A2d, T3d); A3d = T3d;
	A2d = _mm_xor_si128(A2d, A2e);

	// Tag output
	memcpy(tag, &A2d, tag_sz);
}


// =======================================================
// SMAC-1x8
// =======================================================
inline void Smac1x8_compress(__m512i& _A1a, __m512i& _A1b, __m512i& _A2a, __m512i& _A2b, __m512i& _A3a, __m512i& _A3b, __m512i* msg, int blocks)
{
	__m512i A1a = _A1a, A1b = _A1b, A2a = _A2a, A2b = _A2b, A3a = _A3a, A3b = _A3b;
	const __m512i sigma512 = _mm512_broadcast_i32x4(_mm_setr_epi8(0, 7, 14, 11, 4, 13, 10, 1, 8, 15, 6, 3, 12, 5, 2, 9));
	int k;
	for (k = 0; k <= 2 * (blocks - 2); k += 4)
	{
		__m512i T1a, T2a, T3a, T1b, T2b, T3b, M0, M1, M2, M3;
		M0 = _mm512_loadu_si512(msg + k + 0); ROUND512(M0, A1a, A2a, A3a, T1a, T2a, T3a);
		M1 = _mm512_loadu_si512(msg + k + 1); ROUND512(M1, A1b, A2b, A3b, T1b, T2b, T3b);
		M2 = _mm512_loadu_si512(msg + k + 2); ROUND512(M2, T1a, T2a, T3a, A1a, A2a, A3a);
		M3 = _mm512_loadu_si512(msg + k + 3); ROUND512(M3, T1b, T2b, T3b, A1b, A2b, A3b);
	}

	if (blocks & 1)
	{
		__m512i M0, M1;
		M0 = _mm512_loadu_si512(msg + k + 0); ROUND512(M0, A1a, A2a, A3a, _A1a, _A2a, _A3a);
		M1 = _mm512_loadu_si512(msg + k + 1); ROUND512(M1, A1b, A2b, A3b, _A1b, _A2b, _A3b);
	}
	else
	{
		_A1a = A1a, _A1b = A1b, _A2a = A2a, _A2b = A2b, _A3a = A3a, _A3b = A3b;
	}
}


// ct[] and aad[] must reserve 256 bytes
void Smac1x8(uint8_t key[32], uint8_t iv[16], uint8_t* aad, int aad_sz, uint8_t* ct, int ct_sz, uint8_t* tag, int tag_sz)
{
	// Key-IV loading
	__m512i A1a, A1b, A2a, A2b, A3a, A3b, A1r, A2r, A3r;
	iv[15] = (8 - 1) * 16;
	A3r = _mm512_broadcast_i32x4(_mm_loadu_si128((__m128i*)(iv)));
	A1r = A1a = A1b = _mm512_broadcast_i32x4(_mm_loadu_si128((__m128i*)(key + 16)));
	A2r = A2a = A2b = _mm512_broadcast_i32x4(_mm_loadu_si128((__m128i*)(key)));
	A3a = _mm512_xor_si512(A3r, _mm512_setr_epi64(0, 0ULL << 56, 0, 1ULL << 56, 0, 2ULL << 56, 0, 3ULL << 56));
	A3b = _mm512_xor_si512(A3r, _mm512_setr_epi64(0, 4ULL << 56, 0, 5ULL << 56, 0, 6ULL << 56, 0, 7ULL << 56));

	// Initialisation, 9 dummy clocks
	const __m512i sigma512 = _mm512_broadcast_i32x4(_mm_setr_epi8(0, 7, 14, 11, 4, 13, 10, 1, 8, 15, 6, 3, 12, 5, 2, 9));
	const __m512i const512 = _mm512_broadcast_i32x4(_mm_cvtsi32_si128(1));

#if 1
	for (int i = 0; i < 8; i += 2)
	{
		__m512i T1a, T2a, T3a, T1b, T2b, T3b;
		ROUND512(const512, A1a, A2a, A3a, T1a, T2a, T3a);
		ROUND512(const512, A1b, A2b, A3b, T1b, T2b, T3b);
		ROUND512(const512, T1a, T2a, T3a, A1a, A2a, A3a);
		ROUND512(const512, T1b, T2b, T3b, A1b, A2b, A3b);
	}
	__m512i T3a, T3b;
	ROUND512(const512, A1a, A2a, A3a, A1a, A2a, T3a);
	ROUND512(const512, A1b, A2b, A3b, A1b, A2b, T3b);
	A3a = T3a, A3b = T3b;

	// Prepare AAD and CT blocks, add Lengths, padding and align
#endif
	Smac1xn_align_data<8>(aad, aad_sz, ct, ct_sz);

#if 1
	// PRP-PRF switch
	A1a = _mm512_xor_si512(A1a, A1r);
	A1b = _mm512_xor_si512(A1b, A1r);
	A2a = _mm512_xor_si512(A2a, A2r);
	A2b = _mm512_xor_si512(A2b, A2r);
	A3a = _mm512_xor_si512(A3a, A3r);
	A3b = _mm512_xor_si512(A3b, A3r);
	
	// compress AAD & CT in 8x16 blocks
	Smac1x8_compress(A1a, A1b, A2a, A2b, A3a, A3b, (__m512i*)aad, aad_sz);
#endif
	Smac1x8_compress(A1a, A1b, A2a, A2b, A3a, A3b, (__m512i*)ct, ct_sz);

#if 1
	// Finalisation 1
	for (int i = 0; i < 6; i += 2)
	{	__m512i T1a, T2a, T3a, T1b, T2b, T3b;
		ROUND512(const512, A1a, A2a, A3a, T1a, T2a, T3a);
		ROUND512(const512, A1b, A2b, A3b, T1b, T2b, T3b);
		ROUND512(const512, T1a, T2a, T3a, A1a, A2a, A3a);
		ROUND512(const512, T1b, T2b, T3b, A1b, A2b, A3b);
	}

	// XOR
	A1a = _mm512_xor_si512(A1a, A1b);
	A2a = _mm512_xor_si512(A2a, A2b);
	A3a = _mm512_xor_si512(A3a, A3b);
	__m128i A1d = _mm_xor_si128(_mm512_castsi512_si128(A1a), _mm512_extracti64x2_epi64(A1a, 1));
	__m128i A2d = _mm_xor_si128(_mm512_castsi512_si128(A2a), _mm512_extracti64x2_epi64(A2a, 1));
	__m128i A3d = _mm_xor_si128(_mm512_castsi512_si128(A3a), _mm512_extracti64x2_epi64(A3a, 1));
	A1d = _mm_ternarylogic_epi64(A1d, _mm512_extracti64x2_epi64(A1a, 2), _mm512_extracti64x2_epi64(A1a, 3), 0x96);
	A2d = _mm_ternarylogic_epi64(A2d, _mm512_extracti64x2_epi64(A2a, 2), _mm512_extracti64x2_epi64(A2a, 3), 0x96);
	A3d = _mm_ternarylogic_epi64(A3d, _mm512_extracti64x2_epi64(A3a, 2), _mm512_extracti64x2_epi64(A3a, 3), 0x96);

	// Finalisation 2
	__m128i A2e = A2d, T1d, T2d, T3d;
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, T1d, T2d, T3d);
	ROUND128(_mm512_castsi512_si128(const512), T1d, T2d, T3d, A1d, A2d, A3d);
	ROUND128(_mm512_castsi512_si128(const512), A1d, A2d, A3d, A1d, A2d, T3d); A3d = T3d;
	A2d = _mm_xor_si128(A2d, A2e);

	// Tag output
	memcpy(tag, &A2d, tag_sz);
#endif
	memcpy(tag, &A2a, tag_sz);
}

