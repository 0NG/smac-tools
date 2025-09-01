#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <memory.h>
#include <time.h>
//#include <intrin.h>
#include <immintrin.h>


#include "smac_x1.h"
#include "smac_xn.h"


int main(int argn, char * argc[])
{
	uint8_t key[32] = { 1 }, iv[16] = { 2 }, MAC[16], aad[512];

	int msg_sz = (1 << 18);
	if (argn > 1) msg_sz = atoi(argc[1]);
	printf("TEST: VER=%d bytes=%d\n", SMAC_VER, msg_sz);

	//uint8_t* msg = (uint8_t*)_aligned_malloc(msg_sz + 512, 256);
	uint8_t* msg = (uint8_t*)aligned_alloc(256, msg_sz + 512);
	for (int i = 0; i < (msg_sz); i++)
		msg[i] = rand();
	
	long long timeout_sec = 1;
	long long tm = time(NULL);
	
	while (1)
	{
		tm += 1;
		double count = 0;
		while (time(NULL) < tm);
		tm += timeout_sec;
		while (time(NULL) < tm)
		{
			//SMAC(key, iv, aad, 0, msg, msg_sz, MAC, 16);
			Smac1x8(key, iv, aad, 0, msg, msg_sz, MAC, 16);
			//Smac1x4(key, iv, aad, 0, msg, msg_sz, MAC, 16);
			*(uint64_t*)msg ^= *(uint64_t*)MAC;
			count += 1;
		}

		double bits = count * (double)msg_sz * 8.;
		double gbps = bits / 1000. / 1000. / 1000. / (double)timeout_sec;
		printf("dummy = %d, Gbps = %lf (count=%.0lf)\n", (int)MAC[0], gbps, count);
	}
}
