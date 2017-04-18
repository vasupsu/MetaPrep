#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <emmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>
#include <math.h>
#include <assert.h>
	
__m128i mask, zeroReg, both_mask2, msb_mask2, lsb_mask2, n_mask;
__m128i zero32,msb_mask32,lsb_mask32,minus1,mask32,n_mask32;
uint32_t hiMask=0;
int kmerSize=27;
char out_bases[110];

static inline uint32_t quicklog2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}

int getCanonKmer (int startPos, uint32_t *key_h, uint32_t *key_l, uint32_t *rc_h, uint32_t *rc_l, char *in_bases, int *lastNPos)
{
	uint32_t key_l_1=0, key_l_2=0, key_h_1=0, key_h_2=0, key2_l_1=0, key2_l_2=0, key2_h_1=0, key2_h_2=0;
	int pos1=startPos;
	int pos2 = pos1 + kmerSize - 16;
	__m128i *in = (__m128i *)&in_bases[pos2];
	__m128i inReg = _mm_loadu_si128 (in);
	uint32_t lastNVal = _mm_movemask_epi8(_mm_cmpeq_epi8(inReg, n_mask));
	if (lastNVal > 0)
	{
		*lastNPos=pos2+(int)quicklog2(lastNVal);
	}
	else
		*lastNPos=-1;

	__m128i outReg = _mm_shuffle_epi8 (inReg, mask);
	__m128i msb = _mm_cmpgt_epi8(_mm_and_si128(outReg, msb_mask2), zeroReg);
	key_h_1 = _mm_movemask_epi8(msb);
	__m128i lsb = _mm_cmpgt_epi8(_mm_and_si128(outReg, lsb_mask2), zeroReg);
	key_l_1 = _mm_movemask_epi8(lsb);

	outReg = _mm_xor_si128(inReg, msb_mask2);
	msb = _mm_cmpgt_epi8(_mm_and_si128(outReg, msb_mask2), zeroReg);
	key2_h_1 = _mm_movemask_epi8(msb)<<(kmerSize - 16);
	lsb = _mm_cmpgt_epi8(_mm_and_si128(outReg, lsb_mask2), zeroReg);
	key2_l_1 = _mm_movemask_epi8(lsb)<<(kmerSize - 16);

	in = (__m128i *)&in_bases[pos1];
	inReg = _mm_loadu_si128 (in);
	lastNVal = _mm_movemask_epi8(_mm_cmpeq_epi8(inReg, n_mask)) & hiMask;
	if ((*lastNPos == -1) && (lastNVal > 0))
	{
		*lastNPos = pos1+quicklog2(lastNVal);
	}

	outReg = _mm_shuffle_epi8 (inReg, mask);
	msb = _mm_cmpgt_epi8(_mm_and_si128(outReg, msb_mask2), zeroReg);
	key_h_2 = (_mm_movemask_epi8(msb)>>(32-kmerSize))<<16;
	lsb = _mm_cmpgt_epi8(_mm_and_si128(outReg, lsb_mask2), zeroReg);
	key_l_2 = (_mm_movemask_epi8(lsb)>>(32-kmerSize))<<16;

	outReg = _mm_xor_si128(inReg, msb_mask2);
	msb = _mm_cmpgt_epi8(_mm_and_si128(outReg, msb_mask2), zeroReg);
	key2_h_2 = _mm_movemask_epi8(msb) & ((1lu<<(kmerSize - 16))-1);
	lsb = _mm_cmpgt_epi8(_mm_and_si128(outReg, lsb_mask2), zeroReg);
	key2_l_2 = _mm_movemask_epi8(lsb) & ((1lu<<(kmerSize - 16))-1);

	*key_h = key_h_1|key_h_2;
	*key_l = key_l_1|key_l_2;
	*rc_h = key2_h_1|key2_h_2;
	*rc_l = key2_l_1|key2_l_2;
	int canon=-1;
	if (abs(*key_h-*rc_h) > abs(*key_l-*rc_l))//difference lies in hsb
		canon = (*key_h<*rc_h)?0:1;
	else
		canon = (*key_l<*rc_l)?0:1;
	return canon;
}
void getFullStr (uint64_t kmer, FILE *fp)
{
	int i=0;
	for (i=26; i>=0; i--)
	{
		if (fp!=NULL)
			fprintf (fp, "%c", "ACTG"[(((kmer & 1lu<<(32+i)) >> (31+i)) & 2) | ((kmer & 1lu<<i) >> i)]);
		else
			printf ("%c", "ACTG"[(((kmer & 1lu<<(32+i)) >> (31+i)) & 2) | ((kmer & 1lu<<i) >> i)]);
	}
/*	if (fp != NULL)
		fprintf (fp, "\n");
	else
		printf ("\n");*/
}
void init(int kLen)
{
	kmerSize=kLen;
	n_mask = _mm_set_epi8(78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78);
	mask = _mm_set_epi8(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
	zeroReg = _mm_set_epi8(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
	both_mask2 = _mm_set_epi8(6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6);
	msb_mask2 = _mm_set_epi8(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4);
	lsb_mask2 = _mm_set_epi8(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);
	zero32 = _mm_set_epi32(0,0,0,0);
	msb_mask32 = _mm_set_epi32(4,4,4,4);
	lsb_mask32 = _mm_set_epi32(2,2,2,2);
	n_mask32 = _mm_set_epi32(78,78,78,78);
	minus1 = _mm_set_epi32(-1,-1,-1,-1);
	mask32= _mm_set_epi32((1<<kmerSize) -1, (1<<kmerSize) -1, (1<<kmerSize) -1, (1<<kmerSize) -1);
}
int getCanonKmers (char *readStr, int readLen, uint64_t *rangeHist, uint64_t kmerMask, int rotater, uint32_t *kmerFreqCount)
{
	int i,j;
	uint32_t kList_h[4],kList_l[4], rcList_h[4], rcList_l[4], canon[4];
	uint32_t key_h, key_l, rc_h, rc_l;
	uint32_t next4Chars[4];
	int kmerIndex=0;
	char *in_bases=readStr;
	int lastKmerPos, kmersPerLane;
	lastKmerPos = readLen - kmerSize;
	kmersPerLane = lastKmerPos/4;
	hiMask = ~(((1u << (32-kmerSize))-1)<<(kmerSize-16)) & 0xFFFF;
	uint64_t kR=0;

	int lastNPos=0;
	int earliestStartOfs[4]={-1,-1,-1,-1};
	for (i=0; i<4; i++)
	{
		canon[i]=getCanonKmer (i*kmersPerLane, &key_h, &key_l, &rc_h, &rc_l, in_bases, &earliestStartOfs[i]);
		assert (canon[i] != -1);
		kList_h[i]=key_h;
		kList_l[i]=key_l;
		rcList_h[i]=rc_h;
		rcList_l[i]=rc_l;
		int valid = i*kmersPerLane > earliestStartOfs[i];
		if (valid)
		{	
			if (canon[i]==0)
			{
				kR=(uint64_t)key_h<<32 & 0xFFFFFFFF00000000;
				kR |= key_l;
				rangeHist[((kR & kmerMask) >> rotater) + 2]++;
				kmerFreqCount[((kR & kmerMask) >> rotater)]++;
				kmerIndex++;
			}
			else
			{
				kR=(uint64_t)rc_h<<32  & 0xFFFFFFFF00000000;
                        	kR |=rc_l;
				rangeHist[((kR & kmerMask) >> rotater) + 2]++;
				kmerFreqCount[((kR & kmerMask) >> rotater)]++;
				kmerIndex++;
			}
		}
		char *readPos=&in_bases[i*kmersPerLane];
		readPos+=kmerSize;
		next4Chars[i]=*((uint32_t *)readPos);
	}
	__m128i kmerReg_h = _mm_loadu_si128((__m128i *)&kList_h[0]);
	__m128i kmerReg_l = _mm_loadu_si128((__m128i *)&kList_l[0]);
	__m128i rcReg_h = _mm_loadu_si128((__m128i *)&rcList_h[0]);
	__m128i rcReg_l = _mm_loadu_si128((__m128i *)&rcList_l[0]);
	__m128i nextChReg = _mm_loadu_si128((__m128i *)&next4Chars[0]);
	for (i=0; i<kmersPerLane-1; i++)
	{
		__m128i tmp = _mm_and_si128 (nextChReg, msb_mask32);
		__m128i msb = _mm_mullo_epi32 (_mm_cmpgt_epi32 (tmp, zeroReg), minus1);
		tmp = _mm_and_si128 (nextChReg, lsb_mask32);
		__m128i lsb = _mm_mullo_epi32 (_mm_cmpgt_epi32 (tmp, zeroReg), minus1);
		kmerReg_h = _mm_and_si128(_mm_or_si128 (_mm_slli_epi32(kmerReg_h, 1), msb),mask32);
		kmerReg_l = _mm_and_si128(_mm_or_si128 (_mm_slli_epi32(kmerReg_l, 1), lsb),mask32);
		tmp = _mm_and_si128 (_mm_xor_si128(nextChReg,msb_mask32), msb_mask32);
		msb = _mm_mullo_epi32 (_mm_cmpgt_epi32 (tmp, zeroReg), minus1);
		rcReg_h = _mm_or_si128(_mm_srli_epi32(rcReg_h, 1),_mm_slli_epi32(msb,kmerSize-1));
		tmp = _mm_and_si128 (nextChReg, lsb_mask32);
		lsb = _mm_mullo_epi32 (_mm_cmpgt_epi32 (tmp, zeroReg), minus1);
		rcReg_l = _mm_or_si128(_mm_srli_epi32(rcReg_l, 1),_mm_slli_epi32(lsb,kmerSize-1));

		kList_h[0] = _mm_extract_epi32(kmerReg_h, 0);
		kList_l[0] = _mm_extract_epi32(kmerReg_l, 0);
		rcList_h[0] = _mm_extract_epi32(rcReg_h, 0);
		rcList_l[0] = _mm_extract_epi32(rcReg_l, 0);
		kList_h[1] = _mm_extract_epi32(kmerReg_h, 1);
		kList_l[1] = _mm_extract_epi32(kmerReg_l, 1);
		rcList_h[1] = _mm_extract_epi32(rcReg_h, 1);
		rcList_l[1] = _mm_extract_epi32(rcReg_l, 1);
		kList_h[2] = _mm_extract_epi32(kmerReg_h, 2);
		kList_l[2] = _mm_extract_epi32(kmerReg_l, 2);
		rcList_h[2] = _mm_extract_epi32(rcReg_h, 2);
		rcList_l[2] = _mm_extract_epi32(rcReg_l, 2);
		kList_h[3] = _mm_extract_epi32(kmerReg_h, 3);
		kList_l[3] = _mm_extract_epi32(kmerReg_l, 3);
		rcList_h[3] = _mm_extract_epi32(rcReg_h, 3);
		rcList_l[3] = _mm_extract_epi32(rcReg_l, 3);

		for (j=0; j<4; j++)
		{
			key_h = kList_h[j];
			key_l = kList_l[j];
			rc_h = rcList_h[j];
			rc_l = rcList_l[j];
			if (abs(key_h-rc_h) > abs(key_l-rc_l))//difference lies in hsb
                		canon[j] = (key_h<rc_h)?0:1;
		        else
                		canon[j] = (key_l<rc_l)?0:1;
			if (readStr[j*kmersPerLane + i + kmerSize] == 'N')
				earliestStartOfs[j]=j*kmersPerLane + i + kmerSize;
			int valid = ((j*kmersPerLane + i + 1) > earliestStartOfs[j]);
			if (valid)
			{
				if (canon[j]==0)
		                {       
	        	        	        kR = (uint64_t)key_h<<32  & 0xFFFFFFFF00000000;
	                	        	kR |= key_l;
						rangeHist[((kR & kmerMask) >> rotater) + 2]++;
						kmerFreqCount[((kR & kmerMask) >> rotater)]++;
						kmerIndex++;
		                }
	        	        else
	                	{       
	                        		kR = (uint64_t)rc_h<<32 & 0xFFFFFFFF00000000;
		                        	kR |= rc_l;
						rangeHist[((kR & kmerMask) >> rotater) + 2]++;
						kmerFreqCount[((kR & kmerMask) >> rotater)]++;
						kmerIndex++;
	        	        }
			}
		}
		if (i%4 == 3)
		{
			for (j=0; j<4; j++)
			{
				char *readPos=&in_bases[j*kmersPerLane];
		                readPos+=kmerSize+i+1;
				next4Chars[j]=*((uint32_t *)readPos);
			}
			nextChReg = _mm_loadu_si128((__m128i *)&next4Chars[0]);
		}
		else
		{
			nextChReg=_mm_srli_epi32(nextChReg,8);
		}
	}
	for (i=4*kmersPerLane; i<=lastKmerPos; i++)
	{
		int c = getCanonKmer (i, &key_h, &key_l, &rc_h, &rc_l, in_bases, &lastNPos);
		int valid = (i > lastNPos);
		if (valid)
		{
			if (c==0)
        	        {       
        	                kR = (uint64_t)key_h<<32 & 0xFFFFFFFF00000000;
                	        kR |= key_l;
				rangeHist[((kR & kmerMask) >> rotater) + 2]++;
				kmerFreqCount[((kR & kmerMask) >> rotater)]++;
				kmerIndex++;
			}
        	        else
                	{       
                       		kR = (uint64_t)rc_h<<32 & 0xFFFFFFFF00000000;
                        	kR |= rc_l;
				rangeHist[((kR & kmerMask) >> rotater) + 2]++;
				kmerFreqCount[((kR & kmerMask) >> rotater)]++;
				kmerIndex++;
        	        }
                }
	}
	
	return kmerIndex;
}

