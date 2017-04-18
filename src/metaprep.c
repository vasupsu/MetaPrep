#define _GNU_SOURCE
#include <sched.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/stat.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>
#define USE_MPI 1
#if USE_MPI
#include <mpi.h>
#endif
#define USE_OMP 1
//Kmers or Reads to be sent in 1 call
#define ASYNC_SIZE 400000000
//Multi Scan optimization
#define SCAN_OPT  1
#if USE_OMP
#include <omp.h>
#endif
#include "unionfind.h"

typedef struct {
        size_t offset;
        size_t offset2;
        uint32_t readId;
        int fileNo;
        uint32_t kmerFreqCount[1048576+2];
}fileIndex;

typedef struct {
	uint32_t **sendKmers;
	size_t *sendThreadOfs;
	size_t *sendThreadCount;
	size_t *threadGenerateCount;
	size_t *rcvThreadOfs;
	size_t *rcvThreadCount;
}asyncData;

char **fileNames=NULL;
size_t *fileSize=NULL, totalFastqSize=0;
size_t fastqBufferSize=500000000;
char *fastqBuffer1=NULL, *fastqBuffer2=NULL;

int numPEFiles, numSEFiles, numInterleavedFiles, K=0;
int numTasks=1, numThreads, rank, passes;
int totalFastqPartitions=0;
fileIndex *fastqPart = NULL;
uint64_t *merHist;
uint64_t *scanKmerRanges=NULL, *processKmerRanges=NULL, *threadKmerRanges=NULL;
int *fastqRanges=NULL;
size_t *recvThreadCount=NULL;
size_t *recvThreadOfs=NULL;
uint64_t numKmers=0;
uint32_t numReads=0;
int hwThreads=0;
int numCommIters=1;
long *threadKmerGenTime=NULL;
FILE **openFp=NULL;//2*numThreads
int *openFileNo=NULL;

uint32_t *kmerDegree=NULL;
uint32_t *parent=NULL;
uint32_t *compSizes=NULL;
long kmerGenWorkTime=0, kmerGenTime=0, commTime=0, memcpyTime=0, sortTime1=0, sortTime2=0, unionFindTime1=0, unionFindTime2=0, mergeTime=0, outputTime=0;

extern void getCanonKmers (char *readStr, uint32_t **kR, uint64_t *curSendOfs, int readLen, uint64_t startKmer, uint64_t endKmer, uint64_t kmerMask, int rotater, uint32_t curReadId, size_t *lastKmerAddedByThread, uint64_t *scanProcessKmerRange, int tid, int rank);
extern void init(int kLen, int numP, int numThreads);
size_t *sendSizes=NULL;
size_t *sendThreadSizes=NULL;
uint64_t kmerMask=0;
int rotater=0;

char oPrefix[100];

uint32_t *kmerIn=NULL, *kmerOut=NULL, *kmerSorted=NULL;
void getKmerMask(uint32_t k)
{
      kmerMask = (((1ul<<k)-1) >> (k-20)) << (32+k-20);
      rotater = 32+k-20;
}

static int hardware_threads(void)
{
        char name[40];
        struct stat st;
        int cpus = -1;
        do {
                sprintf(name, "/sys/devices/system/cpu/cpu%d", ++cpus);
        } while (stat(name, &st) == 0);
        return cpus;
}


static int parentCompare(const void *p1, const void *p2)
{
        uint32_t f1=*(uint32_t *)p1;
        uint32_t f2=*(uint32_t *)p2;
        if (f1 < f2)
                return -1;
        else if (f1 > f2)
                return 1;
        else
                return 0;
}

void bufSeek (int bufNo, FILE *fp, long ofs, size_t *fqSize, size_t *fqInd, char **fqBuffer)
{
        fqSize[bufNo]=fqInd[bufNo]=0;
        assert (fseek(fp, ofs, SEEK_SET) == 0);

        fqSize[bufNo]=fread(fqBuffer[bufNo], 1, fastqBufferSize, fp);
//	printf ("[%d] fqSize[%d]=%lu\n", rank, bufNo, fqSize[bufNo]);
}

char * getNextLine(int bufNo, FILE *fp, int *len, int print, size_t *fqSize, size_t *fqInd, char **fqBuffer, int *eof, size_t *ofs, int f1)
{
        size_t curInd = fqInd[bufNo];
        if (bufNo == 0)
                *eof=0;
        while (fqBuffer[bufNo][curInd] != '\n')
        {
                if (curInd == fqSize[bufNo])
                {
                        memcpy (fqBuffer[bufNo], &fqBuffer[bufNo][fqInd[bufNo]], curInd-fqInd[bufNo]);
                        fqSize[bufNo] = curInd-fqInd[bufNo];
                        fqInd[bufNo] = 0;
                        curInd = fqSize[bufNo];
                        if (!feof(fp))
                        {
                                size_t nRead = fread(&fqBuffer[bufNo][fqSize[bufNo]], 1, fastqBufferSize-fqSize[bufNo], fp);
                                fqSize[bufNo] += nRead;
                        }
                        else
                        {
                                *eof=1;
                        }
                }
                if (*eof) break;
                curInd++;
        }
        if ((*eof==0) && (fqInd[bufNo] < fqSize[bufNo]))
        {
                *(fqBuffer[bufNo]+curInd) = '\0';
                *len=curInd-fqInd[bufNo];
                char *retPtr = fqBuffer[bufNo]+fqInd[bufNo];
                if (bufNo == 0)
                {
                        *ofs += (curInd+1-fqInd[bufNo]);
                        if (*ofs == fileSize[f1])
                        {
                                *eof=1;
                        }
                }
                fqInd[bufNo]=curInd+1;
                return retPtr;
        }
        return NULL;
}
//return type - 0 (paired end), 1 - single end , 2 - interleaved, -1 - EO chunk
int getNextFullRead(size_t *fqInd, size_t *fqSize, char **fqBuffer, char **line1, char **line2, int *len1, int *len2, int interleave, size_t interleaveStart)
{
        int ret=1;
        *line1 = NULL; *line2=NULL;
        *len1 = 0; *len2 = 0;
        if (fqInd[0] >= fqSize[0])
                return -1;
        *line1 = &(fqBuffer[0][fqInd[0]]);
	int oldFqInd=fqInd[0];

	if (fqInd[0] < fqSize[0])
        {
		while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                fqInd[0]++;
		fqBuffer[0][fqInd[0]-3]='/';
		fqBuffer[0][fqInd[0]-2]='1';
	        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
	        {
        	        fqInd[0]++;
                	*len1 = *len1 +1;
	        }
		fqInd[0]++;
		while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
        	fqInd[0]++;
		fqInd[0]+=*len1;
        	fqInd[0]++;	
		*len1 = fqInd[0]-oldFqInd;
		if (interleave && (fqInd[0] >= interleaveStart))
                {
			oldFqInd=fqInd[0];
			*line2 = &(fqBuffer[0][fqInd[0]]);
                        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                        fqInd[0]++;
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
                        {
                                fqInd[0]++;
                                *len2 = *len2 +1;
                        }
			fqInd[0]++;
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                        fqInd[0]++;
			fqInd[0]+=*len2;
			fqInd[0]++;
			*len2 = fqInd[0]-oldFqInd;
			ret=2;
			if (fqInd[0] < fqSize[0])
                        {
			}
                }
        }
	if (ret==1)
        {
                if (fqInd[1] < fqSize[1])
                {
			oldFqInd=fqInd[1];
                        *line2 = &(fqBuffer[1][fqInd[1]]);
                        while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                        fqInd[1]++;
			fqBuffer[1][fqInd[1]-3]='/';
			fqBuffer[1][fqInd[1]-2]='2';
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n'))
                        {
                                fqInd[1]++;
                                *len2 = *len2 +1;
                        }
			ret=0;
			fqInd[1]++;
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                        fqInd[1]++;
			fqInd[1]+=*len2;
			fqInd[1]++;
			*len2 = fqInd[1]-oldFqInd;
			if (fqInd[1] < fqSize[1])
                        {
			}
                }
        }
        return ret;
}
//return type - 0 (paired end), 1 - single end , 2 - interleaved, -1 - EO chunk
int getNextRead(size_t *fqInd, size_t *fqSize, char **fqBuffer, char **line1, char **line2, int *len1, int *len2, int interleave, size_t interleaveStart)
{
        int ret=1;
        *line1 = NULL; *line2=NULL;
        *len1 = 0; *len2 = 0;
        if (fqInd[0] >= fqSize[0])
                return -1;
        *line1 = &(fqBuffer[0][fqInd[0]]);

        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
        {
                fqInd[0]++;
                *len1 = *len1 +1;
        }
	fqInd[0]++;
	while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
        fqInd[0]++;
	fqInd[0]+=*len1;
        fqInd[0]++;	
	if (fqInd[0] < fqSize[0])
        {
		while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                fqInd[0]++;
		if (interleave && (fqInd[0] >= interleaveStart))
                {
			*line2 = &(fqBuffer[0][fqInd[0]]);
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n'))
                        {
                                fqInd[0]++;
                                *len2 = *len2 +1;
                        }
			fqInd[0]++;
			while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                        fqInd[0]++;
			fqInd[0]+=*len2;
			fqInd[0]++;
			ret=2;
			if (fqInd[0] < fqSize[0])
                        {
                                while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
                                fqInd[0]++;
			}
                }
        }
	if (ret==1)
        {
                if (fqInd[1] < fqSize[1])
                {
                        *line2 = &(fqBuffer[1][fqInd[1]]);
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n'))
                        {
                                fqInd[1]++;
                                *len2 = *len2 +1;
                        }
			ret=0;
			fqInd[1]++;
			while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                        fqInd[1]++;
			fqInd[1]+=*len2;
			fqInd[1]++;
			if (fqInd[1] < fqSize[1])
                        {
                                while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                                fqInd[1]++;
			}
                }
        }
        return ret;
}
	
void loadChunk (int chunk, size_t *fqInd, size_t *fqSize, char **fqBuffer, int *interleave, size_t *interleaveStart, int tid)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	FILE *fp1=NULL, *fp2=NULL;
        fqInd[0] = fqInd[1]=0;
        fqSize[0] = fqSize[1]=0;
        int curFile=fastqPart[chunk].fileNo;
        size_t curOfs1 = fastqPart[chunk].offset;
        size_t curOfs2 = fastqPart[chunk].offset2;
        int endFile=fastqPart[chunk+1].fileNo;
        int startFile = curFile;
        size_t startOfs1 = curOfs1;
        size_t startOfs2 = curOfs2;

        size_t endOfs1 = fastqPart[chunk+1].offset;
        size_t endOfs2 = fastqPart[chunk+1].offset2;

        size_t bytesToRead1=0;
        size_t bufferBytesRemaining1 = fastqBufferSize;
        size_t bufferOfs1=0;

        size_t bytesToRead2=0;
        size_t bufferBytesRemaining2 = fastqBufferSize;
        size_t bufferOfs2=0;
        while ((curFile < endFile) || ((curFile == endFile) && (curOfs1 < endOfs1)))
        {
		if (curFile >= (numPEFiles + numSEFiles))//Interleaved
                {
	                if (*interleave != 1)
        	        {
        	        	*interleave=1;
	        	        *interleaveStart = bufferOfs1;
        	        }
                }
                fp1 = fopen(fileNames[curFile], "r");
		if (fp1==NULL) 
		{
			printf ("Cant open file %s\n", fileNames[curFile]);
			int pid = getpid();
			char line[2048];
			sprintf (line, "/proc/%d/status", pid);
			FILE *statFile = fopen(line, "r");
			assert (statFile != NULL);
			fgets (line, 2048, statFile);
			while (!feof (statFile))
			{
				if (strstr(line,"VmPeak") || strstr(line,"VmHWM"))
				{
					printf ("[%d] %s", rank, line);
				}
				fgets (line, 2048, statFile);
			}
			fclose (statFile);
			printf ("[%d] tid %d fp1 NULL, fileNo %d\n", rank, tid, curFile);
			sleep(60);
		}
                assert (fp1 != NULL);
                assert (fseek (fp1, curOfs1, SEEK_SET) == 0);
                fp2 = NULL;
                if (curFile < numPEFiles)
                {
	                fp2 = fopen(fileNames[curFile+1], "r");
                        assert (fp2 != NULL);
                        assert (fseek (fp2, curOfs2, SEEK_SET) == 0);
                }
                if (curFile < endFile)
                {
                        bytesToRead1 = fileSize[curFile]-curOfs1;
                        if (curFile < numPEFiles)
                                bytesToRead2 = fileSize[curFile+1]-curOfs2;
                }
                else
                {
                        bytesToRead1 = endOfs1-curOfs1;
                        if (curFile < numPEFiles)
                                bytesToRead2 = endOfs2-curOfs2;
                }
                assert (bufferBytesRemaining1 > bytesToRead1);
                assert (fread (&(fqBuffer[0][bufferOfs1]), 1, bytesToRead1, fp1) == bytesToRead1);
                bufferBytesRemaining1 -= bytesToRead1;
                bufferOfs1 += bytesToRead1;
		if (curFile < numPEFiles)
                {
                        assert (bufferBytesRemaining2 > bytesToRead2);
                        assert (fread (&(fqBuffer[1][bufferOfs2]), 1, bytesToRead2, fp2) == bytesToRead2);
                        bufferBytesRemaining2 -= bytesToRead2;
                        bufferOfs2 += bytesToRead2;
                }

                fclose (fp1);
                if (fp2 != NULL)
                        fclose (fp2);
                if (curFile < numPEFiles) curFile++;
                curFile++;
                curOfs1=curOfs2=0;
        }
        fqSize[0]=bufferOfs1;
        fqSize[1]=bufferOfs2;
	gettimeofday (&eTime, NULL);
}
uint32_t getReadAndKmerCount (int chunk, size_t *fqInd, size_t *fqSize, char **fqBuffer, int interleave, size_t interleaveStart, uint32_t **sendKmers, uint64_t *curSendOfs, uint64_t startKmer, uint64_t endKmer, uint64_t kmerMask, int rotater, uint32_t curReadId, size_t *lastKmerAddedByThread, uint64_t *scanProcessKmerRange, int tid)
{
        int len1=0, len2=0;
        char *read1=NULL, *read2=NULL;
        uint32_t numPEReads=0, numSEReads=0, numInterleavedReads=0;
        int ret = getNextRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
#if SCAN_OPT
	if (parent==NULL)
	{
#endif
	        while (ret != -1)
        	{
        	       	if (read1 != NULL)
    				getCanonKmers (read1, sendKmers, curSendOfs, len1, startKmer, endKmer, kmerMask, rotater, curReadId, lastKmerAddedByThread, scanProcessKmerRange, tid, rank);
                	if (read2 != NULL)
	    			getCanonKmers (read2, sendKmers, curSendOfs, len2, startKmer, endKmer, kmerMask, rotater, curReadId, lastKmerAddedByThread, scanProcessKmerRange, tid, rank);
                	if (ret == 0)
	                {
                	      numPEReads++;
	                }
        	        if (ret == 1) numSEReads++;
                	if (ret == 2) numInterleavedReads++;
			curReadId++;
        	        ret = getNextRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
	        }
#if SCAN_OPT
	}
	else
	{
		int numRandAccesses=0;
	        while (ret != -1)
        	{
                	if (read1 != NULL)
				getCanonKmers (read1, sendKmers, curSendOfs, len1, startKmer, endKmer, kmerMask, rotater, findRoot(parent, curReadId, &numRandAccesses), lastKmerAddedByThread, scanProcessKmerRange, tid, rank);
	                if (read2 != NULL)
				getCanonKmers (read2, sendKmers, curSendOfs, len2, startKmer, endKmer, kmerMask, rotater, findRoot(parent, curReadId, &numRandAccesses), lastKmerAddedByThread, scanProcessKmerRange, tid, rank);
                	if (ret == 0)
	                {
                	      numPEReads++;
	                }
        	        if (ret == 1) numSEReads++;
                	if (ret == 2) numInterleavedReads++;
			curReadId++;
        	        ret = getNextRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
	        }
	}
#endif
        return curReadId;

}

#if USE_MPI
void Alltoallv_p2p (uint32_t **sendKmers, size_t *sendCount, uint64_t *curSendOfs, size_t *rcvCount, uint64_t *rcvOfs)
{
	int maxSendSize = 999999999*2/numTasks;

	while (maxSendSize % 3 != 0)
		maxSendSize++;
	int commStep=0,i=0;
	MPI_Request reqs[1000];
	MPI_Status stats[1000];
	struct timeval sTime1, eTime1;
	for (commStep=0; commStep < numTasks; commStep++)
	{
		gettimeofday (&sTime1, NULL);
		int fromProc = rank-commStep;
		if (fromProc < 0)
			fromProc += numTasks;
		int toProc = (rank+commStep) % numTasks;
		size_t totalSendSize = sendCount[toProc];
		size_t totalRcvSize = rcvCount[fromProc];
		int sendIters = totalSendSize/maxSendSize;
		if ((totalSendSize % maxSendSize) != 0)
			sendIters++;
		int recvIters = totalRcvSize/maxSendSize;
		if ((totalRcvSize % maxSendSize) != 0)
			recvIters++;
		assert (recvIters <= 1000);
		size_t rcvSoFar=0;
		uint32_t *curRcvPtr=&kmerOut[rcvOfs[fromProc]];
		uint32_t *curSendPtr=sendKmers[numThreads*toProc];
		if (commStep == 0)
		{
			assert ((fromProc == toProc) && (totalSendSize == totalRcvSize));
			memcpy (curRcvPtr, curSendPtr, totalSendSize*sizeof(uint32_t));
	      		gettimeofday (&eTime1, NULL);
			long elapsed =(eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
			continue;
		}
		for (i=0; i<recvIters; i++)
		{
			int curRecv=0;
			size_t remToRcv = totalRcvSize-rcvSoFar;
			if (remToRcv > (size_t)maxSendSize)
				curRecv=maxSendSize;
			else
				curRecv=(int)remToRcv;
			MPI_Irecv(curRcvPtr, curRecv, MPI_UNSIGNED, fromProc, i, MPI_COMM_WORLD, &reqs[i]);
			rcvSoFar += curRecv;
			curRcvPtr+= curRecv;
		}
		assert (rcvSoFar == totalRcvSize);
		size_t sendSoFar=0;
		for (i=0; i<sendIters; i++)
		{
			int curSend=0;
			size_t remToSend = totalSendSize-sendSoFar;
			if (remToSend > (size_t)maxSendSize)
				curSend = maxSendSize;
			else
				curSend = (int)remToSend;
			MPI_Send (curSendPtr, curSend, MPI_UNSIGNED, toProc, i, MPI_COMM_WORLD);
			curSendPtr += curSend;
			sendSoFar+= curSend;
		}
		assert (sendSoFar == totalSendSize);
		MPI_Waitall (recvIters, reqs, stats);
  		gettimeofday (&eTime1, NULL);
		long elapsed =(eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
	}
}
#endif

int getNumRecvIters (size_t *rcvThreadCount, int maxSendSize)
{
	int numIters = 0;
	int i=0;
	for (i=0; i<numTasks*numThreads; i++)
	{
//		printf ("[%d] rcvThreadCount[%d]= %lu, numIters %d\n", rank, i, rcvThreadCount[i], numIters);
		if ((i/numThreads) == rank)
			continue;
		numIters += rcvThreadCount[i]/maxSendSize;
		if ((rcvThreadCount[i] % maxSendSize) != 0)
			numIters ++;
	}
	return numIters;
}
void writeOutput(uint32_t maxCompId)
{
        struct timeval sTime, eTime;
        gettimeofday (&sTime, NULL);
        int i=0, j=0;
#define OUTBUF_SIZE 100000000
#if USE_OMP
        #pragma omp parallel private(i) num_threads (numThreads)
        {
                int tid = omp_get_thread_num();
                char fName[100];
                sprintf (fName, "%s/Fk%d_Rank%d_Tid%d.paired_LC.fastq", oPrefix, K, rank, tid);
                FILE *fp_paired_lc = fopen(fName, "w");
                assert (fp_paired_lc != NULL);
                sprintf (fName, "%s/Fk%d_Rank%d_Tid%d.paired_other.fastq", oPrefix, K, rank, tid);
                FILE *fp_paired_other = fopen(fName, "w");
                assert (fp_paired_other != NULL);

                char *outP = (char *)malloc(OUTBUF_SIZE);
                char *outU = (char *)malloc(OUTBUF_SIZE);
                assert ((outP != NULL) && (outU != NULL));
                size_t outP_n=0, outU_n=0;

                char *read1=NULL, *read2=NULL;
    #else
                int tid = 0;
    #endif
                char *line1=NULL,*line2=NULL;
                size_t fqSize[2]={0,0}, fqInd[2]={0,0}, interleaveStart=0;
                int len1=0, len2=0;
                char *fqBuffer[2];
                fqBuffer[0]=&fastqBuffer1[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
                fqBuffer[1] = &fastqBuffer2[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
                int interleave=0;
                uint32_t curReadId=0;
                for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
                {
                        curReadId = fastqPart[i].readId;
                        loadChunk (i, fqInd, fqSize, fqBuffer, &interleave, &interleaveStart, tid);
                        int ret = getNextFullRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
                        while (ret != -1)
                        {
                                if (parent[curReadId]==maxCompId)
                                {
                                        if (ret == 0)
                                        {
                                                if ((outP_n+len1+len2) >= OUTBUF_SIZE)
                                                {
                                                        {
                                                                fwrite (outP, 1, outP_n, fp_paired_lc);
                                                        }
                                                        outP_n=0;
                                                }
                                              memcpy (&outP[outP_n], read1, len1);
                                                outP_n+=len1;
                                              memcpy (&outP[outP_n], read2, len2);
                                                outP_n+=len2;
                                        }
                                        if (ret == 1)
                                        {
                                                if ((outU_n+len1) >= OUTBUF_SIZE)
                                                {
                                                        {
//                                                                fwrite (outU, 1, outU_n, fp_unpaired);
                                                        }
                                                        outU_n=0;
                                                }
                                              memcpy (&outU[outU_n], read1, len1);
                                                outU_n+=len1;
                                        }
                                        if (ret == 2)
                                        {
                                                if ((outP_n+len1+len2) >= OUTBUF_SIZE)
                                                {
                                                        {
                                                                fwrite (outP, 1, outP_n, fp_paired_lc);
                                                        }
                                                        outP_n=0;
                                                }
                                              memcpy (&outP[outP_n], read1, len1);
                                                outP_n+=len1;
                                              memcpy (&outP[outP_n], read2, len2);
                                                outP_n+=len2;
                                        }
                                }
                                else
                                {
                                        if (ret != 1)
                                        {
                                                fwrite (read1, 1, len1, fp_paired_other);
                                                fwrite (read2, 1, len2, fp_paired_other);
                                        }
                                        else
                                        {
                                                assert (0);
                                        }
                                }
                                curReadId++;
                                ret = getNextFullRead (fqInd, fqSize, fqBuffer, &read1, &read2, &len1, &len2, interleave, interleaveStart);
                        }
                }
                {
                        if (outP_n > 0)
                                fwrite (outP, 1, outP_n, fp_paired_lc);
//                        if (outU_n > 0)
//                                fwrite (outU, 1, outU_n, fp_unpaired);
                }
                free (outP);
                free (outU);
                fclose (fp_paired_lc);
                fclose (fp_paired_other);
//            fclose (fp_unpaired);
#if USE_OMP
        }
#endif
        gettimeofday (&eTime, NULL);
}
void kmerGen(int scanNo)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	int i=0, j=0;
	size_t *sendThreadCount = (size_t *)malloc(numTasks*numThreads*sizeof(size_t));
	assert (sendThreadCount != NULL);
	recvThreadCount = (size_t *)malloc((numTasks*numThreads + 2)*sizeof(size_t));
	assert (recvThreadCount != NULL);
	recvThreadOfs = (size_t *)malloc((numTasks*numThreads + 2)*sizeof(size_t));
	assert (recvThreadOfs != NULL);
	size_t *sendCount = (size_t *)malloc(numTasks*sizeof(size_t));
	assert (sendCount != NULL);
	size_t *rcvCount = (size_t *)malloc(numTasks * sizeof(size_t));
	assert (rcvCount != NULL);
	uint32_t **sendKmers = (uint32_t **)malloc(numTasks*numThreads * sizeof(uint32_t *));
	assert (sendKmers != NULL);
	uint64_t *curSendOfs = (uint64_t *)calloc (numTasks*numThreads, sizeof(uint64_t));
	assert (curSendOfs != NULL);
	int *curSendIter = (int *)calloc (numTasks*numThreads, sizeof(int));
	assert (curSendIter != NULL);

	size_t totalSendSize=0;
		
	for (i=0; i<numTasks; i++)
	{
		sendCount[i] = sendSizes[scanNo*numTasks + i];
		totalSendSize += sendCount[i];
		for (j=0; j<numThreads; j++)
		{
			sendThreadCount[i*numThreads + j] = sendThreadSizes[scanNo*numTasks*numThreads + i*numThreads + j];
		}
	}
#if USE_MPI
	MPI_Alltoall (sendThreadCount, numThreads, MPI_UNSIGNED_LONG, recvThreadCount, numThreads, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
	memcpy (&recvThreadOfs[2], recvThreadCount, numThreads*numTasks*sizeof(size_t));
#else
	memcpy (&recvThreadOfs[2], sendThreadCount, numThreads*sizeof(size_t));
	memcpy (recvThreadCount, sendThreadCount, numThreads*sizeof(size_t));
#endif
	
	recvThreadOfs[0]=0;
	for (i=1; i<=numTasks * numThreads; i++)
	{
		recvThreadOfs[i]=recvThreadOfs[i-1]+recvThreadOfs[i+1];
	}
	totalSendSize *= 3* sizeof(uint32_t);
	assert (kmerIn != NULL);
	uint64_t curOfs=0;
	for (i=0; i<numTasks*numThreads; i++)
	{
		sendKmers[i] = &kmerIn[curOfs *3];
		curOfs += sendThreadCount[i];
	}

#if USE_MPI
	MPI_Alltoall (sendCount, 1, MPI_UNSIGNED_LONG, rcvCount, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
	rcvCount[0] = sendCount[0];
#endif

	uint64_t *scanProcessKmerRange = &processKmerRanges[scanNo*numTasks]; 
	int curTask=0;
	size_t totalRcv=0;
	for (i=0; i<numTasks; i++)
	{
		totalRcv += rcvCount[i];
	}
	bzero (threadKmerGenTime, numThreads*sizeof(long));
	size_t *lastKmerAddedByThread = (size_t *)calloc(numThreads*numTasks , sizeof(size_t));
	assert (lastKmerAddedByThread != NULL);
	numKmers = totalRcv;
	assert (kmerOut!=NULL);
	int numOut=0;

	int numRcvIters = getNumRecvIters(recvThreadCount, ASYNC_SIZE);
#if USE_OMP
	#pragma omp parallel private(i, j) num_threads (numThreads)
	{
		int tid = omp_get_thread_num();
         
#else 
   		int tid = 0;
#endif
		struct timeval sTime1, eTime1;
    		char *line1=NULL,*line2=NULL;
    		size_t fqSize[2]={0,0}, fqInd[2]={0,0}, interleaveStart=0;
    		int len1=0, len2=0;
    		char *fqBuffer[2];
    		fqBuffer[0]=&fastqBuffer1[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
    		fqBuffer[1] = &fastqBuffer2[(size_t)tid*(fastqBufferSize+10)];//(char *)malloc(fastqBufferSize+10);
    		int interleave=0;
    		uint32_t curReadId=0;
		size_t cpySoFar=0;
    		for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
    		{
    			curReadId = fastqPart[i].readId;
    			loadChunk (i, fqInd, fqSize, fqBuffer, &interleave, &interleaveStart, tid);
		        while ((fqInd[0] < fqSize[0]) && (fqBuffer[0][fqInd[0]] != '\n')) fqInd[0]++;
		        fqInd[0]++;
			assert (fqInd[0] < fqSize[0]);
		        if (fqInd[1] < fqSize[1])
		        {
		                while ((fqInd[1] < fqSize[1]) && (fqBuffer[1][fqInd[1]] != '\n')) fqInd[1]++;
                		fqInd[1]++;
		                assert (fqInd[1] < fqSize[1]);
		        }
    			gettimeofday (&sTime1, NULL);
    			curReadId = getReadAndKmerCount (i, fqInd, fqSize, fqBuffer, interleave, interleaveStart, sendKmers, curSendOfs, scanKmerRanges[scanNo], scanKmerRanges[scanNo+1], kmerMask, rotater, curReadId, &lastKmerAddedByThread[rank*numThreads + tid], scanProcessKmerRange, tid);

    			gettimeofday (&eTime1, NULL);
    			threadKmerGenTime[tid] +=(eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
    		}
#if USE_OMP
	}
#endif
	long maxWorkTime=0;
	for (i=0; i<numThreads; i++)
	{
		if (threadKmerGenTime[i] > maxWorkTime)
			maxWorkTime = threadKmerGenTime[i];
	}
	kmerGenWorkTime += maxWorkTime;
	gettimeofday (&eTime, NULL);
	kmerGenTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	gettimeofday (&sTime, NULL);

#if USE_MPI
	if (numTasks > 1)
	{
		uint64_t *rcvOfs = (uint64_t *)malloc(numTasks * sizeof(uint64_t));
		assert (rcvOfs != NULL);
		rcvOfs[0]=0; curSendOfs[0]=0;
		sendCount[0]*=3; rcvCount[0]*=3;
		for (i=1; i<numTasks; i++)
		{
			rcvOfs[i]=rcvOfs[i-1]+rcvCount[i-1];
			curSendOfs[i]=curSendOfs[i-1]+sendCount[i-1];
			sendCount[i] *= 3;
			rcvCount[i] *= 3;
		}
	
		Alltoallv_p2p (sendKmers, sendCount, curSendOfs, rcvCount, rcvOfs);
		free (rcvOfs);
	}
#endif
	gettimeofday (&eTime, NULL);
	commTime +=(eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);

	free (sendKmers);
	free (lastKmerAddedByThread);
	free (curSendOfs);
	free (rcvCount);
	free (sendCount);
	free (sendThreadCount);
	free (curSendIter);
}

uint64_t radixPass (int pass, uint32_t *in_list, uint32_t *out_list, uint64_t startOfs, uint64_t numKmersToProcess, int scanNo, int tid)
{
	int passOrder[3]={2,0,1};
	uint32_t radixMask=0xFFu;
	int radixRotate=0;
	uint32_t i=0,j=0;
	radixMask = radixMask << (8 * (pass & 3));
	radixRotate = 8 * (pass & 3);
	int uintOfs = passOrder[pass/4];

	uint64_t curKmerMask = 0;

	uint64_t *rangeOfsCopy = (uint64_t *)calloc(256+2 , sizeof(uint64_t));
	assert (rangeOfsCopy != NULL);
	uint64_t maxCompCount=0;
	uint32_t *rangeLastSeenComp=NULL;
	uint64_t *rangeCompCount=NULL;
	uint64_t *suffixCompCount=NULL;
	uint64_t currentSuffix=0;
	uint64_t kmerInd=0;
	uint64_t kmerRemoveCount=0;
	uint64_t kmerCount=0;
	int curM=0;
	uint64_t beginInd=startOfs, beginInd2=startOfs;
	for (kmerInd = startOfs; kmerInd<(startOfs+numKmersToProcess); kmerInd++)
	{
		uint32_t bucket = (in_list[kmerInd*3 + uintOfs] & radixMask) >> radixRotate;
		rangeOfsCopy[bucket+2]++;
	}

#define LMAX 100
#define LSIZE 300
	uint32_t lRadixKmers[256][LSIZE];
	int lCount[256];
	for (i=1; i<=256; i++)
	{
		rangeOfsCopy[i]=rangeOfsCopy[i-1]+rangeOfsCopy[i+1];
		lCount[i-1]=0;
	}
	uint64_t retVal = rangeOfsCopy[256];
	for (kmerInd = startOfs; kmerInd<(startOfs+numKmersToProcess); kmerInd++)
	{
		uint32_t bucket = (in_list[kmerInd*3 + uintOfs] & radixMask) >> radixRotate;
		if (in_list[kmerInd*3 +2] != 0xFFFFFFFF)
		{
			if (lCount[bucket] == LMAX)
			{
				uint64_t curOfs = startOfs+rangeOfsCopy[bucket];
				for (i=0; i<LMAX; i++)
				{
					out_list[curOfs*3]=lRadixKmers[bucket][i*3];
					out_list[curOfs*3 + 1]=lRadixKmers[bucket][i*3 + 1];
					out_list[curOfs*3 + 2]=lRadixKmers[bucket][i*3 + 2];
					curOfs++;
				}
				rangeOfsCopy[bucket]+=LMAX;
				lCount[bucket]=0;
			}
			int curOfs = lCount[bucket];
			lRadixKmers[bucket][curOfs*3]=in_list[kmerInd*3];
			lRadixKmers[bucket][curOfs*3 + 1]=in_list[kmerInd*3 + 1];
			lRadixKmers[bucket][curOfs*3 + 2]=in_list[kmerInd*3 + 2];
			lCount[bucket]++;
		}
	}
	for (i=0; i<256; i++)
	{
		uint64_t curOfs = startOfs+rangeOfsCopy[i];
		for (j=0; j<lCount[i]; j++)
		{
			out_list[curOfs*3]=lRadixKmers[i][j*3];
			out_list[curOfs*3 +1]=lRadixKmers[i][j*3 +1];
			out_list[curOfs*3 +2]=lRadixKmers[i][j*3 +2];
			curOfs++;
		}
		rangeOfsCopy[i]+=lCount[i];
	}
	assert (rangeOfsCopy[255] == rangeOfsCopy[256]);
	free (rangeOfsCopy);
	return retVal;
}

void calcRangeCounts (int tid, size_t *threadRangeCounts, uint64_t *localThreadKmerRanges)
{
	int startThread =  tid * numTasks;
	int endThread =  (tid+1) * numTasks;
	int fastqPartitionsPerThread = totalFastqPartitions/numTasks/numThreads;
	int startFastqPartition = startThread*fastqPartitionsPerThread;
	int endFastqPartition = endThread * fastqPartitionsPerThread;
	size_t startRcvBufInd = recvThreadOfs[startThread];
	size_t endRcvBufInd = recvThreadOfs[endThread];
//	printf ("Rank %d Tid %d startFilePart %d endFilePart %d startRcvBufInd %lu, endRcvBufInd %lu\n", rank, tid, startFastqPartition, endFastqPartition, startRcvBufInd, endRcvBufInd);
/*	if (tid == 0)
	{
		int i=0;
		for (i=0; i<numTasks * numThreads; i++)
			printf ("%lu\t", recvThreadCount[i+1]-recvThreadCount[i]);
		printf ("\n");
	}*/
	int f,t,p;
	for (f=startFastqPartition; f<endFastqPartition; f++)
	{
		for (t=0; t<numThreads; t++)
		{
			for (p=localThreadKmerRanges[t]; p<localThreadKmerRanges[t+1]; p++)
			{
				threadRangeCounts[tid*numThreads + t] += fastqPart[f].kmerFreqCount[p];
			}
		}
	}
}

void localSort (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets)
{
	typedef struct 
	{
		uint32_t hi_kmer;
		uint32_t lo_kmer;
		uint32_t readId;
		int range;
	}kmerTemp;
	size_t *threadRangeCounts =(size_t *)calloc (numThreads * numThreads ,sizeof(size_t));
	assert (threadRangeCounts != NULL);
	size_t *threadRangeOfs =(size_t *)calloc (numThreads * numThreads ,sizeof(size_t));
	assert (threadRangeOfs != NULL);

	struct timeval sTime, eTime;
	int *threadMask = (int *)calloc(1048576, sizeof(int));
	assert (threadMask!=NULL);
	int f=0, i;
	uint64_t i1=0;
	uint64_t processLimitUp = processKmerRanges[scanNo*numTasks+rank + 1];
	uint64_t processLimitDn = processKmerRanges[scanNo*numTasks+rank];
#ifdef DEBUG
	{
		for (i1=0; i1<numKmers; i1++)
		{
			uint64_t curKmer = *((uint64_t *)&kmerOut[i1 * 3]);
			uint64_t bucket=(curKmer & kmerMask) >> rotater;
			if ((bucket >= processLimitUp) || (bucket < processLimitDn))
			{
				printf ("[%d]i1 %lu curKmer %llx bucket %lu(%p) procLimUp %lu, procLimDn %lu", rank, i1, curKmer, bucket, &kmerOut[i1 * 3], processLimitUp, processLimitDn);
				printf ("*\n");
			}
			assert ((bucket < processLimitUp) && (bucket >= processLimitDn));
		}
	}
#endif
	uint32_t prevReadId = kmerOut[2];
	assert (kmerSorted != NULL);
	uint64_t *localThreadKmerRanges = &threadKmerRanges[scanNo*numTasks*numThreads + rank*numThreads];
	uint64_t *localRangeOffsets = (uint64_t *)calloc(numThreads, sizeof(uint64_t));
	assert (localRangeOffsets != NULL);	
	for (f=0; f<totalFastqPartitions; f++)
	{
		int curThread=0;
		uint64_t startRange = (localThreadKmerRanges[0]);
		uint64_t endRange = (localThreadKmerRanges[numThreads]);
		for (i1=startRange; i1<endRange; i1++)
		{
			if (i1 >= localThreadKmerRanges[curThread+1])
			{
				curThread++;
				localRangeOffsets[curThread]=localRangeOffsets[curThread-1]+localRangeCounts[curThread-1];
				rangeOffsets[curThread]=localRangeOffsets[curThread];
			}
			threadMask[i1]=curThread;
			localRangeCounts[curThread]+=fastqPart[f].kmerFreqCount[i1];
		}
	}
	for (i=0; i<numThreads; i++)
	{
		calcRangeCounts (i, threadRangeCounts, localThreadKmerRanges);
	}
//	if (rank == 0)
	{
		uint64_t curCount=0;
		for (i=0; i<numThreads; i++)//column
		{
			for (f=0; f<numThreads; f++)//row
			{
				threadRangeOfs[f*numThreads + i]=curCount;
				curCount += threadRangeCounts[f*numThreads + i];
			}	
		}

	}

	gettimeofday (&sTime, NULL);
#if 1
//Parallel Range Partiton
#if USE_OMP
	#pragma omp parallel num_threads(numThreads)
	{
		int tid=omp_get_thread_num();
#else
	{
		int tid=0;
#endif
		int startThread =  tid * numTasks;
		int endThread = (tid+1) * numTasks;
		size_t startRcvBufInd = recvThreadOfs[startThread];
		size_t endRcvBufInd = recvThreadOfs[endThread];
		int f=0, i;
		uint64_t j1=0;
		uint64_t prevInd=0;
#define BUFSIZE 1000
		kmerTemp lKmers[64][BUFSIZE];
		int rangeWiseCounts[64];
		assert (numThreads <= 64); 
		for (i=0; i<numThreads; i++)
			rangeWiseCounts[i]=0;

		size_t *localPtr = &threadRangeOfs[tid*numThreads];
		for (j1=startRcvBufInd; j1 < endRcvBufInd; j1++)
		{
			uint64_t maskedKmer=*(uint64_t *)(&kmerOut[j1*3]);
			maskedKmer = (maskedKmer & kmerMask) >> rotater;
			int currentBin=0;
			while (maskedKmer >= localThreadKmerRanges[currentBin])
				currentBin++;
			currentBin--;
			if ((currentBin < 0) || (currentBin >= numThreads))
			{
				printf ("[%d] Tid %d numThreads %d j1 %lu curKmer %llx, maskedKmer %lu, localThreadKmerRanges[0]=%llx, currentBin %d kmerMask %llx, cond1 %d, cond2 %d\n",  rank, tid, numThreads, j1, *(uint64_t *)(&kmerOut[j1*3]), maskedKmer, localThreadKmerRanges[0], currentBin, kmerMask, (currentBin < 0), (currentBin >= numThreads));
			}
			assert ((currentBin >= 0) && (currentBin < numThreads));
			uint32_t kmer_l = kmerOut[j1*3];
			uint32_t kmer_h = kmerOut[j1*3 + 1];
			uint32_t readId = kmerOut[j1*3 + 2];
			if (rangeWiseCounts[currentBin]==BUFSIZE)
			{
				uint64_t insertPos = localPtr[currentBin];
				for (i=0; i<BUFSIZE; i++)
				{
					kmerSorted[3*insertPos]=lKmers[currentBin][i].lo_kmer;
					kmerSorted[3*insertPos + 1]=lKmers[currentBin][i].hi_kmer;
					kmerSorted[3*insertPos + 2]=lKmers[currentBin][i].readId;
					insertPos++;
				}
				localPtr[currentBin]+=BUFSIZE;
				lKmers[currentBin][0].hi_kmer = kmer_h;
				lKmers[currentBin][0].lo_kmer = kmer_l;
				lKmers[currentBin][0].readId = readId;
				lKmers[currentBin][0].range = currentBin;
				rangeWiseCounts[currentBin]=1;
			}
			else
			{
				int ofs=rangeWiseCounts[currentBin];
				lKmers[currentBin][ofs].hi_kmer = kmer_h;
				lKmers[currentBin][ofs].lo_kmer = kmer_l;
				lKmers[currentBin][ofs].readId = readId;
				lKmers[currentBin][ofs].range = currentBin;
				rangeWiseCounts[currentBin]++;
			}
		}
		for (i=0; i<numThreads; i++)
		{
			if (rangeWiseCounts[i] > 0)
			{
				uint64_t insertPos = localPtr[i];
				for (f=0; f<rangeWiseCounts[i]; f++)
				{
					kmerSorted[3*insertPos]=lKmers[i][f].lo_kmer;
					kmerSorted[3*insertPos + 1]=lKmers[i][f].hi_kmer;
					kmerSorted[3*insertPos + 2]=lKmers[i][f].readId;
					insertPos++;
				}
				localPtr[i]=insertPos;
			}
		}
	} 
//	if (rank == 0)
#ifdef DEBUG
	{
		for (i1=0; i1<numKmers; i1++)
		{
			uint64_t curKmer = *((uint64_t *)&kmerSorted[i1 * 3]);
			uint64_t bucket=(curKmer & kmerMask) >> rotater;
			if ((bucket >= processLimitUp) || (bucket < processLimitDn))
			{
				printf ("[%d]i1 %lu curKmer %llx bucket %lu", rank, i1, curKmer, bucket);
				printf ("*\n");
			}
			assert ((bucket < processLimitUp) && (bucket >= processLimitDn));
		}
	}
#endif
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
    if (rank==0)
    	printf ("[%d] Range Partition Done (%lf sec)\n", rank, (double)elapsed/1000000);
	sortTime1 += elapsed;
	int curThread=0;
	gettimeofday (&sTime, NULL);
#if USE_OMP
	#pragma omp parallel num_threads(numThreads)
	{
		int curThread=omp_get_thread_num();
		assert (omp_get_num_threads() == numThreads);
#else
	{
		//Serial out-of-place Radix sort
		int curThread=0;
#endif
		localRangeCounts[curThread] = radixPass (4, kmerSorted, kmerOut, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (5, kmerOut, kmerSorted, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (6, kmerSorted, kmerOut, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (7, kmerOut, kmerSorted, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (8, kmerSorted, kmerOut, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (9, kmerOut, kmerSorted, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (10, kmerSorted, kmerOut, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
		localRangeCounts[curThread] = radixPass (11, kmerOut, kmerSorted, rangeOffsets[curThread], localRangeCounts[curThread], scanNo, curThread);
//		break;
	}
	gettimeofday (&eTime, NULL);
	elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
    if (rank == 0)
    	printf ("[%d] Sort Done (%lf sec)\n", rank, (double)elapsed/1000000);
	sortTime2 += elapsed;
	uint64_t *lastAddr = (uint64_t *)(&(kmerSorted[(rangeOffsets[numThreads-1]+localRangeCounts[numThreads-1]-1)*3]));
//	if (rank == 0)
#ifdef DEBUG
	{
		for (i1=0; i1<numKmers; i1++)
		{
			uint64_t curKmer = *((uint64_t *)&kmerSorted[i1 * 3]);
			uint64_t bucket=(curKmer & kmerMask) >> rotater;
			if ((bucket >= processLimitUp) || (bucket < processLimitDn))
			{
				printf ("[%d]i1 %lu curKmer %llx bucket %lu", rank, i1, curKmer, bucket);
				printf ("*\n");
			}
			assert ((bucket < processLimitUp) && (bucket >= processLimitDn));
		}
	}
#endif
#ifdef DEBUG
	for (i=0; i<numThreads; i++)
	{
		uint64_t prevKmer=*((uint64_t *)&kmerSorted[rangeOffsets[i]*3]);
		prevReadId = kmerSorted[rangeOffsets[i]*3 +2];
		for (i1=rangeOffsets[i]+1; i1<(rangeOffsets[i]+localRangeCounts[i]); i1++)
		{		
			uint64_t curKmer = *(uint64_t *)(&(kmerSorted[i1*3]));
			uint32_t curReadId = kmerSorted[i1*3 +2];
			if ((curKmer < prevKmer) )
				printf ("i1 %lu CurKmer %" PRIu64 ", prevKmer %" PRIu64 ", curReadId %u prevReadId %u, numKmers %lu\n", i1, curKmer, prevKmer, curReadId, prevReadId, numKmers);
			assert ((curKmer >= prevKmer) ); 
			prevKmer = curKmer;
			prevReadId = curReadId;
		}
	}
#endif
#endif
	free (localRangeOffsets);
	free (threadMask);
	free (threadRangeCounts);
	free (threadRangeOfs);
}


void localCC (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	#pragma omp parallel num_threads (numThreads)
	{
		uint64_t i=0;
		int numUnions=0;
	        int numRandAccess=0, maxRandAccess=0;
		int tid = omp_get_thread_num();
		uint64_t curKmer = *(uint64_t *)&(kmerSorted[(numKmers-1)*3]);
		uint32_t minReadId = kmerSorted[2];
		uint32_t prevReadId = kmerSorted[2];
		uint32_t curVertexDegree=0;
		int processThisKmer=1;
		int disp=0;
		uint32_t startReadId = tid*numReads/numThreads;
		uint32_t endReadId = (tid+1)*numReads/numThreads;
	        if (tid == (numThreads-1))
			endReadId = numReads;
		uint64_t insPos = rangeOffsets[tid];
		struct timeval sTime1, eTime1;
		gettimeofday (&sTime1, NULL);
		for (i=rangeOffsets[tid]; i<(rangeOffsets[tid]+localRangeCounts[tid]); i++)
		{
			uint64_t newKmer = *((uint64_t *)&(kmerSorted[3*i]));
			uint32_t curReadId = kmerSorted[3*i + 2];
			if (newKmer != curKmer)
			{
				disp=0;
				curKmer = newKmer;
				minReadId = curReadId;
			}	
			else
			{
				if (processThisKmer == 1) 
				{
					uint32_t p_u = minReadId;//findRoot (parent, minReadId, &numRandAccess);
					uint32_t p_v = curReadId;//findRoot (parent, curReadId, &numRandAccess);
					if (p_u != p_v)
					{
						int rAccess=0;
						while (parent[p_u] != p_u)
						{
							if (parent[p_u] != parent[parent[p_u]])
							parent[p_u] = parent[parent[p_u]];
							p_u = parent[p_u];
							rAccess = rAccess + 1;
						}
						numRandAccess = numRandAccess + rAccess;
						if (maxRandAccess < rAccess)
							maxRandAccess = rAccess;
						rAccess=0;
						
						while (parent[p_v] != p_v)
						{
							if (parent[p_v] != parent[parent[p_v]])
								parent[p_v] = parent[parent[p_v]];
							p_v = parent[p_v];
							rAccess = rAccess + 1;
						}
						numRandAccess = numRandAccess + rAccess;
						if (maxRandAccess < rAccess)
							maxRandAccess = rAccess;
					}

					if (p_u != p_v)
					{
						kmerSorted[3*insPos + 2*numUnions] = minReadId;
						kmerSorted[3*insPos + 2*numUnions + 1] = curReadId;
						numUnions++;
						if (p_u > p_v)
							parent[p_v]=p_u;
//							unionOp (parent, compSizes, p_u, p_v);
						else
							parent[p_u]=p_v;
//							unionOp (parent, compSizes, p_v, p_u);
					}
				}
			}
			prevReadId = curReadId;
		}
	    	gettimeofday (&eTime1, NULL);
		long elapsed1 = (eTime1.tv_sec * 1000000 + eTime1.tv_usec) - (sTime1.tv_sec * 1000000 + sTime1.tv_usec);
		localRangeCounts[tid]=numUnions;
	}
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	unionFindTime1 += elapsed;
}
void localCCVerify (int scanNo, uint64_t *localRangeCounts, uint64_t *rangeOffsets)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint64_t curKmer = *(uint64_t *)&(kmerSorted[(numKmers-1)*3]);
	uint32_t prevReadId = kmerSorted[2];
	uint32_t curVertexDegree=0;
	int processThisKmer=1;
	int disp=0;
	#pragma omp parallel num_threads(numThreads)
	{
		int tid=omp_get_thread_num();
        	int numRandAccess = 0;
		while (localRangeCounts[tid] > 0)
		{
			int numUnions=0;
			uint64_t i=0;
			for (i=0; i<localRangeCounts[tid]; i++)
			{
				uint32_t fromReadId = kmerSorted[3*rangeOffsets[tid] + 2*i];
				uint32_t toReadId = kmerSorted[3*rangeOffsets[tid] + 2*i +1];
				uint32_t p_u = findRoot (parent, fromReadId, &numRandAccess);
				uint32_t p_v = findRoot (parent, toReadId, &numRandAccess);
				if (p_u != p_v)
				{
					numUnions++;
					kmerSorted[3*rangeOffsets[tid] + 2*numUnions] = fromReadId;
					kmerSorted[3*rangeOffsets[tid] + 2*numUnions + 1] = toReadId;
					if (p_u > p_v)
						unionOp (parent, compSizes, p_u, p_v);
					else
						unionOp (parent, compSizes, p_v, p_u);
				}
			}
			localRangeCounts[tid] = numUnions;
		}
	}
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	unionFindTime1 += elapsed;
}
void mergeStep (uint32_t *otherRankParent)
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	#pragma omp parallel num_threads(numThreads)
	{
		uint32_t numUnions=1;
		int numRandAccesses = 0;
		int tid = omp_get_thread_num();
	        uint32_t startOfs = (uint32_t)((uint64_t)tid * numReads / numThreads);
	        uint32_t endOfs = (uint32_t)((uint64_t)(tid+1) * numReads / numThreads);
		uint32_t i=0;
	
		while (numUnions > 0)
		{
			numUnions = 0;
	            
			for (i=startOfs; i<endOfs; i++)
			{
				if ((parent[i]==otherRankParent[i]) || (otherRankParent[i] == numReads))
				{
					otherRankParent[i] = numReads;
					continue;
				}
				uint32_t p_u = findRoot (parent, i, &numRandAccesses);
				uint32_t p_v = findRoot (parent, otherRankParent[i], &numRandAccesses);
				if (p_u != p_v)
				{
					if (p_u > p_v)
						unionOp (parent, compSizes, p_u, p_v);
					else
						unionOp (parent, compSizes, p_v, p_u);
					numUnions++;
				}
				else
					otherRankParent[i] = numReads;
			}
		}
	}
	gettimeofday (&eTime, NULL);
	long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
}
uint32_t printStats()
{
	struct timeval sTime, eTime;
	gettimeofday (&sTime, NULL);
	uint32_t readNo = 0, countReads=0, uniqComps=0;
	uint32_t largestCompSize=0, largestCompId=0;
	int numRandAccesses=0;
	for (readNo=0; readNo<numReads; readNo++)
	{
		parent[readNo]=findRoot(parent, readNo, &numRandAccesses);
		if (parent[readNo] == readNo)
		{
			uniqComps++;
		}	
		else
		{
			compSizes[readNo]=0;
			compSizes[parent[readNo]]++;
			if (compSizes[parent[readNo]] > largestCompSize)
			{
				largestCompSize = compSizes[parent[readNo]] ;
				largestCompId = parent[readNo];
			}
		}
	}
	gettimeofday (&eTime, NULL);
	unionFindTime2 += (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	return largestCompId;
}

void printStats2 ()
{
	struct timeval sTime, eTime;
	uint32_t readNo = 0, countReads=0, uniqComps=0;
	gettimeofday (&sTime, NULL);
	qsort (compSizes, numReads, sizeof(uint32_t), parentCompare);
	printf ("[%d]largestComp Size %u\n", rank, compSizes[numReads-1]);
	uint32_t curCompSize=0, curCompSizeCount=0, numComps=0;
	for (readNo=0; readNo<numReads; readNo++)
	{
		if (compSizes[readNo] != curCompSize)
		{
			if (curCompSize > 0)
			{
				printf ("CompSize %u NumComponents %u\n", curCompSize, curCompSizeCount);
				numComps+=curCompSizeCount;
			}
			curCompSize = compSizes[readNo];
			curCompSizeCount=0;
		}
		curCompSizeCount++;
	}
	if (curCompSize > 0)
	{
		printf ("CompSize %u NumComponents %u\n", curCompSize, curCompSizeCount);
		numComps+=curCompSizeCount;
	}
	gettimeofday (&eTime, NULL);
}
int main(int argc, char **argv)
{
	int i=0, j=0, k=0;
	static struct option long_options[] = {
		{"out-prefix", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};
	int c=0;
	while (1)
	{
		int option_index = 0;
		c = getopt_long (argc, argv, "o:", long_options, &option_index);
		if (c==-1) 
		{
			break;
		}
		switch (c)
		{
			case 'o':
				strcpy (oPrefix, optarg);
//				printf ("optarg %s\n", optarg);
				break;
		}
	}
    int providedThreadSupport;
#if USE_MPI
	MPI_Init_thread(&argc,&argv, MPI_THREAD_MULTIPLE, &providedThreadSupport);
	MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0)
        	printf ("[%d] providedThreadSupport %d\n", rank, providedThreadSupport);
#else
	numTasks = 1;
	rank = 0;
#endif
	assert (numTasks < 1024);
	if (rank==0)
   		printf ("numTasks %d\n", numTasks);
	struct stat buf;
	hwThreads = hardware_threads();
	if (rank==0)
    		printf ("HwThreads %d argc %d optind %d\n", hwThreads, argc, optind);
	if (argc < 11)
	{
		printf ("Usage: %s -o <output-prefix> fastqBufferSize kmerSize numPasses numThreads "
				"<numPEFiles> <numSEFiles> <numInterleavedFiles> <FASTQfiles> "
				"b.fastq ...\n", argv[0]);
		printf ("FASTQfiles - space separated list of unzipped FASTQ files. "
			"numPEFiles*2 files with PE reads example a.1.fq a.2.fq b.1.fq b.2.fq ... "
			"Followed by numSEFiles SE files and numInterleavedFiles interleaved PE files\n");
		exit(1);
	}
	merHist = (uint64_t *)malloc((1048576 +1) * sizeof(uint64_t));
	assert (merHist != NULL);
	
	argv = &argv[optind-1];
	fastqBufferSize = (unsigned long)atol(argv[1]);
	K = atoi(argv[2]);
	getKmerMask (K);
	passes = atoi (argv[3]);
	numThreads = atoi(argv[4]);
	threadKmerGenTime = (long *)malloc (numThreads * sizeof (long));
	assert (threadKmerGenTime != NULL);

	numPEFiles = atoi(argv[5])*2;
	numSEFiles = atoi(argv[6]);
	numInterleavedFiles = atoi(argv[7]);
	argv = &argv[7];
	argc-=(7+optind-1);
	printf ("fastqBufferSize %lu k-mer size %d passes %d numPEFiles %d numSEFiles %d numInterleavedFiles %d\n", fastqBufferSize, K, passes, numPEFiles, numSEFiles, numInterleavedFiles);
	init(K, numTasks, numThreads);
	fileSize = (size_t *)malloc((argc-1) * sizeof(size_t));
	fileNames = (char **)malloc((argc-1) * sizeof(char *));
	assert (fileNames != NULL);
	assert (fileSize != NULL);

	for (i=1; i<argc; i++)
	{
		fileNames[i-1]=(char *)malloc(500);
		assert (fileNames[i-1] != NULL);
		sprintf (fileNames[i-1], "%s", argv[i]);
		int status = stat (fileNames[i-1], &buf);
		assert (status == 0);
		fileSize[i-1]=buf.st_size;
	}
	fastqRanges = (int *)malloc ((numThreads+1) * sizeof(int));
	assert (fastqRanges != NULL);

	char idxFile[200];
	sprintf (idxFile, "%s.KmerIdx", argv[1]);
	if (rank==0)
    		printf ("idxFile %s\n", idxFile);
	int status = stat (idxFile, &buf);
	assert (status == 0);
	totalFastqPartitions = buf.st_size/sizeof(fileIndex) -1;
	int fastqPartitionsPerThread = totalFastqPartitions/numTasks/numThreads;

	fastqPart = (fileIndex *)malloc((totalFastqPartitions + 1) * sizeof(fileIndex));
	assert (fastqPart != NULL);
	FILE *fp_idx = fopen (idxFile, "r");
	assert (fp_idx != NULL);
	assert (fread(fastqPart, sizeof(fileIndex), totalFastqPartitions+1, fp_idx) == (totalFastqPartitions+1));
	fclose (fp_idx);
	numReads = fastqPart[totalFastqPartitions].readId;
	if (rank==0)
        	printf ("[%d] NUMREADS %u\n", rank, numReads);

	for (i=0; i<=numThreads; i++)
	{
		fastqRanges[i]=rank*numThreads*fastqPartitionsPerThread + i*fastqPartitionsPerThread;
	}

	/*Range Hist*/
	sprintf (idxFile, "%s.kmerHist", argv[1]);
	FILE *fp_merHist = fopen(idxFile, "r");
	assert (fp_merHist != NULL);
	assert (fread(merHist, sizeof(uint64_t), 1048576+1, fp_merHist) == (1048576+1));
	fclose (fp_merHist);
	if (rank==0)
		printf ("TotalKmers %lu\n", merHist[1048576]);

	scanKmerRanges = (uint64_t *)calloc((passes+1) , sizeof(uint64_t));
	processKmerRanges = (uint64_t *)calloc(passes*numTasks +1 , sizeof(uint64_t));
	threadKmerRanges = (uint64_t *)calloc((passes*numTasks*numThreads + 1) , sizeof(uint64_t));
	assert ((scanKmerRanges != NULL) && (processKmerRanges != NULL) && (threadKmerRanges != NULL));
	uint64_t avgKmersPerScan = merHist[1048576]/passes;
	int curRange=1;
	scanKmerRanges[0]=0;
	
	for (i=0; i<1048576; i++)
	{
		if (merHist[i] >= avgKmersPerScan*curRange)
		{
			scanKmerRanges[curRange]=i;
			curRange++;
		}
	}
	if (merHist[i] >= avgKmersPerScan*curRange)
	{
		scanKmerRanges[curRange]=i;
	}
	uint64_t avgKmersPerProcessPerScan = avgKmersPerScan/numTasks;
	
	curRange=1;
	for (i=0; i<1048576; i++)
	{
		if (merHist[i] >= curRange*avgKmersPerProcessPerScan)
		{
			if ((curRange % numTasks) == 0)
				processKmerRanges[curRange]=scanKmerRanges[curRange/numTasks];
			else
				processKmerRanges[curRange]=i;
			curRange++;
		}
	}
	if (merHist[i] >= curRange*avgKmersPerProcessPerScan)
	{
		processKmerRanges[curRange]=i;
	}
	uint64_t avgKmersPerThreadPerScan = avgKmersPerProcessPerScan/numThreads;
	curRange=1;
	for (i=0; i<1048576; i++)
	{
		if (merHist[i] >= curRange*avgKmersPerThreadPerScan)
		{
			if ((curRange % numThreads) == 0)
				threadKmerRanges[curRange]=processKmerRanges[curRange/numThreads];
			else
				threadKmerRanges[curRange]=i;
			curRange++;
		}
	}
	if (merHist[i] >= curRange*avgKmersPerThreadPerScan)
	{
		threadKmerRanges[curRange]=i;
	}

	sendSizes = (size_t *)calloc(passes*numTasks , sizeof(size_t));
	assert (sendSizes != NULL);
	sendThreadSizes = (size_t *)calloc(passes*numTasks*numThreads , sizeof(size_t));
	assert (sendThreadSizes != NULL);
	
	int tid=0;
	size_t *sendSizeInScan = (size_t *)calloc(passes, sizeof(size_t));
	assert (sendSizeInScan != NULL);
	for (tid=0; tid<numThreads; tid++)
	{
		for (i=fastqRanges[tid]; i<fastqRanges[tid+1]; i++)
		{
			k=0;
			for (j=0; j<passes*numTasks; j++)
			{
				int curScan=j/numTasks;
				int curTask = j%numTasks;
				while (k < processKmerRanges[curScan*numTasks + curTask + 1])
				{
					sendSizes[curScan*numTasks + curTask]+=fastqPart[i].kmerFreqCount[k];
					sendThreadSizes[curScan*numTasks*numThreads + curTask*numThreads + tid]+=fastqPart[i].kmerFreqCount[k];
					sendSizeInScan[curScan] += fastqPart[i].kmerFreqCount[k];
					k++;
				}
			}
		}
	}
	size_t maxSendSize=0, maxRcvSize=0;
	for (i=0; i<passes; i++)
	{
		if (sendSizeInScan[i] > maxSendSize)
			maxSendSize = sendSizeInScan[i];
	}
	free (sendSizeInScan);
	size_t *rcvSizeInScan = (size_t *)calloc(passes, sizeof(size_t));
	assert (rcvSizeInScan != NULL);
	for (i=0; i<totalFastqPartitions; i++)
	{
		for (j=0; j<passes; j++)
		{
			for (k=processKmerRanges[j*numTasks + rank]; k<processKmerRanges[j*numTasks + rank + 1]; k++)
			{
				rcvSizeInScan[j]+=fastqPart[i].kmerFreqCount[k];
			}
			if (rcvSizeInScan[j] > maxRcvSize)
				maxRcvSize = rcvSizeInScan[j];
		}
	}
	free (rcvSizeInScan);
	if (maxRcvSize > maxSendSize)
		maxSendSize = maxRcvSize;
	kmerIn = (uint32_t *)malloc (maxSendSize * 3 * sizeof(uint32_t));
	assert (kmerIn != NULL);
	if (numTasks > 1)
	{
		kmerOut = (uint32_t *)malloc (maxRcvSize * 3 * sizeof (uint32_t));
		assert (kmerOut != NULL);
		kmerSorted = kmerIn;
	}
	else
	{
		kmerSorted = (uint32_t *)malloc (maxRcvSize * 3 * sizeof (uint32_t));
		assert (kmerSorted != NULL);
		kmerOut = kmerIn;
	}
	fastqBuffer1 = (char *)malloc((size_t)numThreads*(fastqBufferSize+10));
	fastqBuffer2 = (char *)malloc((size_t)numThreads*(fastqBufferSize+10));
	assert ((fastqBuffer1 != NULL) && (fastqBuffer2 != NULL));

	parent = (uint32_t *)malloc(numReads * sizeof(uint32_t));
	compSizes = (uint32_t *)malloc(numReads * sizeof(uint32_t));
	assert ((parent != NULL) && (compSizes != NULL));
	uint32_t readNo=0;
	for (readNo=0; readNo<numReads; readNo++)
	{
		parent[readNo]=readNo;
		compSizes[readNo]=1;
	}
	uint64_t *localRangeCounts = (uint64_t *)calloc((numThreads+1), sizeof(uint64_t));
	assert(localRangeCounts != NULL);
	uint64_t *rangeOffsets = (uint64_t *)calloc(numThreads, sizeof(uint64_t));
	assert (rangeOffsets != NULL);	
	for (i=0; i<passes; i++)
	{
		kmerGen(i);
		bzero (localRangeCounts, (numThreads+1)*sizeof(uint64_t));
		bzero (rangeOffsets, numThreads*sizeof(uint64_t));
		localSort (i, localRangeCounts, rangeOffsets);
		localCC(i, localRangeCounts, rangeOffsets);
		localCCVerify(i, localRangeCounts, rangeOffsets);
		if (rank == 0)
			printf ("[%d] localCC Done\n", rank);
	}
	free (localRangeCounts);
	free (rangeOffsets);
	free (kmerIn);
	if (numTasks > 1)
	{
		free (kmerOut);
	}
	else
	{
		free (kmerSorted);
	}
	free (sendSizes);
	free (sendThreadSizes);
	free (scanKmerRanges);
	free (processKmerRanges);
	free (threadKmerRanges);
	free (merHist);

	//mergeComps
	uint32_t *otherRanksParent = (uint32_t *)calloc(numReads, sizeof(uint32_t));
	assert (otherRanksParent != NULL);

//MergeCC start
	struct timeval sTime, eTime;
#if USE_MPI
	int numMergeIters=0;
	int tasks=numTasks-1;
	while (tasks > 0)
	{
		numMergeIters ++;
		tasks /= 2;
	}
	int commIters = (int)(numReads / ASYNC_SIZE);
	if ((numReads % ASYNC_SIZE) != 0) commIters++;
	for (i=0; i<numMergeIters; i++)
	{
		gettimeofday (&sTime, NULL);
		int commDistance = 1<< i;
		int processDivider = commDistance<<1;
		if (rank % processDivider == commDistance)
		{
			int toRank = rank - commDistance;
        		for (j=0; j<(commIters-1); j++)
    				MPI_Send (&parent[j*ASYNC_SIZE], ASYNC_SIZE, MPI_UNSIGNED, toRank, j, MPI_COMM_WORLD);
    			MPI_Send (&parent[j*ASYNC_SIZE], (int)(numReads%ASYNC_SIZE), MPI_UNSIGNED, toRank, j, MPI_COMM_WORLD);
		}
		else if (rank % processDivider == 0)
		{
			MPI_Status status;
			for (j=0; j<(commIters-1); j++)
			    MPI_Recv (&otherRanksParent[j*ASYNC_SIZE], ASYNC_SIZE, MPI_UNSIGNED, rank+commDistance, j, MPI_COMM_WORLD, &status);
			MPI_Recv (&otherRanksParent[j*ASYNC_SIZE], (int)(numReads%ASYNC_SIZE), MPI_UNSIGNED, rank+commDistance, j, MPI_COMM_WORLD, &status);
			mergeStep(otherRanksParent);
		}
		gettimeofday (&eTime, NULL);
		long elapsed = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
		mergeTime += elapsed;
	}
#endif
//MergeCC end
	long writeComm=0;
	uint32_t maxCompId = 0;
	if (rank == 0)
	{
		maxCompId = printStats();
	}
#if USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	gettimeofday (&sTime, NULL);
	if (numTasks > 1)
	{
		int commIters = (int)(numReads / ASYNC_SIZE);
		if ((numReads % ASYNC_SIZE) != 0) commIters++;
		uint32_t ofs=0;
		for (j=0; j<(commIters-1); j++)
		{
			MPI_Bcast (&parent[ofs], ASYNC_SIZE, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	        	ofs+=ASYNC_SIZE;
		}
		MPI_Bcast (&parent[ofs], (int)(numReads%ASYNC_SIZE), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	}
	gettimeofday (&eTime, NULL);
	writeComm  = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
		MPI_Bcast (&maxCompId, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif
	writeOutput (maxCompId);
	gettimeofday (&eTime, NULL);
	outputTime = (eTime.tv_sec * 1000000 + eTime.tv_usec) - (sTime.tv_sec * 1000000 + sTime.tv_usec);
	if (rank == 0)
	{
		printStats2();
	}

	printf ("[%d] FQRead %lf sec(MaxKmerGen %lf sec), CommTime %lf sec, SortTime %lf sec(%lf + %lf), unionFind %lf sec (%lf + %lf), mergeComps %lf sec, WriteFq %lf sec(%lf sec)\n", rank, (double)(kmerGenTime-kmerGenWorkTime)/1000000, (double)kmerGenWorkTime/1000000, (double)commTime/1000000, (double)sortTime1/1000000 + (double)sortTime2/1000000, (double)sortTime1/1000000, (double)sortTime2/1000000, (double)unionFindTime1/1000000 + (double)unionFindTime2/1000000, (double)unionFindTime1/1000000, (double)unionFindTime2/1000000, (double)mergeTime/1000000, (double)outputTime/1000000 ,(double)writeComm/1000000);
	free (parent);
	free (otherRanksParent);
	free (compSizes);
	free (fastqBuffer1);
	free (fastqBuffer2);
	free (fastqRanges);
	free (fastqPart);
	for (i=0; i<argc-1; i++) {
		free(fileNames[i]);
	}
	free (fileNames);
	free (fileSize);
	if (recvThreadCount != NULL)
	{
		free (recvThreadCount);
		free (recvThreadOfs);
	}
	
	free (openFileNo);
	free (threadKmerGenTime);

	if (rank == 0)
	{
		int pid = getpid();
		char line[2048];
		sprintf (line, "/proc/%d/status", pid);
		FILE *statFile = fopen(line, "r");
	    	assert (statFile != NULL);
		fgets (line, 2048, statFile);
	    	while (!feof (statFile))
		{
			if (strstr(line,"VmPeak") || strstr(line,"VmHWM"))
	    		{
				printf ("[%d] %s", rank, line);
			}
	    		fgets (line, 2048, statFile);
		}
		fclose (statFile);
	}

#if USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
