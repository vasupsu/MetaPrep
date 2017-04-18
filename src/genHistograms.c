#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <sys/time.h>

typedef struct {
        size_t offset;
        size_t offset2;
        uint32_t readId;
        int fileNo;
        uint32_t kmerFreqCount[1048576+2];
}fileIndex;

char **fileNames=NULL;
size_t *fileSize=NULL, totalFastqSize=0;
uint64_t *rangeHist=NULL;
uint32_t curReadId=0;
size_t fastqBufferSize=500000000;
int numPEFiles, numSEFiles, numInterleavedFiles;
uint64_t *kmerList=NULL;

char *fqBuffer[2]={NULL,NULL};
size_t fqInd[2]={0,0};
size_t fqSize[2]={0,0};

int file1=0, file2=1;
struct stat buf;
uint32_t K=0;

uint64_t kmerMask=0;
int rotater=0;
uint64_t totalKmers=0;
void getKmerMask(uint32_t k)
{
	kmerMask = (((1ul<<k)-1) >> (k-20)) << (32+k-20);
	rotater = 32+k-20;
}

extern uint64_t getCanonKmers (char *readStr, int readLen, uint64_t *rangeHist, uint64_t kmerMask, int rotater, uint32_t *kmerFreqCount);
extern void getFullStr (uint64_t kmer, FILE *fp);
extern void init(int kLen);
FILE *fp1_copy=NULL, *fp2_copy;
size_t ofs_copy=0;
int eof=0;
void bufSeek (int bufNo, FILE *fp, long ofs)
{
	fqSize[bufNo]=fqInd[bufNo]=0;
	assert (fseek(fp, ofs, SEEK_SET) == 0);
	
	fqSize[bufNo]=fread(fqBuffer[bufNo], 1, fastqBufferSize, fp);
	if (bufNo == 0)
		ofs_copy=ofs;
}
char * getNextLine(int bufNo, FILE *fp, int *len, int print)
{
	size_t curInd = fqInd[bufNo];
	if (bufNo == 0)
		eof=0;
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
				if (bufNo==0)
					eof=1;
			}
		}
		if (eof) break;
		curInd++;
	}
	if ((eof==0) && (fqInd[bufNo] < fqSize[bufNo]))
	{
		*(fqBuffer[bufNo]+curInd) = '\0';
		*len=curInd-fqInd[bufNo];
		char *retPtr = fqBuffer[bufNo]+fqInd[bufNo];
		if (bufNo == 0)
		{
			ofs_copy += (curInd+1-fqInd[bufNo]);
			if (ofs_copy == fileSize[file1])
			{
				eof=1;
			}
		}
		fqInd[bufNo]=curInd+1;
		return retPtr;
	}
	return NULL;
}

int main(int argc, char **argv)
{
	char *line1_copy=NULL, *line2_copy=NULL;
        int i, k, len;
	if  (argc < 5)
        {
                printf ("Usage: %s kmerSize numPEFiles numSEFiles "
						"numInterleavedFiles a.1.fastq a.2.fastq "
						"b.1.fastq b.2.fastq ...\n", argv[0]);
                exit(1);
        }
	rangeHist = (uint64_t *)calloc (1048576+2, sizeof(uint64_t));
	assert (rangeHist != NULL);
	K=(uint32_t)atoi(argv[1]);
	init(K);
	kmerList = (uint64_t *)malloc(200 * sizeof(uint64_t));
	assert (kmerList != NULL); 
	getKmerMask(K);
	fqBuffer[0] = (char *)malloc(fastqBufferSize+10);
	fqBuffer[1] = (char *)malloc(fastqBufferSize+10);
	assert ((fqBuffer[0] != NULL) && (fqBuffer[1] != NULL));

        numPEFiles = atoi(argv[2])*2;
        numSEFiles = atoi(argv[3]);
        numInterleavedFiles = atoi(argv[4]);
        argv = &argv[4];
        argc-=4;
        fileNames = (char **)malloc((argc-1) * sizeof(char *));
        fileSize = (size_t *)malloc((argc-1) * sizeof(size_t));
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
		totalFastqSize += buf.st_size;
        }
	char idxFile[200];
	char outFile[200];
	sprintf (outFile, "%s.KmerIdx", argv[1]);
	FILE *fp_out = fopen(outFile, "w");
	assert (fp_out != NULL);
	sprintf (idxFile, "%s.idx", argv[1]);
        int status = stat (idxFile, &buf);
        assert (status == 0);
        int numRecs = buf.st_size/sizeof(fileIndex);
        fileIndex *fastqPart = (fileIndex *)malloc(numRecs * sizeof(fileIndex));
        assert (fastqPart != NULL);
        FILE *fp = fopen (idxFile, "r");
        fread(fastqPart, sizeof(fileIndex), numRecs, fp);
        fclose (fp);

	size_t ofs1=0,ofs2=0;
	fp1_copy = fopen (fileNames[file1], "r");
	assert ((fp != NULL) && (fp1_copy != NULL));
	if (file2 < numPEFiles)
	{
		fp2_copy = fopen (fileNames[file2], "r");
		assert (fp2_copy != NULL);
	}
        for (i=0; i<numRecs-1; i++)
	{
		ofs1 = fastqPart[i].offset;
		ofs2 = fastqPart[i].offset2;
		file1 = fastqPart[i].fileNo;
		bufSeek (0, fp1_copy, ofs1);
		if (fp2_copy != NULL)
		{
			file2 = fastqPart[i].fileNo+1;
			bufSeek (1, fp2_copy, ofs2);
		}
		else
			file2 = -1;
		fastqPart[i].readId = curReadId;
		line1_copy=getNextLine (0, fp1_copy, &len, 0);
		if (fp2_copy != NULL)
		{
			line2_copy=getNextLine (1, fp2_copy, &len, 0);
		}
		while ((file1 < fastqPart[i+1].fileNo) || ((file1 == fastqPart[i+1].fileNo) && (ofs_copy < fastqPart[i+1].offset)))
		{
			line1_copy=getNextLine (0, fp1_copy, &len, 0);
			totalKmers += (uint64_t)getCanonKmers(line1_copy, len, rangeHist, kmerMask, rotater, fastqPart[i].kmerFreqCount);
			//genKmers
			if (fp2_copy != NULL)
			{
				line2_copy=getNextLine (1, fp2_copy, &len, 0);
				totalKmers += (uint64_t)getCanonKmers(line2_copy, len, rangeHist, kmerMask, rotater, fastqPart[i].kmerFreqCount);
			}
			//genKmers
			line1_copy=getNextLine (0, fp1_copy, &len, 0);
			line1_copy=getNextLine (0, fp1_copy, &len, 0);
			if (eof==0)
			{
				{
					line1_copy=getNextLine (0, fp1_copy, &len, 0);
				}
			}
			if (fp2_copy != NULL)
			{
				line2_copy=getNextLine (1, fp2_copy, &len, 0);
				line2_copy=getNextLine (1, fp2_copy, &len, 0);
				if (eof==0)
				{
					line2_copy=getNextLine (1, fp2_copy, &len, 0);
				}
			}
			if (file1 >= (numPEFiles + numSEFiles))
			{
				line1_copy=getNextLine (0, fp1_copy, &len, 0);
				totalKmers += (uint64_t)getCanonKmers(line1_copy, len, rangeHist, kmerMask, rotater, fastqPart[i].kmerFreqCount);
				//genKmers
				line1_copy=getNextLine (0, fp1_copy, &len, 0);
				line1_copy=getNextLine (0, fp1_copy, &len, 0);
				if (eof==0)
				{
					line1_copy=getNextLine (0, fp1_copy, &len,0);
				}
			}
			curReadId++;
			if (eof)
			{
				fclose (fp1_copy);
				if (fp2_copy != NULL)
				{
					fclose(fp2_copy);
					file1+=2;
					if (file1 < (argc-1))
						fp1_copy = fopen (fileNames[file1], "r");
					else
						fp1_copy=NULL;
					file2+=2;
					if (file2 < numPEFiles)
					{
						fp2_copy = fopen (fileNames[file2], "r");
						if (fp2_copy != NULL)
						{
							bufSeek (1, fp2_copy, 0);
							line2_copy=getNextLine (1, fp2_copy, &len, 0);
						}
					}
					else
						fp2_copy = NULL;
				}
				else
				{
					file1++;
					fp1_copy = fopen (fileNames[file1], "r");
				}
				if (fp1_copy != NULL)
				{
					bufSeek (0, fp1_copy, 0);
					line1_copy=getNextLine (0, fp1_copy, &len, 0);
				}
			}
		}
	}
	fastqPart[i].readId = curReadId;
	fwrite (fastqPart, sizeof (fileIndex), numRecs, fp_out);
	for (i=1; i<=1048576; i++)
	{
		rangeHist[i]=rangeHist[i-1]+rangeHist[i+1];
	}
	assert (rangeHist[1048576]==totalKmers);
	char histFile[200];
	sprintf (histFile, "%s.kmerHist", argv[1]);
	FILE *fp_hist = fopen(histFile, "w");
	assert (fp_hist != NULL);
	fwrite (rangeHist, sizeof(uint64_t), 1048576+1, fp_hist);
	fclose (fp_hist);

	uint64_t newTotalKmers=0;
	for (i=0; i<numRecs; i++)
	{
		for (k=0; k<1048576; k++)
		{
			newTotalKmers += fastqPart[i].kmerFreqCount[k];
		}
	}
	assert (newTotalKmers == totalKmers);
	free (fastqPart);
	free (fqBuffer[0]);
	free (fqBuffer[1]);
	free (kmerList);
	free (rangeHist);
	fclose (fp_out);
	for (i=0; i<(argc-1); i++)
	{
		free (fileNames[i]);
	}
	free(fileNames);
#if 0
	fclose(fp1_copy);
	fclose(fp2_copy);
#endif
	free(fileSize);

	return 0;
}
