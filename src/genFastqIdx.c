#include <stdio.h>
#include <stdint.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

typedef struct {
	size_t offset;
    size_t offset2;
    uint32_t readId;
	int fileNo;
    uint32_t kmerFreqCount[1048576+2];
} fileIndex;

int readSize=230;
int numPartitions=16;
size_t totalFastqSize = 0;
int numPEFiles, numSEFiles, numInterleavedFiles;
size_t *fileSize=NULL;
char **fileNames = NULL;

void getFileIndex (int chunkNum, size_t sizePerChunk, int numFiles, size_t *fileSize, fileIndex *result)
{
	int currentFileInd=0;
	size_t bytesRemaining = sizePerChunk * chunkNum;
	int fileFactor=1;
	if (currentFileInd < numPEFiles)
		fileFactor=2;
	while (bytesRemaining >= fileSize[currentFileInd]*fileFactor)
	{
		bytesRemaining -= fileSize[currentFileInd]*fileFactor;
		if (currentFileInd < numPEFiles)
		{
			currentFileInd+=2;
		}
		else
			currentFileInd++;
		if (currentFileInd < numPEFiles)
			fileFactor=2;
		else
			fileFactor=1;
	}
	result->fileNo=currentFileInd;
	result->offset=bytesRemaining/fileFactor;
}
int file2No=1;
size_t file2Ofs=0;
size_t findOfsOtherEnd (int fileNo, size_t offset, char *readName)
{
//	printf ("Begin findOfsOtherEnd %s\n", readName);
	size_t delta = fileSize[fileNo]*10/100;
	char *fileBuffer = (char *)malloc(2*delta);
	assert (fileBuffer != NULL);
	char curReadName[100];
	FILE *fp = fopen (fileNames[fileNo+1], "r");
        assert (fp != NULL);
	
	size_t startOfs=0,len=0;
	if (offset >= delta)
		startOfs = offset-delta;
	else
		startOfs = 0;	
	len = fileSize[fileNo+1] - startOfs;

	if ((offset + delta) < fileSize[fileNo+1])
		len = offset + delta - startOfs;
        fseek(fp, startOfs, SEEK_SET);
	fread (fileBuffer, 1, len, fp);
	fclose (fp);
	size_t i=0;
	
	while ((i<len) && ((fileBuffer[i] != '\n') || (fileBuffer[i+1] != '+')))
	{
		i++;
	}
	
	assert (i<len);
	i++;
	while (fileBuffer[i] != '\n') i++;
	i++;
	while (fileBuffer[i] != '\n') i++;
        i++;
	assert (i < len);
	
	int readNameInd=0;
	size_t j=i;
	while (fileBuffer[j] != '\n')
	{
		curReadName[readNameInd++] = fileBuffer[j];
		j++;
	}
	curReadName[readNameInd-2]='\0';

	while (i < len)
	{
		int comp=strcmp(curReadName, readName);
		if (comp == 0) break;
		while (fileBuffer[i] != '\n') {
			i++; 
		}
		i++;
		while (fileBuffer[i] != '\n') {
			i++; 
		}
		i++;
		while (fileBuffer[i] != '\n') {
			i++;
		}	
		i++;
		while (fileBuffer[i] != '\n') {
			i++;
		}	
		i++;
		readNameInd=0;
		j=i;
		while (fileBuffer[j] != '\n')
		{
			curReadName[readNameInd++] = fileBuffer[j];
			j++;
		}
		curReadName[readNameInd-2]='\0';
	}
//	printf ("\t2-%s Pos %lu\n", curReadName, startOfs+i);
	free (fileBuffer);
	return (startOfs+i);
}


int main(int argc, char **argv) {

	int maxReadSize=300;
	//FILE *fp;
	struct stat buf;
	//char readFile[500];

	if (argc < 3) {
		printf("Usage: %s numPartitions numPEFiles numSEFiles "
			"numInterleavedFiles a.1.fastq a.2.fastq b.1.fastq "
			"b.2.fastq ...\n", argv[0]);
		exit(1);
	}

	numPartitions = atoi(argv[1]);
	assert(numPartitions > 0);
	int numChunks=numPartitions;
	//int k;
	int i;
	numPEFiles = atoi(argv[2])*2;
	numSEFiles = atoi(argv[3]);
	numInterleavedFiles = atoi(argv[4]);
	argv = &argv[4];
	argc-=4;
//	printf ("Argc %d\n", argc);

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
//		printf ("File %s\n", fileNames[i-1]);
        	assert (status == 0);
		fileSize[i-1]=buf.st_size;
//		printf ("file %s Size %lu\n", fileNames[i-1], buf.st_size);
//		if (i <= numPEFiles)
//			if ((i % 2) == 0) continue;
		totalFastqSize += buf.st_size;
	}
//	printf ("size %lu\n", totalFastqSize);
	size_t sizePerChunk = totalFastqSize/numChunks;
	size_t j;
	
	char *readChunk = (char *)malloc(5000);

	long writeOfs = 0;
	//FILE *fpOp = NULL;
	char idxFile[200];
	sprintf (idxFile, "%s.idx", fileNames[0]);
	//long startOfs=0;
	FILE *fpOut = fopen(idxFile, "w");
	assert (fpOut != NULL);
	//int eof=0;
	int interleaved=0;
	//size_t curOfs=2;
	fileIndex fI;
	fI.fileNo=0;
	fI.offset=0;
#if 0
	fI.kmerFreqCount = (uint64_t *) malloc((1048576+2)*sizeof(uint64_t));
	assert(fI.kmerFreqCount != NULL);
#endif

	bzero (fI.kmerFreqCount, (1048576+2)*sizeof(uint32_t));
	int flag=0;	
	char curReadName[100];
	for (i=0; i<numChunks; i++)
	{
		flag=0;
		getFileIndex (i, sizePerChunk, argc-1, fileSize, &fI);
		if (fI.fileNo  >= (numPEFiles + numSEFiles))
			interleaved=1;
		if ((fileSize[fI.fileNo]-fI.offset) < maxReadSize)
		{
//			printf ("fileNo+1 %d, argc %d\n", fI.fileNo+1, argc);
			if (fI.fileNo+1 == argc-1) break;
			fI.fileNo++;
			fI.offset=0;
		}
//		printf ("chunk %d: File %d, Offset %lu\n", i, fI.fileNo, fI.offset);

		if (i==0) 
		{
			fI.offset2=0;
			fI.readId=0;
#if 0
			fwrite(&(fI.fileNo), sizeof(int), 1, fpOut);
			fwrite(&(fI.offset), sizeof(size_t), 1, fpOut);
			fwrite(&(fI.offset2), sizeof(size_t), 1, fpOut);
			fwrite(&(fI.readId), sizeof(uint32_t), 1, fpOut);
			fwrite(fI.kmerFreqCount, sizeof(uint64_t), 1048576+2,
					fpOut);
#endif
/*			printf ("size of fileindex is %lu, actual %lu\n", 
					sizeof(fileIndex), 
					sizeof(int)+2*sizeof(size_t)+sizeof(uint32_t)
					+1048578*sizeof(uint32_t));*/
			assert(sizeof(fileIndex) ==
					(sizeof(int)+2*sizeof(size_t)+sizeof(uint32_t)+1048578*sizeof(uint32_t)));

			fwrite (&fI, sizeof(fileIndex), 1, fpOut);
			continue;
		}
		FILE *fp = fopen (fileNames[fI.fileNo], "r");
		assert (fp != NULL);
		fseek(fp, fI.offset, SEEK_SET);
		//curOfs=startOfs;
		size_t nread =fread (readChunk, 1, 3000, fp);
//		printf ("nread %lu\n", nread);
		assert (nread > maxReadSize);
		int k=0;
		for (j=0; j<nread; j++)
		{
			if ((readChunk[j]=='\n') && (readChunk[j+1]=='+')/* && (readChunk[j+2]=='\n')*/)
			{
				flag=1;
				j=j+2;
				while (readChunk[j] != '\n')
					j++;
				j++;
			}
			k=0;
			if (flag && (readChunk[j]=='\n'))
			{
				j++;k++;
//				printf ("Ofs %lu, interleaved %d\n", fI.offset+j, interleaved);
				if (interleaved)
				{	
					for (k=j; readChunk[k]!='/'; k++);
					if (readChunk[k+1]=='2')
					{
						flag=0;
						j=k;
					}
					else
						break;
				}
				else
					break;
			}
		}
		if ((2*k+30) < maxReadSize) maxReadSize = 2*k+30;
		assert(flag==1);
		writeOfs = fI.offset+j;
//		printf ("writeOfs %lu\n", writeOfs);
		fI.offset = writeOfs;
		int readInd=0;
		while (readChunk[j]!='\n')
			curReadName[readInd++]=readChunk[j++];
//		assert (curReadName[readInd-1]=='1');
		curReadName[readInd-2]='\0';
//		printf ("FIle %d curReadName %s\n", fI.fileNo, curReadName);	
		if (fI.fileNo < numPEFiles)
		{	
			size_t ofs2=findOfsOtherEnd (fI.fileNo, fI.offset, curReadName);
			fI.offset2 = ofs2;
		}
		else
		{
			fI.offset2 = 0;
		}
		fI.readId=0;

#if 0
		fwrite(&fI.fileNo, sizeof(int), 1, fpOut);
		fwrite(&fI.offset, sizeof(size_t), 1, fpOut);
		fwrite(&fI.offset2, sizeof(size_t), 1, fpOut);
		fwrite(&fI.readId, sizeof(uint32_t), 1, fpOut);
		fwrite(fI.kmerFreqCount, sizeof(uint64_t), 1048576+2,
					fpOut);
#endif
		fwrite (&fI, sizeof(fileIndex), 1, fpOut);
//		printf ("P%d - %d %ld\n", rank, count, startOfs+i);
//		printf ("Chunk %d: File %d, Offset %lu,%lu\n", i, fI.fileNo, fI.offset, fI.offset2);
		fclose (fp);
	}
	if ((argc-1) == numPEFiles)
	{
		fI.fileNo = argc-3;
		fI.offset = fileSize[argc-3];
		fI.offset2 = fileSize[argc-2];
	}
	else
	{
		fI.fileNo = argc-2;
		fI.offset = fileSize[argc-2];
		fI.offset2 = 0;
	}
	fI.readId=0;

#if 0
	fwrite(&fI.fileNo, sizeof(int), 1, fpOut);
	fwrite(&fI.offset, sizeof(size_t), 1, fpOut);
	fwrite(&fI.offset2, sizeof(size_t), 1, fpOut);
	fwrite(&fI.readId, sizeof(uint32_t), 1, fpOut);
	fwrite(fI.kmerFreqCount, sizeof(uint64_t), 1048576+2,
				fpOut);
#endif

	fwrite (&fI, sizeof(fileIndex), 1, fpOut);
//	printf ("Chunk %d: File %d, Offset %lu,%lu\n", i, fI.fileNo, fI.offset, fI.offset2);
	free (readChunk);
	for (i=0; i<(argc-1); i++)
	{
		free (fileNames[i]);
	}
	free(fileNames);
#if 0
	free(fI.kmerFreqCount);
#endif
	free (fileSize);
	fclose (fpOut);

	return 0;
}
