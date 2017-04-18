echo "Compiling code..."
gcc -o bin/genFastqIdx src/genFastqIdx.c
gcc -o bin/genHistograms src/genHistograms.c src/countMmers.c -mavx
mpicc -O3 -fopenmp -o bin/metaprep src/metaprep.c src/enumerateKmers.c -mavx
echo "Compile Done."

echo "Creating FASTQ index file"
bin/genFastqIdx 2 1 0 0 data/a1.fq data/a2.fq
echo "Created 2 chunks"

echo "Creating merHist and FASTQPart tables"
bin/genHistograms 27 1 0 0 data/a1.fq data/a2.fq

echo "Run MetaPrep using k=27"
bin/metaprep -o ./output 30000 27 1 2 1 0 0 data/a1.fq data/a2.fq
echo "Done. Output FASTQ files present in output directory"

echo "Cleaning index files"
rm data/a1.fq.*
