#!/bin/bash
python kmers.py --fasta ./fastaFiles/mouseA.fa -k 2,3,4,5,6 --name mouseA -n 5
python kmers.py --fasta ./fastaFiles/mouseB.fa -k 2,3,4,5,6 --name mouseB -n 5
python kmers.py --fasta ./fastaFiles/mouseD.fa -k 2,3,4,5,6 --name mouseD -n 5
python kmers.py --fasta ./fastaFiles/mouseE.fa -k 2,3,4,5,6 --name mouseE -n 5
python kmers.py --fasta ../mm10Transcriptome/mm10transcripts.fa -k 2,3,4,5,6 --name mm10Trscpts -n 5
# python kmers.py --fasta ./fastaFiles/noRatAE.fa -k 2,3,4,5,6 --name noRatAE -n 4
# python kmers.py --fasta ./fastaFiles/noRatBCD.fa -k 2,3,4,5,6 --name noRatBCD -n 4
# python kmers.py --fasta ./fastaFiles/noDogBCD.fa -k 2,3,4,5,6 --name noDogBCD -n 4
# python kmers.py --fasta ./fastaFiles/noDogAE.fa -k 2,3,4,5,6 --name noDogAE -n 4
# python kmers.py --fasta ./fastaFiles/noCowAE.fa -k 2,3,4,5,6 --name noCowAE -n 4
# python kmers.py --fasta ./fastaFiles/noCowBCD.fa -k 2,3,4,5,6 --name noCowBCD -n 4
# python kmers.py --fasta ./fastaFiles/noBoarBCD.fa -k 2,3,4,5,6 --name noBoarBCD -n 4
# python kmers.py --fasta ./fastaFiles/noBoarAE.fa -k 2,3,4,5,6 --name noBoarAE -n 4
# python kmers.py --fasta ./fastaFiles/noMarmAE.fa -k 2,3,4,5,6 --name noMarmAE -n 4
# python kmers.py --fasta ./fastaFiles/noMarmBCD.fa -k 2,3,4,5,6 --name noMarmBCD -n 4
# python kmers.py --fasta ./fastaFiles/noPigBCD.fa -k 2,3,4,5,6 --name noPigBCD -n 4
# python kmers.py --fasta ./fastaFiles/noPigAE.fa -k 2,3,4,5,6 --name noPigAE -n 4
# python kmers.py --fasta ./fastaFiles/mamBCD.fa -k 2,3,4,5,6 --name mamBCD -n 4
# python kmers.py --fasta ./fastaFiles/mamAE.fa -k 2,3,4,5,6 --name mamAE -n 4
# python kmers.py --fasta ./fastaFiles/mamAE.fa -k 2,3,4,5,6 --name mamAE -n 4
# python kmers.py --fasta ../hgTranscriptome/hg26Transcripts.fa -k 2,3,4,5,6 --name hg26Transcripts -n 4