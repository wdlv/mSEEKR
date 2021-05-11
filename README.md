# mSEEKR

This is a program for identifying regions of high similarity based on *k*-mer content to some set of query sequences, relative to a null background set of sequences.

## Dependencies

#### Python3.6
check with command 

```
python -V
```

if < 3.6, easiest solution is to download anaconda python3.6/3.7
#### cython

Try 
```
which cython
``` 
if none

```
pip install cython
```

OR

```
conda install -c anaconda cython
```
## Installation

#### Type following commands in desired directory
```
	git clone https://github.com/spragud2/mSEEKR.git
	cd mSEEKR/
	python setup.py build_ext --inplace
```
#### Ignore warnings (unless further steps don't work)
<hr/>

## Counting k-mers 


  1. Curate unique fasta files for queries and null model before hand
  2. Use the following command
```
  python kmers.py --fasta ./fastaFiles/mA.fa -k 2,3,4 --name mouseA --dir ./counts/
  python kmers.py --fasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa -k 2,3,4 --name mm10Trscpts --dir ./counts/
```

#### Parameters:

1. --fasta : Path to fasta file. (Required)
2. -k : Comma delimited string of possible k-mer values. For example, 3,4,5 or just 4. (Required)
3. --name : Desired output name for count file. Default is 'out'
4. --dir : Directory to save output count file. Default is './' 
5. -a : String, Alphabet to generate k-mers (e.g. ATCG). Default is 'ATCG'


  Output:

  Outputs binary .dict files containing count matrices to directory set by --dir

<hr/>

## Training markov models

  0. Count k-mers for queries and null before proceeding  
  1. Run the following command
```
  python train.py --query ./counts/mouseA.dict --null ./counts/mm10Trscpts.dict -k 2,3,4 --qPrefix mouseA --nPrefix mm10Trscpts --qT .9999 --nT .9999 --dir ./markovModels/
```

Parameters:

1. --query : Path to kmer count file for sequences of interest (e.g. functional regions of a ncRNA). (Required)
2. --null : Path to kmer count file that compose null model (e.g. transcriptome, genome, etc). (Required)
3. -k : Comma separated list of values of k to train for, must have been calculated prior. (Required)
4. --qPrefix : prefix file name for query. (Required)
5. --nPrefix : prefix file name for null. (Required)
6. --qT : Query to query transition parameter, query to null is 1 - qT. (Required)
7. --nT : Null to null transition parameter, null to query is 1 - nT. (Required)
8. --dir : Output directory. Default is './'
9. -a : String, Alphabet to generate k-mers (e.g. ATCG). Default is 'ATCG'

  Output:

  Directory containing models for each value of k specified, the directory structure would look like:

    | markovModels
    |
    |--- mouseA_mouseNull
    |------- 2
    |------------- hmm.dict
    |------- 3
    |------------- hmm.dict
    |------- 4
    |------------- hmm.dict
    |--- mouseB_mouseNull
    .
    .
    .

## (OPTIONAL) Optimize transition parameters using the B-W Algorithm
This script will take the output .mkv output file from train.py and find an MLE for the transition parameters. 

  1. Run the following command
```
  python bw.py -k 4 --db ./fastaFiles/xist.fa --prior markovModels/mouseA_mm10Trscpts/4/hmm.dict --its 3 -cf
```

The above command takes the HMM trained on repeat A at k =4, and uses Xist to find an MLE for the transition parameters. More than one sequence can be provided as training, if a multi-entry FASTA file is provided. 

Parameters:

1. -k : Value of k for k-mers. One integer is requried. (Required)
2. --db : Path to fasta file containing training sequences. (Required)
3. --prior : Path to binary .dict file output from train.py(e.g. markovModels/D_null/2/hmm.dict). (Required)
4. --its : Number of iterations of the baum-welch algorithm to run. Default is 20
5. -cf : FLAG, pass this argument to create a new file in the same directory as --prior rather than overwrite

Output:
1. Replaces the file passed as --prior with an updated version that contains the MLE transition matrix
2. A file containing the trajectory of the transition parameters over each loop of the BW algorithm, check this for convergence




## Find HMM state path through sequences of interest

  1. Run the mSEEKR command
```
  python mSEEKR.py --db ./fastaFiles/mm10kcn.fa --prefix kcn_queryMouseA --model ./markovModels/mouseA_mm10Trscpts/4/hmm_MLE.dict -k 4 --fasta
```

Parameters

1. --db : Path to fasta file with sequences to calculate similarity score. (Required)
2. --model : Path to .dict file output from train.py or bw.py. (Required)
3. -k : Value of k to use. Must be the same as the k value used in training. (Required)
4. --prefix : File name for output, useful to include information about the experiment. (Required)




## Calculate background fasta's seekr pearson correlation score's mean and standard deviation

  2. Run the getBackgroundMatrixMeanStd.py command
  ```
  python getBackgroundMatrixMeanStd.py --backgroundFasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa -k 4
```

Parameters

1. --backgroundFasta : Path to lncRNA background sequences fasta file; used to calculate mean and standard deviation of k-mer seekr score matrix. (Required)
2. -k : The same k value used in the python mSEEKR.py step. (Required)
3. --name : Name for the file output background fasta kmer seekr score matrix mean and std. Default is 'bgMatrixMeanStd'
4. --outdir : Directory to save output background fasta kmer seekr score matrix mean and std. Default is './'
5. -a : String, Alphabet to generate k-mers (e.g. ATCG). Default is 'ATCG'





## Calculate seekr pearson correlation score

  3. Run the getSeekrScore command
  ```
  python getSeekrScore.py --queryFasta ./fastaFiles/mA.fa --backgroundFasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa --backgroundMatrixMeanStd ./bgMatrixMeanStd.txt --mSEEKRdataframeDir ./kcn_queryMouseA_4_viterbi.txt -k 4
```

Parameters

1. --queryFasta : Path to query fasta file. (Required)
2. --backgroundFasta : Path to lncRNA background sequences fasta file. (Required)
3. --backgroundMatrixMeanStd : Path to lncRNA background matrix mean and std generated from getBackgroundMatrixMeanStd.py. Note: backgroundMatrixMeanStd and backgroundFasta are using the same fasta file. The reason there are two commands is that computing mean and std of background fasta kmer seekr score requires too much resources. (Required)
4. --mSEEKRdataframeDir : Directory to read in the output dataframe generated from mSEEKR.py. (Required)
5. -k : The same k value used in the mSEEKR.py step. (Required)
6. --dir : Directory to save output dataframe. Default is './'
7. --name : Name for output file. Default is 'output'
8. --minSeqLength : The minimum length of sequence found. Default is 0
9. -a : String, Alphabet to generate k-mers (e.g. ATCG). Default is 'ATCG'

