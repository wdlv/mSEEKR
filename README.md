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
  The mSEEKR package has used functions from seekr so please install the seekr package first.
  It is lcoated at https://github.com/CalabreseLab/seekr
  
  Then lets install mSEEKR package:
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

1. --fasta : path to fasta file
2. -k : comma separated list of values of k
3. --name : name to give output file
4. --dir : output directory 


  Output:

  Outputs binary .skr (seekr) files containing count matrices to --dir

<hr/>

## Training markov models

  0. Count k-mers for queries and null before proceeding  
  1. Run the following command
```
  python train.py --query ./counts/mouseA.dict --null ./counts/mm10Trscpts.dict -k 2,3,4 --qPrefix mouseA --nPrefix mm10Trscpts --qT .9999 --nT .9999 --dir ./markovModels/
```

Parameters:

1. --query : Path to query count file
2. --null : Path to null count file
3. -k : comma separated list of values of k to train for, must have been calculated prior
4. --qPrefix : prefix file name for query
5. --nPrefix : prefix file name for null
6. --qT : Query to query transition parameter, query to null is 1 - qT
7. --nT : Null to null transition parameter, null to query is 1 - nT
8. --dir : output directory

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

1. -k : value of k for k-mers
2. --db : path to FASTA file to optimize parameters on
3. --prior : hmm.dict file output by train.py
4. --its : number of iterations of the baum-welch algorithm to run
5. -cf : FLAG, pass this argument to create a new file in the same directory as --prior 

Output:
1. Replaces the file passed as --prior with an updated version that contains the MLE transition matrix
2. A file containing the trajectory of the transition parameters over each loop of the BW algorithm, check this for convergence




## Find HMM state path through sequences of interest

  1. Run the mSEEKR command
```
  python mSEEKR.py --db ./fastaFiles/mm10kcn.fa --prefix kcn_queryMouseA --model ./markovModels/mouseA_mm10Trscpts/4/hmm_MLE.dict -k 4 --fasta
```

Parameters

1. --db : sequences of interest to run the HMM on
2. --model : Path to .dict file output from train.py or bw.py
3. -k : value of k to use in the analysis (must have been specified in training)
4. -n : Number of processing cores to use. Scales with size of fasta file (# of entries, not sequence length)
5. --prefix : file name for output, useful to include information about the experiment


## Calculate seekr pearson correlation score

  2. Run the getSeekrScore command
  ```
  python getSeekrScore.py --queryFasta ./fastaFiles/mA.fa --backgroundFasta ./fastaFiles/gencode.vM17.lncRNA_transcripts.fa --mSEEKRdataframeDir ./kcn_queryMouseA_4_viterbi.txt -k 4
```

Parameters

1. --queryFasta : query fasta files used at the beginning
2. --backgroundFasta : lncRNA background sequences fasta file; used to calculate mean and standard deviation for each k-mer
3. --mSEEKRdataframeDir : directory to read in the output dataframe generated from command python mSEEKR.py
4. -k : The same k value used in the python mSEEKR.py step
5. --dir : Directory to save output seekr score dataframe

