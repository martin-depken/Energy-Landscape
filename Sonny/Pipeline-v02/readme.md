# Whole-Genome Off-Target Prediction

This is a pipeline that contains the algorithms for whole-genome off-target prediction of CRISPR/Cas systems.

## Getting Started

In order to get the code working, place all necessary files in one folder on the TU Delft cluster hpc05, or on a computer with Python 3. Remember the explicit path to this folder, for example ```\home\sfdejong\11_03_2019\```.

### Prerequisites

The following Python 3 packages are required: 1. [numpy](http://www.numpy.org/); 2. [h5py](https://www.h5py.org/), which deals with HDF5 file formats and their compression; 3. [Biopython](https://biopython.org/), which deals with genomic file formats like FASTA; 4. [pickle](https://docs.python.org/3/library/pickle.html), which can save dictionaries in separate files.

### Python 2
If you really need to use Python 2, make sure to run the ```startup.py``` script in Python 2 too, so that the lookup dictionaries will have supported pickle protocols.

## How the algorithm works

The algorithm scans the entire genome from 5' to 3', applies a sequence-dependent, location-variant model, computes the expected cleavage time and stores it in a dataset.

### Convention for guide RNA

Because the guide is typically published from 5' to 3' in the nonhybridising backbone, there are two implications for our algorithm. To fit with our convention, all nucleotides should be complemented to RNA (A&rarr;U, T&rarr;A etc.), which is carried out using ```.transcribe()``` in Biopython. For Cas9, the on-target is also reversed. These specifics can be modified per Cas type in ```Class CRISPR``` from ```genomicpipeline.py```.

It does not matter whether the guide is provided in RNA or DNA alphabet, as long as it is given from 5' to 3' with respect to the nontarget strand.

### Orientation of target DNA

The algorithm is constructed such, that the variable ```target``` always contains the hybridising sequence from seed to end, regardless whether this is 5' to 3' or not. However, the genomic FASTA file is interpreted as if it is the targeted, hybridising strand from 5' to 3'. Therefore, the  ```partition()``` function depends on the value of ```Cas._5primeseed_wrt_target```, which indicates if the Cas starts hybridising from 5' or 3' (with respect to hybridising strand).

### Orientation of PAM

The variable ```PAM``` is constructed in harmony with the above convention. This implies that its first letter is alwways furthest from the target sequence and its last letter is adjacent.

### Example of the above 

given guide: ```GGGTGGGGGGAGTTTGCTCC``` (from 5' to 3' on nontarget strand)

Cas9: seed is at 3' with respect to target strand

match in FASTA file : ```GGGTGGGGGGAGTTTGCTCCTGG``` (from 5' to 3')

variable ```target``` contains: ```GGAGCAAACTCCCCCCACCC```

variable ```guide``` contains: ```GGAGCAAACUCCCCCCACCC```

variable ```PAM``` contains: ```GGT```


### Something else

Please, observe that the very first and last letters in the FASTA file will be excluded because of the way Python slices. To prevent this, an overly complicated construction using ```or None``` and ```1*(not position)``` was required, which we deemed unnecessary.

### HDF5 format

The data is stored in HDF5 format, which is compressed and has the nice property that Python can read parts of the dataset without too much RAM. 

### FASTA format

The same advantage is to reading out FASTA formats.

## Running on hpc05

### Prework

The ```startup.py``` file is not computationally heavy and creates the dictionary files ```lookuptable_PAM.p``` and ```lookuptable_target.p```. If you do not have these dictionaries already, create them and upload them to the cluster.

### Uploading to hpc05

If the above has been carried out, the following files should be uploaded to the cluster: ```functions.py```, ```genomicpipeline.py```, ```kinetic_model.py```, ```mainpath.txt```, ```jobs.txt```, ```run.sh```. ```chromosomes/```. The latter should contain at least one FASTA file.


### Preparing to run

In the ```jobs.txt``` file, declare all the specifics of possibly multiple runs of the pipeline. Next, write the explicit path (see Getting Started) in the file ```mainpath.txt```. The algorithm will use this to navigate.

### Running
The following command runs your code.

```
qsub run.sh -t 1-1%1 -N tst-sfdejong -e /home/sfdejong/11_03_2019/error.txt -o /home/sfdejong/11_03_2019/terminal.txt
```

## Licence

This project is licensed under Creative Commons BY-NC-SA by S.F. de Jong.

