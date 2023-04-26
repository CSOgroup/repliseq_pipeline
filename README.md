# Repli-seq processing pipeline

> :warning: **For the time being, this pipeline can be run only on the <u>graviton</u> server**: In-progress work to make the pipeline working on other machines!

## About the pipeline
This pipeline is heavily inspired by the [Gilbert pipeline](https://www.nature.com/articles/nprot.2017.148) and the [shart pipeline](https://www.github.com/dvera/shart). For any information about the internal workings of the pipeline, please check these references. 

The main usage of this pipeline is to run the analysis of multiple Repli-seq samples and organize their results in a coherent way.

## Installing the pipeline
Clone the repository and enter in the folder:

```
git clone https://github.com/CSOgroup/repliseq_pipeline.git
cd repliseq_pipeline
```

Install the dependencies using conda/mamba, creating a new environment (`repliseq_pipeline`):
```
conda env create -f environment.yml
```

Test that the pipeline works by running:
```
./run_repliseq_pipeline.sh.sh test_input.csv
```

## Processing repli-seq samples
This is the main step of the pipeline.
```
./run_repliseq_pipeline.sh <input-samples.csv>
```

### Input format
The `run_repliseq_pipeline.sh` script accepts a .csv file (with header) as input having one line for each sample to be processed. The required columns are:
- *sample_path*: path to the sample results
- *raw_path*: path to the sample fastq files. Fastq files are assumed to be paired-ended and located in the same folder. Read1 and Read2 are denoted by `_R1_` and `_R2_` inside the file names.
- *genome_assembly*: genome assembly (`hg19`, `mm10`, etc...)
- *genome_sequence*: path to the fasta file for the reference genome
- *chromsizes*: path to the chromosome sizes file for the reference genome
- *blacklisted_regions* (optional): path to a BED file containing regions that should be excluded from the analysis (like [these](https://github.com/Boyle-Lab/Blacklist/tree/master/lists))

You can check the [test_input.csv](./test_input.csv) file for reference.

### Steps
The pipeline will run the following analyses:
1. Fastq quality control with `fastqc`
2. Read alignment, filtering, sorting and statistics (bwa, samtools)
3. Read deduplication (samtools)
4. Computing repli-seq coverage (number of reads) in equally spaced bins (10Kb), giving the output as a `.bedGraph`