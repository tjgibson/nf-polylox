# nf-polylox

## Introduction
This repository contains a nextflow pipeline for processing PacBio data from lineage tracing experiments using the PolyloxExpress barcoding system.

For more information on the Polylox system, see the following papers:

Pei, W., Feyerabend, T., Rössler, J. et al. Polylox barcoding reveals haematopoietic stem cell fates realized in vivo. Nature 548, 456–460 (2017). <https://doi.org/10.1038/nature23653>

Pei, W., Wang, X., Rössler, J. et al. Using Cre-recombinase-driven Polylox barcoding for in vivo fate mapping in mice. Nat Protoc 14, 1820–1840 (2019). <https://doi.org/10.1038/s41596-019-0163-5>

Pei W., Shang F., Wang X. et al. Resolving Fates and Single-Cell Transcriptomes of Hematopoietic Stem Cell Clones by PolyloxExpress Barcoding. Cell Stem Cell 27, 383-395.e8 (2020). <https://doi.org/10.1016/j.stem.2020.07.018>

## Pipeline
This pipeline executes the following steps:

1. Convert PacBio bam files to fastq format ([`pbtk bam2fastq`](https://github.com/pacificbiosciences/pbtk/))
2. Split fastq files into chunks to speed up downstream processing ([`seqkit split2`](https://bioinf.shenwei.me/seqkit/))
3. Parse the Polylox barcodes from PacBio reads ([`scRPBPBR`](https://github.com/sunlightwang/PolyloxExpress))
4. Merge parsed Polylox barcodes back together for each sample
5. Compute the generation probability (pgen) for each barcode ([`polyloxpgen`](https://github.com/mauricelanghinrichs/polyloxpgen))

## Usage
First, prepare a samplesheet with your input data in the following format:

**samplesheet.csv**:

```csv
sample_id,reads,read_index
test_bam,data/test/test.bam,data/test/test.bam.pbi
test_fastq,data/test/test.fastq.gz,
```

Each row corresponds to a file containing PacBio reads from a PolyloxExpress barcoding experiment. 
Files should contain CCS/hi-fi reads, not raw reads. 
Files can be in BAM or FASTQ format and should be listed in the `reads` column. 
For bam files, provide the corresponding `.pbi` index file. 
For FASTQ files, the `read_index` field can be left blank. 
FASTQ files should be GZIP compressed and end in the file extension `.fastq.gz`.
The `sample_id` column should contain a meaningful sample identifier.

Now, you can run the pipeline using:

```bash
nextflow run https://github.com/tjgibson/nf-polylox \
    --samplesheet samplesheet.csv \
    --results_dir <RESULTSDIR> \
    -profile <docker/singularity/...>
```

You can optionally change how FASTQ files are split into chunks by specifying `--fastq_splitSize <n_reads>`, where `n_reads` is how many reads should be in each chunk.
The default value is set to 5000, which gave good performance in testing, but may need to be modified based on the execution environment. 
Setting this value to a small number will cause many intermediate files to be written. On HPC systems, each intermediate file will be processed as a separate job. 
Consider how you set this parameter and how it will affect the number of files written and number of jobs submitted.


## Pipeline output
For each input sample, this pipeline will produce three results files:
1. **sample.seg_assemble.tsv**

A table containing three columns: PacBio read name, corresponding scRNA-seq cell barcode, Polylox barcode.

2. **sample_merged.stat.tsv**

A Table containing two columns: Polylox barcode, n_reads. This table summarizes number of PacBio reads with a given polylox Identity.

3. **Sample_pgen.txt**

Table containing the generation probability for each Polylox barcode.

Examples of output files can be found under `test_results`.

## To-do
aIn the future, it might be nice to add some of the following features:
- FASTQC/multiQC
- Filtering and QC of polylox barcodes
- Summary metrics, tables, or figures (e.g. PacBio read length distribution, number of unique polylox barcodes, polylox barcode length distribution, barcode pgen, etc.)
- Integration with scRNA-seq data? Likely makes most sense to process scRNA-seq data separately, but this pipeline could include integration of the two data types, possibly by providing a seurat object as an input to this pipeline.