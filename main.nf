#! /usr/bin/env nextflow

log.info """
	PolyLox barcode processing pipeline
	===================================
	bam_file: ${params.bam_file}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()


/*
 * convert a pacbio bam file to fastq
 */

process bam2fastq {
	
	conda "bioconda::pbtk=3.1.1"
	container = "quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0"
	
	input:
	tuple val(sample_id), path(bam)
	
	output:
	path "sample.fastq"
	
	script:
    """
    bam2fastq -u -o sample ${bam[0]}
    """
}

// /*
//  * run published scRPBPBR pipeline on split fastq files
//  */

process parse_polylox_barcodes {

	input:
	path "sample.fastq"
	
	output:
	path "sample.seg_assemble.tsv"
	path "sample.PB_per_BC.summary.tsv"
	path "sample.stat.tsv"
	
	script:
	"""
	touch sample.seg_assemble.tsv
	touch sample.PB_per_BC.summary.tsv
	touch sample.stat.tsv	
	"""
}

workflow {
	bam_ch = Channel.fromFilePairs(params.bam_file , checkIfExists: true)
	
    bam2fastq(bam_ch)
    parse_polylox_barcodes(bam2fastq.out)
}