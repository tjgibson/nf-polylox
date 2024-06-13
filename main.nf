#! /usr/bin/env nextflow

log.info """
	PolyLox barcode processing pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()


/*
 * convert a pacbio bam file to fastq
 */

process bam2fastq {
	tag "$meta.sample_id"
	conda "bioconda::pbtk=3.1.1"
	container = "quay.io/biocontainers/pbtk:3.1.1--h9ee0642_0"
	
	input:
	tuple val(meta), path(bam), path(bam_index)
	
	output:
	tuple val(meta), path("*.fastq")
	
	script:
    """
    bam2fastq -u -o $meta.sample_id $bam
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
	bam_ch = Channel.fromPath(params.samplesheet)
	| splitCsv( header:true )
    | map { row ->
        meta = row.subMap('sample_id')
        [
        	meta, 
        	file(row.bam_file, checkIfExists: true),
            file(row.pbi_index, checkIfExists: true)]
    }
       
    bam2fastq(bam_ch)
//     parse_polylox_barcodes(bam2fastq.out)
}