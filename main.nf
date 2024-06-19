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
	tuple val(meta), path("${meta.sample_id}.fastq")
	
	script:
    """
    bam2fastq -u -o $meta.sample_id $bam
    """
}

/*
 * split fastq files into parts for faster processing
 */

process split_fastq {
	tag "$meta.sample_id"
	conda "bioconda::seqkit=2.8.2"
	container = "quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0"

	input:
	tuple val(meta), path(fastq)
	
	output:
	tuple val(meta), path("*.fastq")
	
	script:
    """
    seqkit split2 $fastq -s $params.fastq_splitSize -o $meta.sample_id -O .
    """

}

// /*
//  * run published scRPBPBR pipeline on split fastq files
//  */

process parse_polylox_barcodes {
	tag "$meta.sample_id"
	conda "bioconda::bowtie2=2.4.4,samtools=1.13"
	container = "quay.io/biocontainers/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:b6524911af823c7c52518f6c886b86916d062940-0"


	input:
	tuple val(meta), path(fastq)
	
	output:
	tuple val(meta), path("*.seg_assemble.tsv"), path("*.stat.tsv")
	
	script:
	"""
	scRPBPBR $fastq sample fastq
	"""
}

process merge_polylox_barcodes {
	tag "$meta.sample_id"
	publishDir params.results_dir, mode: 'copy'
	container = "rocker/tidyverse:4.4.1"
	
	input:
	tuple val(meta),  path("input/*.seg_assemble.tsv"), path("input/*.stat.tsv")
	
	output:
	tuple val(meta), path("${meta.sample_id}_merged.seg_assemble.tsv"), path("${meta.sample_id}_merged.stat.tsv")
	
	script:
	"""
	cat input/*.seg_assemble.tsv > ${meta.sample_id}_merged.seg_assemble.tsv
	merge_stats.R $meta.sample_id
	"""

}

process compute_pgen {
	tag "$meta.sample_id"
	publishDir params.results_dir, mode: 'copy'
	container = "community.wave.seqera.io/library/pip_polyloxpgen:d2ac367a6df9753f"
	
	input:
	tuple val(meta), path("${meta.sample_id}_merged.seg_assemble.tsv"), path("${meta.sample_id}_merged.stat.tsv")
	
	output:
	tuple val(meta), path("${meta.sample_id}_pgen.txt")
	
	script:
	"""
	pgen.py ${meta.sample_id}_merged.stat.tsv $meta.sample_id
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
       
    fq_ch = bam2fastq(bam_ch)
    
    split_fq_ch = split_fastq(fq_ch)
    .transpose()
    
    polylox_ch = parse_polylox_barcodes(split_fq_ch)
    | groupTuple
    | merge_polylox_barcodes
    | compute_pgen
}