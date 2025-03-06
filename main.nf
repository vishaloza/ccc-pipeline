#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { LIANA_ANALYSIS } from './modules/liana'
include { NICHENET_ANALYSIS } from './modules/nichenet'
include { VISUALIZE_RESULTS } from './modules/visualize'

// Default parameters
params.input = null
params.sender_celltype = null
params.receiver_celltype = null
params.outdir = 'results'
params.help = false

// Print help message
def helpMessage() {
    log.info"""
    =======================================
    LIANA + NicheNet Nextflow Pipeline
    =======================================
    
    Usage:
    
    nextflow run main.nf --input [h5ad file] --sender_celltype [cell type] --receiver_celltype [cell type]
    
    Required arguments:
      --input              Path to AnnData object (h5ad file)
      --sender_celltype    Sender cell type name
      --receiver_celltype  Receiver cell type name
    
    Optional arguments:
      --outdir             Output directory (default: 'results')
      --help               Show this message
    """
}

// Show help message
if (params.help || params.input == null || params.sender_celltype == null || params.receiver_celltype == null) {
    helpMessage()
    exit 0
}

// Validate input path
if (params.input) {
    input_file = file(params.input)
    if(!input_file.exists()) {
        exit 1, "Input file not found: ${params.input}"
    }
}

// Log parameters
log.info"""
=============================================
LIANA + NicheNet Nextflow Pipeline
=============================================
Input file    : ${params.input}
Sender        : ${params.sender_celltype}
Receiver      : ${params.receiver_celltype}
Output dir    : ${params.outdir}
=============================================
"""

// Main workflow
workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input)
    
    // Run LIANA analysis
    LIANA_ANALYSIS(
        input_ch,
        params.sender_celltype,
        params.receiver_celltype
    )
    
    // Run NicheNet analysis
    NICHENET_ANALYSIS(
        LIANA_ANALYSIS.out.h5ad,
        LIANA_ANALYSIS.out.liana_results,
        params.sender_celltype,
        params.receiver_celltype
    )
    
    // Visualize integrated results
    VISUALIZE_RESULTS(
        NICHENET_ANALYSIS.out.nichenet_results,
        LIANA_ANALYSIS.out.liana_results
    )
}

workflow.onComplete {
    log.info"""
    ===============================================
    LIANA + NicheNet Pipeline Complete!
    Results directory: ${params.outdir}
    ===============================================
    """
}
