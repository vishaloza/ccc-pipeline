#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { LIANA_ANALYSIS } from './modules/liana'
include { NICHENET_ANALYSIS } from './modules/nichenet'
include { NICHENET_ANALYSIS_STANDALONE } from './modules/nichenet_standalone'
include { MULTINICHENET_ANALYSIS_STANDALONE } from './modules/multinichenet_standalone'
include { VISUALIZE_RESULTS } from './modules/visualize'
include { VISUALIZE_LIANA_RESULTS } from './modules/visualize_liana'
include { VISUALIZE_NICHENET_RESULTS } from './modules/visualize_nichenet'
include { VISUALIZE_MULTINICHENET_RESULTS } from './modules/visualize_multinichenet'

// Default parameters
params.input = null
params.input_type = 'h5ad'  // Default to LIANA-only
params.outdir = 'results'
params.help = false

// Print help message
def helpMessage() {
    log.info"""
    =======================================
    LIANA + NicheNet Nextflow Pipeline
    =======================================
    
    Usage:
    
    nextflow run main.nf --input [file] --input_type [type]
    
    Required arguments:
      --input              Path to input file (.h5ad for LIANA, .rds for NicheNet/MultiNicheNet)
    
    Optional arguments:
      --input_type         Type of input file: 'h5ad' (LIANA only), 'rds' (NicheNet only), 
                           'rds_multi' (MultiNicheNet only), or 'both' (integrated analysis) (default: 'h5ad')
      --outdir             Output directory (default: 'results')
      --help               Show this message
    """
}

// Show help message
if (params.help || params.input == null) {
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
LIANA + NicheNet/MultiNicheNet Nextflow Pipeline
=============================================
Input file    : ${params.input}
Input type    : ${params.input_type}
Output dir    : ${params.outdir}
=============================================
"""

// Main workflow
workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input)
    
    // Choose workflow based on input type
    if(params.input_type == 'h5ad') {
        log.info "Running LIANA-only workflow for AnnData input"
        
        // Run LIANA analysis
        LIANA_ANALYSIS(input_ch)
        
        // Visualize LIANA results only
        VISUALIZE_LIANA_RESULTS(LIANA_ANALYSIS.out.liana_results)
    } 
    else if(params.input_type == 'rds') {
        log.info "Running NicheNet-only workflow for Seurat object input"
        
        // Run NicheNet analysis directly with Seurat object
        NICHENET_ANALYSIS_STANDALONE(input_ch)
        
        // Visualize NicheNet results only
        VISUALIZE_NICHENET_RESULTS(NICHENET_ANALYSIS_STANDALONE.out.nichenet_results)
    }
    else if(params.input_type == 'rds_multi') {
        log.info "Running MultiNicheNet workflow for Seurat object input"
        
        // Run MultiNicheNet analysis directly with Seurat object
        MULTINICHENET_ANALYSIS_STANDALONE(input_ch)
        
        // Visualize MultiNicheNet results
        VISUALIZE_MULTINICHENET_RESULTS(MULTINICHENET_ANALYSIS_STANDALONE.out.multinichenet_results)
    }
    else {
        log.info "Running integrated LIANA + NicheNet workflow"
        
        // Run LIANA analysis
        LIANA_ANALYSIS(input_ch)
        
        // Run NicheNet analysis with LIANA output
        NICHENET_ANALYSIS(
            LIANA_ANALYSIS.out.h5ad,
            LIANA_ANALYSIS.out.raw_expr,
            LIANA_ANALYSIS.out.norm_expr,
            LIANA_ANALYSIS.out.cell_meta,
            LIANA_ANALYSIS.out.gene_meta,
            LIANA_ANALYSIS.out.liana_results,
            LIANA_ANALYSIS.out.top_lr_pairs
        )
        
        // Visualize integrated results
        VISUALIZE_RESULTS(
            NICHENET_ANALYSIS.out.nichenet_results,
            LIANA_ANALYSIS.out.liana_results
        )
    }
}

workflow.onComplete {
    log.info"""
    ===============================================
    LIANA + NicheNet/MultiNicheNet Pipeline Complete!
    Results directory: ${params.outdir}
    ===============================================
    """
}