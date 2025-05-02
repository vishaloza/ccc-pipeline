#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import only the modules needed for MultiNicheNet
include { MULTINICHENET_ANALYSIS_STANDALONE } from './modules/multinichenet_standalone'
include { VISUALIZE_MULTINICHENET_RESULTS } from './modules/visualize_multinichenet'

// Parameters
params.input = null
params.outdir = 'results'
params.celltype_column = 'cell_type'
params.help = false

// Main workflow
workflow {
    input_ch = Channel.fromPath(params.input)
    
    // Run MultiNicheNet analysis
    MULTINICHENET_ANALYSIS_STANDALONE(input_ch)
    
    // Visualize results
    VISUALIZE_MULTINICHENET_RESULTS(MULTINICHENET_ANALYSIS_STANDALONE.out.multinichenet_results)
}