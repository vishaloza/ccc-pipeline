// Visualization module for MultiNicheNet-only results

process VISUALIZE_MULTINICHENET_RESULTS {
    tag "multinichenet_visualization"
    label 'process_medium'
    
    publishDir "${params.outdir}/multinichenet_results", mode: 'copy'
    
    input:
    path multinichenet_results
    
    output:
    path "*.pdf"
    path "multinichenet_summary.csv"
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Better error reporting
    options(warn=1)
    
    # Load required libraries
    tryCatch({
        library(tidyverse)
        library(ggplot2)
        library(pheatmap)
        library(viridis)
        
        # Install and load MultiNicheNet
        if (!requireNamespace("devtools", quietly = TRUE)) {
            install.packages("devtools")
        }
        if (!requireNamespace("multinichenetr", quietly = TRUE)) {
            devtools::install_github("saeyslab/multinichenetr")
        }
        library(multinichenetr)
        
        cat("Libraries loaded successfully\\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\\n")
        stop("Library loading failed")
    })
    
    # Load MultiNicheNet results
    multinichenet_loaded <- tryCatch({
        multinichenet_output <- readRDS("${multinichenet_results}")
        cat("MultiNicheNet results loaded successfully\\n")
        
        if(is.list(multinichenet_output) && "error" %in% names(multinichenet_output)) {
            cat("MultiNicheNet results contain an error:", multinichenet_output\$error, "\\n")
            FALSE
        } else {
            TRUE
        }
    }, error = function(e) {
        cat("Error loading MultiNicheNet results:", conditionMessage(e), "\\n")
        FALSE
    })
    
    if(!multinichenet_loaded) {
        cat("Creating placeholder output due to loading errors\\n")
        write.csv(data.frame(error = "Error loading MultiNicheNet results"), "multinichenet_summary.csv")
        
        pdf("multinichenet_summary_plot.pdf", width = 8, height = 6)
        plot(1, type="n", xlab="", ylab="", main="Error Loading MultiNicheNet Results")
        dev.off()
        
        quit(status = 1)
    }
    
    # Create additional summary visualizations
    cat("Creating additional summary visualizations...\\n")
    
    # Extract components
    prioritization_tables <- multinichenet_output$prioritization_tables
    abundances <- multinichenet_output$abundances
    filtered_abundances <- multinichenet_output$filtered_abundances
    de_results <- multinichenet_output$de_results
    de_summary <- multinichenet_output$de_summary
    top_ligands <- multinichenet_output$top_ligands
    
    # Create a circos plot of top sender-receiver pairs
    tryCatch({
        cat("Creating circos plot...\\n")
        pdf("sender_receiver_circos.pdf", width = 10, height = 10)
        multinichenetr::plot_sender_receiver_circos(
            prioritization_tables,
            top_n = 30
        )
        dev.off()
        cat("Created circos plot\\n")
    }, error = function(e) {
        cat("Error creating circos plot:", conditionMessage(e), "\\n")
        pdf("sender_receiver_circos.pdf", width = 10, height = 10)
        plot(1, type="n", xlab="", ylab="", main="Error Creating Circos Plot")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # Create a dotplot of differentially expressed ligands
    tryCatch({
        cat("Creating DE ligands dotplot...\\n")
        pdf("de_ligands_dotplot.pdf", width = 12, height = 10)
        multinichenetr::plot_de_ligands_dotplot(
            prioritization_tables,
            top_n = 20
        )
        dev.off()
        cat("Created DE ligands dotplot\\n")
    }, error = function(e) {
        cat("Error creating DE ligands dotplot:", conditionMessage(e), "\\n")
        pdf("de_ligands_dotplot.pdf", width = 12, height = 10)
        plot(1, type="n", xlab="", ylab="", main="Error Creating DE Ligands Dotplot")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # Create a summary report
    tryCatch({
        sr_pairs <- names(prioritization_tables)
        top_sr_pairs <- sr_pairs[1:min(10, length(sr_pairs))]
        
        summary_data <- list(
            total_sender_receiver_pairs = length(sr_pairs),
            top_sender_receiver_pairs = top_sr_pairs,
            total_ligands = length(unique(top_ligands$ligand)),
            top_ligands = unique(top_ligands$ligand)[1:min(20, length(unique(top_ligands$ligand)))]
        )
        
        write.csv(as.data.frame(summary_data), "multinichenet_summary.csv")
        cat("Created summary report\\n")
    }, error = function(e) {
        cat("Error creating summary report:", conditionMessage(e), "\\n")
        write.csv(data.frame(error = conditionMessage(e)), "multinichenet_summary.csv")
    })
    
    cat("MultiNicheNet visualization complete!\\n")
    """
}