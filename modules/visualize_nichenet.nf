// Visualization module for NicheNet-only results

process VISUALIZE_NICHENET_RESULTS {
    tag "nichenet_visualization"
    label 'process_medium'
    
    publishDir "${params.outdir}/nichenet_results", mode: 'copy'
    
    input:
    path nichenet_results
    
    output:
    path "*.pdf"
    path "nichenet_summary.csv"
    
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
        cat("Libraries loaded successfully\\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\\n")
        stop("Library loading failed")
    })
    
    # Load NicheNet results
    nichenet_loaded <- tryCatch({
        nichenet_output <- readRDS("${nichenet_results}")
        cat("NicheNet results loaded successfully\\n")
        
        if(is.list(nichenet_output) && "error" %in% names(nichenet_output)) {
            cat("NicheNet results contain an error:", nichenet_output$error, "\\n")
            FALSE
        } else {
            # Check if key components exist
            required_components <- c("ligand_activities", "ligand_target_matrix")
            missing_components <- required_components[!required_components %in% names(nichenet_output)]
            
            if(length(missing_components) > 0) {
                cat("NicheNet results missing key components:", 
                    paste(missing_components, collapse=", "), "\\n")
                cat("Available components:", paste(names(nichenet_output), collapse=", "), "\\n")
                FALSE
            } else {
                TRUE
            }
        }
    }, error = function(e) {
        cat("Error loading NicheNet results:", conditionMessage(e), "\\n")
        FALSE
    })
    
    if(!nichenet_loaded) {
        cat("Creating placeholder output due to loading errors\\n")
        write.csv(data.frame(error = "Error loading NicheNet results"), "nichenet_summary.csv")
        
        pdf("nichenet_summary_plot.pdf", width = 8, height = 6)
        plot(1, type="n", xlab="", ylab="", main="Error Loading NicheNet Results")
        dev.off()
        
        quit(status = 1)
    }
    
    # Create summary plots
    cat("Creating summary plots...\\n")
    
    # 1. Ligand activity summary
    tryCatch({
        if("ligand_activities" %in% names(nichenet_output)) {
            # For standard nichenet_output format
            ligand_activities <- nichenet_output$ligand_activities
            
            # Sort by activity score
            top_ligands <- ligand_activities %>%
                arrange(desc(aupr_corrected)) %>%
                head(20)
            
            # Create barplot
            p1 <- ggplot(top_ligands, aes(x = reorder(test_ligand, aupr_corrected), y = aupr_corrected)) +
                geom_col(fill = "steelblue") +
                coord_flip() +
                labs(title = "Top 20 Ligands by Activity Score",
                     x = "Ligand",
                     y = "Activity Score (aupr_corrected)") +
                theme_minimal()
            
            ggsave("top_ligand_activities.pdf", p1, width = 10, height = 8)
            cat("Created top ligand activities plot\\n")
        } else {
            cat("Cannot create ligand activity plot: missing required data\\n")
        }
    }, error = function(e) {
        cat("Error creating ligand activity plot:", conditionMessage(e), "\\n")
        pdf("top_ligand_activities.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="Error Creating Plot")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # 2. Target genes per ligand
    tryCatch({
        if("ligand_target_matrix" %in% names(nichenet_output)) {
            # Get ligand-target matrix
            ligand_target_matrix <- nichenet_output$ligand_target_matrix
            
            # Count number of targets per ligand (using a threshold)
            threshold <- 0.005  # Adjust as needed
            targets_per_ligand <- rowSums(ligand_target_matrix > threshold)
            
            # Create data frame for plotting
            targets_df <- data.frame(
                ligand = names(targets_per_ligand),
                target_count = as.numeric(targets_per_ligand)
            ) %>%
                arrange(desc(target_count)) %>%
                head(20)
            
            # Create barplot
            p2 <- ggplot(targets_df, aes(x = reorder(ligand, target_count), y = target_count)) +
                geom_col(fill = "darkgreen") +
                coord_flip() +
                labs(title = "Top 20 Ligands by Number of Target Genes",
                     x = "Ligand",
                     y = paste("Number of Target Genes (threshold", threshold, ")")) +
                theme_minimal()
            
            ggsave("ligands_by_target_count.pdf", p2, width = 10, height = 8)
            cat("Created ligands by target count plot\\n")
        } else {
            cat("Cannot create target genes plot: missing required data\\n")
        }
    }, error = function(e) {
        cat("Error creating target genes plot:", conditionMessage(e), "\\n")
        pdf("ligands_by_target_count.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="Error Creating Plot")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # 3. Ligand-target heatmap for top ligands
    tryCatch({
        if(all(c("ligand_activities", "ligand_target_matrix") %in% names(nichenet_output))) {
            # Get top ligands
            top_ligands <- nichenet_output$ligand_activities %>%
                arrange(desc(aupr_corrected)) %>%
                pull(test_ligand) %>%
                head(10)
            
            # Get target genes for these ligands
            ligand_target_matrix <- nichenet_output$ligand_target_matrix
            
            # Subset the matrix
            top_ligand_target_matrix <- ligand_target_matrix[top_ligands, ]
            
            # Find top targets (sum scores across all top ligands)
            target_sums <- colSums(top_ligand_target_matrix)
            top_targets <- names(sort(target_sums, decreasing = TRUE)[1:min(50, length(target_sums))])
            
            # Create matrix for heatmap
            heatmap_matrix <- top_ligand_target_matrix[, top_targets]
            
            # Create heatmap
            pdf("top_ligand_target_heatmap.pdf", width = 12, height = 8)
            pheatmap(heatmap_matrix,
                    color = viridis(100),
                    main = "Top Ligands and Their Target Genes",
                    fontsize_row = 8,
                    fontsize_col = 6,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_colnames = TRUE)
            dev.off()
            cat("Created top ligand-target heatmap\\n")
        } else {
            cat("Cannot create ligand-target heatmap: missing required data\\n")
        }
    }, error = function(e) {
        cat("Error creating ligand-target heatmap:", conditionMessage(e), "\\n")
        pdf("top_ligand_target_heatmap.pdf", width = 12, height = 8)
        plot(1, type="n", xlab="", ylab="", main="Error Creating Heatmap")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # Create a summary report
    tryCatch({
        summary_data <- list(
            total_ligands = if("ligand_activities" %in% names(nichenet_output)) 
                            nrow(nichenet_output$ligand_activities) else 0,
            total_targets = if("ligand_target_matrix" %in% names(nichenet_output)) 
                            ncol(nichenet_output$ligand_target_matrix) else 0
        )
        
        if("ligand_activities" %in% names(nichenet_output)) {
            top_ligands <- nichenet_output$ligand_activities %>%
                arrange(desc(aupr_corrected)) %>%
                head(10)
            
            summary_data$top_ligands <- top_ligands$test_ligand
            summary_data$top_ligand_scores <- top_ligands$aupr_corrected
        }
        
        write.csv(as.data.frame(summary_data), "nichenet_summary.csv")
        cat("Created summary report\\n")
    }, error = function(e) {
        cat("Error creating summary report:", conditionMessage(e), "\\n")
        write.csv(data.frame(error = conditionMessage(e)), "nichenet_summary.csv")
    })
    
    cat("NicheNet visualization complete!\\n")
    """
}