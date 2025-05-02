// Visualization module for LIANA-only results

process VISUALIZE_LIANA_RESULTS {
    tag "liana_visualization"
    label 'process_medium'
    
    publishDir "${params.outdir}/liana_results", mode: 'copy'
    
    input:
    path liana_results
    
    output:
    path "*.pdf"
    path "liana_summary.csv"
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Better error reporting
    options(warn=1)
    
    # Load required libraries
    tryCatch({
        library(tidyverse)
        library(ggplot2)
        library(viridis)
        cat("Libraries loaded successfully\\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\\n")
        stop("Library loading failed")
    })
    
    # Load LIANA results
    liana_loaded <- tryCatch({
        liana_results_df <- read.csv("${liana_results}")
        cat("LIANA results loaded with", nrow(liana_results_df), "rows and", 
            ncol(liana_results_df), "columns\\n")
        cat("LIANA columns:", paste(colnames(liana_results_df), collapse=", "), "\\n")
        TRUE
    }, error = function(e) {
        cat("Error loading LIANA results:", conditionMessage(e), "\\n")
        FALSE
    })
    
    if(!liana_loaded) {
        cat("Creating placeholder output due to loading errors\\n")
        write.csv(data.frame(error = "Error loading LIANA results"), "liana_summary.csv")
        
        pdf("liana_summary_plot.pdf", width = 8, height = 6)
        plot(1, type="n", xlab="", ylab="", main="Error Loading LIANA Results")
        dev.off()
        
        quit(status = 1)
    }
    
    # Determine the main score column(s) in LIANA results
    score_columns <- c("magnitude_rank", "specificity_rank", "magnitude", "specificity_score", 
                       "lr_score", "score")
    liana_score_col <- NULL
    
    for(col in score_columns) {
        if(col %in% colnames(liana_results_df)) {
            liana_score_col <- col
            cat("Using LIANA score column:", liana_score_col, "\\n")
            break
        }
    }
    
    if(is.null(liana_score_col)) {
        # Fall back to first numeric column
        numeric_cols <- sapply(liana_results_df, is.numeric)
        if(any(numeric_cols)) {
            liana_score_col <- names(numeric_cols)[which(numeric_cols)[1]]
            cat("Falling back to first numeric column:", liana_score_col, "\\n")
        } else {
            liana_score_col <- colnames(liana_results_df)[1]
            cat("No numeric columns found, using:", liana_score_col, "\\n")
        }
    }
    
    # Check if 'ascending' should be used for ranking
    ascending <- 'rank' %in% liana_score_col
    cat("Sorting by", liana_score_col, "in", ifelse(ascending, "ascending", "descending"), "order\\n")
    
    # Extract ligand, receptor, source, and target columns if they exist
    ligand_cols <- grep("ligand", colnames(liana_results_df), value = TRUE, ignore.case = TRUE)
    receptor_cols <- grep("receptor", colnames(liana_results_df), value = TRUE, ignore.case = TRUE)
    source_cols <- grep("source", colnames(liana_results_df), value = TRUE, ignore.case = TRUE)
    target_cols <- grep("target", colnames(liana_results_df), value = TRUE, ignore.case = TRUE)
    
    if(length(ligand_cols) == 0) ligand_cols <- "ligand"
    if(length(receptor_cols) == 0) receptor_cols <- "receptor"
    if(length(source_cols) == 0) source_cols <- "source"
    if(length(target_cols) == 0) target_cols <- "target"
    
    ligand_col <- ligand_cols[1]
    receptor_col <- receptor_cols[1]
    source_col <- source_cols[1]
    target_col <- target_cols[1]
    
    # Create summary plots
    cat("Creating summary plots...\\n")
    
    # 1. Top ligand-receptor pairs
    tryCatch({
        # Ensure required columns exist
        required_cols <- c(ligand_col, receptor_col, liana_score_col)
        missing_cols <- required_cols[!required_cols %in% colnames(liana_results_df)]
        
        if(length(missing_cols) > 0) {
            cat("Warning: Missing required columns:", paste(missing_cols, collapse=", "), "\\n")
            # Create dummy data for visualization
            liana_results_df[[ligand_col]] <- if(ligand_col %in% colnames(liana_results_df)) 
                                           liana_results_df[[ligand_col]] else paste0("Ligand_", 1:nrow(liana_results_df))
            liana_results_df[[receptor_col]] <- if(receptor_col %in% colnames(liana_results_df)) 
                                             liana_results_df[[receptor_col]] else paste0("Receptor_", 1:nrow(liana_results_df))
        }
        
        # Create pair labels
        liana_results_df$pair <- paste(liana_results_df[[ligand_col]], liana_results_df[[receptor_col]], sep = " - ")
        
        # Get top pairs
        top_n <- min(20, nrow(liana_results_df))
        top_pairs <- liana_results_df %>%
            arrange(ifelse(ascending, .data[[liana_score_col]], desc(.data[[liana_score_col]]))) %>%
            head(top_n)
        
        # Create plot
        p1 <- ggplot(top_pairs, aes(x = reorder(pair, ifelse(ascending, .data[[liana_score_col]], 
                                                          -.data[[liana_score_col]])), 
                                 y = .data[[liana_score_col]])) +
            geom_col(fill = "steelblue") +
            coord_flip() +
            labs(title = paste("Top", top_n, "Ligand-Receptor Pairs"),
                 x = "Ligand-Receptor Pair",
                 y = liana_score_col) +
            theme_minimal() +
            theme(axis.text.y = element_text(size = 8))
        
        ggsave("top_lr_pairs.pdf", p1, width = 10, height = 8)
        cat("Created top ligand-receptor pairs plot\\n")
    }, error = function(e) {
        cat("Error creating top pairs plot:", conditionMessage(e), "\\n")
        pdf("top_lr_pairs.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="Error Creating Plot")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # 2. Cell type interaction heatmap
    tryCatch({
        if(all(c(source_col, target_col, liana_score_col) %in% colnames(liana_results_df))) {
            # Aggregate scores by source-target pairs
            interaction_scores <- liana_results_df %>%
                group_by(.data[[source_col]], .data[[target_col]]) %>%
                summarize(avg_score = mean(.data[[liana_score_col]], na.rm = TRUE), .groups = 'drop')
            
            # Spread to create a matrix
            interaction_matrix <- interaction_scores %>%
                spread(key = .data[[target_col]], value = avg_score, fill = NA)
            
            # Convert to matrix form
            rownames(interaction_matrix) <- interaction_matrix[[source_col]]
            interaction_matrix <- interaction_matrix %>% select(-all_of(source_col))
            matrix_data <- as.matrix(interaction_matrix)
            
            # Create heatmap
            pdf("cell_type_interaction_heatmap.pdf", width = 10, height = 8)
            image(t(matrix_data), 
                  col = viridis(100, direction = ifelse(ascending, 1, -1)),
                  xlab = "Sender Cell Type", 
                  ylab = "Receiver Cell Type",
                  main = "Cell Type Interaction Scores")
            axis(1, at = seq(0, 1, length.out = ncol(matrix_data)), 
                 labels = colnames(matrix_data), las = 2, cex.axis = 0.7)
            axis(2, at = seq(0, 1, length.out = nrow(matrix_data)), 
                 labels = rownames(matrix_data), las = 2, cex.axis = 0.7)
            dev.off()
            cat("Created cell type interaction heatmap\\n")
        } else {
            cat("Cannot create cell type heatmap: missing required columns\\n")
        }
    }, error = function(e) {
        cat("Error creating cell type heatmap:", conditionMessage(e), "\\n")
        pdf("cell_type_interaction_heatmap.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="Error Creating Heatmap")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # 3. Top ligands and receptors barplots
    tryCatch({
        # Top ligands
        top_ligands <- liana_results_df %>%
            group_by(.data[[ligand_col]]) %>%
            summarize(avg_score = mean(.data[[liana_score_col]], na.rm = TRUE), .groups = 'drop') %>%
            arrange(ifelse(ascending, avg_score, desc(avg_score))) %>%
            head(15)
        
        p2 <- ggplot(top_ligands, aes(x = reorder(.data[[ligand_col]], 
                                             ifelse(ascending, avg_score, -avg_score)), 
                                   y = avg_score)) +
            geom_col(fill = "darkgreen") +
            coord_flip() +
            labs(title = "Top 15 Ligands",
                 x = "Ligand",
                 y = paste("Average", liana_score_col)) +
            theme_minimal()
        
        ggsave("top_ligands.pdf", p2, width = 9, height = 7)
        
        # Top receptors
        top_receptors <- liana_results_df %>%
            group_by(.data[[receptor_col]]) %>%
            summarize(avg_score = mean(.data[[liana_score_col]], na.rm = TRUE), .groups = 'drop') %>%
            arrange(ifelse(ascending, avg_score, desc(avg_score))) %>%
            head(15)
        
        p3 <- ggplot(top_receptors, aes(x = reorder(.data[[receptor_col]], 
                                               ifelse(ascending, avg_score, -avg_score)), 
                                     y = avg_score)) +
            geom_col(fill = "darkred") +
            coord_flip() +
            labs(title = "Top 15 Receptors",
                 x = "Receptor",
                 y = paste("Average", liana_score_col)) +
            theme_minimal()
        
        ggsave("top_receptors.pdf", p3, width = 9, height = 7)
        cat("Created top ligands and receptors plots\\n")
    }, error = function(e) {
        cat("Error creating top ligands/receptors plots:", conditionMessage(e), "\\n")
    })
    
    # Create a summary report
    tryCatch({
        summary_data <- list(
            total_interactions = nrow(liana_results_df),
            unique_ligands = length(unique(liana_results_df[[ligand_col]])),
            unique_receptors = length(unique(liana_results_df[[receptor_col]])),
            top_pairs = top_pairs$pair[1:min(10, nrow(top_pairs))]
        )
        
        write.csv(as.data.frame(summary_data), "liana_summary.csv")
        cat("Created summary report\\n")
    }, error = function(e) {
        cat("Error creating summary report:", conditionMessage(e), "\\n")
        write.csv(data.frame(error = conditionMessage(e)), "liana_summary.csv")
    })
    
    cat("LIANA visualization complete!\\n")
    """
}