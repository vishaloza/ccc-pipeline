// Visualization module for integrated results

process VISUALIZE_RESULTS {
    tag "integrated_visualization"
    label 'process_medium'
    
    publishDir "${params.outdir}/integrated_results", mode: 'copy'
    
    input:
    path nichenet_results
    path liana_results
    
    output:
    path "liana_nichenet_integrated_results.csv" optional true
    path "*.pdf" optional true
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Better error reporting
    options(warn=1)
    
    cat("Working directory:", getwd(), "\n")
    cat("NicheNet results file:", "${nichenet_results}", "\n")
    cat("LIANA results file:", "${liana_results}", "\n")
    
    # Load required libraries with error handling
    tryCatch({
        library(tidyverse)
        library(ggplot2)
        library(pheatmap)
        library(viridis)
        cat("Base libraries loaded successfully\n")
        
        # Optional libraries with fallbacks
        library_loaded <- function(lib_name) {
            if(!require(lib_name, character.only = TRUE, quietly = TRUE)) {
                cat("Optional library", lib_name, "not available\n")
                return(FALSE)
            }
            cat("Loaded optional library", lib_name, "\n")
            return(TRUE)
        }
        
        igraph_available <- library_loaded("igraph")
        ggnetwork_available <- library_loaded("ggnetwork")
        
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\n")
        cat("Available libraries:", paste(.packages(all.available = TRUE), collapse=", "), "\n")
        stop("Library loading failed")
    })
    
    # Function to safely create plots
    safe_plot <- function(plot_name, plot_function) {
        tryCatch({
            cat("Creating", plot_name, "...\n")
            plot_function()
            cat(plot_name, "created successfully\n")
            return(TRUE)
        }, error = function(e) {
            cat("Error creating", plot_name, ":", conditionMessage(e), "\n")
            # Create an error notification plot
            pdf(paste0(plot_name, ".pdf"), width = 8, height = 6)
            plot(1, type="n", xlab="", ylab="", main=paste("Error creating", plot_name))
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
            return(FALSE)
        })
    }
    
    # Load results with error handling
    cat("Loading results...\n")
    
    # Load LIANA results
    liana_loaded <- tryCatch({
        liana_results_df <- read.csv("${liana_results}")
        cat("LIANA results loaded with", nrow(liana_results_df), "rows and", 
            ncol(liana_results_df), "columns\n")
        cat("LIANA columns:", paste(colnames(liana_results_df), collapse=", "), "\n")
        TRUE
    }, error = function(e) {
        cat("Error loading LIANA results:", conditionMessage(e), "\n")
        FALSE
    })
    
    # Load NicheNet results
    nichenet_loaded <- tryCatch({
        nichenet_output <- readRDS("${nichenet_results}")
        cat("NicheNet results loaded successfully\n")
        
        # Check if it's an error object
        if(is.list(nichenet_output) && "error" %in% names(nichenet_output)) {
            cat("NicheNet results contain an error:", nichenet_output$error, "\n")
            FALSE
        } else {
            # Check if key components exist
            required_components <- c("ligand_activities", "ligand_target_matrix")
            missing_components <- required_components[!required_components %in% names(nichenet_output)]
            
            if(length(missing_components) > 0) {
                cat("NicheNet results missing key components:", 
                    paste(missing_components, collapse=", "), "\n")
                cat("Available components:", paste(names(nichenet_output), collapse=", "), "\n")
                FALSE
            } else {
                TRUE
            }
        }
    }, error = function(e) {
        cat("Error loading NicheNet results:", conditionMessage(e), "\n")
        FALSE
    })
    
    # Create placeholder files if loading failed
    if(!liana_loaded || !nichenet_loaded) {
        cat("Creating placeholder output due to loading errors\n")
        
        # Create empty integrated results file
        write.csv(data.frame(
            error = "Error loading input data",
            detail = paste0("LIANA loaded: ", liana_loaded, ", NicheNet loaded: ", nichenet_loaded)
        ), "liana_nichenet_integrated_results.csv")
        
        # Create placeholder plots
        pdf("liana_nichenet_correlation.pdf", width = 8, height = 6)
        plot(1, type="n", xlab="", ylab="", main="Error Loading Input Data")
        text(1, 1, paste0("LIANA loaded: ", liana_loaded, ", NicheNet loaded: ", nichenet_loaded), cex=0.8)
        dev.off()
        
        pdf("top_ligands_summary_heatmap.pdf", width = 10, height = 12)
        plot(1, type="n", xlab="", ylab="", main="Error Loading Input Data")
        text(1, 1, paste0("LIANA loaded: ", liana_loaded, ", NicheNet loaded: ", nichenet_loaded), cex=0.8)
        dev.off()
        
        quit(status = 1)
    }
    
    # Continue with analysis if data loaded successfully
    cat("Processing ligand-target relationships...\n")
    tryCatch({
        # 1. Extract ligand-target relationships
        ligand_target_matrix <- nichenet_output$ligand_target_matrix
        
        # Get the unique ligands from the matrix
        unique_ligands <- rownames(ligand_target_matrix)
        cat("Found", length(unique_ligands), "unique ligands in NicheNet results\n")
        
        # For each ligand, get its top targets
        top_targets_per_ligand <- lapply(unique_ligands, function(lig) {
            if(length(ligand_target_matrix[lig,]) == 0) return(NULL)
            
            targets <- ligand_target_matrix[lig,]
            targets_df <- data.frame(
                ligand = lig,
                target = names(targets),
                target_score = as.numeric(targets)
            )
            targets_df <- targets_df[order(targets_df$target_score, decreasing = TRUE),]
            return(head(targets_df, 10))  # Top 10 targets per ligand
        })
        
        # Remove NULL entries and combine all targets
        all_top_targets <- do.call(rbind, top_targets_per_ligand[!sapply(top_targets_per_ligand, is.null)])
        cat("Processed", nrow(all_top_targets), "ligand-target relationships\n")
        
        # Determine the right ranking score column in LIANA results
        score_columns <- c("magnitude", "specificity_score", "lr_score", "score")
        liana_score_col <- NULL
        for(col in score_columns) {
            if(col %in% colnames(liana_results_df)) {
                liana_score_col <- col
                cat("Using LIANA score column:", liana_score_col, "\n")
                break
            }
        }
        
        if(is.null(liana_score_col)) {
            cat("No recognized score column found in LIANA results. Available columns:", 
                paste(colnames(liana_results_df), collapse=", "), "\n")
            # Fall back to first numeric column
            numeric_cols <- sapply(liana_results_df, is.numeric)
            if(any(numeric_cols)) {
                liana_score_col <- names(numeric_cols)[which(numeric_cols)[1]]
                cat("Falling back to first numeric column:", liana_score_col, "\n")
            } else {
                stop("No numeric score column found in LIANA results")
            }
        }
        
        # Join with LIANA scores
        liana_subset <- liana_results_df %>%
            select(ligand, receptor, source, target, all_of(liana_score_col)) %>%
            distinct()
        
        # Join the datasets
        integrated_results <- all_top_targets %>%
            left_join(liana_subset, by = "ligand")
        
        # Remove NA values
        integrated_results <- na.omit(integrated_results)
        cat("Created integrated results with", nrow(integrated_results), "rows\n")
        
        # Create a comprehensive report
        write.csv(integrated_results, "liana_nichenet_integrated_results.csv", row.names = FALSE)
        cat("Saved integrated results to liana_nichenet_integrated_results.csv\n")
        
        # 2. Create a correlation plot between LIANA scores and NicheNet activities
        safe_plot("liana_nichenet_correlation", function() {
            # Extract ligand activities from NicheNet
            ligand_activities <- nichenet_output$ligand_activities %>%
                as.data.frame() %>% 
                rownames_to_column("ligand") %>%
                rename(nichenet_activity = pearson)
            
            # Merge with LIANA scores
            liana_nichenet_integrated <- liana_results_df %>%
                select(ligand, all_of(liana_score_col)) %>%
                distinct() %>%
                inner_join(ligand_activities, by = "ligand")
            
            # Calculate correlation
            corr <- cor(liana_nichenet_integrated[[liana_score_col]], 
                        liana_nichenet_integrated$nichenet_activity, 
                        method = "spearman",
                        use = "pairwise.complete.obs")
            
            # Plot correlation
            correlation_plot <- ggplot(liana_nichenet_integrated, 
                                   aes(x = .data[[liana_score_col]], y = nichenet_activity)) +
                geom_point(aes(size = nichenet_activity, color = .data[[liana_score_col]]), alpha = 0.7) +
                geom_smooth(method = "lm", color = "red", linetype = "dashed") +
                scale_color_viridis_c() +
                labs(title = paste0("Correlation between LIANA and NicheNet (Spearman rho: ", 
                                   round(corr, 3), ")"),
                     x = paste("LIANA", liana_score_col),
                     y = "NicheNet ligand activity",
                     size = "NicheNet activity",
                     color = paste("LIANA", liana_score_col)) +
                theme_minimal() +
                theme(legend.position = "right")
            
            ggsave("liana_nichenet_correlation.pdf", correlation_plot, width = 8, height = 6)
        })
        
        # 3. Create a summary heatmap
        safe_plot("top_ligands_summary_heatmap", function() {
            # Extract ligand activities from NicheNet
            ligand_activities <- nichenet_output$ligand_activities %>%
                as.data.frame() %>% 
                rownames_to_column("ligand") %>%
                rename(nichenet_activity = pearson)
            
            # Get ranks for each ligand in both analyses
            liana_ranks <- liana_results_df %>%
                select(ligand, all_of(liana_score_col)) %>%
                distinct() %>%
                arrange(desc(.data[[liana_score_col]])) %>%
                mutate(liana_rank = row_number())
            
            nichenet_ranks <- ligand_activities %>%
                arrange(desc(nichenet_activity)) %>%
                mutate(nichenet_rank = row_number())
            
            # Get top ligands from both analyses
            top_liana_ligands <- head(arrange(liana_ranks, liana_rank)$ligand, 15)
            top_nichenet_ligands <- head(arrange(nichenet_ranks, nichenet_rank)$ligand, 15)
            top_ligands <- unique(c(top_liana_ligands, top_nichenet_ligands))
            
            # Create a summary matrix
            summary_df <- data.frame(
                ligand = top_ligands,
                liana_score = NA,
                liana_rank = NA,
                nichenet_score = NA,
                nichenet_rank = NA
            )
            
            # Fill in values
            for (i in 1:nrow(summary_df)) {
                lig <- summary_df$ligand[i]
                
                # LIANA info
                liana_row <- filter(liana_ranks, ligand == lig)
                if (nrow(liana_row) > 0) {
                    summary_df$liana_score[i] <- liana_row[[liana_score_col]][1]
                    summary_df$liana_rank[i] <- liana_row$liana_rank[1]
                }
                
                # NicheNet info
                nichenet_row <- filter(nichenet_ranks, ligand == lig)
                if (nrow(nichenet_row) > 0) {
                    summary_df$nichenet_score[i] <- nichenet_row$nichenet_activity[1]
                    summary_df$nichenet_rank[i] <- nichenet_row$nichenet_rank[1]
                }
            }
            
            # Create a matrix for the heatmap
            summary_matrix <- summary_df %>%
                select(ligand, liana_score, nichenet_score) %>%
                column_to_rownames("ligand") %>%
                as.matrix()
            
            # Scale the data for better visualization
            summary_matrix_scaled <- scale(summary_matrix)
            
            # Create the heatmap
            pdf("top_ligands_summary_heatmap.pdf", width = 10, height = 12)
            pheatmap(summary_matrix_scaled,
                    main = "Top Ligands from LIANA and NicheNet",
                    color = viridis(100),
                    cluster_rows = TRUE,
                    cluster_cols = FALSE,
                    annotation_row = data.frame(
                        LIANA_rank = summary_df$liana_rank,
                        NicheNet_rank = summary_df$nichenet_rank,
                        row.names = summary_df$ligand
                    ),
                    angle_col = 45)
            dev.off()
        })
        
        # Create a ligand-target network plot if igraph and ggnetwork are available
        if(igraph_available && ggnetwork_available) {
            safe_plot("ligand_target_network", function() {
                # Get top interactions
                top_interactions <- all_top_targets %>%
                    arrange(desc(target_score)) %>%
                    head(100)
                
                # Create network
                network_df <- top_interactions %>%
                    select(ligand, target, weight = target_score)
                
                # Create graph
                g <- igraph::graph_from_data_frame(network_df, directed = TRUE)
                
                # Set node type
                V(g)$type <- ifelse(V(g)$name %in% unique(network_df$ligand), "ligand", "target")
                
                # Create plot with ggnetwork
                network_plot <- ggplot(ggnetwork::ggnetwork(g, layout = "fruchtermanreingold"), 
                                    aes(x = x, y = y, xend = xend, yend = yend)) +
                    geom_edges(aes(width = weight), color = "gray70", alpha = 0.5) +
                    geom_nodes(aes(color = type, size = type)) +
                    geom_nodelabel_repel(aes(label = name), size = 3) +
                    scale_color_manual(values = c("ligand" = "darkred", "target" = "darkblue")) +
                    scale_size_manual(values = c("ligand" = 4, "target" = 2)) +
                    scale_edge_width(range = c(0.2, 2)) +
                    theme_void() +
                    labs(title = "Ligand-Target Network") +
                    theme(legend.position = "bottom")
                
                ggsave("ligand_target_network.pdf", network_plot, width = 12, height = 10)
            })
        }
        
        cat("Visualization complete!\n")
    }, error = function(e) {
        cat("Error in visualization process:", conditionMessage(e), "\n")
        # Create minimal output to satisfy pipeline requirements
        write.csv(data.frame(error = conditionMessage(e)), "liana_nichenet_integrated_results.csv")
    })
    """
}
