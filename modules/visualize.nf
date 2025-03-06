// Visualization module for integrated results

process VISUALIZE_RESULTS {
    tag "integrated_visualization"
    label 'process_medium'
    
    publishDir "${params.outdir}/integrated_results", mode: 'copy'
    
    input:
    path nichenet_results
    path liana_results
    
    output:
    path "liana_nichenet_integrated_results.csv"
    path "*.pdf"
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    library(Seurat)
    library(nichenetr)
    library(tidyverse)
    library(ggplot2)
    library(igraph)
    library(ggnetwork)
    library(pheatmap)
    library(viridis)
    
    # Load results
    cat("Loading results...\n")
    nichenet_output <- readRDS("${nichenet_results}")
    liana_results <- read.csv("${liana_results}")
    
    # 1. Extract ligand-target relationships
    cat("Processing ligand-target relationships...\n")
    ligand_target_matrix <- nichenet_output\$ligand_target_matrix
    
    # Get the unique ligands from the matrix
    unique_ligands <- rownames(ligand_target_matrix)
    
    # For each ligand, get its top targets
    top_targets_per_ligand <- lapply(unique_ligands, function(lig) {
      if(length(ligand_target_matrix[lig,]) == 0) return(NULL)
      
      targets <- ligand_target_matrix[lig,]
      targets_df <- data.frame(
        ligand = lig,
        target = names(targets),
        target_score = as.numeric(targets)
      )
      targets_df <- targets_df[order(targets_df\$target_score, decreasing = TRUE),]
      return(head(targets_df, 10))  # Top 10 targets per ligand
    })
    
    # Remove NULL entries and combine all targets
    all_top_targets <- do.call(rbind, top_targets_per_ligand[!sapply(top_targets_per_ligand, is.null)])
    
    # Join with LIANA scores
    liana_subset <- liana_results %>%
      select(ligand, receptor, source, target, specificity_score) %>%
      distinct()
    
    # Join the datasets
    integrated_results <- all_top_targets %>%
      left_join(liana_subset, by = "ligand")
    
    # Remove NA values
    integrated_results <- na.omit(integrated_results)
    
    # Create a comprehensive report
    write.csv(integrated_results, "liana_nichenet_integrated_results.csv")
    cat("Saved integrated results to liana_nichenet_integrated_results.csv\n")
    
    # 2. Create a correlation plot between LIANA scores and NicheNet activities
    cat("Creating correlation analysis...\n")
    
    # Extract ligand activities from NicheNet
    ligand_activities <- nichenet_output\$ligand_activities %>%
      as.data.frame() %>% 
      rownames_to_column("ligand") %>%
      rename(nichenet_activity = pearson)
    
    # Merge with LIANA scores
    liana_nichenet_integrated <- liana_results %>%
      select(ligand, specificity_score) %>%
      distinct() %>%
      inner_join(ligand_activities, by = "ligand")
    
    # Plot correlation
    correlation_plot <- ggplot(liana_nichenet_integrated, aes(x = specificity_score, y = nichenet_activity)) +
      geom_point(aes(size = nichenet_activity, color = specificity_score), alpha = 0.7) +
      geom_smooth(method = "lm", color = "red", linetype = "dashed") +
      scale_color_viridis_c() +
      labs(title = "Correlation between LIANA specificity and NicheNet activity",
           x = "LIANA specificity score",
           y = "NicheNet ligand activity",
           size = "NicheNet activity",
           color = "LIANA score") +
      theme_minimal() +
      theme(legend.position = "right")
    
    ggsave("liana_nichenet_correlation.pdf", correlation_plot, width = 8, height = 6)
    cat("Saved correlation plot to liana_nichenet_correlation.pdf\n")
    
    # 3. Create a summary heatmap
    cat("Creating summary heatmap...\n")
    
    # Get ranks for each ligand in both analyses
    liana_ranks <- liana_results %>%
      select(ligand, specificity_score) %>%
      distinct() %>%
      arrange(desc(specificity_score)) %>%
      mutate(liana_rank = row_number())
    
    nichenet_ranks <- ligand_activities %>%
      arrange(desc(nichenet_activity)) %>%
      mutate(nichenet_rank = row_number())
    
    # Get top ligands from both analyses
    top_liana_ligands <- head(arrange(liana_ranks, liana_rank)\$ligand, 15)
    top_nichenet_ligands <- head(arrange(nichenet_ranks, nichenet_rank)\$ligand, 15)
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
      lig <- summary_df\$ligand[i]
      
      # LIANA info
      liana_row <- filter(liana_ranks, ligand == lig)
      if (nrow(liana_row) > 0) {
        summary_df\$liana_score[i] <- liana_row\$specificity_score[1]
        summary_df\$liana_rank[i] <- liana_row\$liana_rank[1]
      }
      
      # NicheNet info
      nichenet_row <- filter(nichenet_ranks, ligand == lig)
      if (nrow(nichenet_row) > 0) {
        summary_df\$nichenet_score[i] <- nichenet_row\$nichenet_activity[1]
        summary_df\$nichenet_rank[i] <- nichenet_row\$nichenet_rank[1]
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
               LIANA_rank = summary_df\$liana_rank,
               NicheNet_rank = summary_df\$nichenet_rank,
               row.names = summary_df\$ligand
             ),
             angle_col = 45)
    dev.off()
    cat("Saved summary heatmap to top_ligands_summary_heatmap.pdf\n")
    
    cat("Visualization complete!\n")
    """
}
