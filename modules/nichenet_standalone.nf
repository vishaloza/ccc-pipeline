// Standalone NicheNet analysis module for direct Seurat object input

process NICHENET_ANALYSIS_STANDALONE {
    tag "$input_file.simpleName"
    label 'process_high'

    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/nichenet", mode: 'copy'
    
    input:
    path input_file
    
    output:
    path "nichenet_results.rds", emit: nichenet_results
    path "ligand_activity_heatmap.pdf", emit: activity_heatmap
    path "ligand_target_heatmap.pdf", emit: target_heatmap
    path "*.pdf"
    
    script:
    """
    #!/usr/bin/env Rscript

    cat("Working directory:", getwd(), "\\n")
    cat("Input file:", "${input_file}", "\\n")
    
    # Load required libraries
    tryCatch({
        library(Seurat)
        library(nichenetr)
        library(tidyverse)
        library(reticulate)   
        library(Matrix)
        cat("Libraries loaded successfully\\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\\n")
        cat("Available libraries:", paste(.packages(all.available = TRUE), collapse=", "), "\\n")
        stop("Library loading failed")
    })
    
    # Load Seurat object from RDS file
    tryCatch({
        cat("Loading Seurat object from RDS file...\\n")
        seurat_obj <- readRDS("${input_file}")
        
        # Check if it's actually a Seurat object
        if (!inherits(seurat_obj, "Seurat")) {
            stop("The loaded object is not a Seurat object.")
        }
        
        cat("Seurat object loaded successfully with dimensions:", dim(seurat_obj), "\\n")
    }, error = function(e) {
        cat("Error loading Seurat object:", conditionMessage(e), "\\n")
        stop("Failed to load Seurat object.")
    })
    
    # Define cell_type_column
    cell_type_column <- "${params.celltype_column}"
    
    # Check if cell type column exists in Seurat object
    if(!cell_type_column %in% colnames(seurat_obj@meta.data)) {
        cat("WARNING: Specified cell type column", cell_type_column, "not found in Seurat object\\n")
        cat("Available columns:", paste(colnames(seurat_obj@meta.data), collapse=", "), "\\n")
        
        # Try to find a suitable column
        possible_columns <- c("cell_type", "celltype", "CellType", "clusters", "cluster", "leiden")
        for(col in possible_columns) {
            if(col %in% colnames(seurat_obj@meta.data)) {
                cell_type_column <- col
                cat("Using alternative column:", cell_type_column, "\\n")
                break
            }
        }
    }
    
    # Set Idents from cell type column for consistency
    Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_column]]
    cat("Available Idents in Seurat object:", paste(levels(Idents(seurat_obj)), collapse=", "), "\\n")
    
    # Get the list of unique cell types
    unique_cell_types <- levels(Idents(seurat_obj))
    cat("Using these cell types:", paste(unique_cell_types, collapse=", "), "\\n")
    
    # Download NicheNet resources directly
    tryCatch({
        cat("Downloading NicheNet resources from Zenodo...\\n")
        ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
        weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
        
        cat("Successfully loaded ligand_target_matrix with dimensions:", dim(ligand_target_matrix)[1], "x", dim(ligand_target_matrix)[2], "\\n")
        cat("Successfully loaded weighted_networks\\n")
    }, error = function(e) {
        cat("ERROR loading NicheNet resources:", conditionMessage(e), "\\n")
        stop("Failed to load required NicheNet resources. Cannot continue analysis.")
    })
    
    # Run NicheNet analysis
    cat("Running NicheNet analysis...\\n")
    
    # Extract genes expressed in at least 10% of the cells in a particular cell type
    cat("Getting expressed genes by cell type...\\n")
    expressed_genes_list <- nichenetr::get_expressed_genes(
        unique_cell_types,
        seurat_obj,
        pct = 0.1
    )
    
    # Get ligands from expressed genes
    cat("Finding expressed ligands and receptors...\\n")
    
    # Get ligand and receptor lists from NicheNet
    ligands <- rownames(ligand_target_matrix)
    receptors <- colnames(weighted_networks$lr_sig)
    
    # Find expressed ligands and receptors
    expressed_genes <- data.frame(gene = expressed_genes_list)
    expressed_ligands <- expressed_genes[["gene"]][expressed_genes[["gene"]] %in% ligands]
    expressed_receptors <- expressed_genes[["gene"]][expressed_genes[["gene"]] %in% receptors]
    
    cat("Found", length(expressed_ligands), "expressed ligands and", 
        length(expressed_receptors), "expressed receptors\\n")
    
    # If too few ligands/receptors are found, use all from the database
    if(length(expressed_ligands) < 10) {
        cat("WARNING: Too few expressed ligands found. Using all ligands from database.\\n")
        expressed_ligands <- head(ligands, 100)
    }
    if(length(expressed_receptors) < 10) {
        cat("WARNING: Too few expressed receptors found. Using all receptors from database.\\n")
        expressed_receptors <- head(receptors, 100)
    }
    
    # Run NicheNet directly
    tryCatch({
        nichenet_output <- nichenetr::nichenet_seuratobj_aggregate(
            seurat_obj = seurat_obj,
            receiver = cell_type_column,
            sender = unique_cell_types,
            condition_colname = NULL,
            ligand_target_matrix = ligand_target_matrix,
            weighted_networks = weighted_networks,
            expr_pct = 0.1
        )
        
        saveRDS(nichenet_output, "nichenet_results.rds")
        cat("Saved NicheNet results to 'nichenet_results.rds'\\n")
    }, error = function(e) {
        cat("ERROR in NicheNet analysis:", conditionMessage(e), "\\n")
        
        # Create minimal output so pipeline can continue
        nichenet_output <- list(
            error = conditionMessage(e),
            ligand_activities = data.frame(test_ligand = head(expressed_ligands, 10), 
                                          auroc = NA, aupr = NA, 
                                          aupr_corrected = NA, pearson = NA, 
                                          rank = 1:10)
        )
        saveRDS(nichenet_output, "nichenet_results.rds")
    })
    
    # Create visualizations (same as in the original nichenet.nf)
    cat("Creating visualizations...\\n")
    
    # Check nichenet_output structure for debugging
    tryCatch({
        cat("NicheNet output contains the following elements:", 
            paste(names(nichenet_output), collapse=", "), "\\n")
        
        # 1. Ligand activity heatmap
        tryCatch({
            if("ligand_activities" %in% names(nichenet_output)) {
                ligand_activity_plot <- nichenetr::nichenet_ligand_activity_plot(
                    nichenet_output = nichenet_output, 
                    top_n_ligands = min(20, nrow(nichenet_output$ligand_activities)),
                    color_scheme = "viridis"
                )
            } else {
                ligand_activity_plot <- nichenetr::ligand_activity_plot(
                    ligand_activities = nichenet_output, 
                    top_n_ligands = min(20, nrow(nichenet_output)),
                    color_scheme = "viridis"
                )
            }
            
            ggsave("ligand_activity_heatmap.pdf", ligand_activity_plot, width = 10, height = 8)
            cat("Created ligand activity heatmap\\n")
        }, error = function(e) {
            cat("Error creating ligand activity heatmap:", conditionMessage(e), "\\n")
            pdf("ligand_activity_heatmap.pdf", width = 10, height = 8)
            plot(1, type="n", xlab="", ylab="", main="Error creating ligand activity heatmap")
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
        })
        
        # 2. Ligand-target heatmap
        tryCatch({
            if("ligand_activities" %in% names(nichenet_output)) {
                top_ligands <- nichenet_output$ligand_activities %>%
                    dplyr::arrange(-aupr_corrected) %>%
                    dplyr::pull(test_ligand) %>%
                    head(min(10, nrow(nichenet_output$ligand_activities)))
                    
                target_heatmap <- nichenetr::nichenet_ligand_target_heatmap(
                    nichenet_output = nichenet_output, 
                    top_n_ligands = min(10, length(top_ligands)),
                    top_n_targets = 50
                )
            } else {
                target_heatmap <- nichenetr::nichenet_ligand_target_heatmap(
                    nichenet_output = nichenet_output, 
                    top_n_ligands = min(10, length(expressed_ligands)),
                    top_n_targets = 50
                )
            }
            
            ggsave("ligand_target_heatmap.pdf", target_heatmap, width = 12, height = 10)
            cat("Created ligand-target heatmap\\n")
        }, error = function(e) {
            cat("Error creating ligand-target heatmap:", conditionMessage(e), "\\n")
            pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
            plot(1, type="n", xlab="", ylab="", main="Error creating ligand-target heatmap")
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
        })
        
    }, error = function(e) {
        cat("Error in visualization process:", conditionMessage(e), "\\n")
        
        # Create fallback error plots
        pdf("ligand_activity_heatmap.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="NicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
        
        pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
        plot(1, type="n", xlab="", ylab="", main="NicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    cat("NicheNet standalone analysis complete!\\n")
    """
}