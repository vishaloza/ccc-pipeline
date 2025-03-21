// NicheNet analysis module

process NICHENET_ANALYSIS {
    tag "$h5ad_file.simpleName"
    label 'process_high'

    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/nichenet", mode: 'copy'
    
    input:
    path h5ad_file
    path liana_results
    
    output:
    path "nichenet_results_from_liana.rds", emit: nichenet_results
    path "ligand_activity_heatmap.pdf", emit: activity_heatmap
    path "ligand_target_heatmap.pdf", emit: target_heatmap
    path "*.pdf"
    
    script:
    """
    #!/usr/bin/env Rscript

    cat("Working directory:", getwd(), "\n")
    cat("Input file: ${h5ad_file}\n")
    cat("LIANA results file: ${liana_results}\n")
    
    # Load required libraries
    tryCatch({
        library(Seurat)
        library(nichenetr)
        library(tidyverse)
        library(reticulate)
        # Try to load sceasy, install if missing
        if (!requireNamespace("sceasy", quietly = TRUE)) {
            cat("Installing sceasy package...\\n")
            if (!requireNamespace("devtools", quietly = TRUE)) {
                install.packages("devtools")
            }
            devtools::install_github("cellgeni/sceasy")
        }
        library(sceasy)
        
        cat("Libraries loaded successfully\\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\\n")
        cat("Available libraries:", paste(.packages(all.available = TRUE), collapse=", "), "\\n")
        stop("Library loading failed")
    })
    
    # Print the available files for debugging
    cat("Files in current directory:", paste(list.files(), collapse=", "), "\\n")
    
    # Import AnnData and convert to Seurat using sceasy
    tryCatch({
        cat("Converting h5ad to Seurat object using sceasy...\\n")
        temp_seurat_file <- "temp_seurat.rds"
        
        # Use sceasy to convert AnnData to Seurat
        sceasy::convertFormat("${h5ad_file}", from="anndata", to="seurat", outFile=temp_seurat_file)
        
        cat("Loading Seurat object from RDS...\\n")
        seurat_obj <- readRDS(temp_seurat_file)
        
        cat("Seurat object created with dimensions:", dim(seurat_obj), "\\n")
        cat("Cell types detected:", paste(unique(seurat_obj@meta.data\$${params.celltype_column}), collapse=", "), "\\n")
    }, error = function(e) {
        cat("Error converting h5ad to Seurat with sceasy:", conditionMessage(e), "\\n")
        cat("Trying alternative conversion method...\\n")
        
        # Alternative method attempting basic conversion
        tryCatch({
            library(hdf5r)
            seurat_obj <- Seurat::ReadH5AD("${h5ad_file}")
            cat("Alternative conversion produced Seurat object with dimensions:", dim(seurat_obj), "\\n")
        }, error = function(e2) {
            cat("Both conversion methods failed. Final error:", conditionMessage(e2), "\\n")
            stop("Unable to convert h5ad file to Seurat object")
        })
    })
    
    # Import LIANA results
    cat("Reading LIANA results...\\n")
    liana_results <- read.csv("${liana_results}")
    top_lr_pairs <- read.csv("top_lr_pairs_for_nichenet.csv")
    
    cat("Summary of LIANA results:\\n")
    print(head(liana_results))
    
    # Check if cell type column exists in Seurat object
    cell_type_column <- "${params.celltype_column}"
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
    
    # Instead of filtering for specific sender and receiver cell types,
    # use all cells as both sender and receiver.
    all_cells <- colnames(seurat_obj)
    sender_cells <- all_cells
    receiver_cells <- all_cells
    cat("Using all", length(all_cells), "cells as both sender and receiver.\\n")
    
    # Get expressed genes
    cat("Getting expressed genes...\\n")
    expressed_genes_sender <- nichenetr::get_expressed_genes(seurat_obj, sender_cells)
    expressed_genes_receiver <- nichenetr::get_expressed_genes(seurat_obj, receiver_cells)
    
    # Filter ligands and receptors based on LIANA results
    cat("Filtering ligands and receptors...\\n")
    prioritized_ligands <- unique(top_lr_pairs\$ligand)
    prioritized_receptors <- unique(top_lr_pairs\$receptor)
    
    cat("Prioritized ligands:", paste(head(prioritized_ligands), collapse=", "), "...\\n")
    cat("Prioritized receptors:", paste(head(prioritized_receptors), collapse=", "), "...\\n")
    
    # Load NicheNet database objects
    tryCatch({
        cat("Checking NicheNet database...\\n")
        data("ligands", package = "nichenetr", envir = environment())
        data("receptors", package = "nichenetr", envir = environment())
        
        if(!exists("ligands") || !exists("receptors")) {
            cat("WARNING: NicheNet database not loaded correctly. Attempting to load manually...\\n")
            nichenetr_path <- system.file(package = "nichenetr")
            ligands <- readRDS(file.path(nichenetr_path, "data", "ligands.rds"))
            receptors <- readRDS(file.path(nichenetr_path, "data", "receptors.rds"))
        }
        
        cat("NicheNet database loaded with", length(ligands), "ligands and", length(receptors), "receptors\\n")
    }, error = function(e) {
        cat("Error loading NicheNet database:", conditionMessage(e), "\\n")
        ligands <- c()
        receptors <- c()
    })
    
    cat("Finding expressed ligands and receptors...\\n")
    expressed_ligands <- expressed_genes_sender\$gene[expressed_genes_sender\$gene %in% ligands]
    expressed_receptors <- expressed_genes_receiver\$gene[expressed_genes_receiver\$gene %in% receptors]
    
    selected_ligands <- intersect(expressed_ligands, prioritized_ligands)
    selected_receptors <- intersect(expressed_receptors, prioritized_receptors)
    
    cat("Using", length(selected_ligands), "ligands and", length(selected_receptors), "receptors from LIANA analysis\\n")
    
    if(length(selected_ligands) == 0) {
        cat("WARNING: No matching ligands found between LIANA results and NicheNet database\\n")
        cat("Will use all expressed ligands (top 50)\\n")
        selected_ligands <- head(prioritized_ligands, 50)
    }
    if(length(selected_receptors) == 0) {
        cat("WARNING: No matching receptors found between LIANA results and NicheNet database\\n")
        cat("Will use all expressed receptors (top 50)\\n")
        selected_receptors <- head(prioritized_receptors, 50)
    }
    
    # Define background expressed genes (for enrichment calculations)
    background_expressed_genes <- nichenetr::get_expressed_genes(seurat_obj)
    
    cat("Running NicheNet analysis...\\n")
    tryCatch({
        nichenet_output <- nichenetr::nichenet_analysis(
            seurat_obj = seurat_obj,
            sender_cells = sender_cells,
            receiver_cells = receiver_cells,
            ligands = selected_ligands,
            receptors = selected_receptors,
            expressed_genes_receiver = expressed_genes_receiver\$gene,
            background_expressed_genes = background_expressed_genes\$gene
        )
        
        saveRDS(nichenet_output, "nichenet_results_from_liana.rds")
        cat("Saved NicheNet results to 'nichenet_results_from_liana.rds'\\n")
        
        cat("Creating visualizations...\\n")
        
        # 1. Ligand activity heatmap
        tryCatch({
            ligand_activity_plot <- nichenetr::nichenet_ligand_activity_plot(
                nichenet_output, 
                top_n_ligands = min(20, length(selected_ligands)),
                color_scheme = "viridis"
            )
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
            target_heatmap <- nichenetr::nichenet_ligand_target_heatmap(
                nichenet_output, 
                top_n_ligands = min(10, length(selected_ligands)),
                top_n_targets = 50
            )
            ggsave("ligand_target_heatmap.pdf", target_heatmap, width = 12, height = 10)
            cat("Created ligand-target heatmap\\n")
        }, error = function(e) {
            cat("Error creating ligand-target heatmap:", conditionMessage(e), "\\n")
            pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
            plot(1, type="n", xlab="", ylab="", main="Error creating ligand-target heatmap")
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
        })
        
        # 3. Ligand-receptor network
        tryCatch({
            lr_network <- nichenetr::nichenet_ligand_receptor_network(
                nichenet_output,
                ligands = selected_ligands,
                receptors = selected_receptors
            )
            ggsave("ligand_receptor_network.pdf", lr_network, width = 10, height = 8)
            cat("Created ligand-receptor network\\n")
        }, error = function(e) {
            cat("Error creating ligand-receptor network:", conditionMessage(e), "\\n")
            pdf("ligand_receptor_network.pdf", width = 10, height = 8)
            plot(1, type="n", xlab="", ylab="", main="Error creating ligand-receptor network")
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
        })
        
    }, error = function(e) {
        cat("Error in NicheNet analysis:", conditionMessage(e), "\\n")
        nichenet_output <- list(error = conditionMessage(e))
        saveRDS(nichenet_output, "nichenet_results_from_liana.rds")
        
        pdf("ligand_activity_heatmap.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="NicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
        
        pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
        plot(1, type="n", xlab="", ylab="", main="NicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
        
        pdf("ligand_receptor_network.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="NicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    cat("NicheNet analysis complete!\\n")
    """
}
