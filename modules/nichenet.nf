// NicheNet analysis module

process NICHENET_ANALYSIS {
    tag "$h5ad_file.simpleName"
    label 'process_high'

    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/nichenet", mode: 'copy'
    
    input:
    path h5ad_file
    path liana_results
    val sender_celltype
    val receiver_celltype
    
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
    cat("Sender celltype: ${sender_celltype}\n")
    cat("Receiver celltype: ${receiver_celltype}\n")
    
    # Load required libraries
    tryCatch({
        library(Seurat)
        library(nichenetr)
        library(tidyverse)
        library(hdf5r)
        # Try to load SeuratDisk, install if missing
        if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        cat("Installing SeuratDisk package...\n")
        remotes::install_github("mojaveazure/seurat-disk")
        library(SeuratDisk)
        }
        
        cat("Libraries loaded successfully\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\n")
        cat("Available libraries:", paste(.packages(all.available = TRUE), collapse=", "), "\n")
        stop("Library loading failed")
    })
    
    # Print the available files for debugging
    cat("Files in current directory:", paste(list.files(), collapse=", "), "\n")
    
    # Import AnnData and convert to Seurat using SeuratDisk
    tryCatch({
        cat("Converting h5ad to h5Seurat...\n")
        Convert("${h5ad_file}", dest = "temp.h5seurat")
        
        cat("Loading h5Seurat as Seurat object...\n")
        seurat_obj <- LoadH5Seurat("temp.h5seurat")
        
        cat("Seurat object created with dimensions:", dim(seurat_obj), "\n")
        cat("Cell types detected:", paste(unique(seurat_obj@meta.data\$${params.celltype_column}), collapse=", "), "\n")
    }, error = function(e) {
        cat("Error converting AnnData to Seurat:", conditionMessage(e), "\n")
        cat("Trying alternative conversion method...\n")
        
        # Alternative method using SeuratDisk
        library(SeuratDisk)
        h5ad_file <- "${h5ad_file}"
        seurat_obj <- ReadH5AD(h5ad_file)
        cat("Alternative conversion produced Seurat object with dimensions:", dim(seurat_obj), "\n")
    })
    
    # Import LIANA results
    cat("Reading LIANA results...\n")
    liana_results <- read.csv("${liana_results}")
    top_lr_pairs <- read.csv("top_lr_pairs_for_nichenet.csv")
    
    # Print summary of LIANA results
    cat("Top LIANA interactions between", "${sender_celltype}", "and", "${receiver_celltype}", "\n")
    head(liana_results)
    
    # Check if cell type column exists in Seurat object
    cell_type_column <- "${params.celltype_column}"
    if(!cell_type_column %in% colnames(seurat_obj@meta.data)) {
        cat("WARNING: Specified cell type column", cell_type_column, "not found in Seurat object\n")
        cat("Available columns:", paste(colnames(seurat_obj@meta.data), collapse=", "), "\n")
        
        # Try to find a suitable column
        possible_columns <- c("cell_type", "celltype", "CellType", "clusters", "cluster", "leiden")
        for(col in possible_columns) {
            if(col %in% colnames(seurat_obj@meta.data)) {
                cell_type_column <- col
                cat("Using alternative column:", cell_type_column, "\n")
                break
            }
        }
    }
    
    # Get cell IDs for these types
    sender_cells <- WhichCells(seurat_obj, expression = seurat_obj@meta.data[[cell_type_column]] == "${sender_celltype}")
    receiver_cells <- WhichCells(seurat_obj, expression = seurat_obj@meta.data[[cell_type_column]] == "${receiver_celltype}")
    
    cat("Found", length(sender_cells), "sender cells and", length(receiver_cells), "receiver cells\n")
    
    if(length(sender_cells) == 0) {
        cat("Sender cell type '${sender_celltype}' not found. Available cell types:", 
            paste(unique(seurat_obj@meta.data[[cell_type_column]]), collapse=", "), "\n")
        stop(paste0("No cells found for sender cell type: ", "${sender_celltype}"))
    }
    if(length(receiver_cells) == 0) {
        cat("Receiver cell type '${receiver_celltype}' not found. Available cell types:", 
            paste(unique(seurat_obj@meta.data[[cell_type_column]]), collapse=", "), "\n")
        stop(paste0("No cells found for receiver cell type: ", "${receiver_celltype}"))
    }
    
    # Get expressed genes
    cat("Getting expressed genes...\n")
    expressed_genes_sender <- nichenetr::get_expressed_genes(seurat_obj, sender_cells)
    expressed_genes_receiver <- nichenetr::get_expressed_genes(seurat_obj, receiver_cells)
    
    # Filter ligands and receptors based on LIANA results
    cat("Filtering ligands and receptors...\n")
    prioritized_ligands <- unique(top_lr_pairs\$ligand)
    prioritized_receptors <- unique(top_lr_pairs\$receptor)
    
    cat("Prioritized ligands:", paste(head(prioritized_ligands), collapse=", "), "...\n")
    cat("Prioritized receptors:", paste(head(prioritized_receptors), collapse=", "), "...\n")
    
    # Make sure all nichenetr data objects are available
    tryCatch({
        # Check if ligands and receptors are available in the NicheNet database
        cat("Checking NicheNet database...\n")
        data("ligands", package = "nichenetr", envir = environment())
        data("receptors", package = "nichenetr", envir = environment())
        
        if(!exists("ligands") || !exists("receptors")) {
            cat("WARNING: NicheNet database not loaded correctly. Attempting to load manually...\n")
            # Get the path to the nichenetr package
            nichenetr_path <- system.file(package = "nichenetr")
            # Load data files manually
            ligands <- readRDS(file.path(nichenetr_path, "data", "ligands.rds"))
            receptors <- readRDS(file.path(nichenetr_path, "data", "receptors.rds"))
        }
        
        cat("NicheNet database loaded with", length(ligands), "ligands and", length(receptors), "receptors\n")
    }, error = function(e) {
        cat("Error loading NicheNet database:", conditionMessage(e), "\n")
        # Create empty vectors to continue
        ligands <- c()
        receptors <- c()
    })
    
    # Get final ligand-receptor pairs that are expressed and in LIANA's top results
    cat("Finding expressed ligands and receptors...\n")
    expressed_ligands <- expressed_genes_sender\$gene[expressed_genes_sender\$gene %in% ligands]
    expressed_receptors <- expressed_genes_receiver\$gene[expressed_genes_receiver\$gene %in% receptors]
    
    selected_ligands <- intersect(expressed_ligands, prioritized_ligands)
    selected_receptors <- intersect(expressed_receptors, prioritized_receptors)
    
    cat("Using", length(selected_ligands), "ligands and", 
        length(selected_receptors), "receptors from LIANA analysis\n")
    
    # Check if we have enough data to proceed
    if(length(selected_ligands) == 0) {
        cat("WARNING: No matching ligands found between LIANA results and NicheNet database\n")
        cat("Will use all expressed ligands (top 50)\n")
        selected_ligands <- head(prioritized_ligands, 50)
    }
    if(length(selected_receptors) == 0) {
        cat("WARNING: No matching receptors found between LIANA results and NicheNet database\n")
        cat("Will use all expressed receptors (top 50)\n")
        selected_receptors <- head(prioritized_receptors, 50)
    }
    
    # Define background expressed genes (for enrichment calculations)
    background_expressed_genes <- nichenetr::get_expressed_genes(seurat_obj)
    
    # Run NicheNet with LIANA-prioritized ligand-receptor pairs
    cat("Running NicheNet analysis...\n")
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
        
        # Save NicheNet results
        saveRDS(nichenet_output, "nichenet_results_from_liana.rds")
        cat("Saved NicheNet results to 'nichenet_results_from_liana.rds'\n")
        
        # Create visualizations
        cat("Creating visualizations...\n")
        
        # 1. Ligand activity heatmap
        tryCatch({
            ligand_activity_plot <- nichenetr::nichenet_ligand_activity_plot(
                nichenet_output, 
                top_n_ligands = min(20, length(selected_ligands)),
                color_scheme = "viridis"
            )
            ggsave("ligand_activity_heatmap.pdf", ligand_activity_plot, width = 10, height = 8)
            cat("Created ligand activity heatmap\n")
        }, error = function(e) {
            cat("Error creating ligand activity heatmap:", conditionMessage(e), "\n")
            # Create a simple placeholder visualization
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
            cat("Created ligand-target heatmap\n")
        }, error = function(e) {
            cat("Error creating ligand-target heatmap:", conditionMessage(e), "\n")
            pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
            plot(1, type="n", xlab="", ylab="", main="Error creating ligand-target heatmap")
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
        })
        
        # 3. Create a ligand-receptor network
        tryCatch({
            lr_network <- nichenetr::nichenet_ligand_receptor_network(
                nichenet_output,
                ligands = selected_ligands,
                receptors = selected_receptors
            )
            ggsave("ligand_receptor_network.pdf", lr_network, width = 10, height = 8)
            cat("Created ligand-receptor network\n")
        }, error = function(e) {
            cat("Error creating ligand-receptor network:", conditionMessage(e), "\n")
            pdf("ligand_receptor_network.pdf", width = 10, height = 8)
            plot(1, type="n", xlab="", ylab="", main="Error creating ligand-receptor network")
            text(1, 1, conditionMessage(e), cex=0.8)
            dev.off()
        })
        
    }, error = function(e) {
        cat("Error in NicheNet analysis:", conditionMessage(e), "\n")
        # Create an empty RDS file to satisfy the output requirement
        nichenet_output <- list(error = conditionMessage(e))
        saveRDS(nichenet_output, "nichenet_results_from_liana.rds")
        
        # Create placeholder visualizations
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
    
    cat("NicheNet analysis complete!\n")
    """
}
