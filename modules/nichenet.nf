// NicheNet analysis module

process NICHENET_ANALYSIS {
    tag "$h5ad_file.simpleName"
    label 'process_high'

    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/nichenet", mode: 'copy'
    
    input:
    path h5ad_file
    path raw_expr
    path norm_expr
    path cell_meta
    path gene_meta
    path liana_results
    path top_lr_pairs
    
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
        library(Matrix)
        cat("Libraries loaded successfully\\n")
    }, error = function(e) {
        cat("Error loading libraries:", conditionMessage(e), "\\n")
        cat("Available libraries:", paste(.packages(all.available = TRUE), collapse=", "), "\\n")
        stop("Library loading failed")
    })
    
    # Print the available files for debugging
    cat("Files in current directory:", paste(list.files(), collapse=", "), "\\n")
    
    # Import expression output and convert to Seurat 
    tryCatch({
        cat("Loading data from CSV files...\n")
        
        # Check if files exist
        raw_expression_file <- "raw_expression_matrix.csv"
        cell_metadata_file <- "cell_metadata.csv"
        
        files_exist <- file.exists(raw_expression_file) && file.exists(cell_metadata_file)
        
        if (!files_exist) {
            stop("Required CSV files not found. Make sure LIANA exported the necessary data files.")
        }
        
        # Read raw expression matrix
        raw_expr <- read.csv(raw_expression_file, row.names=1, check.names=FALSE)
        raw_expr_sparse <- as(as.matrix(raw_expr), "dgCMatrix")
        # Transpose efficiently
        raw_expr_sparse_t <- t(raw_expr_sparse)
        cat("Transposed raw expr matrix...\n")

        
        # Read cell metadata
        cat("Reading cell metadata...\n")  
        cell_metadata <- read.csv(cell_metadata_file, row.names=1, check.names=FALSE)
        
        # Create Seurat object with raw counts
        cat("Creating Seurat object...\n")
        seurat_obj <- CreateSeuratObject(counts = raw_expr_sparse_t, meta.data = cell_metadata)

        # Define cell_type_column
        cell_type_column <- "${params.celltype_column}"  # Use the Nextflow parameter

        cat("Seurat object created with dimensions:", dim(seurat_obj), "\n")

        # Check if cell type column exists before trying to access it
        if(cell_type_column %in% colnames(seurat_obj@meta.data)) {
            cat("Cell types detected:", paste(unique(seurat_obj@meta.data[[cell_type_column]]), collapse=", "), "\n")
        } else {
            cat("WARNING: Cell type column", cell_type_column, "not found. Available columns:", paste(colnames(seurat_obj@meta.data), collapse=", "), "\n")
        }

        cat("Seurat object created with dimensions:", dim(seurat_obj), "\n")
        cat("Cell types detected:", paste(unique(seurat_obj@meta.data[[cell_type_column]]), collapse=", "), "\n")
    }, error = function(e) {
        cat("Error loading CSV data:", conditionMessage(e), "\n")
        
        # Try to load h5ad as fallback
        cat("Attempting to load h5ad as fallback...\n")
        tryCatch({
            library(hdf5r)
            seurat_obj <- Seurat::ReadH5AD("${h5ad_file}")
            cat("Successfully loaded h5ad file as fallback\n")
        }, error = function(e2) {
            cat("Both CSV and h5ad methods failed. Final error:", conditionMessage(e2), "\n")
            stop("Unable to load input data for NicheNet analysis")
        })
    })
        
    # Import LIANA results
    cat("Reading LIANA results...\\n")
    #liana_results <- read.csv("${liana_results}")
    liana_results <- read.csv("liana_ranked_interactions.csv", header=TRUE)
    cat("Loaded", nrow(liana_results), "ligand-receptor interactions from LIANA\n")
    
    # Extract all unique ligands and receptors
    ligands_from_liana <- unique(liana_results[["ligand_complex"]])  
    receptors_from_liana <- unique(liana_results[["receptor_complex"]]) 
    
    cat("Using", length(ligands_from_liana), "unique ligands and", 
    length(receptors_from_liana), "unique receptors from LIANA\n")
    #top_lr_pairs <- read.csv("top_lr_pairs_for_nichenet.csv")
    
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

# Get Idents from cell type column for consistency
Idents(seurat_obj) <- seurat_obj@meta.data[[cell_type_column]]
cat("Available Idents in Seurat object:", paste(levels(Idents(seurat_obj)), collapse=", "), "\n")

# Get the list of unique cell types
unique_cell_types <- levels(Idents(seurat_obj))
cat("Using these cell types instead:", paste(unique_cell_types, collapse=", "), "\n")

# Get expressed genes using the cell types, not the barcodes
cat("Getting expressed genes by cell type...\n")
expressed_genes_list <- nichenetr::get_expressed_genes(
    unique_cell_types,  # Pass cell types, not barcodes
    seurat_obj,          # The Seurat object
    pct = 0.1                         # The percentage threshold
)

# Convert to the format expected by other NicheNet functions
cat("Converting to expected format...\n")
expressed_genes_sender <- data.frame("gene" = expressed_genes_list)
expressed_genes_receiver <- expressed_genes_sender  # Same for receiver
cat("Got all expressed genes...\n")

# Filter ligands and receptors based on LIANA results
cat("Filtering ligands and receptors...\n")
prioritized_ligands <- ligands_from_liana
prioritized_receptors <- receptors_from_liana

cat("Prioritized ligands:", paste(head(prioritized_ligands), collapse=", "), "...\n")
cat("Prioritized receptors:", paste(head(prioritized_receptors), collapse=", "), "...\n")

# Load NicheNet database objects
tryCatch({
    cat("Checking NicheNet database...\n")
    
    # Try to load from package first
    data("ligands", package = "nichenetr", envir = environment())
    data("receptors", package = "nichenetr", envir = environment())
    
    # If not found, download directly from Zenodo
    if(!exists("ligands") || !exists("receptors")) {
        cat("WARNING: NicheNet database not found in package. Downloading from Zenodo...\n")
        
        # Create a temporary file for downloading
        temp_ligands <- tempfile(fileext = ".rds")
        temp_receptors <- tempfile(fileext = ".rds")
        
        # Download files
        download.file("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds", temp_ligands, mode = "wb")
        download.file("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds", temp_receptors, mode = "wb")
        
        # Load the downloaded files
        ligands <- readRDS(temp_ligands)
        receptors <- readRDS(temp_receptors)
        
        # Clean up temp files
        unlink(temp_ligands)
        unlink(temp_receptors)
    }
    
    cat("NicheNet database loaded with", length(ligands), "ligands and", length(receptors), "receptors\n")
}, error = function(e) {
    cat("Error loading NicheNet database:", conditionMessage(e), "\n")
    
    # Create empty vectors if download fails
    if(!exists("ligands")) ligands <- c()
    if(!exists("receptors")) receptors <- c()
})
cat("Finding expressed ligands and receptors...\n")
expressed_ligands <- expressed_genes_sender[["gene"]][expressed_genes_sender[["gene"]] %in% ligands]
expressed_receptors <- expressed_genes_receiver[["gene"]][expressed_genes_receiver[["gene"]] %in% receptors]

selected_ligands <- intersect(expressed_ligands, prioritized_ligands)
selected_receptors <- intersect(expressed_receptors, prioritized_receptors)

cat("Using", length(selected_ligands), "ligands and", length(selected_receptors), "receptors from LIANA analysis\n")

if(length(selected_ligands) == 0) {
    cat("WARNING: No matching ligands found between LIANA results and NicheNet database\n")
    cat("Will use all expressed ligands (top 100)\n")
    selected_ligands <- head(prioritized_ligands, 100)
}
if(length(selected_receptors) == 0) {
    cat("WARNING: No matching receptors found between LIANA results and NicheNet database\n")
    cat("Will use all expressed receptors (top 100)\n")
    selected_receptors <- head(prioritized_receptors, 100)
}

# Define background expressed genes (for enrichment calculations)
background_genes_list <- nichenetr::get_expressed_genes(
    unique_cell_types,
    seurat_obj
)
background_expressed_genes <- data.frame("gene" = background_genes_list)


# Get sender and receiver cells based on cell types
cat("Getting sender and receiver cell indices...\n")
sender_cell_indices <- WhichCells(seurat_obj, idents = unique_cell_types)
receiver_cell_indices <- sender_cell_indices  # Use same cells for both

# Load NicheNet resources
tryCatch({
    cat("Downloading NicheNet resources from Zenodo...\n")
    ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
    
    cat("Successfully loaded ligand_target_matrix with dimensions:", dim(ligand_target_matrix)[1], "x", dim(ligand_target_matrix)[2], "\n")
    cat("Successfully loaded weighted_networks\n")
}, error = function(e) {
    cat("ERROR loading NicheNet resources:", conditionMessage(e), "\n")
    stop("Failed to load required NicheNet resources. Cannot continue analysis.")
})

cat("Running NicheNet analysis...\n")
tryCatch({
    nichenet_output <- nichenetr::nichenet_seuratobj_aggregate(
    seurat_obj = seurat_obj,
    receiver = cell_type_column,  # The cell type column name
    sender = unique_cell_types,   # All cell types as senders
    condition_colname = NULL,     # No condition comparison (NULL if not doing DE analysis)
    lr_network = data.frame(from = ligands_from_liana, to = receptors_from_liana) %>% 
        distinct(from, to),
    ligand_target_matrix = ligand_target_matrix,  # You'll need to load this
    weighted_networks = weighted_networks
    )
    
    saveRDS(nichenet_output, "nichenet_results_from_liana.rds")
    cat("Saved NicheNet results to 'nichenet_results_from_liana.rds'\n")
}, error = function(e) {
    cat("ERROR in NicheNet analysis:", conditionMessage(e), "\n")

    # Create minimal output so pipeline can continue
    nichenet_output <- list(
        error = conditionMessage(e),
        ligand_activities = data.frame(test_ligand = head(ligands_from_liana, 10), 
                                      auroc = NA, aupr = NA, 
                                      aupr_corrected = NA, pearson = NA, 
                                      rank = 1:10)
    )
    saveRDS(nichenet_output, "nichenet_results_from_liana.rds")
})

    

# Create visualizations
cat("Creating visualizations...\n")

# Check nichenet_output structure for debugging
tryCatch({
    cat("NicheNet output contains the following elements:", paste(names(nichenet_output), collapse=", "), "\n")
    
    if("ligand_activities" %in% names(nichenet_output)) {
        cat("Found ligand_activities with", nrow(nichenet_output$ligand_activities), "rows\n")
    }
    
    # 1. Ligand activity heatmap
    tryCatch({
        # Check if we need to use the new visualization function structure
        if("ligand_activities" %in% names(nichenet_output)) {
            # For newer nichenet_seuratobj_aggregate output
            ligand_activity_plot <- nichenetr::nichenet_ligand_activity_plot(
                nichenet_output = nichenet_output, 
                top_n_ligands = min(20, nrow(nichenet_output$ligand_activities)),
                color_scheme = "viridis"
            )
        } else {
            # Fallback for direct ligand activity matrix
            ligand_activity_plot <- nichenetr::ligand_activity_plot(
                ligand_activities = nichenet_output, 
                top_n_ligands = min(20, nrow(nichenet_output)),
                color_scheme = "viridis"
            )
        }
        
        ggsave("ligand_activity_heatmap.pdf", ligand_activity_plot, width = 10, height = 8)
        cat("Created ligand activity heatmap\n")
    }, error = function(e) {
        cat("Error creating ligand activity heatmap:", conditionMessage(e), "\n")
        pdf("ligand_activity_heatmap.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="Error creating ligand activity heatmap")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # 2. Ligand-target heatmap
    tryCatch({
        # Get top ligands based on output structure
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
            # Fallback for other output structures
            target_heatmap <- nichenetr::nichenet_ligand_target_heatmap(
                nichenet_output = nichenet_output, 
                top_n_ligands = min(10, length(selected_ligands)),
                top_n_targets = 50
            )
        }
        
        ggsave("ligand_target_heatmap.pdf", target_heatmap, width = 12, height = 10)
        cat("Created ligand-target heatmap\n")
    }, error = function(e) {
        cat("Error creating ligand-target heatmap:", conditionMessage(e), "\n")
        cat("Error details:", e$message, "\n")
        cat("Error call:", deparse(e$call), "\n")
        
        pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
        plot(1, type="n", xlab="", ylab="", main="Error creating ligand-target heatmap")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    
    # 3. Ligand-receptor network
    tryCatch({
        # Adjust for different output structures
        if("ligand_activities" %in% names(nichenet_output) && 
           "expressed_receptors" %in% names(nichenet_output)) {
            
            top_ligands <- nichenet_output$ligand_activities %>%
                dplyr::arrange(-aupr_corrected) %>%
                dplyr::pull(test_ligand) %>%
                head(min(10, nrow(nichenet_output$ligand_activities)))
                
            lr_network <- nichenetr::nichenet_ligand_receptor_network(
                nichenet_output = nichenet_output,
                ligands = top_ligands,
                receptors = nichenet_output$expressed_receptors
            )
        } else {
            lr_network <- nichenetr::nichenet_ligand_receptor_network(
                nichenet_output = nichenet_output,
                ligands = selected_ligands,
                receptors = selected_receptors
            )
        }
        
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
    cat("Error in visualization process:", conditionMessage(e), "\n")
    
    # Create fallback error plots
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
}
"""
