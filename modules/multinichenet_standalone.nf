// Standalone MultiNicheNet analysis module for direct Seurat object input

process MULTINICHENET_ANALYSIS_STANDALONE {
    tag "$input_file.simpleName"
    label 'process_high'

    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/multinichenet", mode: 'copy'
    
    input:
    path input_file
    
    output:
    path "multinichenet_results.rds", emit: multinichenet_results
    path "ligand_activity_heatmap.pdf", emit: activity_heatmap
    path "ligand_target_heatmap.pdf", emit: target_heatmap
    path "*.pdf"
    path "*.csv"
    
    script:
    """
    #!/usr/bin/env Rscript

    cat("Working directory:", getwd(), "\\n")
    cat("Input file:", "${input_file}", "\\n")
    
    # Load required libraries
    tryCatch({
        library(Seurat)
        library(tidyverse)
        library(reticulate)   
        library(Matrix)
        # Install and load MultiNicheNet
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
        }
        if (!requireNamespace("devtools", quietly = TRUE)) {
            install.packages("devtools")
        }
        if (!requireNamespace("nichenetr", quietly = TRUE)) {
            devtools::install_github("saeyslab/nichenetr")
        }
        if (!requireNamespace("multinichenetr", quietly = TRUE)) {
            devtools::install_github("saeyslab/multinichenetr")
        }
        if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
            BiocManager::install("SingleCellExperiment")
        }
        if (!requireNamespace("muscat", quietly = TRUE)) {
            BiocManager::install("muscat")
        }
        if (!requireNamespace("LRBaseDb", quietly = TRUE)) {
            BiocManager::install("LRBaseDb")
        }
        library(SingleCellExperiment)
        library(multinichenetr)
        library(nichenetr)
        library(muscat)
        library(LRBaseDb)
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
    
    # Convert Seurat to SingleCellExperiment
    cat("Converting Seurat to SingleCellExperiment...\\n")
    sce <- as.SingleCellExperiment(seurat_obj)
    
    # Extract condition information from metadata
    # Define sample and group IDs based on available columns
    # Default to the sample column if available, otherwise try alternatives
    if ("sample" %in% colnames(sce@colData)) {
        sample_id <- "sample"
    } else if ("orig.ident" %in% colnames(sce@colData)) {
        sample_id <- "orig.ident"
    } else {
        # Try to find a suitable column for sample
        possible_sample_cols <- c("Sample", "sample_id", "donor", "patient", "individual")
        for (col in possible_sample_cols) {
            if (col %in% colnames(sce@colData)) {
                sample_id <- col
                break
            }
        }
        
        # If still not found, create a dummy sample column
        if (!exists("sample_id")) {
            sce$sample <- "sample1"
            sample_id <- "sample"
            cat("WARNING: No sample column found. Created a dummy 'sample' column.\\n")
        }
    }
    
    # Try to find group/condition column
    if ("condition" %in% colnames(sce@colData)) {
        group_id <- "condition"
    } else if ("group" %in% colnames(sce@colData)) {
        group_id <- "group"
    } else {
        # Try to find a suitable column for group
        possible_group_cols <- c("Group", "condition_id", "disease", "treatment", "state")
        for (col in possible_group_cols) {
            if (col %in% colnames(sce@colData)) {
                group_id <- col
                break
            }
        }
        
        # If still not found, create a dummy group column
        if (!exists("group_id")) {
            sce$group <- "groupA"
            group_id <- "group"
            cat("WARNING: No group/condition column found. Created a dummy 'group' column.\\n")
        }
    }
    
    cat("Using sample ID column:", sample_id, "\\n")
    cat("Using group ID column:", group_id, "\\n")
    cat("Using cell type column:", cell_type_column, "\\n")
    
    # Run MultiNicheNet analysis
    tryCatch({
        cat("Running MultiNicheNet analysis...\\n")
        
        # Check if we have multiple conditions
        unique_groups <- unique(sce@colData[[group_id]])
        
        if (length(unique_groups) < 2) {
            cat("WARNING: Only one group/condition found:", paste(unique_groups, collapse=", "), 
                "\\nMultiNicheNet requires at least two groups for comparison.\\n")
            cat("Will proceed with a basic analysis, but results may not be meaningful.\\n")
            
            # Create a dummy second group by splitting the first group
            cell_indices <- 1:ncol(sce)
            group_a_indices <- cell_indices[1:floor(length(cell_indices)/2)]
            group_b_indices <- cell_indices[(floor(length(cell_indices)/2)+1):length(cell_indices)]
            
            # Assign new groups
            sce@colData[[group_id]] <- "groupA"
            sce@colData[[group_id]][group_b_indices] <- "groupB"
            unique_groups <- c("groupA", "groupB")
            
            cat("Created dummy groups 'groupA' and 'groupB' for analysis\\n")
        }
        
        # Get cell abundances
        cat("Calculating cell abundances...\\n")
        abundances <- multinichenetr::get_abundance_data(
            sce, 
            sample_id = sample_id,
            group_id = group_id,
            celltype_id = cell_type_column
        )
        
        # Filter cell types
        cat("Filtering cell types based on abundances...\\n")
        filtered_abundances <- multinichenetr::filter_celltypes_by_abundance(
            abundances, 
            min_cells_per_sample = 5,
            min_frac_samples = 0.5,
            min_frac_groups = 1
        )
        
        # Visualize cell type abundances
        cat("Generating abundance plots...\\n")
        p_abund <- multinichenetr::plot_abundances(
            filtered_abundances,
            filtered = TRUE,
            type = "boxplot"
        )
        
        # Save abundance plot
        pdf("cell_abundances.pdf", width = 10, height = 8)
        print(p_abund)
        dev.off()
        
        # Perform differential expression analysis
        cat("Running differential expression analysis...\\n")
        pb_res <- multinichenetr::calculate_de_genes_pb(
            sce,
            sample_id = sample_id,
            group_id = group_id,
            celltype_id = cell_type_column,
            group_comparison_for_de = unique_groups[1:2],  # Compare the first two groups
            filtering_results = filtered_abundances
        )
        
        # Summarize DE results
        cat("Summarizing differential expression results...\\n")
        summary_de <- multinichenetr::summarize_de_results(pb_res)
        
        # Save DE summary
        write.csv(summary_de, "de_summary.csv")
        
        # Generate sender-receiver combinations
        cat("Generating sender-receiver combinations...\\n")
        sender_receiver_combis <- multinichenetr::generate_sender_receiver_combinations(
            filtered_abundances,
            is_receiver = function(x) TRUE,  # All cell types can be receivers
            is_sender = function(x) TRUE     # All cell types can be senders
        )
        
        # Prepare ligand-receptor interactions
        cat("Preparing ligand-receptor interactions...\\n")
        lr_network_filtered <- LRBaseDb::getLRBaseData()
        
        # Run MultiNicheNet prioritization
        cat("Running MultiNicheNet ligand-target analysis...\\n")
        prioritization_tables <- multinichenetr::multi_nichenet_per_sender_receiver(
            sender_receiver_list = sender_receiver_combis,
            de_results = pb_res,
            lr_network = lr_network_filtered,
            expr_data = sce,
            sample_id = sample_id,
            group_id = group_id,
            celltype_id = cell_type_column
        )
        
        # Extract top ligands
        cat("Extracting top ligands...\\n")
        top_ligands <- multinichenetr::get_top_ligands(
            prioritization_tables,
            n_ligands = 20
        )
        
        # Save MultiNicheNet results
        saveRDS(list(
            prioritization_tables = prioritization_tables,
            abundances = abundances,
            filtered_abundances = filtered_abundances,
            de_results = pb_res,
            de_summary = summary_de,
            top_ligands = top_ligands
        ), "multinichenet_results.rds")
        
        # Create visualizations
        cat("Creating visualizations...\\n")
        
        # Top ligands heatmap
        pdf("ligand_activity_heatmap.pdf", width = 12, height = 10)
        multinichenetr::plot_ligand_heatmap(
            prioritization_tables,
            top_n = 20
        )
        dev.off()
        
        # Create ligand-target heatmap for top ligands
        pdf("ligand_target_heatmap.pdf", width = 14, height = 12)
        multinichenetr::plot_ligand_target_heatmap(
            prioritization_tables,
            top_n_ligands = 10,
            top_n_targets = 50
        )
        dev.off()
        
        # Create network visualization
        pdf("lr_network_visualization.pdf", width = 12, height = 10)
        multinichenetr::plot_lr_network(
            prioritization_tables,
            top_n = 15
        )
        dev.off()
        
        cat("MultiNicheNet analysis complete!\\n")
    }, error = function(e) {
        cat("ERROR in MultiNicheNet analysis:", conditionMessage(e), "\\n")
        cat("Error details:", e$message, "\\n")
        cat("Error call:", deparse(e$call), "\\n")
        
        # Create minimal output so pipeline can continue
        multinichenet_output <- list(
            error = conditionMessage(e)
        )
        saveRDS(multinichenet_output, "multinichenet_results.rds")
        
        # Create fallback error plots
        pdf("ligand_activity_heatmap.pdf", width = 10, height = 8)
        plot(1, type="n", xlab="", ylab="", main="MultiNicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
        
        pdf("ligand_target_heatmap.pdf", width = 12, height = 10)
        plot(1, type="n", xlab="", ylab="", main="MultiNicheNet Analysis Failed")
        text(1, 1, conditionMessage(e), cex=0.8)
        dev.off()
    })
    """
}