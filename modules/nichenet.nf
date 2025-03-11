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
    
    # Load required libraries
    library(Seurat)
    library(nichenetr)
    library(tidyverse)
    library(anndata)
    
    # Import AnnData and convert to Seurat
    print("Importing AnnData object...")
    adata <- read_h5ad("${h5ad_file}")
    seurat_obj <- as.Seurat(adata, counts = "X")

    # Print the available files for debugging
    cat("Files in current directory:", paste(list.files(), collapse=", "), "\n")
    
    # Import LIANA results
    liana_results <- read.csv("${liana_results}")
    top_lr_pairs <- read.csv("top_lr_pairs_for_nichenet.csv")
    
    # Print summary of LIANA results
    print(paste("Top LIANA interactions between", 
                "${sender_celltype}", "and", "${receiver_celltype}"))
    print(head(liana_results[, c("ligand", "receptor", "source", "target", "specificity_score")], 10))
    
    # Get cell IDs for these types
    sender_cells <- WhichCells(seurat_obj, expression = cell_type == "${sender_celltype}")
    receiver_cells <- WhichCells(seurat_obj, expression = cell_type == "${receiver_celltype}")
    
    if(length(sender_cells) == 0) {
        stop(paste0("No cells found for sender cell type: ", "${sender_celltype}"))
    }
    if(length(receiver_cells) == 0) {
        stop(paste0("No cells found for receiver cell type: ", "${receiver_celltype}"))
    }
    
    # Get expressed genes
    expressed_genes_sender <- nichenetr::get_expressed_genes(seurat_obj, sender_cells)
    expressed_genes_receiver <- nichenetr::get_expressed_genes(seurat_obj, receiver_cells)
    
    # Filter ligands and receptors based on LIANA results
    prioritized_ligands <- unique(top_lr_pairs\$ligand)
    prioritized_receptors <- unique(top_lr_pairs\$receptor)
    
    # Get final ligand-receptor pairs that are expressed and in LIANA's top results
    selected_ligands <- intersect(
        expressed_genes_sender\$gene[expressed_genes_sender\$gene %in% nichenetr::ligands],
        prioritized_ligands
    )
    
    selected_receptors <- intersect(
        expressed_genes_receiver\$gene[expressed_genes_receiver\$gene %in% nichenetr::receptors],
        prioritized_receptors
    )
    
    print(paste("Using", length(selected_ligands), "ligands and", 
                length(selected_receptors), "receptors from LIANA analysis"))
    
    if(length(selected_ligands) == 0 || length(selected_receptors) == 0) {
        stop("No matching ligands or receptors found between LIANA results and NicheNet database")
    }
    
    # Define background expressed genes (for enrichment calculations)
    background_expressed_genes <- nichenetr::get_expressed_genes(seurat_obj)
    
    # Run NicheNet with LIANA-prioritized ligand-receptor pairs
    print("Running NicheNet analysis...")
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
    print("Saved NicheNet results to 'nichenet_results_from_liana.rds'")
    
    # Create visualizations
    print("Creating visualizations...")
    
    # 1. Ligand activity heatmap
    ligand_activity_plot <- nichenetr::nichenet_ligand_activity_plot(nichenet_output, 
                                                top_n_ligands = min(20, length(selected_ligands)),
                                                color_scheme = "viridis")
    ggsave("ligand_activity_heatmap.pdf", ligand_activity_plot, width = 10, height = 8)
    
    # 2. Ligand-target heatmap
    target_heatmap <- nichenetr::nichenet_ligand_target_heatmap(nichenet_output, 
                                                top_n_ligands = min(10, length(selected_ligands)),
                                                top_n_targets = 50)
    ggsave("ligand_target_heatmap.pdf", target_heatmap, width = 12, height = 10)
    
    # 3. Create a ligand-receptor network
    lr_network <- nichenetr::nichenet_ligand_receptor_network(nichenet_output,
                                                ligands = selected_ligands,
                                                receptors = selected_receptors)
    ggsave("ligand_receptor_network.pdf", lr_network, width = 10, height = 8)
    
    print("NicheNet analysis complete!")
    """
}
