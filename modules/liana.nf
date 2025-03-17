// LIANA analysis module

process LIANA_ANALYSIS {
    tag "$input_file.simpleName"
    label 'process_medium'
    
    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/liana", mode: 'copy'
    
    input:
    path input_file
    val sender_celltype
    val receiver_celltype
    
    output:
    path "for_nichenet.h5ad", emit: h5ad
    path "liana_ranked_interactions.csv", emit: liana_results
    path "top_lr_pairs_for_nichenet.csv", emit: top_lr_pairs
    path "liana_dotplot.pdf", emit: dotplot
    
    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import liana as li
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    import sys

    print(f"Working directory: {os.getcwd()}", file=sys.stderr)
    
    # Load AnnData object
    adata = sc.read_h5ad("${input_file}")
    
    # Check if the specified cell type column exists
    celltype_column = "${params.celltype_column}"
    if celltype_column not in adata.obs.columns:
        if 'cell_type' in adata.obs.columns:
            print(f"Using 'cell_type' as cell type column instead of {celltype_column}")
            celltype_column = 'cell_type'
        elif 'leiden' in adata.obs.columns:
            print(f"Using 'leiden' as cell type column instead of {celltype_column}")
            celltype_column = 'leiden'
        elif 'clusters' in adata.obs.columns:
            print(f"Using 'clusters' as cell type column instead of {celltype_column}")
            celltype_column = 'clusters'
        else:
            raise ValueError(f"Specified cell type column '{celltype_column}' not found, and no suitable alternatives found")
    
    print(f"Using column '{celltype_column}' for cell type annotations")
    
    # Verify cell types exist
    if "${sender_celltype}" not in adata.obs[celltype_column].unique():
        raise ValueError(f"Sender cell type '${sender_celltype}' not found in data")
    if "${receiver_celltype}" not in adata.obs[celltype_column].unique():
        raise ValueError(f"Receiver cell type '${receiver_celltype}' not found in data")
    
    # Run LIANA
    liana_result = li.liana_pipe(
        adata,
        groupby=celltype_column,
        resource_name=['CellChat', 'CellPhoneDB', 'Consensus'],
        verbose=True
    )
    
    # Filter interactions between sender and receiver
    filtered_interactions = liana_result[
        ((liana_result['source'] == "${sender_celltype}") & (liana_result['target'] == "${receiver_celltype}")) |
        ((liana_result['source'] == "${receiver_celltype}") & (liana_result['target'] == "${sender_celltype}"))
    ]
    
    # Rank by specificity score
    ranked_interactions = filtered_interactions.sort_values(by='specificity_score', ascending=False)
    
    # Save full results
    ranked_interactions.to_csv("liana_ranked_interactions.csv", index=False)
    
    # Extract top ligand-receptor pairs
    top_lr_pairs = ranked_interactions.head(50)
    top_lr_pairs[['ligand', 'receptor']].to_csv("top_lr_pairs_for_nichenet.csv", index=False)
    
    # Create dotplot visualization
    li.pl.dotplot(
        ranked_interactions.head(20),
        source_groups=["${sender_celltype}"],
        target_groups=["${receiver_celltype}"],
        figsize=(12, 10)
    )
    plt.savefig("liana_dotplot.pdf", bbox_inches='tight', dpi=300)
    
    # Export AnnData for NicheNet
    print(f"Writing output to: {os.getcwd()}/for_nichenet.h5ad", file=sys.stderr)
    adata.write_h5ad("for_nichenet.h5ad")
    print(f"Files in current directory: {os.listdir('.')}", file=sys.stderr)
    """
}
