// LIANA analysis module

process LIANA_ANALYSIS {
    tag "$input_file.simpleName"
    label 'process_medium'
    
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
    
    # Load AnnData object
    adata = sc.read_h5ad("${input_file}")
    
    # Check if cell type annotations exist
    if 'cell_type' not in adata.obs.columns:
        if 'leiden' in adata.obs.columns:
            print("Using 'leiden' as cell type column")
            adata.obs['cell_type'] = adata.obs['leiden']
        elif 'clusters' in adata.obs.columns:
            print("Using 'clusters' as cell type column")
            adata.obs['cell_type'] = adata.obs['clusters']
        else:
            raise ValueError("No cell type annotations found in AnnData object")
    
    # Verify cell types exist
    if "${sender_celltype}" not in adata.obs['cell_type'].unique():
        raise ValueError(f"Sender cell type '${sender_celltype}' not found in data")
    if "${receiver_celltype}" not in adata.obs['cell_type'].unique():
        raise ValueError(f"Receiver cell type '${receiver_celltype}' not found in data")
    
    # Run LIANA
    liana_result = li.liana_pipe(
        adata,
        groupby='cell_type',
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
    adata.write_h5ad("for_nichenet.h5ad")
    """
}
