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
    print(f"Input file: ${input_file}", file=sys.stderr)
    print(f"Sender: ${sender_celltype}", file=sys.stderr)
    print(f"Receiver: ${receiver_celltype}", file=sys.stderr)
    print(f"Cell type column: ${params.celltype_column}", file=sys.stderr)
    
    # Print LIANA version for debugging
    print(f"LIANA version: {li.__version__}", file=sys.stderr)
    
    # Load AnnData object
    adata = sc.read_h5ad("${input_file}")
    
    # Print all available columns for debugging
    print(f"Available columns in adata.obs: {list(adata.obs.columns)}", file=sys.stderr)
    
    # Strictly check if the specified column exists
    if "${params.celltype_column}" not in adata.obs.columns:
        raise ValueError(f"ERROR: The specified cell type column '${params.celltype_column}' does not exist in the dataset. Available columns: {list(adata.obs.columns)}")
    
    # Use exactly the specified column
    celltype_column = "${params.celltype_column}"
    
    # Print unique values in the cell type column for debugging
    print(f"Unique values in '{celltype_column}': {adata.obs[celltype_column].unique()}", file=sys.stderr)
    
    # Check if the sender and receiver cell types exist in the specified column
    if "${sender_celltype}" not in adata.obs[celltype_column].unique():
        raise ValueError(f"ERROR: Sender cell type '${sender_celltype}' not found in column '{celltype_column}'. Available cell types: {adata.obs[celltype_column].unique()}")
    
    if "${receiver_celltype}" not in adata.obs[celltype_column].unique():
        raise ValueError(f"ERROR: Receiver cell type '${receiver_celltype}' not found in column '{celltype_column}'. Available cell types: {adata.obs[celltype_column].unique()}")
    
    # Run LIANA using the pipe function (current API)
    print("Running LIANA analysis...", file=sys.stderr)
    liana_result = li.pipe(
        adata,
        groupby=celltype_column,
        resource_name=['cellchat', 'cellphonedb', 'consensus'],
        verbose=True
    )
    
    # Filter interactions between sender and receiver
    print("Filtering interactions between specified cell types...", file=sys.stderr)
    filtered_interactions = liana_result[
        ((liana_result['source'] == "${sender_celltype}") & (liana_result['target'] == "${receiver_celltype}")) |
        ((liana_result['source'] == "${receiver_celltype}") & (liana_result['target'] == "${sender_celltype}"))
    ]
    
    # Check if we have any results
    if len(filtered_interactions) == 0:
        print(f"WARNING: No interactions found between '${sender_celltype}' and '${receiver_celltype}'", file=sys.stderr)
        print(f"Available source-target pairs in results:", file=sys.stderr)
        for pair in liana_result[['source', 'target']].drop_duplicates().values:
            print(f"  {pair[0]} -> {pair[1]}", file=sys.stderr)
    
    # Rank by specificity score
    print("Ranking interactions...", file=sys.stderr)
    ranked_interactions = filtered_interactions.sort_values(by='magnitude', ascending=False)
    
    # Save full results
    ranked_interactions.to_csv("liana_ranked_interactions.csv", index=False)
    
    # Extract top ligand-receptor pairs
    top_lr_pairs = ranked_interactions.head(50)
    top_lr_pairs[['ligand', 'receptor']].to_csv("top_lr_pairs_for_nichenet.csv", index=False)
    
    # Create dotplot visualization
    print("Creating visualization...", file=sys.stderr)
    try:
        li.pl.dotplot(
            ranked_interactions.head(20),
            source_groups=["${sender_celltype}"],
            target_groups=["${receiver_celltype}"],
            figsize=(12, 10)
        )
        plt.savefig("liana_dotplot.pdf", bbox_inches='tight', dpi=300)
    except Exception as e:
        print(f"Error creating dotplot: {str(e)}", file=sys.stderr)
        # Create a simple alternative visualization
        plt.figure(figsize=(12, 10))
        top_n = min(20, len(ranked_interactions))
        if top_n > 0:
            top_pairs = ranked_interactions.head(top_n)
            plt.barh(
                [f"{row['ligand']}-{row['receptor']}" for _, row in top_pairs.iterrows()],
                top_pairs['magnitude'],
                color='skyblue'
            )
            plt.xlabel('Magnitude Score')
            plt.ylabel('Ligand-Receptor Pair')
            plt.title(f'Top {top_n} L-R Interactions between {sender_celltype} and {receiver_celltype}')
        else:
            plt.text(0.5, 0.5, f"No interactions found between {sender_celltype} and {receiver_celltype}",
                    horizontalalignment='center', verticalalignment='center')
        plt.tight_layout()
        plt.savefig("liana_dotplot.pdf", bbox_inches='tight', dpi=300)
    
    # Export AnnData for NicheNet
    print(f"Writing output to: {os.getcwd()}/for_nichenet.h5ad", file=sys.stderr)
    adata.write_h5ad("for_nichenet.h5ad")
    print(f"Files in current directory: {os.listdir('.')}", file=sys.stderr)
    print("LIANA analysis complete!", file=sys.stderr)
    """
}
