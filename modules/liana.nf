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
    path "liana_debug.log"
    
    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import liana as li
    import pandas as pd
    import matplotlib.pyplot as plt
    import os
    import sys
    import importlib
    import time

    # Set up logging
    log_file = open("liana_debug.log", "w")
    def log(msg):
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
        formatted_msg = f"[{timestamp}] {msg}"
        print(formatted_msg, file=log_file, flush=True)
        print(formatted_msg, file=sys.stderr, flush=True)

    log("Starting LIANA analysis")
    log(f"Working directory: {os.getcwd()}")
    log(f"Input file: ${input_file}")
    log(f"Sender: ${sender_celltype}")
    log(f"Receiver: ${receiver_celltype}")
    log(f"Cell type column: ${params.celltype_column}")
    
    # Print LIANA version and available modules
    log(f"LIANA version: {li.__version__}")
    log("Available LIANA modules:")
    for module_name in dir(li):
        if not module_name.startswith('_'):
            log(f"  - {module_name}")
            submodule = getattr(li, module_name)
            if hasattr(submodule, '__dir__'):
                for func_name in dir(submodule):
                    if not func_name.startswith('_'):
                        log(f"    - {module_name}.{func_name}")
    
    # Load AnnData object
    log("Loading AnnData object...")
    adata = sc.read_h5ad("${input_file}")
    log(f"AnnData shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Print all available columns for debugging
    log(f"Available columns in adata.obs: {list(adata.obs.columns)}")
    
    # Strictly check if the specified column exists
    if "${params.celltype_column}" not in adata.obs.columns:
        error_msg = f"ERROR: The specified cell type column '${params.celltype_column}' does not exist in the dataset. Available columns: {list(adata.obs.columns)}"
        log(error_msg)
        raise ValueError(error_msg)
    
    # Use exactly the specified column
    celltype_column = "${params.celltype_column}"
    
    # Print unique values in the cell type column for debugging
    log(f"Unique values in '{celltype_column}': {adata.obs[celltype_column].unique()}")
    
    # Check if the sender and receiver cell types exist in the specified column
    if "${sender_celltype}" not in adata.obs[celltype_column].unique():
        error_msg = f"ERROR: Sender cell type '${sender_celltype}' not found in column '{celltype_column}'. Available cell types: {adata.obs[celltype_column].unique()}"
        log(error_msg)
        raise ValueError(error_msg)
    
    if "${receiver_celltype}" not in adata.obs[celltype_column].unique():
        error_msg = f"ERROR: Receiver cell type '${receiver_celltype}' not found in column '{celltype_column}'. Available cell types: {adata.obs[celltype_column].unique()}"
        log(error_msg)
        raise ValueError(error_msg)
    
    # Run LIANA using the correct API (following the vignette)
    log("Running LIANA analysis...")
    
    # Step 1: Filter the data using pp.filter
    log("Step 1: Filtering data...")
    li.pp.filter(adata, groupby=celltype_column)
    
    # Step 2: Rank using tl.rank
    log("Step 2: Ranking interactions...")
    li.tl.rank(adata, groupby=celltype_column)
    
    # Step 3: Get the results
    log("Step 3: Getting results...")
    
    # Check if the results are stored in adata.uns
    if 'liana_res' not in adata.uns:
        log("WARNING: 'liana_res' not found in adata.uns after running analysis")
        log(f"Available keys in adata.uns: {list(adata.uns.keys())}")
        
        # Try to find any LIANA-related keys
        liana_keys = [key for key in adata.uns.keys() if 'liana' in key.lower()]
        if liana_keys:
            log(f"Found potential LIANA result keys: {liana_keys}")
            liana_result = adata.uns[liana_keys[0]]
        else:
            error_msg = "No LIANA results found in the AnnData object"
            log(error_msg)
            raise KeyError(error_msg)
    else:
        liana_result = adata.uns['liana_res']
    
    # Convert to DataFrame if needed
    if not isinstance(liana_result, pd.DataFrame):
        log("Converting LIANA results to DataFrame...")
        try:
            if isinstance(liana_result, dict):
                # Try to find the aggregated results
                if 'aggregated' in liana_result:
                    liana_result = liana_result['aggregated']
                elif len(liana_result) > 0:
                    # Take the first item if it's a dict of DataFrames
                    first_key = list(liana_result.keys())[0]
                    liana_result = liana_result[first_key]
                    log(f"Using results from key: {first_key}")
                else:
                    raise ValueError("Empty LIANA results dictionary")
            else:
                raise TypeError(f"Unexpected LIANA result type: {type(liana_result)}")
        except Exception as e:
            log(f"Error processing LIANA results: {str(e)}")
            # Create a minimal result set to continue
            liana_result = pd.DataFrame({
                'source': ['unknown'],
                'target': ['unknown'],
                'ligand': ['unknown'],
                'receptor': ['unknown'],
                'magnitude': [0.0]
            })
    
    log(f"LIANA result DataFrame shape: {liana_result.shape}")
    log(f"LIANA result columns: {list(liana_result.columns)}")
    
    # Filter interactions between sender and receiver
    log("Filtering interactions between specified cell types...")
    filtered_interactions = liana_result[
        ((liana_result['source'] == "${sender_celltype}") & (liana_result['target'] == "${receiver_celltype}")) |
        ((liana_result['source'] == "${receiver_celltype}") & (liana_result['target'] == "${sender_celltype}"))
    ]
    
    log(f"Found {len(filtered_interactions)} interactions between the specified cell types")
    
    # Determine the correct scoring column
    score_columns = ['magnitude', 'specificity_score', 'lr_score', 'score']
    score_column = None
    for col in score_columns:
        if col in filtered_interactions.columns:
            score_column = col
            log(f"Using score column: {score_column}")
            break
    
    if score_column is None:
        # Fall back to the first numeric column
        numeric_cols = filtered_interactions.select_dtypes(include=['float64', 'int64']).columns
        if len(numeric_cols) > 0:
            score_column = numeric_cols[0]
            log(f"Falling back to numeric column: {score_column}")
        else:
            score_column = filtered_interactions.columns[0]
            log(f"No numeric columns found, using: {score_column}")
    
    # Rank by score
    log(f"Ranking interactions by {score_column}...")
    ranked_interactions = filtered_interactions.sort_values(by=score_column, ascending=False)
    
    # Save full results
    ranked_interactions.to_csv("liana_ranked_interactions.csv", index=False)
    log("Saved ranked interactions to liana_ranked_interactions.csv")
    
    # Extract top ligand-receptor pairs
    log("Extracting top ligand-receptor pairs...")
    top_n = min(50, len(ranked_interactions))
    top_lr_pairs = ranked_interactions.head(top_n)
    
    if 'ligand' in top_lr_pairs.columns and 'receptor' in top_lr_pairs.columns:
        top_lr_pairs[['ligand', 'receptor']].to_csv("top_lr_pairs_for_nichenet.csv", index=False)
        log(f"Saved top {top_n} L-R pairs to top_lr_pairs_for_nichenet.csv")
    else:
        log(f"WARNING: 'ligand' or 'receptor' columns not found in results. Available columns: {list(top_lr_pairs.columns)}")
        # Create a minimal output file
        pd.DataFrame({'ligand': ['unknown'], 'receptor': ['unknown']}).to_csv("top_lr_pairs_for_nichenet.csv", index=False)
    
    # Create dotplot visualization
    log("Creating visualization...")
    try:
        if hasattr(li, 'pl') and hasattr(li.pl, 'dotplot'):
            log("Using li.pl.dotplot for visualization...")
            li.pl.dotplot(
                ranked_interactions.head(20),
                source_groups=["${sender_celltype}"],
                target_groups=["${receiver_celltype}"],
                figsize=(12, 10)
            )
            plt.savefig("liana_dotplot.pdf", bbox_inches='tight', dpi=300)
            log("Visualization saved to liana_dotplot.pdf")
        else:
            log("LIANA dotplot function not available, creating a basic plot")
            plt.figure(figsize=(12, 10))
            
            if len(ranked_interactions) > 0:
                top_n_viz = min(20, len(ranked_interactions))
                top_pairs = ranked_interactions.head(top_n_viz)
                
                # Create a simple bar chart
                pair_labels = [f"{row['ligand']}-{row['receptor']}" if 'ligand' in row and 'receptor' in row else f"Pair {i+1}" 
                              for i, (_, row) in enumerate(top_pairs.iterrows())]
                
                plt.barh(range(len(pair_labels)), top_pairs[score_column], color='skyblue')
                plt.yticks(range(len(pair_labels)), pair_labels)
                plt.xlabel(score_column)
                plt.ylabel('Ligand-Receptor Pair')
                plt.title(f'Top {top_n_viz} L-R Interactions between {sender_celltype} and {receiver_celltype}')
            else:
                plt.text(0.5, 0.5, f"No interactions found between {sender_celltype} and {receiver_celltype}",
                        horizontalalignment='center', verticalalignment='center')
            
            plt.tight_layout()
            plt.savefig("liana_dotplot.pdf", bbox_inches='tight', dpi=300)
            log("Basic visualization saved to liana_dotplot.pdf")
    except Exception as e:
        log(f"Error creating visualization: {str(e)}")
        # Create a simple error notification plot
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, f"Error creating visualization: {str(e)}",
                horizontalalignment='center', verticalalignment='center')
        plt.savefig("liana_dotplot.pdf", bbox_inches='tight')
        log("Error notification saved to liana_dotplot.pdf")
    
    # Export AnnData for NicheNet
    log(f"Writing AnnData to: {os.getcwd()}/for_nichenet.h5ad")
    adata.write_h5ad("for_nichenet.h5ad")
    log("LIANA analysis complete!")
    
    # Close log file
    log_file.close()
    """
}
