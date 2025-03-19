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
    log("LIANA available modules:")
    for module_name in dir(li):
        if not module_name.startswith('_'):
            log(f"  - {module_name}")
            try:
                submodule = getattr(li, module_name)
                if hasattr(submodule, '__dir__'):
                    for func_name in dir(submodule):
                        if not func_name.startswith('_'):
                            log(f"    - {module_name}.{func_name}")
            except Exception as e:
                log(f"    Error inspecting module: {str(e)}")
    
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
    
    # Run LIANA using the correct API according to the vignette
    log("Running LIANA analysis with mt.rank_aggregate...")
    
    try:
        # First check if mt.rank_aggregate exists
        if hasattr(li, 'mt') and hasattr(li.mt, 'rank_aggregate'):
            log("Using li.mt.rank_aggregate as shown in the vignette...")
            li.mt.rank_aggregate(
                adata,
                groupby=celltype_column,
                resource_name='consensus',
                expr_prop=0.1,
                verbose=True
            )
        # If that fails, try alternative methods
        elif hasattr(li, 'rank_aggregate'):
            log("Using li.rank_aggregate...")
            li.rank_aggregate(
                adata,
                groupby=celltype_column,
                resource_name='consensus',
                expr_prop=0.1,
                verbose=True
            )
        # Last resort - look for any method that might work
        else:
            log("WARNING: Could not find rank_aggregate function. Looking for alternatives...")
            if hasattr(li, 'method') and hasattr(li.method, 'liana_pipe'):
                log("Using li.method.liana_pipe...")
                li.method.liana_pipe(
                    adata,
                    groupby=celltype_column,
                    resource_name=['CellChat', 'CellPhoneDB', 'Consensus'],
                    verbose=True
                )
            elif hasattr(li, 'pipe'):
                log("Using li.pipe...")
                li.pipe(
                    adata,
                    groupby=celltype_column,
                    verbose=True
                )
            else:
                error_msg = "Could not find any suitable LIANA analysis function"
                log(error_msg)
                raise AttributeError(error_msg)
    except Exception as e:
        log(f"Error during LIANA analysis: {str(e)}")
        raise
    
    # Get the results from adata.uns
    log("Extracting LIANA results...")
    
    # Check if results are in adata.uns
    if 'liana_res' in adata.uns:
        log("Found liana_res in adata.uns")
        liana_result = adata.uns['liana_res']
    else:
        log(f"WARNING: liana_res not found in adata.uns")
        log(f"Available keys in adata.uns: {list(adata.uns.keys())}")
        
        # Try to find any LIANA-related results
        liana_keys = [key for key in adata.uns.keys() if 'liana' in key.lower()]
        if liana_keys:
            log(f"Found potential LIANA result keys: {liana_keys}")
            liana_result = adata.uns[liana_keys[0]]
        else:
            # Last resort - check if rank_aggregate results might be named differently
            rank_keys = [key for key in adata.uns.keys() if 'rank' in key.lower()]
            if rank_keys:
                log(f"Found potential rank result keys: {rank_keys}")
                liana_result = adata.uns[rank_keys[0]]
            else:
                error_msg = "No LIANA results found in the AnnData object"
                log(error_msg)
                raise KeyError(error_msg)
    
    log(f"Result type: {type(liana_result)}")
    
    # Convert to DataFrame if needed
    if not isinstance(liana_result, pd.DataFrame):
        log("Converting LIANA results to DataFrame...")
        try:
            if isinstance(liana_result, dict):
                # Try to find the aggregated results
                if 'aggregated' in liana_result:
                    liana_result = liana_result['aggregated']
                    log("Using 'aggregated' results from dictionary")
                elif 'consensus' in liana_result:
                    liana_result = liana_result['consensus']
                    log("Using 'consensus' results from dictionary")
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
                'magnitude_rank': [0.0]
            })
    
    log(f"LIANA result DataFrame shape: {liana_result.shape}")
    log(f"LIANA result columns: {list(liana_result.columns)}")
    
    # Save a sample of the results for debugging
    log("Sample of LIANA results:")
    log(liana_result.head().to_string())
    
    # Filter interactions between sender and receiver
    log("Filtering interactions between specified cell types...")
    filtered_interactions = liana_result[
        ((liana_result['source'] == "${sender_celltype}") & (liana_result['target'] == "${receiver_celltype}")) |
        ((liana_result['source'] == "${receiver_celltype}") & (liana_result['target'] == "${sender_celltype}"))
    ]
    
    log(f"Found {len(filtered_interactions)} interactions between the specified cell types")
    
    # Determine the correct scoring column - according to the vignette, we should use magnitude_rank
    score_columns = ['magnitude_rank', 'specificity_rank', 'magnitude', 'specificity_score']
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
    
    # For rank columns, lower is better, so we need to sort in ascending order
    ascending = 'rank' in score_column
    log(f"Sorting by {score_column} in {'ascending' if ascending else 'descending'} order")
    
    # Rank by score
    ranked_interactions = filtered_interactions.sort_values(by=score_column, ascending=ascending)
    
    # Save full results
    ranked_interactions.to_csv("liana_ranked_interactions.csv", index=False)
    log("Saved ranked interactions to liana_ranked_interactions.csv")
    
    # Extract top ligand-receptor pairs
    log("Extracting top ligand-receptor pairs...")
    top_n = min(50, len(ranked_interactions))
    top_lr_pairs = ranked_interactions.head(top_n)
    
    # Get the actual ligand/receptor column names
    # They might be 'ligand'/'receptor' or 'ligand_complex'/'receptor_complex'
    ligand_cols = [col for col in top_lr_pairs.columns if 'ligand' in col.lower()]
    receptor_cols = [col for col in top_lr_pairs.columns if 'receptor' in col.lower()]
    
    if ligand_cols and receptor_cols:
        ligand_col = ligand_cols[0]
        receptor_col = receptor_cols[0]
        log(f"Using column '{ligand_col}' for ligands and '{receptor_col}' for receptors")
        top_lr_pairs[[ligand_col, receptor_col]].to_csv("top_lr_pairs_for_nichenet.csv", index=False)
        log(f"Saved top {top_n} L-R pairs to top_lr_pairs_for_nichenet.csv")
    else:
        log(f"WARNING: Ligand or receptor columns not found. Available columns: {list(top_lr_pairs.columns)}")
        # Create a minimal output file
        pd.DataFrame({'ligand': ['unknown'], 'receptor': ['unknown']}).to_csv("top_lr_pairs_for_nichenet.csv", index=False)
    
    # Create dotplot visualization using the API shown in the vignette
    log("Creating visualization...")
    try:
        if hasattr(li, 'pl') and hasattr(li.pl, 'dotplot'):
            log("Using li.pl.dotplot for visualization...")
            li.pl.dotplot(
                adata=adata,
                source_labels=["${sender_celltype}"],
                target_labels=["${receiver_celltype}"],
                colour=score_column,
                size='specificity_rank' if 'specificity_rank' in liana_result.columns else score_column,
                inverse_size='rank' in score_column,
                inverse_colour='rank' in score_column,
                top_n=20,
                orderby=score_column,
                orderby_ascending=ascending,
                figure_size=(12, 10)
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
                pair_labels = []
                for _, row in top_pairs.iterrows():
                    # Try to get ligand-receptor pair label
                    ligand = row[ligand_cols[0]] if ligand_cols else 'Unknown_Ligand'
                    receptor = row[receptor_cols[0]] if receptor_cols else 'Unknown_Receptor'
                    pair_labels.append(f"{ligand}-{receptor}")
                
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
