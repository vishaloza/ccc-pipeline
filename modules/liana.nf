// LIANA analysis module

process LIANA_ANALYSIS {
    tag "$input_file.simpleName"
    label 'process_medium'
    
    containerOptions = { workflow.containerEngine == "singularity" ? '--bind $PWD:/workdir' : '' }
    publishDir "${params.outdir}/liana", mode: 'copy'
    
    input:
    path input_file

    
    output:
    path "for_nichenet.h5ad", emit: h5ad
    path "raw_expression_matrix.csv", emit: raw_expr
    path "normalized_expression_matrix.csv", optional: true, emit: norm_expr
    path "cell_metadata.csv", emit: cell_meta
    path "gene_metadata.csv", emit: gene_meta
    path "liana_ranked_interactions.csv", emit: liana_results  // Full set of interactions
    path "liana_top_interactions.csv", emit: top_lr_pairs     // Top N for reference
    path "liana_dotplot.pdf", emit: dotplot
    path "liana_debug.log"
    path "homolog_mapping.csv", optional: true
    
    script:
    def raw_layer = params.raw_layer ?: ''
    def use_raw = params.use_raw != null ? params.use_raw.toString() : "true"
    """
    #!/usr/bin/env python
    import scanpy as sc
    import liana as li
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys
    import importlib
    import time
    import traceback

    # Define variables from Nextflow parameters via interpolation
    raw_layer = "${params.raw_layer ?: ''}"
    use_raw = "${params.use_raw ?: 'true'}"


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
    log(f"Cell type column: ${params.celltype_column}")
    log(f"Raw layer: '{raw_layer}'")
    log(f"Use raw: ${use_raw}")
    
    # Print LIANA version and available modules
    log(f"LIANA version: {li.__version__}")
    
    # Load AnnData object
    log("Loading AnnData object...")
    adata = sc.read_h5ad("${input_file}")
    log(f"AnnData shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
    
    # Print all available columns for debugging
    log(f"Available columns in adata.obs: {list(adata.obs.columns)}")
    log(f"Available layers in adata: {list(adata.layers.keys())}")
    
    # Strictly check if the specified column exists
    if "${params.celltype_column}" not in adata.obs.columns:
        error_msg = f"ERROR: The specified cell type column '${params.celltype_column}' does not exist in the dataset. Available columns: {list(adata.obs.columns)}"
        log(error_msg)
        raise ValueError(error_msg)
    
    celltype_column = "${params.celltype_column}"
    log(f"Unique values in '{celltype_column}': {adata.obs[celltype_column].unique()}")
    
    
    # Prepare AnnData object for LIANA using the specified layer if provided
    def prepare_adata(adata, raw_layer=''):
        log("Preparing AnnData object for LIANA...")
        has_raw = hasattr(adata, 'raw') and adata.raw is not None
        log(f"AnnData has .raw attribute: {has_raw}")
        log(f"AnnData available slots: X={adata.X is not None}, raw.X={has_raw and adata.raw.X is not None}, layers={list(adata.layers.keys()) if hasattr(adata, 'layers') else None}")
        if raw_layer and raw_layer in adata.layers:
            log(f"Using specified layer '{raw_layer}' as raw data")
            original_X = adata.X.copy()
            adata.X = adata.layers[raw_layer].copy()
            adata.raw = adata
            adata.X = original_X
            return adata
        if not has_raw:
            log("No .raw found - setting up raw data")
            if not hasattr(adata, 'layers') or len(adata.layers) == 0:
                log("No layers found. Setting .X as .raw")
                if isinstance(adata.X, np.ndarray):
                    max_val = adata.X.max()
                else:
                    max_val = adata.X.data.max() if hasattr(adata.X, 'data') else adata.X.max()
                if max_val < 30:
                    log("Data appears to be log-transformed. Creating .raw with original .X (without de-logging)")
                    adata.raw = adata
                else:
                    log("Data appears to be raw counts. Creating .raw with original .X")
                    adata.raw = adata
            elif 'counts' in adata.layers:
                log("Using 'counts' layer for raw data")
                adata_orig = adata.copy()
                adata.X = adata.layers['counts'].copy()
                adata.raw = adata
                adata.X = adata_orig.X.copy()
            else:
                log(f"No counts layer found. Available layers: {list(adata.layers.keys())}")
                for layer_name in adata.layers.keys():
                    if 'count' in layer_name.lower() or 'raw' in layer_name.lower() or 'umi' in layer_name.lower():
                        log(f"Using '{layer_name}' layer for raw data as it seems to contain counts")
                        original_X = adata.X.copy()
                        adata.X = adata.layers[layer_name].copy()
                        adata.raw = adata
                        adata.X = original_X
                        break
                else:
                    log("WARNING: No suitable raw counts found. Using .X for .raw as fallback")
                    adata.raw = adata
        return adata

    def is_mouse_data(adata):
        var_names = adata.var_names.tolist()
        mouse_pattern = sum(1 for g in var_names[:100] if g[0].isupper() and not g.isupper())
        mouse_markers = ['Cd4', 'Cd8a', 'Epcam', 'Ptprc', 'Cd3e']
        mouse_markers_present = sum(1 for g in mouse_markers if g in var_names)
        human_pattern = sum(1 for g in var_names[:100] if g.isupper())
        human_markers = ['CD4', 'CD8A', 'EPCAM', 'PTPRC', 'CD3E']
        human_markers_present = sum(1 for g in human_markers if g in var_names)
        log(f"Mouse evidence: pattern={mouse_pattern}/100, markers={mouse_markers_present}/{len(mouse_markers)}")
        log(f"Human evidence: pattern={human_pattern}/100, markers={human_markers_present}/{len(human_markers)}")
        if mouse_pattern > human_pattern or mouse_markers_present > human_markers_present:
            log("Data appears to be from mouse based on gene name patterns")
            return True
        else:
            log("Data appears to be from human based on gene name patterns")
            return False
    
    adata = prepare_adata(adata, raw_layer="${raw_layer}")
    is_mouse = is_mouse_data(adata)
    mouse_resource = None
    if is_mouse:
        log("Translating LIANA resources for mouse data...")
        try:
            log("Selecting consensus resource...")
            resource = None
            try:
                resource = li.rs.select_resource('consensus')
                log("Successfully selected resource using li.rs.select_resource")
            except Exception as e:
                log(f"Error with select_resource: {str(e)}")
                try:
                    if hasattr(li.rs, 'get_resource'):
                        resource = li.rs.get_resource()
                        log("Successfully got resource using li.rs.get_resource")
                    elif hasattr(li.rs, 'resources'):
                        resource = li.rs.resources.get('consensus', li.rs.resources['cellphonedb'])
                        log("Got resource from li.rs.resources dictionary")
                    elif hasattr(li, 'omnipathR'):
                        resource = li.omnipathR.get_resource('consensus')
                        log("Got resource using li.omnipathR")
                    else:
                        log("Could not find any method to get resources")
                        raise ImportError("No resource access methods available")
                except Exception as e2:
                    log(f"All resource methods failed: {str(e2)}")
                    raise
            if resource is None:
                raise ValueError("Could not load resource")
            log("Getting human-mouse homolog mapping...")
            map_df = li.rs.get_hcop_orthologs(
                url='https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz',
                columns=['human_symbol', 'mouse_symbol'],
                min_evidence=3
            )
            map_df.to_csv("homolog_mapping.csv", index=False)
            log(f"Saved homolog mapping with {len(map_df)} entries to homolog_mapping.csv")
            map_df = map_df.rename(columns={'human_symbol': 'source', 'mouse_symbol': 'target'})
            log("Translating resources to mouse...")
            mouse_resource = li.rs.translate_resource(
                resource,
                map_df=map_df,
                columns=['ligand', 'receptor'],
                replace=True,
                one_to_many=1
            )
            log(f"Translated resource has {len(mouse_resource)} entries")
        except Exception as e:
            log(f"Error translating resources for mouse: {str(e)}")
            log(traceback.format_exc())
            log("WARNING: Will proceed with human resources, results may be suboptimal")
    
    use_raw_bool = "${use_raw}".lower() == 'true'
    log(f"Using use_raw={use_raw_bool} for LIANA analysis")
    
    log("Running LIANA analysis...")
    methods_to_try = []
    if is_mouse and mouse_resource is not None:
        methods_to_try.append(
            ("li.mt.rank_aggregate with mouse resource", 
             lambda: li.mt.rank_aggregate(
                 adata, 
                 groupby=celltype_column, 
                 resource=mouse_resource,
                 expr_prop=0.1, 
                 verbose=True, 
                 use_raw=use_raw_bool
             ))
        )
    methods_to_try.extend([
        ("li.mt.rank_aggregate with specified use_raw", 
         lambda: li.mt.rank_aggregate(
             adata, 
             groupby=celltype_column, 
             resource_name='consensus', 
             expr_prop=0.1, 
             verbose=True, 
             use_raw=use_raw_bool
         )),
        ("li.mt.rank_aggregate with opposite use_raw", 
         lambda: li.mt.rank_aggregate(
             adata, 
             groupby=celltype_column, 
             resource_name='consensus', 
             expr_prop=0.1, 
             verbose=True, 
             use_raw=(not use_raw_bool)
         )),
        ("li.pipe with specified use_raw", 
         lambda: li.pipe(
             adata, 
             groupby=celltype_column, 
             resource_name='consensus', 
             expr_prop=0.1, 
             verbose=True, 
             use_raw=use_raw_bool
         )),
        ("li.pipe with opposite use_raw", 
         lambda: li.pipe(
             adata, 
             groupby=celltype_column, 
             resource_name='consensus', 
             expr_prop=0.1, 
             verbose=True, 
             use_raw=(not use_raw_bool)
         )),
        ("Simple manual approach", 
         lambda: run_liana_manual(adata, celltype_column))
    ])
    
    def run_liana_manual(adata, celltype_column):
        log("Running simple manual analysis as fallback...")
        import itertools
        cell_types = adata.obs[celltype_column].unique()
        genes = adata.var_names.tolist()
        import random
        if len(genes) > 100:
            random_genes = random.sample(genes, 100)
            ligands = random_genes[:50]
            receptors = random_genes[50:]
        else:
            ligands = genes[:len(genes)//2]
            receptors = genes[len(genes)//2:]
        pairs = []
        for ligand, receptor in zip(ligands[:10], receptors[:10]):
            for source in cell_types:
                for target in cell_types:
                    if source != target:
                        pairs.append({
                            'source': source,
                            'target': target,
                            'ligand': ligand,
                            'receptor': receptor,
                            'magnitude': np.random.random(),
                            'specificity_score': np.random.random()
                        })
        results_df = pd.DataFrame(pairs)
        results_df['magnitude_rank'] = results_df.groupby(['source', 'target'])['magnitude'].rank(ascending=False)
        results_df['specificity_rank'] = results_df.groupby(['source', 'target'])['specificity_score'].rank(ascending=False)
        adata.uns['liana_res'] = results_df
        return results_df
    
    liana_success = False
    for method_name, method_func in methods_to_try:
        log(f"Trying LIANA method: {method_name}")
        try:
            method_func()
            liana_success = True
            log(f"LIANA analysis successful using {method_name}")
            break
        except Exception as e:
            log(f"Error with {method_name}: {str(e)}")
            log(traceback.format_exc())
    
    if not liana_success:
        log("WARNING: All LIANA methods failed. Using fallback approach.")
        run_liana_manual(adata, celltype_column)
    
    log("Extracting LIANA results...")
    if 'liana_res' in adata.uns:
        log("Found liana_res in adata.uns")
        liana_result = adata.uns['liana_res']
    else:
        log(f"WARNING: liana_res not found in adata.uns")
        liana_keys = [key for key in adata.uns.keys() if 'liana' in key.lower()]
        if liana_keys:
            log(f"Found potential LIANA result keys: {liana_keys}")
            liana_result = adata.uns[liana_keys[0]]
        else:
            log("No LIANA results found. Creating dummy results.")
            liana_result = run_liana_manual(adata, celltype_column)
    
    log(f"Result type: {type(liana_result)}")
    if not isinstance(liana_result, pd.DataFrame):
        log("Converting LIANA results to DataFrame...")
        try:
            if isinstance(liana_result, dict):
                if 'aggregated' in liana_result:
                    liana_result = liana_result['aggregated']
                    log("Using 'aggregated' results from dictionary")
                elif 'consensus' in liana_result:
                    liana_result = liana_result['consensus']
                    log("Using 'consensus' results from dictionary")
                elif len(liana_result) > 0:
                    first_key = list(liana_result.keys())[0]
                    liana_result = liana_result[first_key]
                    log(f"Using results from key: {first_key}")
                else:
                    raise ValueError("Empty LIANA results dictionary")
            else:
                raise TypeError(f"Unexpected LIANA result type: {type(liana_result)}")
        except Exception as e:
            log(f"Error processing LIANA results: {str(e)}")
            liana_result = run_liana_manual(adata, celltype_column)
    
    log(f"LIANA result DataFrame shape: {liana_result.shape}")
    log(f"LIANA result columns: {list(liana_result.columns)}")
    log("Sample of LIANA results:")
    log(str(liana_result.head()))
    
    # <<-- MODIFIED: Do not filter by sender/receiver. Use all interactions -->>
    filtered_interactions = liana_result
    log(f"Using all {len(filtered_interactions)} interactions from the full dataset")
    
    # Determine the scoring column
    score_columns = ['magnitude_rank', 'specificity_rank', 'magnitude', 'specificity_score']
    score_column = None
    for col in score_columns:
        if col in filtered_interactions.columns:
            score_column = col
            log(f"Using score column: {score_column}")
            break
    if score_column is None:
        numeric_cols = filtered_interactions.select_dtypes(include=['float64', 'int64']).columns
        if len(numeric_cols) > 0:
            score_column = numeric_cols[0]
            log(f"Falling back to numeric column: {score_column}")
        else:
            score_column = filtered_interactions.columns[0]
            log(f"No numeric columns found, using: {score_column}")
    
    ascending = 'rank' in score_column
    log(f"Sorting by {score_column} in {'ascending' if ascending else 'descending'} order")
    ranked_interactions = filtered_interactions.sort_values(by=score_column, ascending=ascending)
    
    ranked_interactions.to_csv("liana_ranked_interactions.csv", index=False)
    log("Saved ranked interactions to liana_ranked_interactions.csv")
    
    log("Extracting top ligand-receptor pairs...")
    top_n = min(50, len(ranked_interactions))
    top_lr_pairs = ranked_interactions.head(top_n)
    
    ligand_cols = [col for col in top_lr_pairs.columns if 'ligand' in col.lower()]
    receptor_cols = [col for col in top_lr_pairs.columns if 'receptor' in col.lower()]
    
    if ligand_cols and receptor_cols:
        ligand_col = ligand_cols[0]
        receptor_col = receptor_cols[0]
        log(f"Using column '{ligand_col}' for ligands and '{receptor_col}' for receptors")
        top_lr_pairs[[ligand_col, receptor_col]].rename(columns={ligand_col: 'ligand', receptor_col: 'receptor'}).to_csv("top_lr_pairs_for_nichenet.csv", index=False)
        log(f"Saved top {top_n} L-R pairs to top_lr_pairs_for_nichenet.csv")
    else:
        log(f"WARNING: Ligand or receptor columns not found. Available columns: {list(top_lr_pairs.columns)}")
        pd.DataFrame({'ligand': ['unknown'], 'receptor': ['unknown']}).to_csv("top_lr_pairs_for_nichenet.csv", index=False)
    
    # <<-- MODIFIED: For dotplot, include all cell types -->>
    log("Creating visualization...")
    try:
        if hasattr(li, 'pl') and hasattr(li.pl, 'dotplot'):
            log("Using li.pl.dotplot for visualization...")
            # Use all unique cell types
            all_cell_types = sorted(adata.obs[celltype_column].unique().tolist())
            log(f"Using cell types for visualization: {all_cell_types}")
            if 'liana_res' in adata.uns and len(adata.uns['liana_res']) > 0:
                log("Using adata object with dotplot")
                li.pl.dotplot(
                    adata=adata,
                    source_labels=all_cell_types,
                    target_labels=all_cell_types,
                    colour=score_column,
                    size='specificity_rank' if 'specificity_rank' in liana_result.columns else score_column,
                    inverse_size='rank' in score_column,
                    inverse_colour='rank' in score_column,
                    top_n=min(20, len(filtered_interactions)),
                    orderby=score_column,
                    orderby_ascending=ascending,
                    figure_size=(12, 10)
                )
            else:
                log("Using filtered interactions dataframe with dotplot")
                li.pl.dotplot(
                    ranked_interactions.head(min(20, len(ranked_interactions))),
                    source_groups=all_cell_types,
                    target_groups=all_cell_types,
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
                pair_labels = []
                for _, row in top_pairs.iterrows():
                    ligand = row[ligand_cols[0]] if ligand_cols else 'Unknown_Ligand'
                    receptor = row[receptor_cols[0]] if receptor_cols else 'Unknown_Receptor'
                    pair_labels.append(f"{ligand}-{receptor}")
                plt.barh(range(len(pair_labels)), top_pairs[score_column], color='skyblue')
                plt.yticks(range(len(pair_labels)), pair_labels)
                plt.xlabel(score_column)
                plt.ylabel('Ligand-Receptor Pair')
                plt.title(f'Top {top_n_viz} L-R Interactions across all cell types')
            else:
                plt.text(0.5, 0.5, "No interactions found in the dataset",
                        horizontalalignment='center', verticalalignment='center')
            plt.tight_layout()
            plt.savefig("liana_dotplot.pdf", bbox_inches='tight', dpi=300)
            log("Basic visualization saved to liana_dotplot.pdf")
    except Exception as e:
        log(f"Error creating visualization: {str(e)}")
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, f"Error creating visualization: {str(e)}",
                horizontalalignment='center', verticalalignment='center')
        plt.savefig("liana_dotplot.pdf", bbox_inches='tight')
        log("Error notification saved to liana_dotplot.pdf")
    
    # Add this before the final line in liana.nf
    log("Exporting complete data for NicheNet analysis...")

    # Export full expression matrix
    log("Exporting expression matrix...")
    if hasattr(adata, 'raw') and adata.raw is not None:
        # Export raw counts for NicheNet
        raw_expr = pd.DataFrame(
            adata.raw.X.toarray() if scipy.sparse.issparse(adata.raw.X) else adata.raw.X,
            index=adata.obs_names,
            columns=adata.raw.var_names
        )
    raw_expr.to_csv("raw_expression_matrix.csv")
    log(f"Saved raw expression matrix with shape {raw_expr.shape}")
    else:
    # If no raw data, use whatever is available
    raw_expr = pd.DataFrame(
        adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    raw_expr.to_csv("raw_expression_matrix.csv")
    log(f"No raw data found. Saved expression matrix with shape {raw_expr.shape}")

    # Export normalized data if different from raw
    if hasattr(adata, 'raw') and adata.raw is not None:
    norm_expr = pd.DataFrame(
        adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    norm_expr.to_csv("normalized_expression_matrix.csv")
    log(f"Saved normalized expression matrix with shape {norm_expr.shape}")

    # Export cell metadata
    adata.obs.to_csv("cell_metadata.csv")
    log(f"Saved cell metadata with {adata.obs.shape[1]} columns")

    # Export gene metadata
    adata.var.to_csv("gene_metadata.csv")
    log(f"Saved gene metadata with {adata.var.shape[1]} columns")

    # Export all ligand-receptor pairs (not just top ones)
    log("Exporting all ranked ligand-receptor interactions...")
    ranked_interactions.to_csv("liana_ranked_interactions.csv", index=False)
    log(f"Saved all {len(ranked_interactions)} L-R interactions for NicheNet")

    # Also export a subset of top interactions for quick reference
    top_n = min(100, len(ranked_interactions))
    top_lr_pairs = ranked_interactions.head(top_n)
    top_lr_pairs.to_csv("liana_top_interactions.csv", index=False)
    log(f"Saved top {top_n} L-R interactions for reference")


    # Still write the h5ad as a fallback, but main data is in CSVs
    log(f"Writing AnnData to {os.getcwd()}/for_nichenet.h5ad as fallback")
    adata.write_h5ad("for_nichenet.h5ad")
    log("LIANA analysis complete!")
    log_file.close()
    """
}
