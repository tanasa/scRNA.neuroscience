#!/usr/bin/env python3
"""
Read h5ad file, check for velocity matrices, and save them in STARsolo format
"""

import scanpy as sc
import scipy.io
import pandas as pd
import os

# ============================================================================
# STEP 1: Read h5ad file and check for velocity layers
# ============================================================================

print("=" * 70)
print("READING H5AD FILE")
print("=" * 70)

# Read the file
adata = sc.read_h5ad("JK4488_01_starsolo_filtered_matrix.h5ad")

print(f"\n✓ File loaded successfully")
print(f"  Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")

# Check for velocity layers
print("\nAvailable layers:", list(adata.layers.keys()))

# ============================================================================
# STEP 2: Check if spliced/unspliced exist
# ============================================================================

print("\n" + "=" * 70)
print("CHECKING FOR VELOCITY MATRICES")
print("=" * 70)

has_spliced = 'spliced' in adata.layers
has_unspliced = 'unspliced' in adata.layers
has_ambiguous = 'ambiguous' in adata.layers

print(f"\n{'✅' if has_spliced else '❌'} spliced: {has_spliced}")
print(f"{'✅' if has_unspliced else '❌'} unspliced: {has_unspliced}")
print(f"{'✅' if has_ambiguous else '❌'} ambiguous: {has_ambiguous}")

# Verdict
if has_spliced and has_unspliced:
    print("\n" + "=" * 70)
    print("✅✅✅ YES - Has spliced/unspliced matrices!")
    print("=" * 70)
    print(f"\n  Spliced shape: {adata.layers['spliced'].shape}")
    print(f"  Unspliced shape: {adata.layers['unspliced'].shape}")
    if has_ambiguous:
        print(f"  Ambiguous shape: {adata.layers['ambiguous'].shape}")
    
    # ========================================================================
    # STEP 3: Save matrices to MTX format (STARsolo/velocyto format)
    # ========================================================================
    
    print("\n" + "=" * 70)
    print("SAVING MATRICES TO MTX FORMAT")
    print("=" * 70)
    
    # Create output directory
    output_dir = "velocity_matrices_output"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}/")
    
    # Save spliced matrix
    print("\nSaving matrices...")
    scipy.io.mmwrite(
        os.path.join(output_dir, "spliced.mtx"),
        adata.layers['spliced'].T  # Transpose to genes × cells (MTX format)
    )
    print("  ✓ spliced.mtx")
    
    # Save unspliced matrix
    scipy.io.mmwrite(
        os.path.join(output_dir, "unspliced.mtx"),
        adata.layers['unspliced'].T
    )
    print("  ✓ unspliced.mtx")
    
    # Save ambiguous if exists
    if has_ambiguous:
        scipy.io.mmwrite(
            os.path.join(output_dir, "ambiguous.mtx"),
            adata.layers['ambiguous'].T
        )
        print("  ✓ ambiguous.mtx")
    
    # Save barcodes (cell names)
    barcodes_df = pd.DataFrame(adata.obs_names)
    barcodes_df.to_csv(
        os.path.join(output_dir, "barcodes.tsv"),
        sep='\t',
        header=False,
        index=False
    )
    print("  ✓ barcodes.tsv")
    
    # Save features (gene names)
    # Create features.tsv with gene_id and gene_name
    if hasattr(adata.var, 'gene_ids'):
        # If gene IDs exist in var
        features_df = pd.DataFrame({
            'gene_id': adata.var['gene_ids'],
            'gene_name': adata.var_names,
            'feature_type': 'Gene Expression'
        })
    else:
        # Otherwise use gene names as IDs
        features_df = pd.DataFrame({
            'gene_id': adata.var_names,
            'gene_name': adata.var_names,
            'feature_type': 'Gene Expression'
        })
    
    features_df.to_csv(
        os.path.join(output_dir, "features.tsv"),
        sep='\t',
        header=False,
        index=False
    )
    print("  ✓ features.tsv")
    
    # ========================================================================
    # STEP 4: Verify saved files
    # ========================================================================
    
    print("\n" + "=" * 70)
    print("VERIFICATION - Files Created:")
    print("=" * 70)
    print(f"\n{output_dir}/")
    
    expected_files = ['spliced.mtx', 'unspliced.mtx', 'features.tsv', 'barcodes.tsv']
    if has_ambiguous:
        expected_files.insert(2, 'ambiguous.mtx')
    
    for i, filename in enumerate(expected_files):
        filepath = os.path.join(output_dir, filename)
        exists = os.path.exists(filepath)
        
        # Tree symbol
        symbol = "└──" if i == len(expected_files) - 1 else "├──"
        
        # Status
        status = "✅" if exists else "❌"
        
        # File size
        if exists:
            size = os.path.getsize(filepath)
            if size < 1024:
                size_str = f"{size} B"
            elif size < 1024**2:
                size_str = f"{size/1024:.1f} KB"
            elif size < 1024**3:
                size_str = f"{size/(1024**2):.1f} MB"
            else:
                size_str = f"{size/(1024**3):.1f} GB"
            
            print(f"{symbol} {filename:<20} {status}  ({size_str})")
        else:
            print(f"{symbol} {filename:<20} {status}")
    
    print("\n" + "=" * 70)
    print("✅ EXPORT COMPLETE!")
    print("=" * 70)
    print(f"\nYou can now use these files with velocyto.R:")
    print(f"  setwd('{output_dir}')")
    print(f"  spliced <- readMM('spliced.mtx')")
    print(f"  unspliced <- readMM('unspliced.mtx')")
    print("=" * 70)

else:
    print("\n" + "=" * 70)
    print("❌ NO - This file does NOT have velocity matrices")
    print("=" * 70)
    print("\nThis is just a standard filtered count matrix")
    print("Available data:")
    print(f"  Main matrix (X): {adata.X.shape}")
    if adata.layers:
        print(f"  Layers: {list(adata.layers.keys())}")
    else:
        print("  Layers: None")
    
    print("\nTo get velocity matrices, you need to:")
    print("  1. Run STARsolo with --soloFeatures Velocyto")
    print("  2. Or use velocyto command-line tool")
    print("=" * 70)
