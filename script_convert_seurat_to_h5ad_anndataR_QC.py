#!/usr/bin/env python3
# check_h5ad_full.py

import sys
import anndata
import numpy as np

def check_h5ad_full(filename):
    """Check X, layers, raw, PCA, UMAP, and t-SNE in h5ad file."""
    
    print(f"Reading: {filename}\n")
    adata = anndata.read_h5ad(filename)
    
    # ========== MAIN MATRIX (X) ==========
    print("=" * 70)
    print("MAIN MATRIX (X)")
    print("=" * 70)
    
    if adata.X is not None:
        print(f"✓ X exists")
        print(f"  Type: {type(adata.X).__name__}")
        print(f"  Shape: {adata.X.shape}")
        print(f"  Data type: {adata.X.dtype}")
        
        # Get sample
        sample = adata.X[0:5, 0:5]
        if hasattr(sample, 'toarray'):
            sample = sample.toarray()
        
        print(f"\n  Sample values (first 5x5):")
        print(sample)
        
        # Check if integer or float
        if np.allclose(sample, np.round(sample)):
            print(f"\n  → INTEGER values (likely raw counts)")
        else:
            print(f"\n  → FLOAT values (likely normalized/log-transformed)")
    else:
        print("✗ X is None (MISSING!)")
    
    # ========== LAYERS ==========
    print("\n" + "=" * 70)
    print("LAYERS")
    print("=" * 70)
    
    if adata.layers:
        print(f"\nFound {len(adata.layers)} layer(s):\n")
        
        for layer_name in adata.layers.keys():
            layer_data = adata.layers[layer_name]
            
            print(f"Layer: '{layer_name}'")
            print(f"  Type: {type(layer_data).__name__}")
            print(f"  Shape: {layer_data.shape}")
            print(f"  Data type: {layer_data.dtype}")
            
            # Get sample
            sample = layer_data[0:5, 0:5]
            if hasattr(sample, 'toarray'):
                sample = sample.toarray()
            
            print(f"  Sample values (first 5x5):")
            print(f"    {sample}")
            
            # Check if integer or float
            if np.allclose(sample, np.round(sample)):
                print(f"  → INTEGER values (likely raw counts)")
            else:
                print(f"  → FLOAT values (likely normalized/log-transformed)")
            
            print()
    else:
        print("\nNo layers found")
    
    # ========== RAW COUNTS ==========
    print("=" * 70)
    print("RAW COUNTS (adata.raw)")
    print("=" * 70)
    
    if adata.raw is not None:
        print(f"✓ Raw counts exist")
        print(f"  Type: {type(adata.raw.X).__name__}")
        print(f"  Shape: {adata.raw.X.shape}")
        print(f"  Data type: {adata.raw.X.dtype}")
        
        # Get sample
        sample = adata.raw.X[0:5, 0:5]
        if hasattr(sample, 'toarray'):
            sample = sample.toarray()
        
        print(f"\n  Sample values (first 5x5):")
        print(sample)
        
        # Check if integer
        if np.allclose(sample, np.round(sample)):
            print(f"\n  → INTEGER values (correct for raw counts)")
        else:
            print(f"\n  → FLOAT values (unexpected for raw counts)")
    else:
        print("✗ No raw counts")
    
    # ========== DIMENSIONAL REDUCTIONS ==========
    print("\n" + "=" * 70)
    print("DIMENSIONAL REDUCTIONS (obsm)")
    print("=" * 70)
    
    if adata.obsm:
        print(f"\n✓ Found {len(adata.obsm)} embedding(s):\n")
        
        for key, value in adata.obsm.items():
            print(f"  {key}:")
            print(f"    Type: {type(value).__name__}")
            print(f"    Shape: {value.shape}")
            print(f"    Sample (first 3 cells, first 5 dims):")
            print(f"      {value[0:3, 0:min(5, value.shape[1])]}")
            print()
    else:
        print("\n✗ No embeddings found")
    
    # ========== SPECIFIC EMBEDDING CHECKS ==========
    print("=" * 70)
    print("SPECIFIC EMBEDDING CHECKS")
    print("=" * 70)
    
    has_pca = 'X_pca' in adata.obsm
    has_umap = 'X_umap' in adata.obsm
    has_tsne = 'X_tsne' in adata.obsm
    
    if has_pca:
        print(f"✓ PCA (X_pca):")
        print(f"    Shape: {adata.obsm['X_pca'].shape}")
        print(f"    Number of PCs: {adata.obsm['X_pca'].shape[1]}")
    else:
        print("✗ PCA (X_pca): Not found")
    
    if has_umap:
        print(f"✓ UMAP (X_umap):")
        print(f"    Shape: {adata.obsm['X_umap'].shape}")
    else:
        print("✗ UMAP (X_umap): Not found")
    
    if has_tsne:
        print(f"✓ t-SNE (X_tsne):")
        print(f"    Shape: {adata.obsm['X_tsne'].shape}")
    else:
        print("✗ t-SNE (X_tsne): Not found")
    
    # ========== SUMMARY ==========
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  {'✓' if adata.X is not None else '✗'} Main matrix (X): {'Present' if adata.X is not None else 'MISSING'}")
    print(f"  {'✓' if adata.raw is not None else '✗'} Raw counts: {'Present' if adata.raw is not None else 'MISSING'}")
    print(f"  {'✓' if adata.layers else '○'} Layers: {len(adata.layers) if adata.layers else 0} found")
    print(f"  {'✓' if has_pca else '✗'} PCA: {'Present' if has_pca else 'Missing'}")
    print(f"  {'✓' if has_umap else '✗'} UMAP: {'Present' if has_umap else 'Missing'}")
    print(f"  {'✓' if has_tsne else '✗'} t-SNE: {'Present' if has_tsne else 'Missing'}")
    
    # Recommendations
    if adata.X is None and adata.layers:
        print("\n⚠️  ISSUE DETECTED:")
        print("   X is None but layers exist.")
        print("   Recommendation: Move 'data' layer → X and 'counts' layer → raw")
    
    if adata.X is not None and adata.raw is not None:
        print("\n✓ File structure looks good!")
    elif adata.X is not None and adata.raw is None:
        print("\n⚠️  X exists but raw counts are missing")
    
    print("=" * 70)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python check_h5ad_full.py <filename.h5ad>")
        sys.exit(1)
    
    check_h5ad_full(sys.argv[1])
