import scanpy as sc
import scanorama
import sys
from pathlib import Path

input_dir = sys.argv[1]
output_dir = sys.argv[2]
dataset_names = sys.argv[3:]

print('='*50)
print('Loading QC-filtered datasets...')
print('='*50)

adatas = []
h5_files = sorted(Path(input_dir).glob('*_qc_filtered.h5'))

if len(h5_files) != len(dataset_names):
    print(f'Error: Found {len(h5_files)} files but expected {len(dataset_names)}')
    sys.exit(1)

for i, h5_file in enumerate(h5_files):
    print(f'Loading {h5_file.name}...')
    adata = sc.read_10x_h5(str(h5_file))
    adata.obs['dataset'] = dataset_names[i]
    adata.obs['orig.ident'] = dataset_names[i]
    adatas.append(adata)
    print(f'  {dataset_names[i]}: {adata.n_obs} cells x {adata.n_vars} genes')

print('\n' + '='*50)
print('Running Scanorama integration...')
print('='*50)

scanorama.integrate_scanpy(adatas, dimred=50)

print('\nMerging datasets...')
integrated = sc.concat(adatas, label='batch', keys=dataset_names)
print(f'Integrated: {integrated.n_obs} cells x {integrated.n_vars} genes')

print('\nRunning UMAP...')
sc.pp.neighbors(integrated, use_rep='X_scanorama', n_neighbors=30)
sc.tl.umap(integrated)
sc.tl.leiden(integrated, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)

Path(output_dir).mkdir(parents=True, exist_ok=True)
h5ad_file = Path(output_dir) / 'integrated_scanorama.h5ad'
integrated.write_h5ad(h5ad_file)
print(f'\nSaved: {h5ad_file}')
print('\nIntegration complete!')
