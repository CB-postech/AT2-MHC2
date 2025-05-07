### Use mouse alveolar epithelium reference to training Celltypist and predict inhouse dataset

### Training
# select project_lung_exercise
import celltypist
from celltypist import models

import scanpy as sc

adata = sc.read_h5ad('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/merged_w_GSE262927_full.h5ad')

adata.X.expm1().sum(axis = 1)

mouse_alveolar_regeneration = celltypist.train(adata, labels = 'annotation', n_jobs = 32, feature_selection = True)
mouse_alveolar_regeneration.write('/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/alveolar_with_prolif_interferon_Celltypist_model.pkl')

### prediction
import scanpy as sc
import celltypist
from celltypist import models
import pyreadr
import scipy
import pandas as pd
import numpy as np

data_set = 'Alveolar_full'
adata = sc.read_h5ad('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/w_cycling_IFN.h5ad')

from scipy import sparse
adata.X = sparse.csr_matrix(adata.X)
adata.X.expm1().sum(axis = 1)

alveolar_regeneration_reference_prediction = celltypist.annotate(adata, majority_voting = True, over_clustering = 'RNA_snn_res.1', model = '/home/sjcho/datas/reference_atlas/mouse_lung_alveolar/outs/20250217_public_figures/merged_w_GSE262927_full_celltypist_model.pkl')

celltypist.dotplot(alveolar_regeneration_reference_prediction, use_as_reference = 'RNA_snn_res.1', use_as_prediction = 'majority_voting')
celltypist.dotplot(alveolar_regeneration_reference_prediction, use_as_reference = 'RNA_snn_res.1', use_as_prediction = 'predicted_labels')
adata_predicted = alveolar_regeneration_reference_prediction.to_adata()
adata_predicted.obs[['predicted_labels', 'over_clustering', 'majority_voting', 'conf_score']].to_csv('/home/sjcho/projects/AT2_MHC2/20241113_after_cellbender/manuscript_figures/figure4/figure4_and_exd6_20250225/' + 'alv_full_celltypist_result_overcluster1.csv')
