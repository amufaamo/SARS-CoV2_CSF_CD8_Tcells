#env= cpdb

import pandas as pd
import glob
import os

cpdb_version = 'v5.0.0'
cpdb_target_dir = os.path.join('cpdbdata', cpdb_version)

from cellphonedb.utils import db_utils

db_utils.download_database(cpdb_target_dir, cpdb_version)

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

#patient
cpdb_resultsp = cpdb_statistical_analysis_method.call(
         cpdb_file_path = 'cpdbdata/v5.0.0/cellphonedb.zip',
         meta_file_path = 'patient/CellphoneDB_patient_meta.csv',
         counts_file_path = 'patient/CellphoneDB_patient_count.txt',
         counts_data = 'hgnc_symbol',
         score_interactions = True,
         threshold = 0.1,
         threads = 8,
         output_path = 'results_patient')

#hc
cpdb_resultshc = cpdb_statistical_analysis_method.call(
         cpdb_file_path = 'cpdbdata/v5.0.0/cellphonedb.zip',
         meta_file_path = 'hc/CellphoneDB_hc_meta.csv',
         counts_file_path = 'hc/CellphoneDB_hc_count.txt',
         counts_data = 'hgnc_symbol',
         score_interactions = True,
         threshold = 0.1,
         threads = 8,
         output_path = 'results_hc')


import anndata as ad
import ktplotspy as kpy
import matplotlib.pyplot as plt

#heatmap https://ktplotspy.readthedocs.io/en/latest/notebooks/tutorial.html

adata = ad.read_text("patient/CellphoneDB_patient_count.txt")
means = pd.read_csv("results_patient/statistical_analysis_means_06_21_2024_115853.txt", sep="\t")
pvals = pd.read_csv("results_patient/statistical_analysis_pvalues_06_21_2024_115853.txt", sep="\t")
decon = pd.read_csv("results_patient/statistical_analysis_deconvoluted_06_21_2024_115853.txt", sep="\t")

#bidirection
kpy.plot_cpdb_heatmap(pvals=pvals, figsize=(5, 5), title="Sum of significant interactions")
plt.show()

#monodirection
kpy.plot_cpdb_heatmap(
    pvals=pvals,
    figsize=(5, 5),
    title="Sum of significant interactions",
    symmetrical=False,
)
plt.show()

