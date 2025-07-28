#!/bin/bash



#cd /data/shoreline/Simulator_datasets/

# Fit intensity versus variability from metaSPARSim's R1 preset:
Rscript scripts/Fit_I-V_NB.R fitting_variability

# Fit three distributions to the number of passes for Kinnex (16S) and Pfizer (Titan) datasets:
bash scripts/00_Fit_np_distributions.sh > NP_fitting_log.txt

# Run initial analysis (just to get edit distance distributions -- v. lazy syntax)
bash scripts/01_run_initial_analysis.sh 

# Run subread accuracy optimization:
bash scripts/02_run_all_linesearches.sh > linesearch_log.txt

# Run, time, and profile memory usage for final analyses, with optimized values:
/usr/bin/time -v bash ./scripts/03_final_analyses.sh > final_log.txt

# Generate publication figures, collate to Publication_figures directory
python scripts/Generate_publication_figures.py > figs_log.txt
python scripts/Plot_amplitypes.py

cp NP_analysis/np_distribution_comparison.png Publication_figures/
cp fitting_variability/nb_model_fit.png Publication_figures/fitting_variability.png






