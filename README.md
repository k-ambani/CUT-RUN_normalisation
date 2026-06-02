# CUT-RUN_normalisation
Code for the paper _"Benchmarking normalisation methods for differential binding analysis in CUT&amp;RUN"_. 

The code that processed the raw files, conducted the differential analysis on the BRG1 and ER experiments, simulation experiments for the CTCF dataset and the code for figure generation are available here.

To view the html files in this repo (particular the figure generation files), prepend to the following URL including the https:// in the link  _https://html-preview.github.io/?url=_

BRG1_dataset 

  - 01_raw_files_processing - workflow to process raw fastQ files to BAM and peak files. 
  - 02_normalisation_method_testing - code to generate scale factors and conduct differential analysis using the normalisation methods  
  - 03_figure_generation - code to generate figures displayed in paper

CTCF_dataset  
  - 01_raw_files_processing - workflow used to process raw fastQ files to BAM and peak files.
  - 02_negative_experiment_analysis -  code to test normalisation methods and code for figure generation
  - 03_redistribution_experiment_analysis -  code to test normalisation methods and code for figure generation
  - 04_global_shift_experiment_analysis -  code to test normalisation methods and code for figure generation

ER_dataset  
  - 01_raw_files_processing - workflow used to process raw fastQ files to BAM and peak files.
  - 02_normalisation_methods_testing - code to test normalsiation methods
  - 03_figure_generation - code for figure generation
