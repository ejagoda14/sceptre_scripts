*To run sceptre itself*

These are the scripts I use to run sceptre for dc-tap-seq data on OAK

To set up the environmnet:
1. create a conda environment from the r_env.yml:
conda-env create -f envs/r_env.yml -p envs/R_env
2. with environment active, open R and install sceptre following the developers instructions
install.packages("devtools")
devtools::install_github("katsevich-lab/sceptre")


Run sceptre using the R script: new_sceptre_for_server.R ags1-13

The inputs to this script are as follows [some of these are dumby variables for future commits, that should not be changed]:

Arg 1 = path to folder with untarred 10x matrix

Arg2 = “x” [dumby]

Arg3 = gene x guide targets to test, for example see: gene_gRNA_group_pairs_wtc11_encode_for_sceptre.txt 

Arg4 = gRNA guides x guide targets file, for example see: gRNA_groups_table_wtc11_encode.txt

Arg5 = sample_name

Arg6 = all [dumby]

Arg7 = no [dumby]

Arg8 = moi 
either  “high” or “low”

Arg9 side = “both” or “left”

Arg10 = grna_integration_strategy either “union” or "bonferroni"

Arg11 = guide assignment strategy either “thresholding” or mixture

Arg12 = output directory path

Arg13 = are there positive controls? “yes” or “no”


*To run power analysis*

Create enviromnent from analyze_crispr_screen_new.yml
with enviromnent active, open R and install sceptre following the developers instructions
install.packages("devtools")
devtools::install_github("katsevich-lab/sceptre")

Run power analysis:
Rscript for_submit_cleaner_and_w_reps_power_simulation_workflow_EJ_new_w_sceptre_fixed_new_1.23.24_w_no_mast_input.R arg1 arg2

arg 1 is the args file, for example see: inputs_tab_wtc11_16_lanes.txt. 

arg 2 is a file with the list of guides to simluate the power analysis, for example see: guide_file.txt (note, you can split a guide file into many smaller files to run multiple jobs in parallel)




