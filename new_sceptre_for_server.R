#!/usr/bin/env Rscript
#for encode wtc11 tap random screeen
#Evelyn Jagoda

#in this version pooled all reps and lanes together
args = commandArgs(trailingOnly=TRUE)

#devtools::install_github("katsevich-lab/sceptre")

library(dplyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(stringr)
library(Matrix)
library(sceptre)

setwd("/oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq")
mtx_name = "/oak/stanford/groups/engreitz/Users/jray/230925-WTC11-ENCODE-DC-TAP-MOI5-Frozen-8-Lanes/Outputs_cellranger/Aggregate_MOI5_A-H_8lanes"
perturb_status_name = "x"
gene_gRNA_group_pairs_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230901_WTC11_ENCODE_DC_TAP_MOI5_FreshvsFrozen/gene_gRNA_group_pairs_wtc11_encode_for_sceptre.txt"
gRNA_groups_table_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230901_WTC11_ENCODE_DC_TAP_MOI5_FreshvsFrozen/gRNA_groups_table_wtc11_encode.txt"
sample_name = "Aggregate_MOI5_A-H_8lanes_testing_with_new_sceptre"
sample_number = "all"
guides_pre_assigned = "no"
moi = "high"
side = "both"
grna_integration_strategy = "union" #or "bonferroni", start with union but try bonferroni
guide_assignment_method = "mixture" #"thresholding" #"mixture" #that's the best but slow, or "thresholding", default is 5
formula_object_mod = "no_plus_1"
outdir = getwd()

mtx_name = args[1] #moi5_sample_1_filtered_feature_bc_matrix"
perturb_status_name = args[2] #moi5_sample_1_get_perturbation_status_transposed_from_sh.txt fixed_moi5_test_pertubation_matrix_u.txt
gene_gRNA_group_pairs_name = args[3] #gene_gRNA_group_pairs_use.txt
gRNA_groups_table_name = args[4] #"gRNA_groups_table_ej_use.txt"
sample_name = args[5] #"moi5_sample1_only_testing_for_his_threshold_server"
sample_number = args[6] #should be either 1 or all, will do the different options based on that"
guides_pre_assigned = args[7] #should be yes or no, if yes use perturb, if no use sceptre
moi = args[8]
side = args[9]
grna_integration_strategy = args[10]
guide_assignment_method = args[11]
outdir = args[12]
do_pos = args[13]
formula_object = args[14]

setwd(outdir)

sample_name = paste0(sample_name,"_side_",side,"_grna_integration_",grna_integration_strategy,"_guide_assignment_",guide_assignment_method,"_formula_mod_",formula_object_mod)

mtx = Read10X(mtx_name)
gene_gRNA_group_pairs = read.table(gene_gRNA_group_pairs_name,header=T,sep = '\t')

gRNA_groups_table_ej = read.table(gRNA_groups_table_name,header=T,sep = '\t')

gRNA_matrix <- mtx[["CRISPR Guide Capture"]]
gene_matrix <- mtx[["Gene Expression"]]

#need to add to server version
saveRDS(gRNA_matrix,paste0(outdir,sample_name,"_gRNA_matrix.RDS"))
saveRDS(gene_matrix,paste0(outdir,sample_name,"_gene_matrix.RDS"))

batch_tab = data.frame(cbind(colnames(gene_matrix)))

for (i in 1:nrow(batch_tab)){
  batch = str_split(batch_tab[i,1],"-")[[1]][2] ##or - depending
  batch_tab$batch[i] = batch
}
row.names(batch_tab) = batch_tab$cbind.colnames.gene_matrix..
batch_tab$cbind.colnames.gene_matrix.. = NULL

extra_covariates = batch_tab


if (guides_pre_assigned == "yes"){
  gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(perturbation_matrix),]
}

if (guides_pre_assigned == "no"){
  gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(gRNA_matrix),]
}

gene_gRNA_group_pairs1 = gene_gRNA_group_pairs[gene_gRNA_group_pairs$response_id %in% row.names(gene_matrix),]
  
response_names = row.names(gene_matrix)

colnames(gRNA_groups_table_ej_use)[2] = "grna_target"

#takeout extra covariates if it's 1 sample - but then have to take out any zeros ahead of time

sceptre_object <- import_data(
  response_matrix = gene_matrix,
  grna_matrix = gRNA_matrix,
  grna_target_data_frame =gRNA_groups_table_ej_use,
  moi = moi,
  extra_covariates = extra_covariates,
  response_names = response_names
)

if (formula_object_mod == "yes"){
  covariate_tab =  sceptre_object@covariate_data_frame
  bad_cells = c()
  for (i in 1:nrow(covariate_tab)){
    row = covariate_tab[i,]
    if (any(row == 0) | any(is.na(row))){
      bad_cells = c(bad_cells,rownames(row))
    }
  }
  good_cells = setdiff(row.names(covariate_tab),bad_cells)
  gene_matrix = gene_matrix[,good_cells]
  gRNA_matrix = gRNA_matrix[,good_cells]
  extra_covariates = data.frame(extra_covariates[good_cells,])
  rownames(extra_covariates) = good_cells
  colnames(extra_covariates)[1] = "batch"
  gRNA_groups_table_ej_use = gRNA_groups_table_ej[gRNA_groups_table_ej$grna_id %in% row.names(gRNA_matrix_filtered),]
  gene_gRNA_group_pairs1 = gene_gRNA_group_pairs[gene_gRNA_group_pairs$response_id %in% row.names(gene_matrix_filtered),]
  colnames(gRNA_groups_table_ej_use)[2] = "grna_target"
  
  response_names = row.names(gene_matrix_filtered)
  sceptre_object <- import_data(
    response_matrix = gene_matrix,
    grna_matrix = gRNA_matrix,
    grna_target_data_frame =gRNA_groups_table_ej_use,
    moi = moi,
    extra_covariates = extra_covariates,
    response_names = response_names
  )
}
#check for zeros and remove

###new thing to update, make positive controls tab something you premake
#make like this


#need to have this be more comlicated, there's a few that are like TSS_1, TSS_2
colnames(gene_gRNA_group_pairs1)[2] = "grna_target"
gene_gRNA_group_pairs1$type = "x"

#for (i in 1:nrow(gene_gRNA_group_pairs1)){
#  if (gene_gRNA_group_pairs1$target_type[i] == "TSS" & (gene_gRNA_group_pairs1$grna_target[i] == paste0(gene_gRNA_group_pairs1$response_id[i],"_TSS"))){
#   gene_gRNA_group_pairs1$type[i] = "pos"
#  }
#}

for (i in 1:nrow(gene_gRNA_group_pairs1)){
  if (gene_gRNA_group_pairs1$target_type[i] == "TSS" & (grepl(gene_gRNA_group_pairs1$grna_target[i],pattern = paste0(gene_gRNA_group_pairs1$response_id[i],"_TSS")))){
    gene_gRNA_group_pairs1$type[i] = "pos"
  }
}


positive_control_pairs =  gene_gRNA_group_pairs1[gene_gRNA_group_pairs1$type == "pos",c("grna_target", "response_id")]
discovery_pairs = gene_gRNA_group_pairs1[gene_gRNA_group_pairs1$type != "pos",c("grna_target", "response_id")]



if (formula_object_mod == "yes"){
  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = unique(discovery_pairs),
    positive_control_pairs = positive_control_pairs,
    side = side,
    grna_integration_strategy = grna_integration_strategy,
    formula(~ log(response_n_nonzero) + log(response_n_umis) + log(grna_n_nonzero) + log(grna_n_umis) + batch)
  )
} else{
  sceptre_object <- set_analysis_parameters(
    sceptre_object = sceptre_object,
    discovery_pairs = unique(discovery_pairs),
    positive_control_pairs = positive_control_pairs,
    side = side,
    grna_integration_strategy = grna_integration_strategy
  )
}



print(sceptre_object)
#log(response_n_nonzero) + log(response_n_umis) + log(grna_n_nonzero + 1) + log(grna_n_umis + 1) + batch is the default? wanna put it back to what I had


png(paste0(sample_name,"_plot_grna_count_distributions.png"))
p = plot_grna_count_distributions(
  sceptre_object = sceptre_object,
  n_grnas_to_plot = 9)
print(p)
dev.off()

#use the default
sceptre_object_guides_assigned <- assign_grnas(
  sceptre_object = sceptre_object,
  method = guide_assignment_method, parallel = F
)

png(paste0(sample_name,"_plot_grna_cell_assignment_method",guide_assignment_method,".png"))
p = plot(sceptre_object_guides_assigned,n_grnas_to_plot =5)
print(p)
dev.off()

sceptre_object = sceptre_object_guides_assigned

print(sceptre_object)

png(paste0(sample_name,"_plot_covariates.png"))
p = plot_covariates(sceptre_object)
print(p)
dev.off()

sceptre_object = run_qc(sceptre_object = sceptre_object)

png(paste0(sample_name,"_plot_qc.png"))
p = plot(sceptre_object)
print(p)
dev.off()

print(sceptre_object)

sceptre_object <- run_calibration_check(
  sceptre_object = sceptre_object,
  parallel = FALSE #set to T on comp
)

png(paste0(sample_name,"_calibration_check.png"))
p = plot(sceptre_object)
print(p)
dev.off()

calibration_result <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_calibration_check"
)

write.table(calibration_result,paste0(sample_name,"_calibration_result.txt"),quote = F,sep = '\t',row.names = F)

print(sceptre_object)

if (do_pos == "yes"){
  sceptre_object <- run_power_check(
    sceptre_object = sceptre_object,
    parallel = FALSE
  )
  png(paste0(sample_name,"_power_check.png"))
  p = plot(sceptre_object)
  print(p)
  dev.off()
}



sceptre_object <- run_discovery_analysis(
  sceptre_object = sceptre_object,
  parallel = FALSE
)

png(paste0(sample_name,"_discovery_result.png"))
p = plot(sceptre_object)
print(p)
dev.off()

discovery_result <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)
head(discovery_result)

write.table(discovery_result,paste0(sample_name,"_discovery_result.txt"),quote = F,sep = '\t',row.names = F)

print(sceptre_object)

write_outputs_to_directory(sceptre_object = sceptre_object,directory = paste0(getwd(),"/",sample_name,"secptre_output_all"))

saveRDS(sceptre_object,paste0(sample_name,"final_sceptre_object.rds"))
