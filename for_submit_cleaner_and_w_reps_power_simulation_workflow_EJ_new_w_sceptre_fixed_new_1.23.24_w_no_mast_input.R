#to set things up, start the flow like simulate_diff_expr which is the function you actually call in the power_simulations.R function

#tmp_logfile <- tools::file_path_sans_ext(snakemake@log[[1]])
#log <- file(tmp_logfile, open = "w")
#sink(log)
#sink(log, type = "message")


#To do for k562:
#1. missing a few genes, still need to fix that to get all (for now pull those from mast if look similar overall)
#2. parallel-ize this with a wrapper script, since it goes 1 perturb at a time, just can divide it up
#3. add reps, since getting the pert object is the thing that's slow and because we're not sampling cells, 
#   just get the object and repeat the other steps x times, default is 20
#statistical power for each pair was defined by the proportion of simulations in which the injected perturbation effect was detected as statistically significant.
#before you go ham, do some comparison


#in /oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example
args = commandArgs(trailingOnly=TRUE)

#pwr_sim_name = args[1] #"pwr_sim.rda"
#load(pwr_sim_name)
#åsceptre_object_name = args[1]
#full_grna_matrix_name =args[2]
#full_response_matrix_name = args[3]
#sce_object_name = args[4]
#n_pert = args[5]
#outprefix = args[6]
#effect = args[7]
#reps = args[8]
#center = args[9]
#guide_file_name = args[10]


#wtc11

#sceptre_object_name = "/oak/stanford/groups/engreitz/Users/ejagoda/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes_new_sceptre_side_both_grna_integration_union_guide_assignment_thresholdingfinal_sceptre_object.rds"
#full_grna_matrix_name = "/oak/stanford/groups/engreitz/Users/ejagoda/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes_new_sceptre_side_both_grna_integration_union_guide_assignment_thresholding_gRNA_matrix.RDS"
#full_response_matrix_name = "/oak/stanford/groups/engreitz/Users/ejagoda/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes_new_sceptre_side_both_grna_integration_union_guide_assignment_thresholding_gene_matrix.RDS"
#sce_object_name = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/results/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/sce_object_from_sceptre_workflow_w_disp.rds"
#n_pert = "all"
#outprefix =  "231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes"
#effect = 0.25
#reps = 3
#center = "FALSE"
#guide_file_name = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/sceptre_power_analysis_logs/temp_dir/wtc11_16_lanes_power/guides00"

#args_file_name = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/resources/231114-WTC11-ENCODE-DC-TAP-MOI5-Frozen-16-Lanes/inputs_tab_wtc11_16_lanes.txt"
#guide_file_name = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/sceptre_power_analysis_logs/temp_dir/wtc11_16_lanes_power/guides00"


args_file_name = args[1]
guide_file_name = args[2]






#to parallelize, just need to have the discovery file be a subset
#also add a step to get the inputs as one clean thing

#for the full thing
#pwr_sim_name = "pwr_sim.rda"
#sceptre_object_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq/testing_new_sceptre/all_moi5_w_fresh_aggregate_no_norm_from_batch_8.29.31_w_new_covar_encode_filters_side_both_grna_integration_union_guide_assignment_thresholdingfinal_sceptre_object.rds"
#full_grna_matrix_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq/testing_new_sceptre/all_moi5_w_fresh_aggregate_no_norm_from_batch_8.29.31_w_new_covar_encode_filters_side_both_grna_integration_union_guide_assignment_thresholding_gRNA_matrix.RDS"
#sce_object_name = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/resources/230327_Encode_Tap_MOI5_AGGREGATED_w_freshlane_5k_controls_w_batches_w_n_gudies_w_lane_7.14.23/perturb_sce_disp.rds"
#n_pert = "all"
#outprefix = "testing_submission"
#effect = 0.25
#over_ride_sample_name = "yes"
#rep = "2"
#full_response_matrix_name = "/oak/stanford/groups/engreitz/Users/ejagoda/230327_Encode_K562_Tap_seq_full_seq/all_moi5_w_fresh_aggregate_no_norm_from_batch_8.29.31_w_new_covar_new_sceptre_side_both_grna_integration_union_guide_assignment_thresholding_gene_matrix.RDS"
#center = FALSE
#reps = 5
#guide_file_name = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/sceptre_power_analysis_logs/temp_dir/guides07"

#load(pwr_sim_name)

snakemake_scriptdir = "/oak/stanford/groups/engreitz/Users/ejagoda/Encode_K562_Random_TAP_screen/CRISPRiScreen-example/workflow/scripts/"

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(sceptre)
  library(stringr)
  #library(MatrixExtra)
  source(file.path(snakemake_scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake_scriptdir, "R_functions/power_simulations_fun.R"))
})


args_file = read.table(args_file_name,header=T)

for (i in 1:nrow(args_file)){
  assign(x = args_file$variable[i],value = args_file$value[i])
}



#it's gotta work if you pull from the response matrix instead of the sce object, try this
sceptre_object = readRDS(sceptre_object_name)
full_grna_matrix = readRDS(full_grna_matrix_name)
full_response_matrix = readRDS(full_response_matrix_name)


outdir = paste0("results/", outprefix,"/sceptre_power_sim/")
outname = paste0(outprefix,"_sceptre_power_sim_effect_size_",effect,"_guide_subset_",basename(guide_file_name), "_n_perts_",n_pert,"_rep_",reps,"_center_effect_size_",center,".txt")

if (!file.exists(paste0("results/", outprefix,"/"))){
  dir.create(file.path(paste0("results/", outprefix,"/")))
}

if (!file.exists(outdir)){
  dir.create(file.path(outdir))
}

message(paste0("output directory: ",outdir,"\nfinal output file: ",outname,"\ntemporary output file: ",paste0("tempfile_",outname)))


# prepare data =====================================================================================

# load prepared input data stored in SingleCellExperiment object
message("Loading input data.")
sce <- readRDS(sce_object_name)

# infer perturbation level based on strategy
pert_level <- switch("perCRE", "perGRNA" = "grna_perts", "perCRE" = "cre_perts",
                     stop("incorrect strategy argument"))

# filter for minimum number of cells per perturbation

# perform power simulations ========================================================================
message("Not normalizing transcript counts because this is sceptre.")


# convert 'percentage decrease' effect size to 'relative expression level'
effect_size <- 1 - as.numeric(effect)

# simulate Perturb-seq data and perform differential gene expression tests
message("Performing power simulations.")

#Step 0, do all this (pulling from simluate_diff_expr) get inputs to get effect size matrix

#set up
effect_size = effect_size
pert_level = pert_level
genes_iter = FALSE
guide_sd = as.numeric(0.13)
#center = FALSE
#rep = 1
norm = "real"


col_names <- colnames(colData(sce)) 
if ("pert" %in% col_names) stop("'pert' cannot be a colData name, please rename.", call. = FALSE)

# check that pert_level is valid
if (!pert_level %in% altExpNames(sce)) {
  stop("pert_level must be one of the altExp names: ", paste(altExpNames(sce), collapse = ", "),
       call. = FALSE)
}

# check that guide_sd is valid
if (!is.numeric(guide_sd) | guide_sd < 0) {
  stop("Invalid 'guide_sd' value. Must be numeric >= 0.", call. = FALSE)
}

# get required functions -------------------------------------------------------------------------
#trying without sampling so that we can just use straight


n_ctrl = FALSE
# get function to generate input data for one perturbation
if (is.numeric(n_ctrl)) {
  pert_input_function <- pert_input_sampled
  n_ctrl <- as.integer(n_ctrl)
} else if (n_ctrl == FALSE) {
  pert_input_function <- pert_input
} else {
  stop("Invalid 'n_ctrl' argument.", call. = FALSE)
}


message("Following sceptre, not sampling control cells")

if (norm == "real") {
  sim_function <- simulate_diff_expr_pert_real
} else if (norm == "sim_nonpert") {
  sim_function <- simulate_diff_expr_pert_sim
} else {
  stop("Invalid 'norm' argument.", call. = FALSE)
}


#
#function to get the simulate count matrix
sim_counts_submit <- function(sce, effect_size_mat) {
  
  gene_means = rowData(sce)[, "mean"]
  cell_size_factors = colData(sce)[, "size_factors"]
  gene_ids = names(sce)
  cell_ids = names(cell_size_factors)
  # simulate Perturb-seq count data with parameters from SCE object
  sim_counts <- simulate_tapseq_counts(gene_means = rowData(sce)[, "mean"],
                                       gene_dispersions = rowData(sce)[, "dispersion"],
                                       cell_size_factors = colData(sce)[, "size_factors"],
                                       effect_size_mat = effect_size_mat, gene_ids = gene_ids,cell_ids = cell_ids)
  row.names(sim_counts) = gene_ids

  return(sim_counts)
  
}

# perform simulated DE tests ---------------------------------------------------------------------

#TO DO NOTE: for sceptre, want to make these control cells and stuff exactly the same - come back to that
#get the discovery_pairs to test
discovery_relevant_pairs_all = sceptre_object@discovery_pairs_with_info

guide_file = read.table(guide_file_name,header = F,sep = '\t')

discovery_relevant_pairs_in_guide_subset = discovery_relevant_pairs_all[discovery_relevant_pairs_all$grna_group %in% guide_file$V1,]

guide_targets_sce <- create_guide_targets(sce, pert_level = pert_level)

discovery_relevant_pairs = discovery_relevant_pairs_in_guide_subset[discovery_relevant_pairs_in_guide_subset$response_id %in% rownames(counts(sce)) & discovery_relevant_pairs_in_guide_subset$grna_group %in% rownames(altExp(sce,pert_level)),]

perts = unique(discovery_relevant_pairs$grna_group)


guide_targets = sceptre_object@grna_target_data_frame
colnames(guide_targets)[2] = "target_id"

discovery_results = data.frame()

if (n_pert == "all"){
  n_pert = length(perts)
}


cat(paste0(c("response_id","grna_target","n_nonzero_trt","n_nonzero_cntrl","pass_qc","p_value","log_2_fold_change","significant","rep\n"),collapse = "\t"),
    file = paste0(outdir,"tempfile_",outname))

for (pert in perts[1:n_pert]){
  discovery_relevant_pairs_pert = discovery_relevant_pairs[discovery_relevant_pairs$grna_group == pert & discovery_relevant_pairs$pass_qc == "TRUE",]
  
  if (nrow(discovery_relevant_pairs_pert) > 0 ){
    #pert_object <- pert_input_function(pert, sce = sce, pert_level = pert_level,
    #                                   cell_batches = cell_batches, n_ctrl = n_ctrl)
    
    #TO DO: should make sure the guide assignments match the sceptre guide assignments
    #can do it like this sceptre_object@grna_assignments$grna_group_idxs
    #this gets the cell indexes sceptre_object@grna_assignments$grna_group_idxs[pert]
    
    pert_guides <- guide_targets[guide_targets$target_id == pert, "grna_id"]
    
    #what if i did that via sceptre like
    #pert_guides_sceptre = sceptre_object@grna_target_data_frame[sceptre_object@grna_target_data_frame$grna_group == pert,"grna_id"]
    # sceptre_object@grna_assignments$grna_group_idxs[pert][[1]] this doesn't work, but the mast ones look okay, gotta figure this out
    #this gets it but haven't figure out how to map the indexes to the right cells  
    
    pert_object = pert_input_function(pert,sce = sce, pert_level = pert_level)
    
    # get perturbation status and gRNA perturbations for all cells
    pert_status <- colData(pert_object)$pert
    
    grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
    
    #2.14.24 edit
    pert_guides_use =  pert_guides[pert_guides %in% rownames(grna_perts)]
    
    pert_guides = pert_guides_use
    
    grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts,
                                                 pert_guides = pert_guides)
    
    pert_genes =discovery_relevant_pairs_pert$response_id
    
    pert_object =  pert_object[pert_genes,]
    
    effect_sizes <- structure(rep(1, nrow(pert_object)), names = rownames(pert_object))
    
    effect_sizes[pert_genes] <- effect_size #does this only for the genes you care about will need to adjust to match sceptre
    
    #for testing 
    pert_cells = rownames(colData(pert_object)[colData(pert_object)$pert == 1,])
    
    
    
    for (i in 1:as.numeric(reps)){
      
    message("creating simulation for rep,",i)
    #okay the problem is with this function, fix that, fix the probablm!! when you get this working it works
    # create effect size matrix 
    es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                        gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
    
    #es_mat[, pert_status == 0] <- 1
    #es_mat[!rownames(es_mat) %in% pert_genes, ] <- 1
    
    # center effect sizes on specified gene-level effect sizes
    if (center == TRUE) {
      es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status,
                                          gene_effect_sizes = effect_sizes)
    }
    
    ###okay figured it out, the es_matrix is ordered with the pert cells in the front, and then the others
    #the order is off, need to fix it
    
    
    es_mat_use = es_mat[,colnames(counts(pert_object))]
    
    #sim_object <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat_use) #don't normalize for sceptre presumably, now you've got a simulated_object
    
  
  #so now we have the simulated object for the pertrubation with every gene it should be tested with (might need to change with sceptre flow))
  sim_counts = sim_counts_submit(pert_object, effect_size_mat = es_mat_use)
  
    
  #trying this
    
    #sim_counts = as.matrix(sim_counts_submit(pert_object, effect_size_mat = es_mat_use))
  if (length(pert_genes) > 1){
    full_response_matrix_sim = full_response_matrix[pert_genes,]
    shared_cells = intersect(colnames(sim_counts),colnames(full_response_matrix))
    full_response_matrix_sim[rownames(sim_counts),shared_cells] = sim_counts[,shared_cells]
  }else{
    sim_counts_mx = as.matrix(sim_counts)
    #need to add another row to keep the structure
    if (pert_genes[1] != rownames(full_response_matrix)[1]){
      alt_gene = rownames(full_response_matrix)[1]
    }else{
      alt_gene = rownames(full_response_matrix)[2]
    }
        
    full_response_matrix_sim = full_response_matrix[c(pert_genes,alt_gene),]
    shared_cells = intersect(colnames(sim_counts_mx),colnames(full_response_matrix))
    full_response_matrix_sim[pert_genes,shared_cells] = sim_counts_mx[pert_genes,shared_cells]
    
  }

  full_response_matrix_sim_sparse = as(as.matrix(full_response_matrix_sim),"RsparseMatrix")
  sceptre_object_use = sceptre_object
  sceptre_object_use@response_matrix =  full_response_matrix_sim_sparse
  sceptre_object_use@discovery_pairs_with_info = discovery_relevant_pairs_pert
  
  ###
  
  
  sceptre_object2 <- run_discovery_analysis(
    sceptre_object = sceptre_object_use,
    parallel = FALSE
  )
  #okay now check, should have ~25% effect ideally
  
  discovery_result <- get_result(
    sceptre_object = sceptre_object2,
    analysis = "run_discovery_analysis"
  )
  discovery_result$rep = i
  write.table(discovery_result, col.names=F, append = T, row.names = F,quote = F,sep = '\t',file = paste0(outdir,"tempfile_",outname))
  discovery_results = data.frame(rbind(discovery_results,discovery_result))
    }
  }
}

message("Processing output.")
disp_outlier <- data.frame(gene = rownames(rowData(sce)),
                           disp_outlier_deseq2 = rowData(sce)[, "disp_outlier_deseq2"],
                           stringsAsFactors = FALSE)

av_expr = data.frame(rowMeans(counts(sce)))
colnames(av_expr)[1] = "average_expression_all_cells"
av_expr$response_id = rownames(av_expr)

colnames(disp_outlier)[1] = "response_id"


# add to output
discovery_results <- left_join(discovery_results, disp_outlier, by = "response_id")
discovery_results <- left_join(discovery_results, av_expr, by = "response_id")

# save simulation output
message("Saving output to file.")
#write_tsv(discovery_results, file = paste0("testing_sceptre_power_analysis",snakemake@wildcards$sample,snakemake@wildcards$rep.tsv))
  
write_tsv(discovery_results, file = paste0(outdir,outname))
