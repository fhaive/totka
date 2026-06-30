### LOAD LIBRARIES
library(TOTKA)
library(ggplot2)
library(PfasDosim)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(htmlwidgets)
library(visNetwork)
library(openxlsx)

### CONFIGURE PANDOC LOCATION
# Use this if running the script in the terminal. To be configured to where the pandoc folder is
if (!nzchar(Sys.getenv("RSTUDIO"))) {
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
}

### SET SEED
set.seed(123)

### CONFIGURE OUTPUT PATH
path_to_folder_output = "../totka_output/"
experiment_path = "high_130/"

### CREATE LOG
sink(file = paste(path_to_folder_output, experiment_path,"log.txt",sep = ""),split = T)


#### LOAD REQUIRED INPUTS
### LOAD BMDx analysis results
optimal_models_stats =  readxl::read_excel(path = "sample_data/liver_spheroids/01_BMD_results.xlsx")
optimal_models_stats = as.data.frame(optimal_models_stats)

### LOAD AOPfingeprint KE enrichment starting from BMDx results
Enrichment_data = bmdx::read_excel_allsheets(filename = "sample_data/liver_spheroids/02_BMD_KE_annotated.xlsx",tibble = F,first_col_as_rownames = F,check_numeric = F)
Enrichment_data = rbind(Enrichment_data$PFOS, Enrichment_data$PFOA)
Enrichment_data$BMD = as.numeric(Enrichment_data$BMD)
Enrichment_data$BMDL = as.numeric(Enrichment_data$BMDL)
Enrichment_data$BMDU = as.numeric(Enrichment_data$BMDU)
Enrichment_data$time = as.numeric(Enrichment_data$time)
unique_exp = unique(Enrichment_data$Experiment)

### SET PBK parameters
body_weight = 70
exposure_time = 40
substance = "BOTH"
predefined_ingestion_rate = c("Vestergreen High PFOS (130 ngkgday)" = 130)

### RUN PBK FORWARD DOSIMETRY
out130 <- PfasDosim::forward_dosimetry(chemical = substance,
                                       ingestion = as.numeric(predefined_ingestion_rate),
                                       BW = body_weight,
                                       duration = exposure_time,
                                       stat_sum = F)
### SAVE PBK RESULTS
save(out130, file = paste(path_to_folder_output,experiment_path,"forward_dosimetry_high_exposure_scenario_130_ngkgday.RData", sep = ""))

## PLOT PBK SIMULATION
p = TOTKA::plot_PBK_simulations(out130,
                                scenario_labels = "Vestergreen High PFOS (130 ngkgday)",
                                chemicals = c("PFOS", "PFOA"),
                                y_label = "Concentrations Liver (µg/L)")

ggsave(p, file = paste(path_to_folder_output,experiment_path,"forward_dosimetry_plot.png",sep = ""),width = 10, height = 7)
ggsave(p, file = paste(path_to_folder_output,experiment_path,"forward_dosimetry_plot.pdf",sep = ""),width = 10, height = 7)

### PLOT PBK AND OMICS CONTEXTUALIZAITON
p = TOTKA::plot_PBK_POD(Enrichment_data, 
                         out130, 
                         exposure_time=40, 
                         exposure_level="X ng/kg/day",
                         y_label = "Concentration log10(µg/L)",
                         x_label = "Time (day)",
                         log_y_scale = T,
                         free_y_param = F,
                         theme_element_size = 20)
ggsave(p, file=paste(path_to_folder_output,experiment_path,"PBK_omics_contextualization.pdf",sep = ""),width = 14, height = 6)


### FORMAT OUTPUT FOR DOWNSTREAM ANALYSES
Simulation_results = TOTKA::compiling_forward_dosimetry_simulation_dataframes_all_samples(out130)

### SET INTEGRATIVE ANALYSIS PARAMETERS

# Number of permutations for BMD random selection in BMDL - BMDU range
NPerm = 100
# Percentile of the BMD distribution across all the BMD of the genes that belong to a KE
BMD_summarization_percentile = 0.05
tissues_of_interest = c("liver","serum")
exposure_time = 40
## Variable that will be samples. It must be BMD
pod_variable = "BMD"
### Upper limit of sampling, 
upper_limit = "BMDU"

### RUN INTEGRATIVE ANALYSIS 
probability_matrices = TOTKA::compute_PBK_BMD_probability_matrix_activation(tissues_of_interest = tissues_of_interest,
                                                                            unique_exp = unique_exp,
                                                                            optimal_models_stats = optimal_models_stats,
                                                                            Enrichment_data = Enrichment_data,
                                                                            NPerm =NPerm,
                                                                            BMD_summarization_percentile = BMD_summarization_percentile,
                                                                            exposure_time =  exposure_time,
                                                                            pod_variable =pod_variable,
                                                                            Simulation_results = Simulation_results,
                                                                            bmd_upper_level = upper_limit)  


probability_matrices_liver = probability_matrices$liver
probability_matrices_serum = probability_matrices$serum

### SAVE INTEGRATIVE ANALYSIS 
# results when the BMD is sampled  in BMDL - BMDU range
save(probability_matrices_liver, file = paste(path_to_folder_output,experiment_path,"probability_matrices_liver_Vestergreen_High (130 ngkgday) 40 years BMDU upper.RData",sep=""))
save(probability_matrices_serum, file = paste(path_to_folder_output,experiment_path,"probability_matrices_serum_Vestergreen_High (130 ngkgday) 40 years BMDU upper.RData",sep=""))


load(file = paste(path_to_folder_output,experiment_path,"probability_matrices_liver_Vestergreen_High (130 ngkgday) 40 years BMDU upper.RData",sep=""))
load(file = paste(path_to_folder_output,experiment_path,"probability_matrices_serum_Vestergreen_High (130 ngkgday) 40 years BMDU upper.RData",sep=""))

## PLOT PROBABILITY MAPPING FOR ALL KEs IN LIVER
df_tissue_specific = TOTKA::convert_tensor_to_dataframes(probability_matrices = probability_matrices_liver,
                                                         time_points = c(1,4,10,14), 
                                                         time_unit = "Day")

#### WRITE PROBABILITIES TO EXCEL
write.xlsx(
  df_tissue_specific,
  file = paste(path_to_folder_output,experiment_path,"KE_probability_liver_compartment_130ngkgday.xlsx",sep = "")
)

gp = TOTKA::plot_KE_probabilities_over_PBK_time(
  df_tissue_specific,
  level_order = c("Molecular", "Cellular", "Tissue", "Organ", "Individual"),
  pfoa_label = "",
  pfos_label = " ",
  both_label = "   ",
  fill_colours = c("#BDBDBD","darkblue", "lightblue", "lightgreen", "yellow", "pink", "red"),
  prob_values = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 1),
  prob_breaks = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 1),
  prob_labels = c("0","0.1", "0.2", "0.4", "0.6", "0.8", "1"),
  x_lab = "PBK Years Simulation",
  y_lab = "Key Events",
  labs_fill_text = "BMD",
  guide_fill_title = "Probability of activation",
  y_wrap_n = 10,
  axis_text_x_angle = 45,
  axis_text_x_hjust = 1,
  axis_text_x_size = 12,
  axis_text_y_size = 14,
  axis_title_size  = 14,
  strip_text_size  = 14,
  legend_title_size = 14,
  legend_text_size  = 14,
  strip_text_y_angle = 0,
  strip_text_y_hjust = 0,
  panel_spacing_lines = 1
)

ggsave(gp, file = paste(path_to_folder_output,experiment_path,"probability_over_KE_vs_liver_compartment_predictions 40 years.pdf",sep = ""),width = 25,height = 30)


### PLOT PROBABILITY MAPPING FOR ALL KEs IN SERUM
df_tissue_specific = TOTKA::convert_tensor_to_dataframes(probability_matrices_serum,
                                                         time_points = c(1,4,10,14), 
                                                         time_unit = "Day")
#### WRITE PROBABILITIES TO EXCEL
write.xlsx(
  df_tissue_specific,
  file = paste(path_to_folder_output,experiment_path,"KE_probability_serum_compartment_130ngkgday.xlsx",sep = "")
)

gp =  TOTKA:::plot_KE_probabilities_over_PBK_time(
  df_tissue_specific,
  level_order = c("Molecular", "Cellular", "Tissue", "Organ", "Individual"),
  pfoa_label = "",
  pfos_label = " ",
  both_label = "   ",
  fill_colours = c("#BDBDBD","darkblue", "lightblue", "lightgreen", "yellow", "pink", "red"),
  prob_values = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 1),
  prob_breaks = c(0,0.1, 0.2, 0.4, 0.6, 0.8, 1),
  prob_labels = c("0","0.1", "0.2", "0.4", "0.6", "0.8", "1"),
  x_lab = "PBK Years Simulation",
  y_lab = "Key Events",
  labs_fill_text = "BMD",
  guide_fill_title = "Probability of activation",
  y_wrap_n = 10,
  axis_text_x_angle = 45,
  axis_text_x_hjust = 1,
  axis_text_x_size = 12,
  axis_text_y_size = 14,
  axis_title_size  = 14,
  strip_text_size  = 14,
  legend_title_size = 14,
  legend_text_size  = 14,
  strip_text_y_angle = 0,
  strip_text_y_hjust = 0,
  panel_spacing_lines = 1
)

ggsave(gp, file = paste(path_to_folder_output,experiment_path,"probability_over_KE_vs_serum_compartment_predictions 40 years.pdf",sep = ""),width = 25,height = 30)


#### PERFORM NETWORK ANALYSIS

### LOAD KE ENRICHMENT OVER DD GENES - should be run by the APP
Enrichment_data_BMD = bmdx::read_excel_allsheets(filename = "sample_data/liver_spheroids/02_BMD_KE_annotated.xlsx",
                                                 tibble = F,
                                                 first_col_as_rownames = F,
                                                 check_numeric = F)
Enrichment_data_BMD = rbind(Enrichment_data_BMD$PFOS, Enrichment_data_BMD$PFOA)

### LOAD KE ENRICHMENT OVER DEGs GENES - should be run by the APP
Enrichment_data_DEG = bmdx::read_excel_allsheets(filename = "sample_data/liver_spheroids/04_DEG_KE_annotated.xlsx",
                                                 tibble = F,
                                                 first_col_as_rownames = F,
                                                 check_numeric = F)
Enrichment_data_DEG = Enrichment_data_DEG$Sheet1
Enrichment_data_DEG$Day = as.numeric(Enrichment_data_DEG$Timepoint)

#### 

load(paste(path_to_folder_output, experiment_path,"probability_matrices_liver_Vestergreen_High (130 ngkgday) 40 years BMDU upper.RData",sep = ""))

probability_matrices = probability_matrices_liver

KE_annotated = TOTKA:::get_docking_ke_mapping()

net_list_res = TOTKA::build_network_lists(
  probability_matrices,
  Enrichment_data_DEG,
  Enrichment_data_BMD,
  KE_annotated,
  experiments = names(probability_matrices),
  PBK_time = 40,
  max_path_length = 3,
  n_AOs = 1,
  n_MIEs = 3,
  mode = "out",
  enlarge_ke_selection = T,
  ke_id = "ke",
  numerical_variables = c("prob"),
  pval_variable = "padj",
  gene_variable = "Genes",
  convert_to_gene_symbols = FALSE,
  group_by = "rank",
  target_node = "Increased, Liver Steatosis",
  normalize = NULL
)

# List of networks to save
nets_to_save <- list(
  PFOA_1  = net_list_res$visnet_list$PFOA_1,
  PFOA_4  = net_list_res$visnet_list$PFOA_4,
  PFOA_10 = net_list_res$visnet_list$PFOA_10,
  PFOA_14 = net_list_res$visnet_list$PFOA_14,
  PFOS_1  = net_list_res$visnet_list$PFOS_1,
  PFOS_4  = net_list_res$visnet_list$PFOS_4,
  PFOS_10 = net_list_res$visnet_list$PFOS_10,
  PFOS_14 = net_list_res$visnet_list$PFOS_14
)



save(net_list_res, nets_to_save, file = paste(path_to_folder_output,experiment_path,"KEKE_prob_analysis_liver.RData",sep=""))

# Save all networks
for (nm in names(nets_to_save)) {

  w <- nets_to_save[[nm]]
  
  w <- htmlwidgets::onRender(
    w,
    '
  function(el) {
    el.style.width = "100%";
    el.style.height = "100vh";
  }
  '
  )
  saveWidget(
    w,
    paste0(path_to_folder_output, experiment_path, "html_files_liver/",nm, "_network.html"),
    selfcontained = TRUE
  )
  
  ######
  
  nodes= net_list_res$net_list[[nm]]$nodes
  edges= net_list_res$net_list[[nm]]$edges
  
  nodes <- nodes[
    !is.na(nodes$rank) &
      grepl("(^|,\\s*)(1|2|3)(,|$)", nodes$rank),
  ]
  
  
  edges <- edges[
    edges$from %in% nodes$id &
      edges$to   %in% nodes$id,
  ]
  
  
  vn <- get_visnet(nodes=nodes, edges=edges, group_by="aop", col_pal = net_list_res$pal_node_fill_list[[nm]])
  
  vn <- htmlwidgets::onRender(
    vn,
    '
  function(el) {
    el.style.width = "100%";
    el.style.height = "100vh";
  }
  '
  )
  saveWidget(
    vn,
    paste0(path_to_folder_output, experiment_path, "html_files_liver/",nm, "_network_top_3_paths.html"),
    selfcontained = TRUE
  )
  
}


sink()
