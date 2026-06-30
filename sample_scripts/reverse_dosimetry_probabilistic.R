### LOAD LIBRARIES
library(TOTKA)
library(ggplot2)
library(PfasDosim)
library(dplyr)
library(tidyr)
library(readxl)
library(openxlsx)
library(ggnewscale)
library(ggrepel)
library(stringr)


### CONFIGURE PANDOC LOCATION
# Use this if running the script in the terminal. To be configured to where the pandoc folder is
if (!nzchar(Sys.getenv("RSTUDIO"))) {
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
}

### SET SEED
set.seed(123)

### CONFIGURE OUTPUT PATH
path_to_folder_output = "../totka_output/"
experiment_path = "reverse_dosimetry_probabilistic/"
sink(file = paste(path_to_folder_output, experiment_path,"log.txt",sep = ""),split = T)


PFAS = c("PFOS","PFOA")
compartment = "liver"
n_cores = 8

library(PfasDosim)
library(TOTKA)

Enrichment_data_BMD = bmdx::read_excel_allsheets(filename = "sample_data/liver_spheroids/02_BMD_KE_annotated.xlsx",tibble = F,first_col_as_rownames = F,check_numeric = F)
Enrichment_data_BMD = rbind(Enrichment_data_BMD$PFOS, Enrichment_data_BMD$PFOA)

results_reverse_dosimetry = reverse_dosimetry_wrapping(Enrichment_data_BMD,
                                                       duration = 40,
                                                       BW = 70,
                                                       PFAS_mw = c("PFOA"=414.07,"PFOS"=500.13),
                                                       pod = "BMDL",
                                                       PFAS = PFAS,
                                                       time = c(1,4,10,14),
                                                       compartments = compartment,
                                                       isStochastic = TRUE,
                                                       isParallel = T, 
                                                       Nsamples = 1000,
                                                       n_cores = n_cores,
                                                       free = FALSE)


reverse_dosimetry_list = results_reverse_dosimetry$reverse_dosimetry_list
reverse_dosimetry_list_summarized = results_reverse_dosimetry$reverse_dosimetry_list_summarized

TOTKA:::plot_reverse_dosimetry_heatmap(reverse_dosimetry_list_summarized = reverse_dosimetry_list_summarized)

save(reverse_dosimetry_list, file = paste(path_to_folder_output,experiment_path,"reverse_dosimetry_probabilistic_list_BMDL_",compartment,".RData",sep =""))
save(reverse_dosimetry_list_summarized, file = paste(path_to_folder_output,experiment_path,"reverse_dosimetry_probabilistic_list_summarized_BMDL_",compartment,".RData",sep =""))
