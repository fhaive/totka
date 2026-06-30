
# BMDx data input
# User input UI
# card BMDx stats input
card_bmd_stats_input <- \(){bslib::card(fill=FALSE,
  bslib::card_header("BMD Stats"),
  bslib::card_body(
    shiny::fileInput(inputId="bmd_stats_file", label="Upload POD file (.xlsx). Expected POD units µM", accept=c(".xlsx", ".xls"))
  ),
  bslib::card_body(
    shiny::uiOutput("bmd_load_str")
  )
)}

# card BMDx KE enrichment input
card_bmd_ke_enrichment_input <- \(){bslib::card(fill=FALSE,
  bslib::card_header("Key Event Enrichment"),
  bslib::card_body(
    shiny::fileInput(inputId="bmd_ke_enrichment_file", label="Upload BMD Key Event Enrichment file (.xlsx)", accept=c(".xlsx", ".xls"))
  ),
  bslib::card_body(
    shiny::uiOutput("bmd_ke_load_str")
  )
)}

# sidebar BMDx input
sidebar_bmdx_input <- \(){bslib::sidebar(width="20%",
  title="Import",
  card_bmd_stats_input(),
  card_bmd_ke_enrichment_input(),
  shiny::actionButton(inputId="upload_bmd_submit", label="Import", width="100%")
)}

# navset BMDx stats
navset_card_bmdx_stats_input <- \(){bslib::navset_card_tab(id="navset_card_bmdx_stats_input",
  bslib::nav_panel("Table",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      DT::DTOutput("bmd_stats_data_dt")
    )
  ),
  bslib::nav_panel("Heatmap",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(heatmap_id="bmdx_stats_heatmap", layout="1-(2|3)", width1=1000, height1=1400, output_ui=NULL)
    )
  )
)}

# navset BMDx KE enrichment
navset_card_bmdx_ke_enrichment_input <- \(){bslib::navset_card_tab(id="navset_card_bmdx_ke_enrichment_input",
  bslib::nav_panel("Table",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      DT::DTOutput("bmd_ke_enrichment_data_dt")
    )
  ),
  bslib::nav_panel("Heatmap",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(heatmap_id="bmdx_ke_heatmap", layout="1-(2|3)", width1=1400, height1=1000, output_ui=NULL)
    )
  ),
  bslib::nav_panel("Distribution",
    bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
      bslib::layout_columns(col_widths=c(2,2,-8),
        shiny::selectInput("score_dist", "Score", choices=c("BMDL", "BMD", "BMDU"), selected="BMDL"),
        shiny::selectInput("level_dist", "Key Event Type", choices=c("Overall"="all", "Per Level"="level"), selected="all")
      ),
      bslib::card(full_screen=TRUE,
        bslib::card_body(
          shiny::plotOutput("bmdx_ke_dist")
        )
      )
    )
  )
)}

# navset BMDx 
navset_card_bmdx_input <- \(){bslib::navset_card_tab(id="navset_card_bmdx_input",
  bslib::nav_panel("BMDx Stats",
    shiny::uiOutput("select_experiment_bmd_stats_data"),
    navset_card_bmdx_stats_input()   
  ),
  bslib::nav_panel("BMDx KE Enrichment",
    shiny::uiOutput("select_experiment_bmd_ke_enrichment_data"),
    navset_card_bmdx_ke_enrichment_input()
  )
)}

# nav_panel BMDx
nav_panel_bmdx_input <- \(){bslib::nav_panel("BMDx Input",
  bslib::card(
    bslib::card_header("BMDx Data"),
    bslib::layout_sidebar(
      sidebar=sidebar_bmdx_input(),
      bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
        navset_card_bmdx_input()
      )
    )
  )
)}

# DE data input
# User input UI
# card DE input
card_deg_input <- \(){bslib::card(fill=FALSE,
  bslib::card_header("Differentially Expressed Genes"),
  bslib::card_body(
    shiny::fileInput(inputId="deg_file", label="Upload Differentially Expressed Genes (.xlsx)", accept=c(".xlsx", ".xls"))
  ),
  bslib::card_body(
    shiny::uiOutput("deg_load_str")
  )
)}

# card DE KE enrichment input
card_deg_ke_enrichment_input <- \(){bslib::card(fill=FALSE,
  bslib::card_header("Key Event Enrichment"),
  bslib::card_body(
    shiny::fileInput(inputId="deg_ke_enrichment_file", label="Upload Differential Analysis Key Event Enrichment file (.xlsx)", accept=c(".xlsx", ".xls"))
  ),
  bslib::card_body(
    shiny::uiOutput("deg_ke_load_str")
  )
)}

# sidebar DE input
sidebar_de_input <- \(){bslib::sidebar(width="20%",
  title="Import",
  card_deg_input(),
  card_deg_ke_enrichment_input(),
  shiny::actionButton(inputId="upload_deg_submit", label="Import", width="100%")
)}

# navset DE
navset_card_deg_input <- \(){bslib::navset_card_tab(id="navset_card_deg_input",
  bslib::nav_panel("Table",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      DT::DTOutput("deg_data_dt")
    )
  ),
  bslib::nav_panel("Heatmap",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(heatmap_id="deg_heatmap", layout="1-(2|3)", width1=1000, height1=1400, output_ui=NULL)
    )
  )
)}

# navset DE KE enrichment
navset_card_de_ke_enrichment_input <- \(){bslib::navset_card_tab(id="navset_card_de_ke_enrichment_input",
  bslib::nav_panel("Table",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      DT::DTOutput("deg_ke_enrichment_data_dt")
    )
  ),
  bslib::nav_panel("Heatmap",
    #bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
    bslib::layout_columns(col_widths=c(12),
      InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(heatmap_id="deg_ke_heatmap", layout="1-(2|3)", width1=1400, height1=1000, output_ui=NULL)
    )
  )
)}

# navset DE
navset_card_de_input <- \(){bslib::navset_card_tab(id="navset_card_de_input",
  bslib::nav_panel("Differential Analysis Stats",
    #shiny::uiOutput("select_experiment_deg_data"),
    navset_card_deg_input()
  ),
  bslib::nav_panel("Differential Analysis KE Enrichment",
    shiny::uiOutput("select_experiment_deg_ke_enrichment_data"),
    navset_card_de_ke_enrichment_input()
  )
)}

# nav_panel DE
nav_panel_de_input <- \(){bslib::nav_panel("Differential Analysis Input",
  bslib::card(
    bslib::card_header("Differential Analysis Data"),
    bslib::layout_sidebar(
      sidebar=sidebar_de_input(),
      bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
        navset_card_de_input()
      )
    )
  )
)}

# nav_panel input
nav_panel_input <- \(){bslib::nav_panel("Input",
  bslib::navset_card_tab(
    nav_panel_bmdx_input(),
    nav_panel_de_input()
  )
)}

# Forward Dosimetry UI
# well_panel models FAQ
well_panel_fw_dosimetry <- \(){shiny::wellPanel(
  shiny::h4("Available Models"),
  shiny::tags$p(shiny::tags$em('All models return concentrations in "µg/L".')),
  shiny::tags$h5("Modelled compartments: Liver & serum"),
  shiny::tags$ul(
    shiny::tags$li("Worley et al. 2017 ", shiny::tags$span("PFOA", class = "label label-primary")),
    shiny::tags$li("Loccisano et al. 2012 ", shiny::tags$span("PFOA & PFOS", class = "label label-info")),
    shiny::tags$li("Fabrega et al. 2016 ", shiny::tags$span("PFOA & PFOS", class = "label label-info")),
    shiny::tags$li("Chou et al. 2019 ", shiny::tags$span("PFOS", class = "label label-success"))
  ),
  shiny::tags$h5("Modelled compartments: serum"),
  shiny::tags$ul(shiny::tags$li("Chiu et al. 2022 ", shiny::tags$span("PFOA & PFOS", class = "label label-info"))),
  shiny::tags$h5("Consensus (Bayesian Ensemble Model, BEM)"),
  shiny::tags$ul(shiny::tags$li("Aggregates the above models where available; output aligns with available compartments (Liver & serum when supported, otherwise serum-only)."))
)}

# sidebar forward dosimetery
sidebar_fw_dosimetry <- \(){bslib::sidebar(width="20%",
  title="Configuration",
  shiny::selectizeInput(inputId="substance", label="Select Chemical(s)", choices=c("PFOS", "PFOA"), selected=c("PFOS", "PFOA"), multiple=TRUE),
  shiny::numericInput("body_weight", "Body Weight (kg)", value = 70, min = 0.1),
  shiny::numericInput("exposure_time", "Exposure Time - Years", value = 40, min = 1, max = 100),
  shiny::numericInput("ingestion_rate", "Ingestion Rate (ng/kg/day)", value = 130.0, min = 1.0, max = 10000.0),
  shiny::actionButton(inputId="set_pbpk_submit", label="Run", width="100%")
)}

# download forward dosimetry
download_fw_dosimetry <- \(){bslib::popover(
  bsicons::bs_icon("download", title="Download", size="2em"),
  bslib::card(
    bslib::card_body(
      shiny::downloadButton("download_fw_dosimetry_single_tables", "Download Single Model Results"),
      shiny::downloadButton("download_fw_dosimetry_consensus_tables", "Download Consensus Results"),
      shiny::downloadButton("download_fw_dosimetry_consensus_summary_tables", "Download Consensus Summary Results")
    )
  ),
  title = "Download"
)}

# nav_panel forward dosimetery
nav_panel_fw_dosimetry <- \(){bslib::nav_panel("Forward Dosimetry",
  bslib::card(
    bslib::card_header("Forward Dosimetry Analysis", download_fw_dosimetry(), class="d-flex justify-content-between"),
    bslib::layout_sidebar(
      sidebar=sidebar_fw_dosimetry(),
      bslib::navset_card_tab(id="navset_card_fw_dosimetry_result",
        bslib::nav_panel("Result Table",
          bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
            bslib::layout_columns(col_widths=c(3,3,3,-3),
              shiny::uiOutput("render_select_fw_dosimetry_scenario"),
              shiny::uiOutput("render_select_fw_dosimetry_chemical"),
              shiny::selectInput("compartment_forward", "Compartment", choices = c("liver", "serum"), selected = "liver"),
            ),
            bslib::navset_card_tab(id="navset_card_fw_dosimetry_result_dt",
              bslib::nav_panel("Single Model",
                DT::DTOutput("fw_dosimetry_simulation_results_single_model_DT")
              ),
              bslib::nav_panel("Consensus Summary",
                DT::DTOutput("fw_dosimetry_simulation_results_consensus_summary_DT")
              )
            )
          )
        ),
        bslib::nav_panel("Result Plot",
          bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
            bslib::layout_columns(col_widths=c(3,3,3,-3),
              shiny::uiOutput("render_selectize_fw_dosimetry_scenario"),
              shiny::uiOutput("render_selectize_fw_dosimetry_chemical"),
              shiny::selectizeInput(inputId="compartment_forward_plot", label="Compartment", choices=c("liver","serum"), selected=c("liver","serum"), multiple=TRUE)
            ),
            bslib::card(full_screen=TRUE,
              bslib::card_body(
                shiny::plotOutput("fw_dosimetry_simulation_results_plot")
              )
            )
          )
        )
      )
    )
  )
)}

# Integrated analysis UI
# sidebar integrated analysis
sidebar_integrated <- \(){bslib::sidebar(width="20%",
  title="Configuration",
  shiny::uiOutput("render_select_fw_dosimetry_scenario_combined"),
  shiny::selectInput(inputId="compartment", label="Compartment", choices=c("liver","serum"), selected="liver"),
  bslib::input_switch(id="filter_keyevents", 
    label=bslib::tooltip(trigger=list("Filter Key Events", bsicons::bs_icon("info-circle", size="1em")), "Select Key Events by Organ/Tissue and Customize."), 
    value=FALSE, 
  width=NULL),
  shiny::conditionalPanel(
    condition = "input.filter_keyevents",
    bslib::card(
      shiny::uiOutput("render_selectize_organ"),
      shiny::uiOutput("render_selectize_organ_ke")
    )
  ),
  shiny::numericInput(inputId="n_permutations", label="Number of Permutations", value=10, min=2, max=1000),
  shiny::uiOutput("render_selectize_time_point"),
  shiny::numericInput(inputId="bmd_summarization_percentile", label="BMD Summarization Percentile", value=0.05, min=0.01, max=1, step=0.01),
  shiny::actionButton(inputId="combine_models_submit", label="Run", width="100%")
)}

# sidebar integrated analysis network
sidebar_network <- \(){bslib::sidebar(width="15%",position="right",
  title="Network Configuration",
  shiny::uiOutput("render_select_exposure_time_point"),
  bslib::input_switch(id="enlarge_keyevents", 
    label=bslib::tooltip(trigger=list("Enlarge Key Events", bsicons::bs_icon("info-circle", size="1em")), "Enlarge the Set Key Events by Selecting Closest Adverse Outcome and Molecular Initiating Events using Key Event to Key Event Relationships."), 
    value=FALSE,
  width=NULL),
  shiny::conditionalPanel(
    condition = "input.enlarge_keyevents",
    bslib::card(
      shiny::numericInput(inputId="max_path_length", label="Maximum Path Length", value=3, min=2, max=100),
      shiny::numericInput(inputId="n_AO", label="Maximum Number of Adverse Outcomes", value=1, min=1, max=100),
      shiny::numericInput(inputId="n_MIE", label="Maximum Number of Molecular Initiating Events", value=3, min=1, max=100)
    )
  ),
  
  shiny::actionButton(inputId="net_submit", label="Run", width="100%")
)}

# download integrated analysis
download_integrated <- \(){bslib::popover(
  bsicons::bs_icon("download", title="Download", size="2em"),
  bslib::card(
    bslib::card_body(
      shiny::downloadButton("download_integrated_tables", "Download Summarized Table"),
      shiny::downloadButton("download_integrated_matrices", "Download Probability Matrices"),
      shiny::downloadButton("download_integrated_visnet_data_all", "Download Networks")
    )
  ),
  title = "Download"
)}

# download integrated network
download_integrated_visnet <- \(){bslib::popover(
  bsicons::bs_icon("download", title="Download", size="2em"),
  bslib::card(
    bslib::card_body(
      shiny::downloadButton("download_integrated_visnet_display", "Download HTML"),
      shiny::downloadButton("download_integrated_visnet_data_display", "Download Data")
    )
  ),
  title = "Download"
)}

# nav_panel integrated analysis
nav_panel_integrated <- \(){bslib::nav_panel("Integrated Analysis",
  bslib::card(
    bslib::card_header("Integrated Analysis", download_integrated(), class="d-flex justify-content-between"),
    bslib::layout_sidebar(
      sidebar=sidebar_integrated(),
      bslib::navset_card_tab(id="navset_card_itegrated_result",
        bslib::nav_panel("Result Table",
          DT::DTOutput("integrated_probability_results_DT")
        ),
        bslib::nav_panel("Result Plot",
          shiny::uiOutput("render_selectize_integrated_ke"),
          bslib::navset_card_tab(id="navset_card_itegrated_result_probability",
            bslib::nav_panel("Activation Probability",
              bslib::card(full_screen=TRUE,
                bslib::card_body(
                  shiny::plotOutput("integrated_probability_results_plot", height="100%")
                )
              )
            ),
            bslib::nav_panel("Activation Bubble Plot",
              bslib::card(full_screen=TRUE,
                bslib::card_body(
                  shiny::plotOutput("integrated_bubble_plot", height="100%")
                )
              )
            ),
            bslib::nav_panel("Key Event PODs and PBK Estimate",
              bslib::card(full_screen=TRUE,
                bslib::card_body(
                  shiny::plotOutput("pbk_pod_plot", height="100%")
                )
              )
            )
          )
        ),
        bslib::nav_panel("Network Plot",
          bslib::layout_sidebar(
            sidebar=sidebar_network(),
            bslib::card(full_screen=TRUE,
              bslib::card_header("Key Event Network", 
                class="d-flex justify-content-between",                  
                download_integrated_visnet()
              ),
              bslib::card_body(fillable=TRUE,
                bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
                  bslib::layout_columns(col_widths=c(3,3,3,-3),
                    shiny::uiOutput("render_select_experiment_net"),
                    shiny::selectInput(inputId="group_visnet_by", label="Group By", choices=c("Adverse Outcomes"="aop","Path to Adverse Outcome of Interest"="rank"), selected="aop"),
                    shiny::conditionalPanel(
                      condition = "input.group_visnet_by == 'rank'",
                      shiny::uiOutput("render_select_target_ke")
                    )
                  ),
                  visNetwork::visNetworkOutput(outputId="ke_visnetwork", height="1000px")
                )
              )
            )
          )
        )
      )
    )
  )
)}

# Reverse Dosimetry UI
# sidebar reverse dosimetry
sidebar_rev_dosimetry <- \(){bslib::sidebar(width="20%",
  title="Configuration",
  shiny::uiOutput("render_selectize_time_point_rev_dosimetry"),
  shiny::selectizeInput(inputId="pfas_reverse_dosimetry", label="Select Chemical(s)",
    choices=c("PFOS","PFOA"), 
    selected=c("PFOS","PFOA"),
  multiple=TRUE),
  # render select input key event for reverse dosimetry",
  bslib::input_switch(id="use_top_keyevents", 
    label=bslib::tooltip(trigger=list("Only Top Key Events", bsicons::bs_icon("info-circle", size="1em")), "Limit this computationally intensive analysis to only top Key Events per Experiment and Time Point by score."), 
    value=FALSE, 
  width=NULL),
  shiny::conditionalPanel(
    condition = "input.use_top_keyevents",
    bslib::card(
      shiny::numericInput(inputId="n_key_events", label="Top (x) Key Events", value=10, min=1, max=100),
      shiny::selectInput(inputId="filter_type_reverse_dosimetry", label="By", choices=c("Count"="n", "Percentage"="perc"), selected="n"),
      shiny::selectInput(inputId="score_type_reverse_dosimetry", label="Score Type", choices=c("BMD"="BMD", "BMDL"="BMDL", "Adjusted P-Value"="padj"), selected="BMD")
    )
  ),
  shiny::numericInput(inputId="body_weight_reverse_dosimetry", label="Body Weight (kg)", value=70, min=0.1),
  shiny::selectizeInput(inputId="compartment_reverse_dosimetry", label="Compartment", choices=c("liver","serum"), selected=c("liver","serum"), multiple=TRUE),
  shiny::numericInput(inputId="exposure_duration_reverse_dosimetry", label="Exposure Time - Years", value=40, min=1, max=100),
  shiny::selectInput(inputId="pod_reverse_dosimetry", label="POD Variable", choices=c("BMD", "BMDL", "BMDU"), selected="BMD"),
  bslib::input_switch(id="is_stochastic", 
    label=bslib::tooltip(trigger=list("Stochastic", bsicons::bs_icon("info-circle", size="1em")), "involving chance or probability : probabilistic"), 
    value=FALSE, 
  width=NULL),
  shiny::conditionalPanel(
    condition = "input.is_stochastic",
    bslib::card(
      shiny::numericInput(inputId="n_samples", label="N Samples for Stochastic", value=10, min=2),
    )
  ),
  bslib::input_switch(id="is_parallel", 
    label=bslib::tooltip(trigger=list("Parallel Processing", bsicons::bs_icon("info-circle", size="1em")), "Use available cpu threads for multithreadded processing."), 
    value=FALSE, 
  width=NULL),
  shiny::conditionalPanel(
    condition = "input.is_parallel",
    shiny::uiOutput("render_select_n_cores")
  ),
  shiny::actionButton(inputId="run_reverse_dosimetry", label="Run", width="100%")
)}

# sidebar reverse dosimetry network
sidebar_rev_dosimetry_network <- \(){bslib::sidebar(width="15%",position="right",
  title="Network Configuration",
  bslib::input_switch(id="enlarge_rev_dosimetry_keyevents", 
    label=bslib::tooltip(trigger=list("Enlarge Key Events", bsicons::bs_icon("info-circle", size="1em")), "Enlarge the Set Key Events by Selecting Closest Adverse Outcome and Molecular Initiating Events using Key Event to Key Event Relationships."), 
    value=FALSE,
  width=NULL),
  shiny::conditionalPanel(
    condition = "input.enlarge_rev_dosimetry_keyevents",
    bslib::card(
      shiny::numericInput(inputId="max_path_length_rev_dosimetry", label="Maximum Path Length", value=3, min=2, max=100),
      shiny::numericInput(inputId="n_AO_rev_dosimetry", label="Maximum Number of Adverse Outcomes", value=1, min=1, max=100),
      shiny::numericInput(inputId="n_MIE_rev_dosimetry", label="Maximum Number of Molecular Initiating Events", value=3, min=1, max=100)
    )
  ),
  shiny::actionButton(inputId="net_submit_rev_dosimetry", label="Run", width="100%")
)}

# download reverse dosimetry
download_rev_dosimetry <- \(){bslib::popover(
  bsicons::bs_icon("download", title="Download", size="2em"),
  bslib::card(
    bslib::card_body(
      shiny::downloadButton("download_rev_dosimetry_tables", "Download Results"),
      shiny::downloadButton("download_rev_dosimetry_summarized_tables", "Download Summarized Results"),
      shiny::downloadButton("download_rev_dosimetry_visnet_data_all", "Download Network")
    )
  ),
  title = "Download"
)}

# download reverse dosimetry network
download_rev_dosimetry_visnet <- \(){bslib::popover(
  bsicons::bs_icon("download", title="Download", size="2em"),
  bslib::card(
    bslib::card_body(
      shiny::downloadButton("download_rev_dosimetry_visnet_display", "Download HTML"),
      shiny::downloadButton("download_rev_dosimetry_visnet_data_display", "Download Data")
    )
  ),
  title = "Download"
)}

# nav_panel reverse dosimetry
nav_panel_rev_dosimetry <- \(){bslib::nav_panel("Reverse Dosimetry",
  bslib::card(
    bslib::card_header("Reverse Dosimetry", download_rev_dosimetry(), class="d-flex justify-content-between"),
    bslib::layout_sidebar(
      sidebar=sidebar_rev_dosimetry(),
      bslib::navset_card_tab(id="navset_card_rev_dosimetry_result",
        bslib::nav_panel("Result Table",
          DT::DTOutput("reverse_dosimetry_results_DT")#,
        ),
        bslib::nav_panel("Density Plot",
          bslib::card(full_screen=TRUE,
            bslib::card_body(
              shiny::plotOutput("reverse_dosimetry_ridgeplot", height="100%")
            )
          )
        ),
        bslib::nav_panel("Heatmap Plot",
          InteractiveComplexHeatmap::InteractiveComplexHeatmapOutput(heatmap_id="reverse_dosimetry_heatmap", layout="1-(2|3)", width1=1400, height1=1000, output_ui=NULL)
        ),
        bslib::nav_panel("Network Plot",
          bslib::layout_sidebar(
            sidebar=sidebar_rev_dosimetry_network(),
            bslib::card(full_screen=TRUE,
              bslib::card_header("Key Event Network", 
                class="d-flex justify-content-between",                  
                download_rev_dosimetry_visnet()
              ),
              bslib::card_body(fillable=TRUE,
                bslib::layout_columns(col_widths=c(12,12), row_heights=c(1,11),
                  bslib::layout_columns(col_widths=c(2,1,-9),
                    shiny::uiOutput("render_select_experiment_reverse_dosimetry_net")
                  ),
                  visNetwork::visNetworkOutput(outputId="reverse_dosimetry_ke_visnetwork", height="1000px")
                )
              )
            )
          )
        )
      )
    )
  )
)}


#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import bslib
#' @importFrom plotly renderPlotly
#' @importFrom DT DTOutput
#' @importFrom waiter useWaiter
#' @importFrom InteractiveComplexHeatmap InteractiveComplexHeatmapOutput
#' @importFrom visNetwork visNetworkOutput
#' @importFrom bsicons bs_icon
#' @noRd
app_ui <- function(request) {
  bslib::page_navbar(title="TOTKA - Toxicogenomics & Toxicokinetics mapping to AOP",
    theme=app_theme(),
    navbar_options = bslib::navbar_options(
      bg = "#5200ff",
      underline=TRUE
    ),
    header={golem_add_external_resources()
    waiter::useWaiter()},
    nav_panel_input(),
    nav_panel_fw_dosimetry(),
    nav_panel_integrated(),
    nav_panel_rev_dosimetry(),
    bslib::nav_spacer(),
    #bslib::nav_item(
    #  shiny::actionLink("options", "Options", icon=shiny::icon("sliders"))
    #),
    bslib::nav_item(
      shiny::actionLink("reset", "Reset", icon=shiny::icon("redo"))
    ),
    bslib::nav_menu(title="Links", align="right",
      #bslib::nav_panel("Help",
      #  "Help/FAQ"
      #),
      bslib::nav_item(align="left",
        a("FHAIVE GitHub", href="https://github.com/fhaive", target="_blank"),
        a("FHAIVE Docker Hub", href="https://hub.docker.com/u/fhaive", target="_blank")
      )
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @import bslib
#' @import InteractiveComplexHeatmap
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @importFrom plotly renderPlotly
#' @importFrom DT DTOutput
#' @importFrom waiter useWaiter
#' @noRd
golem_add_external_resources <- function(tmpDIR=tempdir()) {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  add_resource_path(
    "tmpDIR",
    tmpDIR
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "TOTKA"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
