
## Sanitizing errors to hide error details
#options(shiny.sanitize.errors = TRUE)
#
## shiny reactiveValues
#r <- shiny::reactiveValues()
#r$file_based <- FALSE
#r$tempFiles <- c()

#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyvalidate
#' @import waiter
#' @import PfasDosim
#' @import bmdx
#' @import plotly
#' @import ggridges
#' @import ggplot2
#' @import patchwork
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import visNetwork
#' @importFrom readxl read_excel
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic
  options(warn=-1)
  suppressPackageStartupMessages(library(InteractiveComplexHeatmap))

  # Sanitizing errors to hide error details
  options(shiny.sanitize.errors = TRUE)

  # shiny reactiveValues
  r <- shiny::reactiveValues()
  r$file_based <- FALSE
  r$tempFiles <- c()

  # input validation for forward dosimetry analysis
  # input validation object
  iv.fw_dosimetry <- shinyvalidate::InputValidator$new()

  # input validation rules
  iv.fw_dosimetry$add_rule("body_weight", shinyvalidate::sv_gte(0.1))
  iv.fw_dosimetry$add_rule("exposure_time", shinyvalidate::sv_between(1, 100))
  iv.fw_dosimetry$add_rule("ingestion_rate", shinyvalidate::sv_between(1.0, 10000.0))

  # displaying errors in the UI
  iv.fw_dosimetry$enable()

  # input validation for integrated analysis
  # input validation object
  iv.integrated <- shinyvalidate::InputValidator$new()

  # input validation rules
  iv.integrated$add_rule("n_permutations", shinyvalidate::sv_between(2, 1000))
  iv.integrated$add_rule("bmd_summarization_percentile", shinyvalidate::sv_between(0.01, 1))

  # displaying errors in the UI
  iv.integrated$enable()

  # input validation for integrated analysis key event network
  shiny::observe({
    # input validation object
    iv.integrated.net <- shinyvalidate::InputValidator$new()

    # input validation rules
    iv.integrated.net$add_rule("select_exposure_time_point", shinyvalidate::sv_required())

    if(input$enlarge_keyevents){
      iv.integrated.net$add_rule("max_path_length", shinyvalidate::sv_between(2, 100))
      iv.integrated.net$add_rule("n_AO", shinyvalidate::sv_between(1, 100))
      iv.integrated.net$add_rule("n_MIE", shinyvalidate::sv_between(1, 100))
    }

    # displaying errors in the UI
    iv.integrated.net$enable()
    r$iv.integrated.net <- iv.integrated.net
  }) |> shiny::bindEvent(input$enlarge_keyevents)
  
  # input validation for reverse dosimetry analysis
  shiny::observe({
    # input validation object
    iv.rev_dosimetry <- shinyvalidate::InputValidator$new()

    # input validation rules
    iv.rev_dosimetry$add_rule("body_weight_reverse_dosimetry", shinyvalidate::sv_gte(0.1))
    iv.rev_dosimetry$add_rule("exposure_duration_reverse_dosimetry", shinyvalidate::sv_between(1, 100))

    if(input$use_top_keyevents){
      iv.rev_dosimetry$add_rule("n_key_events", shinyvalidate::sv_between(1, 100))
    }

    if(input$is_stochastic){
      iv.rev_dosimetry$add_rule("n_samples", shinyvalidate::sv_gte(2))
    }

    # displaying errors in the UI
    iv.rev_dosimetry$enable()
    r$iv.rev_dosimetry <- iv.rev_dosimetry
  }) |> shiny::bindEvent(c(input$use_top_keyevents, input$is_stochastic))

  # input validation for reverse dosimetry analysis key event network
  shiny::observe({
    # input validation object
    iv.rev_dosimetry.net <- shinyvalidate::InputValidator$new()

    # input validation rules
    if(input$enlarge_rev_dosimetry_keyevents){
      iv.rev_dosimetry.net$add_rule("max_path_length_rev_dosimetry", shinyvalidate::sv_between(2, 100))
      iv.rev_dosimetry.net$add_rule("n_AO_rev_dosimetry", shinyvalidate::sv_between(1, 100))
      iv.rev_dosimetry.net$add_rule("n_MIE_rev_dosimetry", shinyvalidate::sv_between(1, 100))
    }

    # displaying errors in the UI
    iv.rev_dosimetry.net$enable()
    r$iv.rev_dosimetry.net <- iv.rev_dosimetry.net
  }) |> shiny::bindEvent(input$enlarge_rev_dosimetry_keyevents)

  r$iv.rev_dosimetry.net.is_valid <- shiny::reactive({
    validity_check <- NULL
    if(input$enlarge_rev_dosimetry_keyevents){
     validity_check <- r$iv.rev_dosimetry.net$is_valid() 
    } else{
     validity_check <- TRUE 
    }
    return(validity_check)
  })

  # loading screen transparent~ish background
  loaderGif <- shiny::tagList(img(src="www/fhaive_seamless_smooth_gears.gif"))
  w <- waiter::Waiter$new(
    #id = "plot",
    #html = waiter::spin_3(),
    html = loaderGif,
    color = waiter::transparent(.5)
  )
 
  # Preview BMD file
  output$bmd_preview_table <- renderTable({
    shiny::validate(need(!is.null(input$bmd_stats_file$datapath), "No BMD POD file!"))
    shiny::req(input$bmd_stats_file$datapath)
    bmd_stats_data_filename <- input$bmd_stats_file$datapath
    DF <- readxl::read_excel(path=bmd_stats_data_filename) |> as.data.frame() |> head(5)
    DF
  }, spacing="s", width="50px")

  r$bmd_colnames <- shiny::reactive({
    shiny::req(input$bmd_stats_file$datapath)
    bmd_stats_data_filename <- input$bmd_stats_file$datapath
    DF <- readxl::read_excel(path=bmd_stats_data_filename) |> as.data.frame() |> head(5)
    cNames <- DF |> colnames()
    cNames
  })

  # Import BMD dialog
  observe({
    shiny::req(input$bmd_stats_file$datapath, r$bmd_colnames())
    choices <- r$bmd_colnames()
    bmd_stats_data_filename <- input$bmd_stats_file$datapath
    DF <- readxl::read_excel(path=bmd_stats_data_filename) |> as.data.frame()
    column_types <- c("Experiment"="Experiment|Chem","Time"="Time","Feature"="Feature|Gene","BMDL"="BMDL","BMD"="BMD","BMDU"="BMDU","Adverse_Direction"="Adverse|Direction")
    column_sel <- column_types |> purrr::map(\(val){get_best_match(choices, val)})
    print("str(column_sel)")
    print(str(column_sel))
    shiny::showModal(myModalDialog(size="full",
      title = "Import POD File",
      bslib::layout_columns(col_widths=c(12), row_heights=c(3),
        bslib::card(
          bslib::card_header("BMD PODs Preview"),
          bslib::card_body(
            shiny::tableOutput("bmd_preview_table")
          )
        )
      ),
      shiny::tags$hr(),
      bslib::layout_columns(col_widths=c(12,12,12,12,12), row_heights=c(1,1,1,1,1),
        bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_experiment_col", 
            label=bslib::tooltip(trigger=list("Experiment Column", bsicons::bs_icon("info-circle", size="1em")), "Experiment information, same as chemical, e.g. PFOS"), 
            choices=choices, selected=column_sel$Experiment),
          shiny::selectInput("bmd_time_col", 
            label=bslib::tooltip(trigger=list("Time Column", bsicons::bs_icon("info-circle", size="1em")), "Time information, e.g. 10"), 
            choices=choices, selected=column_sel$Time),
          shiny::selectInput("bmd_feature_col", 
            label=bslib::tooltip(trigger=list("Feature Column", bsicons::bs_icon("info-circle", size="1em")), "Feature information, like gene name, e.g. CNN1"), 
            choices=choices, selected=column_sel$Feature),
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_bmdl_col", 
            label=bslib::tooltip(trigger=list("BMDL Column", bsicons::bs_icon("info-circle", size="1em")), "Benchmark Dose lower confidence limit, e.g. 0.6813788"), 
            choices=choices, selected=column_sel$BMDL),
          shiny::selectInput("bmd_bmd_col", 
            label=bslib::tooltip(trigger=list("BMD Column", bsicons::bs_icon("info-circle", size="1em")), "Benchmark Dose, e.g. 0.9759008"), 
            choices=choices, selected=column_sel$BMD),
          shiny::selectInput("bmd_bmdu_col", 
            label=bslib::tooltip(trigger=list("BMDU Column", bsicons::bs_icon("info-circle", size="1em")), "Benchmark Dose upper confidence limit, e.g. 1.339539"), 
            choices=choices, selected=column_sel$BMDU)
        ),bslib::layout_columns(col_widths=c(4,-8),
          shiny::selectInput("bmd_adverse_direction_col", 
            label=bslib::tooltip(trigger=list("Adverse Direction Column", bsicons::bs_icon("info-circle", size="1em")), "Adverse direction in BMD results indicates whether the overall dose–response trend shows an increase or a decrease in the measured effect. It is determined by fitting a model across all dose levels and evaluating the sign of the slope (positive or negative), e.g. -1"), 
            choices=choices, selected=column_sel$Adverse_Direction)
        ),
        shiny::tags$hr(),
        bslib::layout_columns(col_widths=c(-4,-4,4),
          shiny::actionButton(inputId="import_bmd_submit", label="Import", width="100%"),
          shiny::uiOutput("bmd_load_str_modal")
        )
      )
    ))
  }) |> shiny::bindEvent(input$bmd_stats_file$datapath)

  # Import BMD submit
  observeEvent(input$import_bmd_submit, {
    shiny::validate(
      need(!is.null(input$bmd_stats_file$datapath), "No BMD POD file!"),
    )
    shiny::req(
      input$bmd_stats_file$datapath, 
      input$bmd_experiment_col,
      input$bmd_time_col,
      input$bmd_feature_col,
      input$bmd_bmdl_col,
      input$bmd_bmd_col,
      input$bmd_bmdu_col,
      input$bmd_adverse_direction_col
    )
    bmd_stats_data_filename <- input$bmd_stats_file$datapath
    cNames <- c(input$bmd_experiment_col, input$bmd_time_col, input$bmd_feature_col, input$bmd_bmdl_col, input$bmd_bmd_col, input$bmd_bmdu_col, input$bmd_adverse_direction_col) 
    dup_check <- cNames |> duplicated() |> sum()
    if(dup_check > 0){
      shiny::showNotification("Duplicated column names. Please check and update.", duration=NULL)
    }
    shiny::req(dup_check==0)
    r$bmd_import_file <- bmd_stats_data_filename
    shiny::removeModal()
    bslib::nav_select(id="navset_card_bmdx_input", selected="BMDx Stats Table", session=session)
  })

  # bmd import status
  output$bmd_load_str <- renderUI({
    if(is.null(r$bmd_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  output$bmd_load_str_modal <- renderUI({
    if(is.null(r$bmd_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })
  
  # Preview BMD KE file
  output$bmd_ke_preview_table <- renderTable({
    shiny::validate(need(!is.null(input$bmd_ke_enrichment_file$datapath), "No BMD Key Event Enrichment file!"))
    shiny::req(input$bmd_ke_enrichment_file$datapath)
    bmd_ke_enrichment_filename <- input$bmd_ke_enrichment_file$datapath
    DF <- readxl::read_excel(path=bmd_ke_enrichment_filename) |> as.data.frame() |> head(2)
    DF
  }, spacing="s", width="50px")

  r$bmd_ke_colnames <- shiny::reactive({
    shiny::req(input$bmd_ke_enrichment_file$datapath)
    bmd_ke_enrichment_filename <- input$bmd_ke_enrichment_file$datapath
    DF <- readxl::read_excel(path=bmd_ke_enrichment_filename) |> as.data.frame() |> head(5)
    cNames <- DF |> colnames()
    cNames
  })

  # Import BMD KE dialog
  observe({
    shiny::req(input$bmd_ke_enrichment_file$datapath, r$bmd_ke_colnames())
    choices <- r$bmd_ke_colnames()
    bmd_ke_enrichment_filename <- input$bmd_ke_enrichment_file$datapath
    DF <- readxl::read_excel(path=bmd_ke_enrichment_filename) |> as.data.frame()
    column_types <- c("Experiment"="Experiment","Time"="Time","PFAS"="PFAS|Chem","Gene"="Gene|Feature","Gene_Count"="Gene|Feature|Count|Num","padj"="padj|adjusted|p.value|p-value|p value","Organ_Tissue"="Organ|Tissue","BMDL"="BMDL","BMD"="BMD","BMDU"="BMDU","AOP_ID"="AOP|ID","KE_ID"="Key Event|Key_Event|KE|ID","KE_Name"="Key Event|Key_Event|KE|Name","KE_Description"="Key Event|Key_Event|KE|Desc","KE_Type"="Key Event|Key_Event|KE|Type","KE_Level"="Key Event|Key_Event|KE|Level")
    column_sel <- column_types |> purrr::map(\(val){get_best_match(choices, val)})
    shiny::showModal(myModalDialog(size="full",
      title = "Import Key Event Enrichment File",
      bslib::layout_columns(col_widths=c(12), row_heights=c(3),
        bslib::card(
          bslib::card_header("BMD Key Event Enrichment Preview"),
          bslib::card_body(
            shiny::tableOutput("bmd_ke_preview_table")
          )
        )
      ),
      shiny::tags$hr(),
      bslib::layout_columns(col_widths=c(12,12,12,12,12,12,12,12), row_heights=c(1,1,1,1,1,1,1,1),
        bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_ke_experiment_col", 
            label=bslib::tooltip(trigger=list("Experiment Column", bsicons::bs_icon("info-circle", size="1em")), "Experiment and time information as <exp>_<time>, e.g. PFOS_10"), 
            choices=choices, selected=column_sel$Experiment),
          shiny::selectInput("bmd_ke_time_col", 
            label=bslib::tooltip(trigger=list("Time Column", bsicons::bs_icon("info-circle", size="1em")), "Time information, e.g. 10"), 
            choices=choices, selected=column_sel$Time),
          shiny::selectInput("bmd_ke_pfas_col", 
            label=bslib::tooltip(trigger=list("PFAS Column", bsicons::bs_icon("info-circle", size="1em")), "PFAS chemical information, e.g. PFOS"), 
            choices=choices, selected=column_sel$PFAS)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_ke_bmdl_col", 
            label=bslib::tooltip(trigger=list("BMDL Column", bsicons::bs_icon("info-circle", size="1em")), "Benchmark Dose lower confidence limit, e.g. 0.6813788"), 
            choices=choices, selected=column_sel$BMDL),
          shiny::selectInput("bmd_ke_bmd_col", 
            label=bslib::tooltip(trigger=list("BMD Column", bsicons::bs_icon("info-circle", size="1em")), "Benchmark Dose, e.g. 0.9759008"), 
            choices=choices, selected=column_sel$BMD),
          shiny::selectInput("bmd_ke_bmdu_col", 
            label=bslib::tooltip(trigger=list("BMDU Column", bsicons::bs_icon("info-circle", size="1em")), "Benchmark Dose upper confidence limit, e.g. 1.339539"), 
            choices=choices, selected=column_sel$BMDU)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_ke_genes_col", 
            label=bslib::tooltip(trigger=list("Genes Column", bsicons::bs_icon("info-circle", size="1em")), "Genes information, semicolon separated, e.g. ALDH2;ETFA;IDO1;PHGDH;POR;STEAP4"), 
            choices=choices, selected=column_sel$Gene),
          shiny::selectInput("bmd_ke_relevant_genes_col", 
            label=bslib::tooltip(trigger=list("Genes Count Column", bsicons::bs_icon("info-circle", size="1em")), "Number of genes, e.g. 6"), 
            choices=choices, selected=column_sel$Gene_Count),
          shiny::selectInput("bmd_ke_padj_col", 
            label=bslib::tooltip(trigger=list("Adjusted P-value Column", bsicons::bs_icon("info-circle", size="1em")), "Adjusted P-value information, e.g. 0.01999476"), 
            choices=choices, selected=column_sel$padj)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_ke_organ_tissue_col", 
            label=bslib::tooltip(trigger=list("Organ Tissue Column", bsicons::bs_icon("info-circle", size="1em")), "Organ or tissue information, e.g. liver"), 
            choices=choices, selected=column_sel$Organ_Tissue),
          shiny::selectInput("bmd_ke_aop_col", 
            label=bslib::tooltip(trigger=list("AOP ID Column", bsicons::bs_icon("info-circle", size="1em")), "Adverse Outcome Pathway ID, e.g. Aop:437"), 
            choices=choices, selected=column_sel$AOP_ID),
          shiny::selectInput("bmd_ke_termid_col", 
            label=bslib::tooltip(trigger=list("Key Event ID Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event ID, e.g. Event:105"), 
            choices=choices, selected=column_sel$KE_ID)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("bmd_ke_key_event_name_col", 
            label=bslib::tooltip(trigger=list("Key Event Name Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event name, e.g. Inhibition, Mitochondrial Electron Transport Chain Complexes"), 
            choices=choices, selected=column_sel$KE_Name),
          shiny::selectInput("bmd_ke_key_event_description_col", 
            label=bslib::tooltip(trigger=list("Key Event Description Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event description, e.g. Inhibition, Mitochondrial Electron Transport Chain Complexes"), 
            choices=choices, selected=column_sel$KE_Description),
          shiny::selectInput("bmd_ke_key_event_type_col", 
            label=bslib::tooltip(trigger=list("Key Event Type Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event type, e.g. MolecularInitiatingEvent"), 
            choices=choices, selected=column_sel$KE_Type)
        ),bslib::layout_columns(col_widths=c(4,-8),
          shiny::selectInput("bmd_ke_level_col", 
            label=bslib::tooltip(trigger=list("Key Event Level Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event level, e.g. Molecular"), 
            choices=choices, selected=column_sel$KE_Level)
        ),
        shiny::tags$hr(),
        bslib::layout_columns(col_widths=c(-4,-4,4),
          shiny::actionButton(inputId="import_bmd_ke_submit", label="Import", width="100%"),
          shiny::uiOutput("bmd_ke_load_str_modal")
        )
      )
    ))
  }) |> shiny::bindEvent(input$bmd_ke_enrichment_file$datapath)

  # Import BMD KE submit
  observeEvent(input$import_bmd_ke_submit, {
    shiny::validate(
      need(!is.null(input$bmd_ke_enrichment_file$datapath), "No BMD POD file!"),
    )
    shiny::req(
      input$bmd_ke_enrichment_file$datapath, 
      input$bmd_ke_termid_col,
      input$bmd_ke_experiment_col,
      input$bmd_ke_time_col,
      input$bmd_ke_pfas_col,
      input$bmd_ke_bmdl_col,
      input$bmd_ke_bmd_col,
      input$bmd_ke_bmdu_col,
      input$bmd_ke_genes_col,
      input$bmd_ke_relevant_genes_col,
      input$bmd_ke_padj_col,
      input$bmd_ke_aop_col,
      input$bmd_ke_key_event_name_col,
      input$bmd_ke_key_event_description_col,
      input$bmd_ke_key_event_type_col,
      input$bmd_ke_level_col,
      input$bmd_ke_organ_tissue_col
    )
    bmd_ke_enrichment_filename <- input$bmd_ke_enrichment_file$datapath
    cNames <- c(input$bmd_ke_termid_col, input$bmd_ke_experiment_col, input$bmd_ke_time_col, input$bmd_ke_pfas_col, input$bmd_ke_bmdl_col, input$bmd_ke_bmd_col, input$bmd_ke_bmdu_col, input$bmd_ke_genes_col, input$bmd_ke_relevant_genes_col, input$bmd_ke_padj_col, input$bmd_ke_aop_col, input$bmd_ke_key_event_name_col, input$bmd_ke_key_event_description_col, input$bmd_ke_key_event_type_col, input$bmd_ke_level_col, input$bmd_ke_organ_tissue_col) 
    dup_check <- cNames |> duplicated() |> sum()
    if(dup_check > 0){
      shiny::showNotification("Duplicated column names. Please check and update.", duration=NULL)
    }
    shiny::req(dup_check==0)
    r$bmd_ke_import_file <- bmd_ke_enrichment_filename
    shiny::removeModal()
    bslib::nav_select(id="navset_card_bmdx_input", selected="BMDx Stats Table", session=session)
  })

  # bmd ke import status
  output$bmd_ke_load_str <- renderUI({
    if(is.null(r$bmd_ke_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  output$bmd_ke_load_str_modal <- renderUI({
    if(is.null(r$bmd_ke_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  # Upload BMD PODs and BMD KE
  observeEvent(input$upload_bmd_submit, {
    shiny::validate(
      need(!is.null(r$bmd_import_file), "No BMD POD file!"),
      need(!is.null(r$bmd_ke_import_file), "No BMD Key Event enrichment file!")
    )
    shiny::req(r$bmd_import_file, r$bmd_ke_import_file)

    w$show()
    on.exit({
      w$hide()
    })

    tryCatch({
      # read BMDx stats data from file
      bmd_stats_data_filename <- r$bmd_import_file
      bmd_stats_data <- readxl::read_excel(path=bmd_stats_data_filename) |> as.data.frame()
      #bmd_stats_data <- bmd_stats_data |> dplyr::rename(Experiment={input$bmd_experiment_col}, time={input$bmd_time_col}, Feature={input$bmd_feature_col}, BMDL={input$bmd_bmdl_col}, BMD={input$bmd_bmd_col}, BMDU={input$bmd_bmdu_col}, Adverse_direction={input$bmd_adverse_direction_col}) 
      bmd_stats_data <- bmd_stats_data |> dplyr::select(Experiment={input$bmd_experiment_col}, time={input$bmd_time_col}, Feature={input$bmd_feature_col}, BMDL={input$bmd_bmdl_col}, BMD={input$bmd_bmd_col}, BMDU={input$bmd_bmdu_col}, Adverse_direction={input$bmd_adverse_direction_col}) 
      bmd_stats_data <- bmd_stats_data |> dplyr::mutate(dplyr::across(c(time, BMDL, BMD, BMDU, Adverse_direction), as.numeric)) 
      r$bmd_stats_data <- bmd_stats_data

      # read BMDx KE enrichment data from file
      bmd_ke_enrichment_data_filename <- r$bmd_ke_import_file
      bmd_ke_enrichment_data <- bmdx::read_excel_allsheets(filename=bmd_ke_enrichment_data_filename, tibble=FALSE, first_col_as_rownames=FALSE, check_numeric=FALSE)
      #bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::rename(TermID=input$bmd_ke_termid_col, Experiment={input$bmd_ke_experiment_col}, time={input$bmd_ke_time_col}, PFAS={input$bmd_ke_pfas_col}, BMDL={input$bmd_ke_bmdl_col}, BMD={input$bmd_ke_bmd_col}, BMDU={input$bmd_ke_bmdu_col}, Genes={input$bmd_ke_genes_col}, relevantGenesInGeneSet={input$bmd_ke_relevant_genes_col}, padj={input$bmd_ke_padj_col}, Aop={input$bmd_ke_aop_col}, key_event_name={input$bmd_ke_key_event_name_col}, Ke_description={input$bmd_ke_key_event_description_col}, Ke_type={input$bmd_ke_key_event_type_col}, level={input$bmd_ke_level_col}, organ_tissue={input$bmd_ke_organ_tissue_col})})
      bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::select(TermID=input$bmd_ke_termid_col, Experiment={input$bmd_ke_experiment_col}, time={input$bmd_ke_time_col}, PFAS={input$bmd_ke_pfas_col}, BMDL={input$bmd_ke_bmdl_col}, BMD={input$bmd_ke_bmd_col}, BMDU={input$bmd_ke_bmdu_col}, Genes={input$bmd_ke_genes_col}, relevantGenesInGeneSet={input$bmd_ke_relevant_genes_col}, padj={input$bmd_ke_padj_col}, Aop={input$bmd_ke_aop_col}, key_event_name={input$bmd_ke_key_event_name_col}, Ke_description={input$bmd_ke_key_event_description_col}, Ke_type={input$bmd_ke_key_event_type_col}, level={input$bmd_ke_level_col}, organ_tissue={input$bmd_ke_organ_tissue_col})})
      bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::mutate(dplyr::across(c(time, BMDL, BMD, BMDU, relevantGenesInGeneSet, padj), as.numeric))}) 
      bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
      bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> dplyr::mutate(chem=Experiment |> stringr::str_replace(pattern="_[:alnum:]*$", replacement=""), time=Experiment |> stringr::str_replace(pattern=".*_", replacement="") |> as.numeric()) |> dplyr::group_by(chem) |> dplyr::group_split() |> as.list()
      nm <- bmd_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::pull(chem) |> unique()}) |> purrr::flatten_chr()
      names(bmd_ke_enrichment_data) <- nm
      r$bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> purrr::map(as.data.frame)
    }, error=\(e){
      str <- paste0("Encountered Error: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })
    bslib::nav_select(id="navset_card_bmdx_input", selected="BMDx Stats Table", session=session)
  })

  # Preview DEG file
  output$deg_preview_table <- renderTable({
    shiny::validate(need(!is.null(input$deg_file$datapath), "No DEG file!"))
    shiny::req(input$deg_file$datapath)
    deg_data_filename <- input$deg_file$datapath
    DF <- readxl::read_excel(path=deg_data_filename) |> as.data.frame() |> head(5)
    DF
  }, spacing="s", width="50px")

  r$deg_colnames <- shiny::reactive({
    shiny::req(input$deg_file$datapath)
    deg_data_filename <- input$deg_file$datapath
    DF <- readxl::read_excel(path=deg_data_filename) |> as.data.frame() |> head(5)
    cNames <- DF |> colnames()
    cNames
  })

  # Import DEG dialog
  observe({
    shiny::req(input$deg_file$datapath, r$deg_colnames())
    choices <- r$deg_colnames()
    deg_data_filename <- input$deg_file$datapath
    DF <- readxl::read_excel(path=deg_data_filename) |> as.data.frame()
    column_types <- c("Experiment"="Experiment|Chem","Time"="Time","Feature"="Feature|Gene","Dose"="Dose|Dosage","Log2FC"="Log2FC|Fold_Change|FoldChange|Fold Change")
    column_sel <- column_types |> purrr::map(\(val){get_best_match(choices, val)})
    shiny::showModal(myModalDialog(size="full",
      title = "Import DEG File",
      bslib::layout_columns(col_widths=c(12), row_heights=c(3),
        bslib::card(
          bslib::card_header("DEG Preview"),
          bslib::card_body(
            shiny::tableOutput("deg_preview_table")
          )
        )
      ),
      shiny::tags$hr(),
      bslib::layout_columns(col_widths=c(12,12,12,12), row_heights=c(1,1,1,1),
        bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("deg_experiment_col", 
            label=bslib::tooltip(trigger=list("Experiment Column", bsicons::bs_icon("info-circle", size="1em")), "Experiment information, same as chemical, e.g. PFOS"), 
            choices=choices, selected=column_sel$Experiment),
          shiny::selectInput("deg_time_col", 
            label=bslib::tooltip(trigger=list("Time Column", bsicons::bs_icon("info-circle", size="1em")), "Time information, e.g. 10"), 
            choices=choices, selected=column_sel$Time),
          shiny::selectInput("deg_feature_col", 
            label=bslib::tooltip(trigger=list("Feature Column", bsicons::bs_icon("info-circle", size="1em")), "Feature information, like gene name, e.g. CNN1"), 
            choices=choices, selected=column_sel$Feature),
        ),bslib::layout_columns(col_widths=c(4,4,-4),
          shiny::selectInput("deg_dose_col", 
            label=bslib::tooltip(trigger=list("Dose Column", bsicons::bs_icon("info-circle", size="1em")), "Dose information, e.g. 10uM"), 
            choices=choices, selected=column_sel$Dose),
          shiny::selectInput("deg_log2fc_col", 
            label=bslib::tooltip(trigger=list("Log2FC Column", bsicons::bs_icon("info-circle", size="1em")), "Log2 fold change information, e.g. 0.614899"), 
            choices=choices, selected=column_sel$Log2FC)
        ),
        shiny::tags$hr(),
        bslib::layout_columns(col_widths=c(-4,-4,4),
          shiny::actionButton(inputId="import_deg_submit", label="Import", width="100%"),
          shiny::uiOutput("deg_load_str_modal")
        )
      )
    ))
  }) |> shiny::bindEvent(input$deg_file$datapath)

  # Import DEG submit
  observeEvent(input$import_deg_submit, {
    shiny::validate(
      need(!is.null(input$deg_file$datapath), "No DEG file!"),
    )
    shiny::req(
      input$deg_file$datapath, 
      input$deg_experiment_col,
      input$deg_time_col,
      input$deg_feature_col,
      input$deg_dose_col,
      input$deg_log2fc_col
    )
    deg_data_filename <- input$deg_file$datapath
    cNames <- c(input$deg_experiment_col, input$deg_time_col, input$deg_feature_col, input$deg_dose_col, input$deg_log2fc_col) 
    dup_check <- cNames |> duplicated() |> sum()
    if(dup_check > 0){
      shiny::showNotification("Duplicated column names. Please check and update.", duration=NULL)
    }
    shiny::req(dup_check==0)
    r$deg_import_file <- deg_data_filename
    shiny::removeModal()
    bslib::nav_select(id="navset_card_de_input", selected="Differential Analysis Stats", session=session)
  })

  # deg import status
  output$deg_load_str <- renderUI({
    if(is.null(r$deg_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  output$deg_load_str_modal <- renderUI({
    if(is.null(r$deg_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  # Preview DEG KE file
  output$deg_ke_preview_table <- renderTable({
    shiny::validate(need(!is.null(input$deg_ke_enrichment_file$datapath), "No DEG Key Event Enrichment file!"))
    shiny::req(input$deg_ke_enrichment_file$datapath)
    deg_ke_enrichment_filename <- input$deg_ke_enrichment_file$datapath
    DF <- readxl::read_excel(path=deg_ke_enrichment_filename) |> as.data.frame() |> head(2)
    DF
  }, spacing="s", width="50px")

  r$deg_ke_colnames <- shiny::reactive({
    shiny::req(input$deg_ke_enrichment_file$datapath)
    deg_ke_enrichment_filename <- input$deg_ke_enrichment_file$datapath
    DF <- readxl::read_excel(path=deg_ke_enrichment_filename) |> as.data.frame() |> head(5)
    cNames <- DF |> colnames()
    cNames
  })

  # Import DEG KE dialog
  observe({
    shiny::req(input$deg_ke_enrichment_file$datapath, r$deg_ke_colnames())
    choices <- r$deg_ke_colnames()
    deg_ke_enrichment_filename <- input$deg_ke_enrichment_file$datapath
    DF <- readxl::read_excel(path=deg_ke_enrichment_filename) |> as.data.frame()
    column_types <- c("Experiment"="Experiment","Time"="Time","PFAS"="PFAS|Chem","Dose"="Dose|Dosage","Log2FC"="Log2FC|Fold_Change|FoldChange|Fold Change","Gene"="Gene|Feature","Gene_Count"="Gene|Feature|Count|Num","padj"="padj|adjusted|p.value|p-value|p value","Organ_Tissue"="Organ|Tissue","AOP_ID"="AOP|ID","KE_ID"="Key Event|Key_Event|KE|ID","KE_Name"="Key Event|Key_Event|KE|Name","KE_Description"="Key Event|Key_Event|KE|Desc","KE_Type"="Key Event|Key_Event|KE|Type","KE_Level"="Key Event|Key_Event|KE|Level")
    column_sel <- column_types |> purrr::map(\(val){get_best_match(choices, val)})
    shiny::showModal(myModalDialog(size="full",
      title = "Import Key Event Enrichment File",
      bslib::layout_columns(col_widths=c(12), row_heights=c(3),
        bslib::card(
          bslib::card_header("DEG Key Event Enrichment Preview"),
          bslib::card_body(
            shiny::tableOutput("deg_ke_preview_table")
          )
        )
      ),
      shiny::tags$hr(),
      bslib::layout_columns(col_widths=c(12,12,12,12,12,12,12), row_heights=c(1,1,1,1,1,1,1),
        bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput(inputId="deg_ke_experiment_col", 
            label=bslib::tooltip(trigger=list("Experiment Column", bsicons::bs_icon("info-circle", size="1em")), "Experiment and time information as <exp>_<time>, e.g. PFOS_10"), 
            choices=choices, selected=column_sel$Experiment),
          shiny::selectInput(inputId="deg_ke_time_col", 
            label=bslib::tooltip(trigger=list("Time Column", bsicons::bs_icon("info-circle", size="1em")), "Time information, e.g. 10"), 
            choices=choices, selected=column_sel$Time),
          shiny::selectInput("deg_ke_pfas_col", 
            label=bslib::tooltip(trigger=list("PFAS Column", bsicons::bs_icon("info-circle", size="1em")), "PFAS chemical information, e.g. PFOS"), 
            choices=choices, selected=column_sel$PFAS)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("deg_ke_dose_col", 
            label=bslib::tooltip(trigger=list("Dose Column", bsicons::bs_icon("info-circle", size="1em")), "Dose information, e.g. 10uM"), 
            choices=choices, selected=column_sel$Dose),
          shiny::selectInput("deg_ke_genes_col", 
            label=bslib::tooltip(trigger=list("Genes Column", bsicons::bs_icon("info-circle", size="1em")), "Genes information, semicolon separated, e.g. ALDH2;ETFA;IDO1;PHGDH;POR;STEAP4"), 
            choices=choices, selected=column_sel$Gene),
          shiny::selectInput("deg_ke_relevant_genes_col", 
            label=bslib::tooltip(trigger=list("Genes Count Column", bsicons::bs_icon("info-circle", size="1em")), "Number of genes, e.g. 6"), 
            choices=choices, selected=column_sel$Gene_Count)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("deg_ke_log2fc_col", 
            label=bslib::tooltip(trigger=list("Log2FC Column", bsicons::bs_icon("info-circle", size="1em")), "Log2 fold change information, e.g. 0.614899"), 
            choices=choices, selected=column_sel$Log2FC),
          shiny::selectInput("deg_ke_padj_col", 
            label=bslib::tooltip(trigger=list("Adjusted P-value Column", bsicons::bs_icon("info-circle", size="1em")), "Adjusted P-value information, e.g. 0.01999476"), 
            choices=choices, selected=column_sel$padj),
          shiny::selectInput("deg_ke_organ_tissue_col", 
            label=bslib::tooltip(trigger=list("Organ Tissue Column", bsicons::bs_icon("info-circle", size="1em")), "Organ or tissue information, e.g. liver"), 
            choices=choices, selected=column_sel$Organ_Tissue)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("deg_ke_aop_col", 
            label=bslib::tooltip(trigger=list("AOP ID Column", bsicons::bs_icon("info-circle", size="1em")), "Adverse Outcome Pathway ID, e.g. Aop:437"), 
            choices=choices, selected=column_sel$AOP_ID),
          shiny::selectInput("deg_ke_termid_col", 
            label=bslib::tooltip(trigger=list("Key Event ID Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event ID, e.g. Event:105"), 
            choices=choices, selected=column_sel$KE_ID),
          shiny::selectInput("deg_ke_key_event_name_col", 
            label=bslib::tooltip(trigger=list("Key Event Name Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event name, e.g. Inhibition, Mitochondrial Electron Transport Chain Complexes"), 
            choices=choices, selected=column_sel$KE_Name)
        ),bslib::layout_columns(col_widths=c(4,4,4),
          shiny::selectInput("deg_ke_key_event_description_col", 
            label=bslib::tooltip(trigger=list("Key Event Description Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event description, e.g. Inhibition, Mitochondrial Electron Transport Chain Complexes"), 
            choices=choices, selected=column_sel$KE_Description),
          shiny::selectInput("deg_ke_key_event_type_col", 
            label=bslib::tooltip(trigger=list("Key Event Type Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event type, e.g. MolecularInitiatingEvent"), 
            choices=choices, selected=column_sel$KE_Type),
          shiny::selectInput("deg_ke_level_col", 
            label=bslib::tooltip(trigger=list("Key Event Level Column", bsicons::bs_icon("info-circle", size="1em")), "Key Event level, e.g. Molecular"), 
            choices=choices, selected=column_sel$KE_Level)
        ),
        shiny::tags$hr(),
        bslib::layout_columns(col_widths=c(-4,-4,4),
          shiny::actionButton(inputId="import_deg_ke_submit", label="Import", width="100%"),
          shiny::uiOutput("deg_ke_load_str_modal")
        )
      )
    ))
  }) |> shiny::bindEvent(input$deg_ke_enrichment_file$datapath)

  # Import DEG KE submit
  observeEvent(input$import_deg_ke_submit, {
    shiny::validate(
      need(!is.null(input$deg_ke_enrichment_file$datapath), "No DEG Key Event Enrichment file!"),
    )
    shiny::req(
      input$deg_ke_enrichment_file$datapath, 
      input$deg_ke_termid_col,
      input$deg_ke_experiment_col,
      input$deg_ke_time_col,
      input$deg_ke_pfas_col,
      input$deg_ke_dose_col,
      input$deg_ke_log2fc_col,
      input$deg_ke_genes_col,
      input$deg_ke_relevant_genes_col,
      input$deg_ke_padj_col,
      input$deg_ke_aop_col,
      input$deg_ke_key_event_name_col,
      input$deg_ke_key_event_description_col,
      input$deg_ke_key_event_type_col,
      input$deg_ke_level_col,
      input$deg_ke_organ_tissue_col
    )
    deg_ke_enrichment_filename <- input$deg_ke_enrichment_file$datapath
    cNames <- c(input$deg_ke_termid_col, input$deg_ke_experiment_col, input$deg_ke_time_col, input$deg_ke_pfas_col, input$deg_ke_dose_col, input$deg_ke_log2fc_col, input$deg_ke_genes_col, input$deg_ke_relevant_genes_col, input$deg_ke_padj_col, input$deg_ke_aop_col, input$deg_ke_key_event_name_col, input$deg_ke_key_event_description_col, input$deg_ke_key_event_type_col, input$deg_ke_level_col, input$deg_ke_organ_tissue_col) 
    dup_check <- cNames |> duplicated() |> sum()
    if(dup_check > 0){
      shiny::showNotification("Duplicated column names. Please check and update.", duration=NULL)
    }
    shiny::req(dup_check==0)
    r$deg_ke_import_file <- deg_ke_enrichment_filename
    shiny::removeModal()
    bslib::nav_select(id="navset_card_de_input", selected="Differential Analysis Stats", session=session)
  })

  # deg ke import status
  output$deg_ke_load_str <- renderUI({
    if(is.null(r$deg_ke_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  output$deg_ke_load_str_modal <- renderUI({
    if(is.null(r$deg_ke_import_file)){
      shiny::tags$span(style="color:red", "Waiting for Import!")
    }else{
      shiny::tags$span(style="color:green", "File Ready!")
    }
  })

  # Upload Differential Expression 
  observeEvent(input$upload_deg_submit, {
    shiny::validate(
      need(!is.null(r$deg_import_file), "No DEG file!"),
      need(!is.null(r$deg_ke_import_file), "No DEG Key Event enrichment file!")
    )
    shiny::req(r$deg_import_file, r$deg_ke_import_file)

    w$show()
    on.exit({
      w$hide()
    })

    tryCatch({
      # read DE data from file
      deg_data_filename <- r$deg_import_file
      #deg_data <- bmdx::read_excel_allsheets(filename=deg_data_filename, tibble=FALSE, first_col_as_rownames=FALSE, check_numeric=FALSE)
      deg_data <- readxl::read_excel(path=deg_data_filename) |> as.data.frame()
      ##deg_data <- deg_data |> purrr::map(\(x){x |> dplyr::rename(Experiment={input$deg_experiment_col}, time={input$deg_time_col}, Feature={input$deg_feature_col}, dose={input$deg_dose_col}, log2FoldChange={input$deg_log2fc_col})})
      #deg_data <- deg_data |> purrr::map(\(x){x |> dplyr::select(Experiment={input$deg_experiment_col}, time={input$deg_time_col}, Feature={input$deg_feature_col}, dose={input$deg_dose_col}, log2FoldChange={input$deg_log2fc_col})})
      #deg_data <- deg_data |> purrr::map(\(x){x |> dplyr::mutate(dplyr::across(c(time, log2FoldChange), as.numeric))})
      deg_data <- deg_data |> dplyr::select(Experiment={input$deg_experiment_col}, time={input$deg_time_col}, Feature={input$deg_feature_col}, dose={input$deg_dose_col}, log2FoldChange={input$deg_log2fc_col})
      deg_data <- deg_data |> dplyr::mutate(dplyr::across(c(time, log2FoldChange), as.numeric))
      r$deg_data <- deg_data

      # read DE KE enrichment data from file
      deg_ke_enrichment_data_filename <- r$deg_ke_import_file
      deg_ke_enrichment_data <- bmdx::read_excel_allsheets(filename=deg_ke_enrichment_data_filename, tibble=FALSE, first_col_as_rownames=FALSE, check_numeric=FALSE)
      #deg_ke_enrichment_data <- deg_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::rename(TermID=input$deg_ke_termid_col, Experiment={input$deg_ke_experiment_col}, time={input$deg_ke_time_col}, PFAS={input$deg_ke_pfas_col}, dose={input$deg_ke_dose_col}, log2FoldChange={input$deg_ke_log2fc_col}, Genes={input$deg_ke_genes_col}, relevantGenesInGeneSet={input$deg_ke_relevant_genes_col}, padj={input$deg_ke_padj_col}, Aop={input$deg_ke_aop_col}, key_event_name={input$deg_ke_key_event_name_col}, Ke_description={input$deg_ke_key_event_description_col}, Ke_type={input$deg_ke_key_event_type_col}, level={input$deg_ke_level_col}, organ_tissue={input$deg_ke_organ_tissue_col})})
      deg_ke_enrichment_data <- deg_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::select(TermID=input$deg_ke_termid_col, Experiment={input$deg_ke_experiment_col}, time={input$deg_ke_time_col}, PFAS={input$deg_ke_pfas_col}, dose={input$deg_ke_dose_col}, log2FoldChange={input$deg_ke_log2fc_col}, Genes={input$deg_ke_genes_col}, relevantGenesInGeneSet={input$deg_ke_relevant_genes_col}, padj={input$deg_ke_padj_col}, Aop={input$deg_ke_aop_col}, key_event_name={input$deg_ke_key_event_name_col}, Ke_description={input$deg_ke_key_event_description_col}, Ke_type={input$deg_ke_key_event_type_col}, level={input$deg_ke_level_col}, organ_tissue={input$deg_ke_organ_tissue_col})})
      deg_ke_enrichment_data <- deg_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::mutate(dplyr::across(c(time, log2FoldChange, relevantGenesInGeneSet, padj), as.numeric))})
      deg_ke_enrichment_data <- deg_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
      deg_ke_enrichment_data <- deg_ke_enrichment_data |> dplyr::mutate(chem=Experiment |> stringr::str_replace(pattern="_[:alnum:]*$", replacement=""), time=Experiment |> stringr::str_replace(pattern=".*_", replacement="") |> as.numeric()) |> dplyr::group_by(chem) |> dplyr::group_split() |> as.list()
      nm <- deg_ke_enrichment_data |> purrr::map(\(x){x |> dplyr::pull(chem) |> unique()}) |> purrr::flatten_chr()
      names(deg_ke_enrichment_data) <- nm
      r$deg_ke_enrichment_data <- deg_ke_enrichment_data |> purrr::map(as.data.frame)
    }, error=\(e){
      str <- paste0("Encountered Error: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })
    bslib::nav_select(id="navset_card_de_input", selected="DEG Table", session=session)
  })

  # render DT BMDx stats
  output$bmd_stats_data_dt <- DT::renderDT({
    shiny::validate(need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!"))

    DF <- r$bmd_stats_data

    DT::datatable(DF, filter="top", selection='single',
      options=list(
        search=list(regex=TRUE, caseInsensitive=FALSE),
        scrollX=TRUE,
        ordering=TRUE,
        rownames=FALSE)
    )
  }, server=TRUE)

  # render BMDx stats heatmap
  observe({
    shiny::validate(
      need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!")
    )
    DF <- r$bmd_stats_data
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, 
      ComplexHeatmap::draw(
        plot_heatmap_preprocess(DF, score_col="BMD", name="BMD", toDraw=FALSE), 
        heatmap_legend_side = "left", 
        annotation_legend_side = "top", 
        padding = unit(c(2, 2, 2, 20),"mm")
      ), 
      heatmap_id="bmdx_stats_heatmap")
  })

  # render UI BMDx KE enrichment
  # render select experiment BMD KE enrichment
  output$select_experiment_bmd_ke_enrichment_data <- renderUI({
    shiny::validate(need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"))
    choices <- r$bmd_ke_enrichment_data |> names()
    selectInput(inputId="select_experiment_bmd_ke_enrichment_data_id", label="Select Experiment", choices=choices, selected=choices, multiple=TRUE)
  })

  # render DT BMDx KE enrihcment
  output$bmd_ke_enrichment_data_dt <- DT::renderDT({
    shiny::validate(
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_bmd_ke_enrichment_data_id), "No BMD KE enrichment experiment selected!")
    )
    shiny::req(input$select_experiment_bmd_ke_enrichment_data_id)

    exp_id <- input$select_experiment_bmd_ke_enrichment_data_id
    DF <- r$bmd_ke_enrichment_data[exp_id] |> purrr::reduce(dplyr::bind_rows)

    DT::datatable(DF, filter="top", selection='single',
      options=list(
        search=list(regex=TRUE, caseInsensitive=FALSE),
        scrollX=TRUE,
        ordering=TRUE,
        rownames=FALSE)
    )
  }, server=TRUE)

  # render BMDx KE Heatmap
  observeEvent(input$select_experiment_bmd_ke_enrichment_data_id, {
    shiny::validate(
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_bmd_ke_enrichment_data_id), "No BMD KE enrichment experiment selected!")
    )
    shiny::req(input$select_experiment_bmd_ke_enrichment_data_id)
    exp_id <- input$select_experiment_bmd_ke_enrichment_data_id
    DF <- r$bmd_ke_enrichment_data[exp_id] |> purrr::reduce(dplyr::bind_rows)
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, 
      ComplexHeatmap::draw(
        plot_heatmap_ke(DF, score_col="BMD", name="BMD", toDraw=FALSE), 
        heatmap_legend_side = "left", 
        annotation_legend_side = "top", 
        padding = unit(c(2, 2, 2, 290),"mm")
      ), 
      heatmap_id="bmdx_ke_heatmap")
  })

  # render BMDx KE distribution plot
  output$bmdx_ke_dist <- renderPlot({
    shiny::validate(
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_bmd_ke_enrichment_data_id), "No BMD KE enrichment experiment selected!")
    )
    shiny::req(input$select_experiment_bmd_ke_enrichment_data_id, input$score_dist, input$level_dist)
    score <- input$score_dist
    level <- input$level_dist
    if(level=="All"){
      level=NULL
    }
    exp_id <- input$select_experiment_bmd_ke_enrichment_data_id
    DF <- r$bmd_ke_enrichment_data[exp_id] |> purrr::reduce(dplyr::bind_rows)

    p <- dist_plot(KE_annotated=DF, level=level, score=score)
    print(p)
  })

  # render UI DE stats
  # render select experiment DE
  #output$select_experiment_deg_data <- renderUI({
  #  shiny::validate(need(!is.null(r$deg_data), "No DEG data loaded!"))
  #  choices <- r$deg_data |> names()
  #  selectInput(inputId="select_experiment_deg_data_id", label="Select Experiment", choices=choices, selected=choices, multiple=TRUE)
  #})

  # render DT DE
  output$deg_data_dt <- DT::renderDT({
    #shiny::validate(
    #  need(!is.null(r$deg_data), "No DEG data loaded!"),
    #  need(!is.null(input$select_experiment_deg_data_id), "No DEG experiment selected!")
    #)
    #shiny::req(input$select_experiment_deg_data_id)
    shiny::validate(need(!is.null(r$deg_data), "No DEG data loaded!"))
    shiny::req(r$deg_data)

    #exp_id <- input$select_experiment_deg_data_id
    #DF <- r$deg_data[exp_id] |> purrr::reduce(dplyr::bind_rows)
    DF <- r$deg_data

    DT::datatable(DF, filter="top", selection='single',
      options=list(
        search=list(regex=TRUE, caseInsensitive=FALSE),
        scrollX=TRUE,
        ordering=TRUE,
        rownames=FALSE)
    )
  }, server=TRUE)
  
  # render DE Heatmap
  observeEvent(r$deg_data, {
  #observeEvent(input$select_experiment_deg_data_id, {
    #shiny::validate(
    #  need(!is.null(r$deg_data), "No DEG data loaded!"),
    #  need(!is.null(input$select_experiment_deg_data_id), "No DEG experiment selected!")
    #)
    #shiny::req(input$select_experiment_deg_data_id)
    shiny::validate(need(!is.null(r$deg_data), "No DEG data loaded!"))
    shiny::req(r$deg_data)

    #exp_id <- input$select_experiment_deg_data_id
    #DF <- r$deg_data[exp_id] |> purrr::reduce(bind_rows)
    DF <- r$deg_data
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, 
      ComplexHeatmap::draw(
        plot_heatmap_preprocess(DF, score_col="log2FoldChange", name="log2FoldChange", toDraw=FALSE), 
        heatmap_legend_side = "left", 
        annotation_legend_side = "top", 
        padding = unit(c(2, 2, 2, 20),"mm")
      ), 
      heatmap_id="deg_heatmap")
  })

  # render UI DE KE enrichment
  # render select DE KE enrichment
  output$select_experiment_deg_ke_enrichment_data <- renderUI({
    shiny::validate(need(!is.null(r$deg_ke_enrichment_data), "No DEG KE enrichment data loaded!"))
    choices <- r$deg_ke_enrichment_data |> names()
    selectInput(inputId="select_experiment_deg_ke_enrichment_data_id", label="Select Experiment", choices=choices, selected=choices, multiple=TRUE)
  })

  # render DT KE enrichment
  output$deg_ke_enrichment_data_dt <- DT::renderDT({
    shiny::validate(
      need(!is.null(r$deg_ke_enrichment_data), "No DEG KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_deg_ke_enrichment_data_id), "No DEG KE enrichment experiment selected!")
    )
    shiny::req(input$select_experiment_deg_ke_enrichment_data_id)

    exp_id <- input$select_experiment_deg_ke_enrichment_data_id
    DF <- r$deg_ke_enrichment_data[exp_id] |> purrr::reduce(dplyr::bind_rows)

    DT::datatable(DF, filter="top", selection='single',
      options=list(
        search=list(regex=TRUE, caseInsensitive=FALSE),
        scrollX=TRUE,
        ordering=TRUE,
        rownames=FALSE)
    )
  }, server=TRUE)

  # render DE KE Heatmap
  observeEvent(input$select_experiment_deg_ke_enrichment_data_id, {
    shiny::validate(
      need(!is.null(r$deg_ke_enrichment_data), "No DEG KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_deg_ke_enrichment_data_id), "No DEG KE enrichment experiment selected!"),
    )
    shiny::req(input$select_experiment_deg_ke_enrichment_data_id)
    exp_id <- input$select_experiment_deg_ke_enrichment_data_id
    DF <- r$deg_ke_enrichment_data[exp_id] |> purrr::reduce(dplyr::bind_rows)
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, 
      ComplexHeatmap::draw(
        plot_heatmap_ke(DF, score_col="log2FoldChange", name="log2FoldChange", toDraw=FALSE), 
        heatmap_legend_side = "left", 
        annotation_legend_side = "top", 
        padding = unit(c(2, 2, 2, 290),"mm")
      ), 
      heatmap_id="deg_ke_heatmap")
  })

  # perform Forward Dosimetry analysis
  observeEvent(input$set_pbpk_submit, {
    shiny::validate(
      need(!is.null(input$substance), "No Chemical selected!"),
      need(!is.null(input$body_weight), "No Body Weight value!"),
      need(!is.null(input$ingestion_rate), "No Ingestion value!"),
      need(!is.null(input$exposure_time), "No Exposure Time value!")
    )
    shiny::req(input$substance, input$body_weight, input$ingestion_rate, input$exposure_time, iv.fw_dosimetry$is_valid())

    w$show()
    on.exit({
      w$hide()
    })
    
    substance <- input$substance
    body_weight <- input$body_weight
    exposure_time <- input$exposure_time
    stat_summary <- input$stat_summary
    
    if(length(substance)==2){
      substance <- "BOTH"
    }
    
    ingestion_rates <- input$ingestion_rate |> unique()

    resFWD <- list()
    resFWD.summary <- list()
    simulation_results <- list()
    
    i=1
    scenarios <- list()

    tryCatch({
      for(ingestion_rate in ingestion_rates){
        x <- paste0("scenario_",i)
        
        resFWD[[x]] <- PfasDosim::forward_dosimetry(chemical=substance, ingestion=as.numeric(ingestion_rate), BW=body_weight, duration=exposure_time, stat_sum=FALSE)
        resFWD[[x]] <- resFWD[[x]] |> purrr::discard(is.null) 
        resFWD.summary[[x]] <- PfasDosim::forward_dosimetry(chemical=substance, ingestion=as.numeric(ingestion_rate), BW=body_weight, duration=exposure_time, stat_sum=TRUE)
        resFWD.summary[[x]] <- resFWD.summary[[x]] |> purrr::discard(is.null) 
        simulation_results[[x]] <- compiling_forward_dosimetry_simulation_dataframes_all_samples(resFWD[[x]])
        scenarios[[x]] <- ingestion_rate
        i <- i+1
      }  
      r$resFWD <- resFWD

      r$resFWD.summary <- resFWD.summary
      r$simulation_results <- simulation_results
      r$scenarios <- scenarios

      if(!is.null(r$resFWD)){
        print("done")
      }
    }, error=\(e){
      str <- paste0("Encountered Error: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })
    
    bslib::nav_select(id="navset_card_fw_dosimetry_result", selected="Result Table", session=session)
  }) 

  # render UI Forward Dosimetry
  # render select Forward Dosimetry scenario
  output$render_select_fw_dosimetry_scenario <- renderUI({
    shiny::validate(need(!is.null(r$resFWD), "No Forward Dosimetry results!"))
    scenarios <- r$scenarios
    choices <- names(scenarios) |> as.list() |> setNames(paste0("Ingestion Rate: ", unname(scenarios), " ng/kg/day"))
    shiny::selectInput(inputId="select_fw_dosimetry_scenario", label="Select Ingestion Scenario",
      choices=choices,
      selected=choices[1],
      multiple=FALSE,
      width="100%"
    )
  })

  # render selectize Forward Dosimetry scenario
  output$render_selectize_fw_dosimetry_scenario <- renderUI({
    shiny::validate(need(!is.null(r$resFWD), "No Forward Dosimetry results!"))
    scenarios <- r$scenarios
    choices <- names(scenarios) |> as.list() |> setNames(paste0("Ingestion Rate: ", unname(scenarios), " ng/kg/day"))
    shiny::selectizeInput(inputId="selectize_fw_dosimetry_scenario", label="Select Ingestion Scenario",
      choices=choices,
      selected=choices[1],
      multiple=TRUE,
      width="100%"
    )
  })

  # render select Forward Dosimetry chemical
  output$render_select_fw_dosimetry_chemical <- renderUI({
    shiny::validate(need(!is.null(input$substance), "No Chemical selected!"))
    choices <- input$substance
    shiny::selectInput(inputId="pfas_forward", label="Select Chemical",
      choices=choices,
      selected=choices[1],
      multiple=FALSE,
      width="100%"
    )
  })

  # render selectize Forward Dosimetry chemical
  output$render_selectize_fw_dosimetry_chemical <- renderUI({
    shiny::validate(need(!is.null(input$substance), "No Chemical selected!"))
    choices <- input$substance
    shiny::selectizeInput(inputId="pfas_forward_plot", label="Select Chemical(s)",
      choices=choices,
      selected=choices,
      multiple=TRUE,
      width="100%"
    )
  })

  # render DT Forward Dosimetry 
  output$fw_dosimetry_simulation_results_single_model_DT <- DT::renderDT({
    shiny::validate(
      need(!is.null(r$resFWD), "No Forward Dosimetry results!"),
      need(!is.null(input$pfas_forward), "No Chemical selected!"),
      need(!is.null(input$compartment_forward), "No Compartment selected!")
    )
    shiny::req(input$pfas_forward, input$compartment_forward, input$select_fw_dosimetry_scenario)
    
    scenario <- input$select_fw_dosimetry_scenario
    pfas_forward <- input$pfas_forward
    compartment_forward <- input$compartment_forward |> stringr::str_c("_df")
    DF <- r$resFWD[[scenario]][[pfas_forward]][["single_models"]][[compartment_forward]]
    
    DT::datatable(DF, filter="top",
                  selection='single',
                  options=list(
                    search=list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=TRUE,
                    rownames=FALSE
                  ))
  },server=TRUE)

  # render DT Forward Dosimetry consensus 
  output$fw_dosimetry_simulation_results_consensus_summary_DT <- DT::renderDT({
    shiny::validate(
      need(!is.null(r$resFWD.summary), "No Forward Dosimetry results!"),
      need(!is.null(input$pfas_forward), "No Chemical selected!"),
      need(!is.null(input$compartment_forward), "No Compartment selected!")
    )
    shiny::req(input$pfas_forward, input$compartment_forward, input$select_fw_dosimetry_scenario)
    
    scenario <- input$select_fw_dosimetry_scenario
    pfas_forward <- input$pfas_forward
    compartment_forward <- input$compartment_forward |> stringr::str_c("_df")
    DF <- r$resFWD.summary[[scenario]][[pfas_forward]][["bma"]][[compartment_forward]]
    
    DT::datatable(DF, filter="top",
                  selection='single',
                  options=list(
                    search=list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=TRUE,
                    rownames=FALSE
                  ))
  },server=TRUE)

  # download Fowrard Dosimetry single model results xlsx
  output$download_fw_dosimetry_single_tables <- downloadHandler(
    filename = "forward_dosimetry_result_tables.xlsx",
    content = function(file) {
      shiny::req(r$resFWD)
      resFWD <- r$resFWD
      nm.sc <- resFWD |> names()
      nm.mod <- c()
      for(nm in nm.sc){
        nm.mod <- append(nm.mod, paste0("IR-", r$scenarios[[nm]], "ngPkgPday"))
      }
      names(resFWD) <- nm.mod

      DF.list <- resFWD |> purrr::list_flatten() |> purrr::map("single_models") |> purrr::list_flatten()

      WriteXLS::WriteXLS(DF.list, ExcelFileName=file)
    }
  )

  # download Fowrard Dosimetry consensus results xlsx
  output$download_fw_dosimetry_consensus_tables <- downloadHandler(
    filename = "forward_dosimetry_consensus_result_tables.xlsx",
    content = function(file) {
      shiny::req(r$resFWD)
      resFWD <- r$resFWD
      nm.sc <- resFWD |> names()
      nm.mod <- c()
      for(nm in nm.sc){
        nm.mod <- append(nm.mod, paste0("IR-", r$scenarios[[nm]], "ngPkgPday"))
      }
      names(resFWD) <- nm.mod

      DF.list <- resFWD |> purrr::list_flatten() |> purrr::map("bma") |> purrr::list_flatten()

      WriteXLS::WriteXLS(DF.list, ExcelFileName=file)
    }
  )

  # download Fowrard Dosimetry consensus summary results xlsx
  output$download_fw_dosimetry_consensus_summary_tables <- downloadHandler(
    filename = "forward_dosimetry_consensus_summary_result_tables.xlsx",
    content = function(file) {
      shiny::req(r$resFWD.summary)
      resFWD <- r$resFWD.summary
      nm.sc <- resFWD |> names()
      nm.mod <- c()
      for(nm in nm.sc){
        nm.mod <- append(nm.mod, paste0("IR-", r$scenarios[[nm]], "ngPkgPday"))
      }
      names(resFWD) <- nm.mod

      DF.list <- resFWD |> purrr::list_flatten() |> purrr::map("bma") |> purrr::list_flatten()

      WriteXLS::WriteXLS(DF.list, ExcelFileName=file)
    }
  )

  # render Forward Dosimetry simulation plot
  output$fw_dosimetry_simulation_results_plot <- renderPlot({
    shiny::validate(
      need(!is.null(r$resFWD), "No Forward Dosimetry results!"),
      need(!is.null(input$pfas_forward_plot), "No Chemical selected!"),
      need(!is.null(input$compartment_forward_plot), "No Compartment selected!"),
      need(!is.null(input$selectize_fw_dosimetry_scenario), "No Scenario selected!")
    )
    shiny::req(input$pfas_forward_plot, input$compartment_forward_plot, input$selectize_fw_dosimetry_scenario)

    scenario <- input$selectize_fw_dosimetry_scenario
    ingestion_rate <- paste0("Ingestion Rate: ", r$scenarios[[scenario]], " ng/kg/day")
    pfas_forward <- input$pfas_forward_plot
    compartment_forward <- input$compartment_forward_plot
    resFWD <- r$resFWD

    plot_PBK_simulations_func <- function(...){
      if(length(compartment_forward)>1){
        PBK_plot <- plot_PBK_simulations(...)
        PBK_plot_serum <- plot_PBK_simulations_serum(...)
        combined <- cowplot::plot_grid(
          #PBK_plot + scale_y_continuous(limits = c(0, 1500)),
          #PBK_plot_serum + scale_y_continuous(limits = c(0, 1500)),
          PBK_plot,
          PBK_plot_serum,
          ncol = 1,
          align = "h"
        )
        return(combined)
      }else if(compartment_forward=="serum"){
        plot_PBK_simulations_serum(...)
      }else if(compartment_forward=="liver"){
        plot_PBK_simulations(...)
      }
    }

    p <- do.call(plot_PBK_simulations_func, resFWD |> append(list(scenario_labels=ingestion_rate, chemicals=pfas_forward)))
    print(p)
  })
  
  # render select Forward Dosimetry scenario
  output$render_select_fw_dosimetry_scenario_combined <- renderUI({
    shiny::validate(need(!is.null(r$resFWD), "No Forward Dosimetry results!"))
    scenarios <- r$scenarios
    choices <- names(scenarios) |> as.list() |> setNames(paste0("Ingestion Rate: ", unname(scenarios), " ng/kg/day"))
    shiny::selectInput(inputId="select_fw_dosimetry_scenario_combined", label="Select Ingestion Scenario",
      choices=choices,
      selected=choices[1],
      multiple=FALSE
    )
  })

  # render selectize time point from BMDx stats
  output$render_selectize_time_point <- renderUI({
    shiny::validate(need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!"))
    bmd_stats_data <- r$bmd_stats_data
    choices <- bmd_stats_data |> dplyr::pull(time) |> unique()
    shiny::selectizeInput(inputId="select_time_point", label="Select Experiment Time Point",
      choices=choices,
      selected=choices,
      multiple=TRUE
    )
  })

  # render selectize organ from BMDx KE enrichment
  output$render_selectize_organ <- renderUI({
    shiny::validate(need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"))
    bmd_ke_enrichment_data <- r$bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
    choices <- bmd_ke_enrichment_data |> dplyr::pull(organ_tissue) |> unique()
    shiny::selectizeInput(inputId="select_organ_tissue", label="Select Organ/Tissue",
      #choices=c("NA", choices),
      #selected="NA",
      choices=choices,
      selected=choices[1],
      multiple=TRUE
    )
  })

  # render selectize KE for organ from BMDx KE enrichment
  output$render_selectize_organ_ke <- renderUI({
    shiny::validate(need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"))
    if(is.null(input$select_organ_tissue)) {
      choices <- "NA"
      shiny::selectizeInput(inputId="select_organ_tissue_ke", label="Select Key Events", choices=choices, selected=choices, multiple=TRUE)
    }
    organs <- input$select_organ_tissue
    bmd_ke_enrichment_data <- r$bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
    choices <- bmd_ke_enrichment_data |> dplyr::filter(organ_tissue %in% organs) |> dplyr::select(key_event_name, TermID) |> dplyr::distinct() |> tibble::deframe()
    shiny::selectizeInput(inputId="select_organ_tissue_ke", label="Select Key Events",
      choices=choices,
      selected=choices,
      multiple=TRUE
    )
  })

  # perform Integrated analysis
  observeEvent(input$combine_models_submit, {
    shiny::validate(
      need(!is.null(input$compartment), "No Compartment selected!"),
      need(!is.null(input$substance), "No Chemical selected!"),
      need(!is.null(r$simulation_results), "No Forward Dosimetry simulation performed!"),
      need(!is.null(r$resFWD), "No Forward Dosimetry results!"),
      need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!"),
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!")
    )
    shiny::req(input$compartment, input$substance, input$select_fw_dosimetry_scenario_combined, iv.integrated$is_valid())

    w$show()
    on.exit({
      w$hide()
    })

    scenario <- input$select_fw_dosimetry_scenario_combined
    PFAS <- input$substance
    tissues_of_interest <- input$compartment
    exposure_time <- input$exposure_time |> as.numeric()
    time_points <- as.numeric(input$select_time_point)
    n_permutations <- as.numeric(input$n_permutations)
    bmd_summarization_percentile <- as.numeric(input$bmd_summarization_percentile)
    pod_variable <- "BMD" # static for prob
    time_unit <- "Day"

    Simulation_results <- r$simulation_results[[scenario]]
    substance <- Simulation_results |> names()
    Enrichment_data <- r$bmd_ke_enrichment_data[substance] |> purrr::reduce(dplyr::bind_rows)
    optimal_models_stats <- r$bmd_stats_data

    # Filter by Key Event
    if(input$filter_keyevents){
      if(all(!is.na(input$select_organ_tissue_ke)) | all(input$select_organ_tissue_ke!="NA")){
        organ_tissue_ke <- input$select_organ_tissue_ke
        Enrichment_data <- Enrichment_data |> dplyr::filter(TermID %in% organ_tissue_ke)
      }
    }

    Enrichment_data$BMDL <- as.numeric(Enrichment_data$BMDL)
    Enrichment_data$BMDU <- as.numeric(Enrichment_data$BMDU)
    Enrichment_data$time <- as.numeric(Enrichment_data$time)
    unique_exp <- unique(Enrichment_data$Experiment)

    tryCatch({
      probability_matrices <- compute_PBK_BMD_probability_matrix_activation(
        tissues_of_interest = tissues_of_interest,
        unique_exp = unique_exp,
        optimal_models_stats = optimal_models_stats,
        Enrichment_data = Enrichment_data,
        NPerm = n_permutations,
        BMD_summarization_percentile = bmd_summarization_percentile,
        exposure_time = exposure_time,
        pod_variable=pod_variable,
        Simulation_results = Simulation_results
      )

      probability_matrices <- probability_matrices |> purrr::discard(is.null)
      if(length(probability_matrices)==0){
        shiny::showNotification("No valid results for chosen scenario. Adjust and try again", duration=NULL)
      }else{
        df_tissue_specific_list <- probability_matrices |> purrr::map(\(x){convert_tensor_to_dataframes(probability_matrices=x, time_points=time_points, time_unit=time_unit)})

        r$probability_matrices <- probability_matrices
        r$df_tissue_specific_list <- df_tissue_specific_list
      }
    }, error=\(e){
      str <- paste0("Encountered Error: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })

    #r$probability_matrices <- readRDS("data/integrated_result_probability_matrices.Rds")
    #r$df_tissue_specific_list <- readRDS("data/integrated_result_table.Rds")

    bslib::nav_select(id="navset_card_itegrated_result", selected="Result Table", session=session)
  })

  # render DT Integrated 
  output$integrated_probability_results_DT <- DT::renderDT({
    shiny::validate(need(!is.null(r$df_tissue_specific_list), "No Integrated Analysis results!"))
    
    DF <- r$df_tissue_specific_list[[1]]
    
    DT::datatable(DF , filter="top",
      selection='single',
      options=list(
        search=list(regex=TRUE, caseInsensitive=FALSE),
        scrollX=TRUE,
        ordering=TRUE,
        rownames=FALSE
    ))
  },server=TRUE)

  # render selectize Integrated KE
  output$render_selectize_integrated_ke <- renderUI({
    shiny::validate(need(!is.null(r$df_tissue_specific_list), "No Integrated Analysis results!"))
    df_tissue_specific <- r$df_tissue_specific_list[[1]]
    choices <- df_tissue_specific |> dplyr::pull(KeyEvent) |> unique()

    shiny::selectizeInput(inputId="select_integrated_ke", label="Select Key Events",
      choices=choices,
      selected=choices,
      width="100%",
      multiple=TRUE
    )
  })

  # render activation plot Integrated
  output$integrated_probability_results_plot <- renderPlot({
    shiny::validate(need(!is.null(r$df_tissue_specific_list), "No Integrated Analysis results!"))
    shiny::req(input$select_integrated_ke) 

    df_tissue_specific <- r$df_tissue_specific_list[[1]]

    integrated_ke <- input$select_integrated_ke
    df_tissue_specific <- df_tissue_specific |> dplyr::filter(KeyEvent %in% integrated_ke)

    gp <- plot_KE_probabilities_over_PBK_time(
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

    print(gp)
  })

  # render KE bubble plot Integrated
  output$integrated_bubble_plot <- renderPlot({
    shiny::validate(
      need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!"),
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(r$deg_data), "No DEG data loaded!"),
      need(!is.null(r$df_tissue_specific_list), "No Integarted Analysis results!")
    )
    shiny::req(input$select_integrated_ke) 

    optimal_models_stats <- r$bmd_stats_data
    Enrichment_data <- r$bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
    #deg_statistics <- r$deg_data |> purrr::reduce(dplyr::bind_rows)
    deg_statistics <- r$deg_data
    df_tissue_specific <- r$df_tissue_specific_list[[1]]
    rel_ke <- input$select_integrated_ke
    PFASs <- df_tissue_specific |> dplyr::pull(PFAS) |> unique() 

    plotPatch <- NULL
    for(PFAS in PFASs){
      bbPlot <- bubble_plot(optimal_models_stats, Enrichment_data, deg_statistics, df_tissue_specific, rel_ke, PFAS)
      plotPatch <- plotPatch | bbPlot
    }
    print(plotPatch)
  })

  # render Forward Dosimetry PBK estimate with PODs BMD Boxplot plot
  output$pbk_pod_plot <- renderPlot({
    shiny::validate(
      need(!is.null(r$resFWD), "No Forward Dosimetry results!"),
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(input$exposure_time), "No Exposure Time value!"),
      need(!is.null(input$select_fw_dosimetry_scenario_combined), "No Scenario selected!")
    )
    shiny::req(input$select_fw_dosimetry_scenario_combined, input$exposure_time)

    scenario <- input$select_fw_dosimetry_scenario_combined
    ingestion_rate <- paste0("Ingestion Rate: ", r$scenarios[[scenario]], " ng/kg/day")
    exposure_time <- input$exposure_time
    resFWD <- r$resFWD[[scenario]]
    Enrichment_data <- r$bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)

    p <- plot_PBK_POD(Enrichment_data_BMD=Enrichment_data, output_forward_dosimetry=resFWD, exposure_time=exposure_time, exposure_level=ingestion_rate)
    print(p)
  })
    
  # download Integrated results 2D data.frame
  output$download_integrated_tables <- downloadHandler(
    filename = "integrated_summarized_results_table.xlsx",
    content = function(file) {
      shiny::req(r$df_tissue_specific_list)
      DF <- r$df_tissue_specific_list[[1]]
      WriteXLS::WriteXLS(DF, ExcelFileName=file, SheetNames="Summarized_Probability")
    }
  )

  # download Integrated results 3D probability matrices dataset
  output$download_integrated_matrices <- downloadHandler(
    filename = "integrated_probability_matrices.json",
    content = function(file) {
      shiny::req(r$probability_matrices)
      mat.json <- jsonlite::toJSON(r$probability_matrices)
      jsonlite::write_json(mat.json, path=file)
    }
  )

  # render select exposure time point 
  output$render_select_exposure_time_point <- renderUI({
    shiny::validate(need(!is.null(r$probability_matrices), "No Integrated Analysis results!"))
    probability_matrix <- r$probability_matrices[[input$compartment]][[1]]
    exposure_time <- probability_matrix |> dim() |> purrr::pluck(3)
    choices <- seq(exposure_time)
    shiny::selectInput(inputId="select_exposure_time_point", label="Select Exposure Time Point",
      choices=choices,
      selected=choices[1]
    )
  })

  # infer KE network Integrated
  observeEvent(input$net_submit, {
    shiny::validate(
      need(!is.null(input$compartment), "No Compartment selected!"),
      need(!is.null(r$df_tissue_specific_list), "No Integrated Analysis results!"),
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      #need(!is.null(r$deg_ke_enrichment_data), "No DEG KE enrichment data loaded!"),
      need(!is.null(r$probability_matrices), "No Integrated Analysis results!")
    )
    shiny::req(input$compartment, r$iv.integrated.net$is_valid())

    w$show()
    on.exit({
      w$hide()
    })

    PBK_time <- input$select_exposure_time_point |> as.numeric()
    # initial grouping by AOP, because user does not know the target 
    group_by <- "aop"

    enlarge_ke_selection <- input$enlarge_keyevents |> as.logical()
    if(enlarge_ke_selection){
      max_path_length <- input$max_path_length |> as.numeric()
      n_AOs <- input$max_path_length |> as.numeric()
      n_MIEs <- input$max_path_length |> as.numeric()
    }else{
      max_path_length <- NULL
      n_AOs <- NULL
      n_MIEs <- NULL
    }

    KE_annotated <- get_docking_ke_mapping()
    probability_matrices <- r$probability_matrices[[input$compartment]]

    Enrichment_data_BMD <- r$bmd_ke_enrichment_data
    Enrichment_data_BMD <- Enrichment_data_BMD |> purrr::reduce(dplyr::bind_rows)
    Enrichment_data_BMD <- Enrichment_data_BMD |> dplyr::mutate(expi=Experiment) |> dplyr::distinct()

    if(!is.null(r$deg_ke_enrichment_data)){
      Enrichment_data_DEG <- r$deg_ke_enrichment_data
      Enrichment_data_DEG <- Enrichment_data_DEG |> purrr::reduce(dplyr::bind_rows)
      Enrichment_data_DEG <- Enrichment_data_DEG |> dplyr::mutate(expi=Experiment) |> dplyr::distinct()
    }

    mode <- "out"
    ke_id <- "ke"
    numerical_variables <- c("prob")
    pval_variable <- "padj"
    gene_variable <- "Genes"
    convert_to_gene_symbols <- FALSE

    tryCatch({
      if(!is.null(r$deg_ke_enrichment_data)){
        net_list_res <- build_network_lists(
          probability_matrices = probability_matrices,
          Enrichment_data_DEG = Enrichment_data_DEG,
          Enrichment_data_BMD = Enrichment_data_BMD,
          KE_annotated = KE_annotated,
          experiments = names(probability_matrices),
          PBK_time = PBK_time,
          max_path_length = max_path_length,
          n_AOs = n_AOs,
          n_MIEs = n_MIEs,
          mode = mode,
          enlarge_ke_selection = enlarge_ke_selection,
          ke_id = ke_id,
          numerical_variables = numerical_variables,
          pval_variable = pval_variable,
          gene_variable = gene_variable,
          convert_to_gene_symbols = convert_to_gene_symbols,
          group_by = group_by
        )
        shiny::updateSelectInput(inputId="group_visnet_by", label="Group By", choices=c("Adverse Outcomes"="aop","Path to Adverse Outcome of Interest"="rank"), selected="aop")
      } else{
        net_list_res <- build_network_lists_only_PBK_BMD(
          probability_matrices = probability_matrices,
          Enrichment_data_BMD = Enrichment_data_BMD,
          experiments = names(probability_matrices),
          PBK_time = PBK_time,
          max_path_length = max_path_length,
          n_AOs = n_AOs,
          n_MIEs = n_MIEs,
          mode = mode,
          enlarge_ke_selection = enlarge_ke_selection,
          ke_id = ke_id,
          numerical_variables = numerical_variables,
          pval_variable = pval_variable,
          gene_variable = gene_variable,
          convert_to_gene_symbols = convert_to_gene_symbols,
          group_by = group_by
        )
      shiny::updateSelectInput(session, "group_visnet_by", choices=c("Adverse Outcomes"="aop"), selected="aop")
      }
      r$net_list_res <- net_list_res
    }, error=\(e){
      str <- paste0("Encountered Error: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })
  })

  # render select experiment KE net Integrated
  output$render_select_experiment_net <- renderUI({
    shiny::validate(need(!is.null(r$net_list_res), "No Network data!"))
    net_list_res <- r$net_list_res
    choices <- net_list_res$visnet_list |> names()
    selectInput(inputId="select_experiment_net", label="Select Experiment", choices=choices, selected=choices[1])
  })

  # render select target KE net Integrated
  r$aop_ke_df <- AOPfingerprintR::aop_ke_table_hure |> dplyr::select(Ke, Ke_type) |> dplyr::add_count(Ke_type) |> dplyr::distinct() |> tidyr::pivot_wider(names_from=Ke_type, values_from=n) |> dplyr::mutate(Ke_type=ifelse(is.na(AdverseOutcome),ifelse(is.na(MolecularInitiatingEvent),"KeyEvent","MolecularInitiatingEvent"),ifelse(is.na(MolecularInitiatingEvent),"AdverseOutcome",ifelse(AdverseOutcome>MolecularInitiatingEvent,"AdverseOutcome","MolecularInitiatingEvent"))))
  output$render_select_target_ke <- renderUI({
    shiny::validate(
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_net), "No Experiment selected!")
    )
    shiny::req(input$select_experiment_net)

    exp_id <- input$select_experiment_net
    net_list_res <- r$net_list_res
    inferred_visNet_data <- net_list_res$net_list[[exp_id]]

    map_df <- inferred_visNet_data$nodes |> dplyr::select(id, ke_description) |> dplyr::left_join(r$aop_ke_df |> dplyr::select(Ke, ke_type=Ke_type), by=c("id"="Ke"))
    choices <- map_df |> dplyr::filter(ke_type=="AdverseOutcome") |> dplyr::pull(ke_description) |> unique()
    #choices <- inferred_visNet_data$nodes |> dplyr::pull(ke_description) |> unique()
    sel <- choices[1]
    if(!is.null(r$target_node_val) && !is.na(r$target_node_val)){
      idx <- which(choices %in% r$target_node_val)
      if(length(idx)>0){
        sel <- choices[idx]
      }
    }

    shiny::selectInput(inputId="select_target_ke", label="Select Adverse Outcome of Interest",
      choices=choices,
      selected=sel
    )
  })

  # render KE network Integrated
  output$ke_visnetwork <- visNetwork::renderVisNetwork({
    shiny::validate(
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!"),
      need(!is.null(input$select_experiment_net), "No Experiment selected!"),
      need(!is.null(r$net_list_res), "No Network data!")
    )
    shiny::req(input$select_experiment_net, input$group_visnet_by)
    Enrichment_data_BMD <- r$bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
    exp_id <- input$select_experiment_net
    net_list_res <- r$net_list_res
    group_by <- input$group_visnet_by

    inferred_visNet_data <- net_list_res$net_list[[exp_id]]
    if(group_by=="rank" & !is.null(input$select_target_ke)){
      shiny::req(input$select_target_ke)
      target_node <- input$select_target_ke
      r$target_node_val <- target_node
      inferred_igraph <- net_list_res$igraph_list[[exp_id]]
      pal_node_fill <- net_list_res$pal_node_fill_list[[exp_id]]
      ke_list_exp <- net_list_res$ke_list_exp
      ke_all_deg_list <- net_list_res$ke_all_deg_list
      docking_list <- net_list_res$docking_list

      res <- get_paths_to_plot(
        inferred_igraph=inferred_igraph, 
        inferred_visNet_data=inferred_visNet_data, 
        ke_list_exp=ke_list_exp, 
        ke_all_deg_list=ke_all_deg_list, 
        docking_list=docking_list, 
        Enrichment_data_BMD=Enrichment_data_BMD, 
        expi=exp_id, 
        target_node=target_node)

      inferred_visNet_data <- res$inferred_visNet_data
      visnet <- get_visnet(nodes=inferred_visNet_data$nodes, edges=inferred_visNet_data$edges, group_by=group_by, col_pal=pal_node_fill)
    }else{
      visnet <- net_list_res[["visnet_list"]][[exp_id]]
    }
    r$visnet_data_integrated <- inferred_visNet_data
    r$visnet_integrated <- visnet
    visnet
  })

  # download integrated visNetwork data all experiments
  output$download_integrated_visnet_data_all <- downloadHandler(
    filename = "integrated_visnetwork_data_all_experiments.xlsx",
    content = function(file) {
      shiny::req(r$net_list_res)
      visNet_data.list <- r$net_list_res$net_list |> purrr::list_flatten() 
      WriteXLS::WriteXLS(visNet_data.list, ExcelFileName=file)
    }
  )

  # download integrated visNetwork on display
  output$download_integrated_visnet_display <- downloadHandler(
    filename = "integrated_visnetwork.html",
    content = function(file) {
      shiny::req(r$visnet_integrated)
      visNetwork::visSave(r$visnet_integrated, file, selfcontained=TRUE, background="white")
    }
  )
  
  # download integrated visNetwork data on display
  output$download_integrated_visnet_data_display <- downloadHandler(
    filename = "integrated_visnetwork_data.xlsx",
    content = function(file) {
      shiny::req(r$visnet_data_integrated)
      WriteXLS::WriteXLS(r$visnet_data_integrated, ExcelFileName=file)
    }
  )

  # render selectize time point for Reverse Dosimetry from BMDx stats
  output$render_selectize_time_point_rev_dosimetry <- renderUI({
    shiny::validate(need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!"))
    bmd_stats_data <- r$bmd_stats_data
    choices <- bmd_stats_data |> dplyr::pull(time) |> unique()
    shiny::selectizeInput(inputId="select_time_point_rev_dosimetry", label="Select Experiment Time Point",
      choices=choices,
      selected=choices,
      multiple=TRUE
    )
  })
  
  # render select number of cores for Reverse Dosimetry
  output$render_select_n_cores <- renderUI({
    choices <- parallelly::availableCores() |> as.numeric() |> seq(from=2)
    shiny::selectInput(inputId="n_cores", label="N Threads",
      choices=choices,
      selected=choices[1],
      multiple=FALSE
    )
  })

  # perform Reverse Dosimetry analysis
  observeEvent(input$run_reverse_dosimetry,{
    shiny::validate(
      need(!is.null(input$compartment_reverse_dosimetry), "No Compartment selected!"),
      need(!is.null(input$pfas_reverse_dosimetry), "No Chemical selected!"),
      need(!is.null(r$bmd_stats_data), "No BMD stats data loaded!"),
      need(!is.null(r$bmd_ke_enrichment_data), "No BMD KE enrichment data loaded!")
    )
    shiny::req(
      input$compartment_reverse_dosimetry, 
      input$pfas_reverse_dosimetry, 
      input$select_time_point_rev_dosimetry, 
      input$pod_reverse_dosimetry, 
      input$body_weight_reverse_dosimetry,
      r$iv.rev_dosimetry$is_valid()
    )

    w$show()
    on.exit({
      w$hide()
    })
    
    time_point <- as.numeric(input$select_time_point_rev_dosimetry)
    pfas <- input$pfas_reverse_dosimetry
    compartment <-  input$compartment_reverse_dosimetry
    pod <- input$pod_reverse_dosimetry
    duration <- input$exposure_duration_reverse_dosimetry
    BW <- input$body_weight_reverse_dosimetry
    isStochastic <- input$is_stochastic |> as.logical()
    n_samples <- input$n_samples |> as.numeric()
    isParallel <- input$is_parallel |> as.logical()
    n_cores <- input$n_cores |> as.numeric()

    bmd_ke_enrichment_data <- r$bmd_ke_enrichment_data |> purrr::reduce(dplyr::bind_rows)
    if(input$use_top_keyevents){
      shiny::req(input$filter_type_reverse_dosimetry, input$score_type_reverse_dosimetry, input$n_key_events)
      filter_type <- input$filter_type_reverse_dosimetry
      score_type <- input$score_type_reverse_dosimetry
      bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> dplyr::mutate(SCORE=get(score_type))
      if(filter_type=="n"){
        n <- input$n_key_events |> as.numeric()
        bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> dplyr::group_by(Experiment) |> dplyr::slice_min(order_by=get(score_type), n=n) |> dplyr::ungroup()
      }else if(filter_type=="perc"){
        prop <- input$n_key_events |> as.numeric() |> magrittr::divide_by(100)
        bmd_ke_enrichment_data <- bmd_ke_enrichment_data |> dplyr::group_by(Experiment) |> dplyr::slice_min(order_by=get(score_type), prop=prop) |> dplyr::ungroup()
      }
    }

    pfas_mw <- c("PFOA"=414.07,"PFOS"=500.13)

    tryCatch({
      results_reverse_dosimetry <- reverse_dosimetry_wrapping(bmd_ke_enrichment_data,
        duration = duration,
        BW = BW,
        PFAS_mw = pfas_mw,
        pod = pod,
        PFAS = pfas,
        time = time_point,
        compartments = compartment,
        isStochastic = isStochastic,
        isParallel = isParallel, 
        Nsamples = n_samples,
        n_cores = n_cores)

      r$results_reverse_dosimetry <- results_reverse_dosimetry
    }, error=\(e){
      str <- paste0("Encountered Error: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })
    #r$results_reverse_dosimetry <- readRDS("data/results_reverse_dosimetry.Rds")
    InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(input, output, session, 
      ComplexHeatmap::draw(
        plot_reverse_dosimetry_heatmap(r$results_reverse_dosimetry$reverse_dosimetry_list_summarized, toDraw=FALSE), 
        heatmap_legend_side = "left", 
        annotation_legend_side = "top", 
        padding = unit(c(2, 2, 2, 100),"mm")
      ), 
      heatmap_id="reverse_dosimetry_heatmap")
    
    bslib::nav_select(id="navset_card_rev_dosimetry_result", selected="Result Table", session=session)
  })

  # download reverse dosimetry results
  output$download_rev_dosimetry_tables <- downloadHandler(
    filename = "reverse_dosimetry_results_table.xlsx",
    content = function(file) {
      shiny::req(r$results_reverse_dosimetry)
      reverse_dosimetry_list <- r$results_reverse_dosimetry$reverse_dosimetry_list
      DF <- reverse_dosimetry_list |> reshape2::melt(id.vars=reverse_dosimetry_list[[1]] |> colnames()) |> dplyr::mutate(Experiment=L1 |> stringr::str_replace(pattern="_[:alpha:]*$", replacement="")) |> dplyr::mutate(Experiment=factor(Experiment, levels=Experiment |> sort_alphnum()))

      Biological_system_annotations <- AOPfingerprintR::Biological_system_annotations |> as.data.frame()
      Biological_system_annotations$organ_tissue[is.na(Biological_system_annotations$organ_tissue)] <- "General"
      rownames(Biological_system_annotations) <- Biological_system_annotations$key_event_name

      DF <- DF |> dplyr::inner_join(Biological_system_annotations, by=c("key_event"="key_event_name"))
      DF <- DF |> dplyr::mutate(level=level |> factor(levels=c("Molecular","Cellular","Tissue","Organ","Individual")), pfas=Experiment |> as.vector() |> stringr::str_replace(pattern="_.*", replacement=""), time=time |> as.numeric())
      WriteXLS::WriteXLS(DF, ExcelFileName=file)
    }
  )

  # download reverse dosimetry summarized reults
  output$download_rev_dosimetry_summarized_tables <- downloadHandler(
    filename = "reverse_dosimetry_summarized_results_table.xlsx",
    content = function(file) {
      shiny::req(r$results_reverse_dosimetry)
      reverse_dosimetry_list_summarized <- r$results_reverse_dosimetry$reverse_dosimetry_list_summarized
      DF <- reverse_dosimetry_list_summarized |> reshape2::melt(id.vars=reverse_dosimetry_list_summarized[[1]] |> colnames()) |> dplyr::mutate(Experiment=L1 |> stringr::str_replace(pattern="_[:alpha:]*$", replacement="")) |> dplyr::mutate(Experiment=factor(Experiment, levels=Experiment |> sort_alphnum()))

      Biological_system_annotations <- AOPfingerprintR::Biological_system_annotations |> as.data.frame()
      Biological_system_annotations$organ_tissue[is.na(Biological_system_annotations$organ_tissue)] <- "General"
      rownames(Biological_system_annotations) <- Biological_system_annotations$key_event_name

      DF <- DF |> dplyr::inner_join(Biological_system_annotations, by=c("key_event"="key_event_name"))
      DF <- DF |> dplyr::mutate(level=level |> factor(levels=c("Molecular","Cellular","Tissue","Organ","Individual")), pfas=Experiment |> as.vector() |> stringr::str_replace(pattern="_.*", replacement=""), time=time |> as.numeric())
      WriteXLS::WriteXLS(DF, ExcelFileName=file)
    }
  )

  # render DT Reverse Dosimetry
  output$reverse_dosimetry_results_DT <- DT::renderDT({
    shiny::validate(need(!is.null(r$results_reverse_dosimetry), "No Reverse Dosimetry results!"))
   
    reverse_dosimetry_list_summarized <- r$results_reverse_dosimetry$reverse_dosimetry_list_summarized
    DF <- reverse_dosimetry_list_summarized |> reshape2::melt(id.vars=reverse_dosimetry_list_summarized[[1]] |> colnames()) |> dplyr::mutate(Experiment=L1 |> stringr::str_replace(pattern="_[:alpha:]*$", replacement="")) |> dplyr::mutate(Experiment=factor(Experiment, levels=Experiment |> sort_alphnum()))

    Biological_system_annotations <- AOPfingerprintR::Biological_system_annotations |> as.data.frame()
    Biological_system_annotations$organ_tissue[is.na(Biological_system_annotations$organ_tissue)] <- "General"
    rownames(Biological_system_annotations) <- Biological_system_annotations$key_event_name

    DF <- DF |> dplyr::inner_join(Biological_system_annotations, by=c("key_event"="key_event_name"))
    DF <- DF |> dplyr::mutate(level=level |> factor(levels=c("Molecular","Cellular","Tissue","Organ","Individual")), pfas=Experiment |> as.vector() |> stringr::str_replace(pattern="_.*", replacement=""), time=time |> as.numeric())
    
    DT::datatable(DF , filter="top",
      selection='single',
      options=list(
        search=list(regex=TRUE, caseInsensitive=FALSE),
        scrollX=TRUE,
        ordering=TRUE,
        rownames=FALSE
    ))
  },server=TRUE)

  # render Reverse Dosimetry density plot
  output$reverse_dosimetry_ridgeplot <- shiny::renderPlot({
    shiny::validate(need(!is.null(r$results_reverse_dosimetry), "No Reverse Dosimetry results!"))
    reverse_dosimetry_list <- r$results_reverse_dosimetry$reverse_dosimetry_list
    gp <- plot_reverse_dosimetry_density(reverse_dosimetry_list)

    print(gp)
  })

  # infer KE network Reverse Dosimetry
  observeEvent(input$net_submit_rev_dosimetry, {
    shiny::validate(
      need(!is.null(r$results_reverse_dosimetry), "No Reverse Dosimetry results!")
    )
    group_by <- "aop"
    reverse_dosimetry_list_summarized <- r$results_reverse_dosimetry$reverse_dosimetry_list_summarized

    enlarge_ke_selection <- input$enlarge_rev_dosimetry_keyevents |> as.logical()
    if(enlarge_ke_selection){
      shiny::req(input$max_path_length, input$n_AO, input$n_MIE, r$iv.rev_dosimetry.net$is_valid())
      max_path_length <- input$max_path_length |> as.numeric()
      n_AOs <- input$n_AO |> as.numeric()
      n_MIEs <- input$n_MIE |> as.numeric()
    }else{
      max_path_length <- NULL
      n_AOs <- NULL
      n_MIEs <- NULL
    }

    w$show()
    on.exit({
      w$hide()
    })

    mode <- "out"
    ke_id <- "ke"
    convert_to_gene_symbols <- FALSE

    tryCatch({
      net_list_res <- reverse_dosimetry_list_summarized |> purrr::map(\(x){
        get_ke_net_reverse_dosimetry(
          reverse_dosimetry_summarized = x,
          max_path_length = max_path_length,
          n_AOs = n_AOs,
          n_MIEs = n_MIEs,
          mode = mode,
          enlarge_ke_selection = enlarge_ke_selection,
          ke_id = ke_id,
          group_by = group_by
        )
      })

      r$net_list_rev_dosimetry <- net_list_res
    }, error=\(e){
      str <- paste0("Error Encountered: ", e)
      print(str)
      shiny::showNotification("Encountered Error: Please contact application maintainers. Check the console for the error information", duration=NULL)
    })
  })

  # render select experiment Reverse Dosimetry from net
  output$render_select_experiment_reverse_dosimetry_net <- renderUI({
    shiny::validate(need(!is.null(r$net_list_rev_dosimetry), "No Network data!"))
    net_list_res <- r$net_list_rev_dosimetry
    choices <- net_list_res |> names()
    selectInput(inputId="select_experiment_reverse_dosimetry_net", label="Select Experiment", choices=choices, selected=choices[1])
  })

  # render reverse dosimetry network plot
  output$reverse_dosimetry_ke_visnetwork <- visNetwork::renderVisNetwork({
    shiny::validate(
      need(!is.null(r$results_reverse_dosimetry), "No Reverse Dosimetry results!"),
      need(!is.null(input$select_experiment_reverse_dosimetry_net), "No Experiment selected!"),
      need(!is.null(r$net_list_rev_dosimetry), "No Network data!")
    )
    shiny::req(input$select_experiment_reverse_dosimetry_net)
    reverse_dosimetry_list_summarized <- r$results_reverse_dosimetry$reverse_dosimetry_list_summarized

    exp_id <- input$select_experiment_reverse_dosimetry_net
    net_list_res <- r$net_list_rev_dosimetry
    visnet_data <- net_list_res[[exp_id]][["net"]]
    visnet <- net_list_res[[exp_id]][["visnet"]]

    r$visnet_data_rev_dosimetry <- visnet_data
    r$visnet_rev_dosimetry <- visnet
    visnet
  })

  # download reverse dosimetry visNetwork data all experiments
  output$download_rev_dosimetry_visnet_data_all <- downloadHandler(
    filename = "reverse_dosimetry_visnetwork_data_all_experiments.xlsx",
    content = function(file) {
      shiny::req(r$net_list_rev_dosimetry)
      visNet_data.list <- r$net_list_rev_dosimetry |> purrr::map("net") |> purrr::list_flatten() 
      WriteXLS::WriteXLS(visNet_data.list, ExcelFileName=file)
    }
  )

  # download reverse dosimetry visNetwork on display
  output$download_rev_dosimetry_visnet_display <- downloadHandler(
    filename = "reverse_dosimetry_visnetwork.html",
    content = function(file) {
      shiny::req(r$visnet_rev_dosimetry)
      visNetwork::visSave(r$visnet_rev_dosimetry, file, selfcontained=TRUE, background="white")
    }
  )
  
  # download reverse dosimetry visNetwork data on display
  output$download_rev_dosimetry_visnet_data_display <- downloadHandler(
    filename = "reverse_dosimetry_visnetwork_data.xlsx",
    content = function(file) {
      shiny::req(r$visnet_data_rev_dosimetry)
      WriteXLS::WriteXLS(r$visnet_data_rev_dosimetry, ExcelFileName=file)
    }
  )

  observe({
    session$reload()
  }) |> shiny::bindEvent(input$reset)

  # Settings dialog
  observe({
    shiny::showModal(shiny::modalDialog(size="m",
      title = "Settings",
      bslib::layout_columns(col_widths=c(12), row_heights=c(3),
        bslib::card(
          bslib::card_header("Font Size"),
          bslib::card_body(fillable=FALSE,
            shiny::actionButton("decrease", "-", icon=NULL),
            shiny::actionButton("increase", "+", icon=NULL)
          ),
          bslib::card_body(
            "Settings"
          )
        )
      )
    ))
  }) |> shiny::bindEvent(input$options)

  # Scale font
  r$font_scale <- 1 
  the_theme <- bs_current_theme()
  #print("str(the_theme)")
  #print(str(the_theme))
  observe({
    ui_theme <- app_theme()
    #print("class(ui_theme)")
    #print(class(ui_theme))
    the_theme <- bs_current_theme()
    #the_theme <- session$getCurrentTheme()
    res <- getShinyOption("bootstrapTheme", default = NULL)
    #print("res")
    #print(res)
    #print("str(the_theme)")
    #print(str(the_theme))
    ##print("str(ui_theme)")
    #print(str(ui_theme))
    font_scale <- r$font_scale 
    font_scale <- min(font_scale-0.2,0.2) 
    ui_theme <- bslib::bs_theme_update(ui_theme, font_scale = font_scale)
    session$setCurrentTheme(ui_theme)
    r$font_scale <- font_scale
  }) |> shiny::bindEvent(input$decrease)

  observe({
    ui_theme <- bslib::bs_current_theme()
    font_scale <- r$font_scale 
    font_scale <- min(font_scale+0.2,2) 
    ui_theme <- bslib::bs_theme_update(ui_theme, font_scale = font_scale)
    session$setCurrentTheme(ui_theme)
    r$font_scale <- font_scale
  }) |> shiny::bindEvent(input$increase)

  # remove temp files on session end
  session$onSessionEnded(function() {
    tempFiles <- isolate(r$tempFiles)
    if(!is.null(tempFiles)){
      #print("str(tempFiles)")
      #print(str(tempFiles))
      unlink(tempFiles)
    }
  })
}

