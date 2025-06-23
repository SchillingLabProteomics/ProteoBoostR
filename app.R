library(shiny)
library(shinydashboard)
library(shinyjs)
library(DT)
library(dplyr)
library(tibble)
library(caret)
library(xgboost)
library(pROC)
library(rBayesianOptimization)
library(ggplot2)
library(fontawesome)

preprocessData <- function(df, annotationColumn, neg_label, pos_label) {
  # filter out rows with NA in annotationColumn
  df <- df[!is.na(df[[annotationColumn]]), ]
  
  # keep only rows where annotationColumn is either pos_label or neg_label
  df <- df[df[[annotationColumn]] %in% c(neg_label, pos_label), ]
  
  # convert annotationColumn to a factor with specified levels
  df[[annotationColumn]] <- factor(df[[annotationColumn]], levels = c(neg_label, pos_label))
  
  return(df)
}

resolveOutputDir <- function(userPath) {
  configFile <- "ProteoBoostR.filesystem"
  if (file.exists(configFile)) {
    # Docker mode: enforce folder name (no slashes)
    if (grepl("[/\\\\]", userPath)) {
      return(NA)
    }
    systemRoot <- trimws(readLines(configFile, n = 1))
    return(file.path(systemRoot, userPath))
  } else {
    # local mode: use the user input as provided
    return(userPath)
  }
}

ui <- dashboardPage(
    dashboardHeader(
        title = tags$span(
            tags$img(src = "ProteoBoostR_logo.png",
                   height = 40),
            "ProteoBoostR"
            )
    ),
    dashboardSidebar(
        useShinyjs(),
        tags$style(HTML("
            .shiny-input-container {margin-bottom: 5px;}
        ")),
        tags$head(tags$style(HTML("
            .tab-disabled { pointer-events: none; opacity: 0.5; }
            .valid-path { color: green; font-weight: bold; }
            .invalid-path { color: red; font-weight: bold; }

            /* Font body content */
            body { font-family: inherit !important; }

            /* Font header title */
            .skin-blue .main-header .logo {
              background-color: #0f3041 !important;
              color: #fbfbf8;
              border-bottom: 0;
              font-family: 'Britannic Bold', serif !important;
            }

            .skin-blue .main-header .navbar {
              background-color: #0f3041 !important;
            }

            /* Sidebar styling */
            .skin-blue .main-sidebar {
              background-color: #fbfbf8 !important;
            }

            .skin-blue .main-sidebar .sidebar .sidebar-menu a {
              color: black !important;
            }

            /* Active tab styling */
            .skin-blue .main-sidebar .sidebar .sidebar-menu .active a {
              background-color: #f05c42 !important;
            }

            /* Hover effect for sidebar menu items */
            .skin-blue .main-sidebar .sidebar .sidebar-menu .li a:hover {
              background-color: #f05c42 !important;
              border-left-color: #ffc133 !important;
            }

            .control-sidebar-light

            /* Content background */
            .content-wrapper, .right-side {
              background-color: #fbfbf8;
            }

            /* Change tab highlighting */
            .skin-blue .main-sidebar .sidebar .sidebar-menu .active > a {
              border-left-color: #ffc133 !important;
            }

            /* Change action button hover state */
            .btn-primary:hover, .btn-primary:active, .btn-primary:focus {
              background-color: #f05c42 !important;
              border-color: #ffc133 !important;
            }

            .btn-info:hover, .btn-info:active, .btn-info:focus {
              background-color: #f05c42 !important;
              border-color: #ffc133 !important;
            }
        "))),
        sidebarMenu(id = "tabs",
            menuItem("Landing", tabName = "landing_tab", icon = icon("home")),
            menuItem("Input", tabName = "input_tab", icon = icon("upload")),
            menuItem("Model Training", tabName = "train_tab", icon = icon("cogs")),
            menuItem("Model Testing", tabName = "test_tab", icon = icon("play")),
            menuItem("Performance Visualization", tabName = "perf_tab", icon = icon("chart-line")),
            menuItem("Log", tabName = "log_tab", icon = icon("clipboard-list"))
        )
    ),
    dashboardBody(
        tabItems(
            # landing Page
            tabItem(tabName = "landing_tab",
                fluidRow(
                    box(title = "Welcome to ProteoBoostR (v1.0.0)", status = "primary", width = 12,
                        br(),
                        p("ProteoBoostR is a Shiny-based tool for supervised classification (e.g. biomarker identification) in proteomics data. It leverages the powerful XGBoost algorithm combined with Bayesian optimization to train and evaluate predictive models. The tool automatically merges proteomics expression data with sample annotations, performs data preprocessing, and outputs key files for further analysis."),
                        br(),
                        p("Quick Start (Defaults)"),
                        tags$ol(
                          tags$li("Upload your annotation and protein matrix files."),
                          tags$li("Set the annotation column and class labels."),
                          tags$li("(Optionally) Provide a protein subset."),
                          tags$li("Specify your output directory."),
                          tags$li("Click", tags$b("Continue to Model Training"), "to process data and train the model."),
                          tags$li("Go to", tags$b("Continue to Model Testing"), "to evaluate the model and automatically save all test outputs."),
                          tags$li("Review the ROC curve in", tags$b("Performance Visualization"), "and check logs in", tags$b("Log"), "."),
                        ),
                        br(),
                        actionButton("goLandingContinue", "Let's go!", icon = icon("arrow-right"))
                    )
                )
            ),
            # input Tab
            tabItem(tabName = "input_tab",
                fluidRow(
                    box(title = "Upload Files (TSV Only)", status = "primary", width = 6,
                        fileInput("annotationFile", "Annotation File (.tsv only)"),
                        helpText("Annotation must be a TSV file (first column is sample_id).", style = "margin-top: -20px"),
                        br(),
                        fileInput("proteinFile", "Protein Matrix (.tsv only)"),
                        helpText("Protein matrix with rows = protein IDs and columns = sample IDs.", style = "margin-top: -20px"),
                        br(),
                        textInput("outputDir", "Output Directory:", value = ""),
                        helpText("Enter folder name only when running in Docker."),
                        uiOutput("outputDirStatus")
                    ),
                    box(title = "Annotation Column Selection & Class Labels", status = "primary", width = 6,
                        uiOutput("annotationColUI"),
                        helpText("Pick the column for class label (sample_id excludfed)."),
                        br(),
                        uiOutput("classSelectionUI"),
                        helpText("Define negative (0) and positive (1) classes.")
                    )
                ),
                fluidRow(
                    box(title = "Protein Subset (Optional)", status = "info", width = 6,
                        fileInput("subsetFile", "Upload Protein Subset IDs (.txt or .tsv, no header)"),
                        helpText("If uploaded, the text area is disabled. One protein ID per line.", style = "margin-top: -20px"),
                        br(),
                        textAreaInput("proteinSubset", "List of Proteins (one per line):", rows = 5),
                        helpText("Leave empty to use all proteins.")
                    ),
                    box(title = "Train/Test Split", status = "info", width = 6,
                        sliderInput("trainSplit", "Train/Test Split (p):", min = 0, max = 1, step = 0.05, value = 0.70),
                        helpText("These settings determine how the data is split before training.")
                    )
                ),
                fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("inputContinueButton", "Continue to Model Training", icon = icon("arrow-right")),
                        actionButton("trainToTestButton", "Jump to Model Testing", icon = icon("arrow-right"))
                    )
                )
            ),
            # model Training Tab
            tabItem(tabName = "train_tab",
                fluidRow(
                    box(title = "Bayesian Optimization Ranges (General)", status = "primary",
                        helpText("General parameters: controls learning rate, tree depth, sampling."),
                        sliderInput("eta_range", "eta Range:", min = 0, max = 1, value = c(0.001, 0.2), step = 0.001),
                        helpText("Controls how much the model updates at each boosting step. A lower eta slows learning but allows for better convergence."),
                        helpText("Use a higher eta (0.1 - 0.3) for small datasets to avoid underfitting. For large datasets, use a lower eta (0.01 - 0.1) with more boosting rounds to improve performance."),
                        sliderInput("md_range", "max_depth Range:", min = 0, max = 20, value = c(2,4), step = 1),
                        helpText("Determines the depth of each decision tree. A larger depth captures more complexity but increases overfitting risk."),
                        helpText("For small datasets, keep it low (2-4) to avoid overfitting. For large datasets, experiment with higher values (4-10) if needed."),
                        sliderInput("subsample_range", "subsample Range:", min = 0, max = 1, value = c(0.8, 1.0),  step = 0.1),
                        helpText("Specifies the fraction of training samples used to grow each tree, introducing randomness and reducing overfitting."),
                        helpText("Use lower values (0.6 - 0.9) for small datasets to prevent overfitting. For very large datasets, subsample ~0.8 can speed up training without losing performance."),
                        sliderInput("colsample_range", "colsample_bytree Range:", min = 0, max = 1, value = c(0.8, 1.0),  step = 0.1),
                        helpText("Controls how many features are randomly sampled when growing a tree. Helps prevent reliance on specific features."),
                        helpText("If you have few features (<10), keep it high (0.8 - 1.0). If you have many features (>50), tune it lower (0.5 - 0.8) to avoid overfitting.")
                    ),
                    box(title = "Bayesian Optimization Ranges (Regularization)", status = "primary",
                        helpText("Regularization parameters: regulate model complexity & reduce overfitting."),
                        sliderInput("child_weight_range", "min_child_weight Range:", min = 0, max = 20, value = c(2,10), step = 1),
                        helpText("Sets the minimum sum of instance weights needed to make a split. Higher values make the model conservative."),
                        helpText("Use higher values (3-10) for small datasets to avoid overfitting. For large datasets, lower values (1-3) allow more splits for better learning."),
                        sliderInput("gamma_range", "gamma Range:", min = 0, max = 20, value = c(0, 2), step = 0.1),
                        helpText("Forces tree splits to have a significant loss reduction, controlling overfitting."),
                        helpText("Use higher values (≥1) for small datasets to avoid unnecessary splits. In large datasets, tune it from 0 to 5 based on performance."),
                        numericInput("alpha_min", "alpha (min)", 0),
                        numericInput("alpha_max", "alpha (max)", 1),
                        helpText("Encourages sparsity in feature selection by penalizing the absolute values of leaf weights, helping with feature selection."),
                        helpText("If you have many features, use higher values (1-5) to eliminate irrelevant ones. If you have few features, keep it low (0-1)."),
                        numericInput("lambda_min", "lambda (min)", 1),
                        numericInput("lambda_max", "lambda (max)", 10),
                        helpText("Penalizes large weight values, helping to smooth the model and reduce overfitting."),
                        helpText("Use higher values (1-10) for small datasets to control overfitting. For large datasets, set it moderate (1-5) unless overfitting is observed.")
                    )
                ),
                fluidRow(
                    box(title = "Best Hyperparameters", status = "primary", width = 12,
                         dataTableOutput("bestParamsTable")
                    )
                ),
                fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("trainContinueButton", "Train the Model", icon = icon("gears")),
                        actionButton("testContinueButton", "Continue to Model Testing", icon = icon("arrow-right"))
                    )
                )
            ),
            # model testing tab
            tabItem(tabName = "test_tab",
                fluidRow(
                    box(title = "Test Model", status = "primary", width = 12,
                        fileInput("pretrainedModel", "Pretrained Model (.rds)"),
                        helpText("Disabled if an in‑session model exists.", style = "margin-top: -20px"),
                        br(),
                        actionButton("evaluateContinueButton", "Continue to Evaluate", icon = icon("play"))
                    )
                ),
                fluidRow(
                    box(title = "Confusion Matrix", status = "info", width = 6,
                        verbatimTextOutput("confusionMatrix")
                    ),
                    box(title = "Evaluation Metrics", status = "info", width = 6,
                        verbatimTextOutput("testResults")
                    )
                ),
                fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("testToPerfButton", "Continue to Visualization", icon = icon("arrow-right"))
                    )
                )
            ),
            # performance visualization tab
            tabItem(tabName = "perf_tab",
            fluidRow(
                    box(title = "ROC Curve", status = "primary", width = 8,
                        plotOutput("rocPlot", height = "400px")
                    )
                ),
                fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("perfToLogButton", "Continue to Log", icon = icon("arrow-right"))
                    )
                )
            ),
            # log tab
            tabItem(tabName = "log_tab",
                fluidRow(
                    box(title = "Log Output", status = "primary", width = 12,
                        verbatimTextOutput("logOutput")
                    )
                ),
                fluidRow(
                    box(title = "Session Info", status = "primary", width = 12,
                        verbatimTextOutput("sessionInfo")
                    )
                )
            )
        )
    ),
    title = "ProteoBoostR"
)

server <- function(input, output, session) {

    rv <- reactiveValues(
        annot = NULL,
        proteins = NULL,
        mergedData = NULL,
        trainData = NULL,
        testData = NULL,
        xgbModel = NULL,
        bestParams = NULL,
        confMat = NULL,
        rocObj = NULL,
        bestThresh = NULL,
        pred_probs = NULL,
        logs = "",
        trainingInProgress = FALSE,
        subsetFileIDs = character(0),
        session_ts = format(Sys.time(), "%Y%m%d%H%M%S"),
        logFile = NULL,
        logSet = FALSE
    )

    # appends messages to the reactive log and to the log file if set
    appendLog <- function(msg) {
        timestamped <- paste0("[", Sys.time(), "] ", msg, "\n")
        rv$logs <- paste0(rv$logs, timestamped)
        if (!is.null(rv$logFile)) {
            cat(timestamped, file = rv$logFile, append = TRUE)
        }
    }

    # set log file when outputDir is valid and user presses Continue on Input tab
    observeEvent(input$outputDir, {
        req(input$outputDir)
        resolvedPath <- resolveOutputDir(input$outputDir)
        if (!is.na(resolvedPath) && dir.exists(resolvedPath) && !rv$logSet) {
            rv$logFile <- file.path(resolvedPath, paste0("ProteoBoostR_", rv$session_ts, ".log"))
            writeLines(rv$logs, con = rv$logFile)
            rv$logSet <- TRUE
        }
    })

    # landing page -> go to input tab
    observeEvent(input$goLandingContinue, {
        updateTabItems(session, "tabs", "input_tab")
    })

    # live check of outputDir
    # set log file only once when valid and user presses Continue on input tab
    output$outputDirStatus <- renderUI({
        req(input$outputDir)
        resolvedPath <- resolveOutputDir(input$outputDir)
        if (is.na(resolvedPath)) {
            span(icon("times"), "Invalid Output Path: enter folder name only.", class = "invalid-path")
        } else if (dir.exists(resolvedPath)) {
            span(icon("check"), paste("Valid Output Path:", resolvedPath), class = "valid-path")
        } else {
            span(icon("times"), paste("Invalid Output Path:", resolvedPath), class = "invalid-path")
        }
    })

    # subset file upload - disable text area
    observeEvent(input$subsetFile, {
        if (!is.null(input$subsetFile)) {
            lines <- readLines(input$subsetFile$datapath, warn = FALSE)
            lines <- trimws(lines)
            # remove quotes and empty lines
            lines <- gsub("\"", "", lines)
            lines <- gsub("'", "", lines)
            lines <- lines[lines != ""]
            rv$subsetFileIDs <- lines
            appendLog(paste("Subset file uploaded with", length(lines), "IDs."))
            shinyjs::disable("proteinSubset")
        } else {
            shinyjs::enable("proteinSubset")
        }
    })

    # disable pretrained model input if an in‑session model exists
    observe({
        if (!is.null(rv$xgbModel)) {
            shinyjs::disable("pretrainedModel")
        } else {
            shinyjs::enable("pretrainedModel")
        }
    })

    # tab enabling logic
    disableAllExceptLog <- function() {
        shinyjs::addClass(selector = "a[data-value='input_tab']", class = "tab-disabled")
        shinyjs::addClass(selector = "a[data-value='train_tab']", class = "tab-disabled")
        shinyjs::addClass(selector = "a[data-value='test_tab']", class = "tab-disabled")
        shinyjs::addClass(selector = "a[data-value='perf_tab']", class = "tab-disabled")
    }

    enableTabsBasedOnInputs <- function() {
        validInput <- !is.null(rv$annot) && !is.null(rv$proteins) && nchar(input$outputDir) > 0 && dir.exists(input$outputDir)
        if (validInput) {
            shinyjs::removeClass(selector = "a[data-value='train_tab']", class = "tab-disabled")
            shinyjs::removeClass(selector = "a[data-value='test_tab']", class = "tab-disabled")
            shinyjs::removeClass(selector = "a[data-value='perf_tab']", class = "tab-disabled")
        } else {
            shinyjs::addClass(selector = "a[data-value='train_tab']", class = "tab-disabled")
            shinyjs::addClass(selector = "a[data-value='test_tab']", class = "tab-disabled")
            shinyjs::addClass(selector = "a[data-value='perf_tab']", class = "tab-disabled")
        }
    }
    observe({
        if (rv$trainingInProgress) {
            disableAllExceptLog()
        } else {
            enableTabsBasedOnInputs()
        }
    })

    doMerge <- function() {
        appendLog("Auto-merging annotation + protein data...")
        prot_t <- t(rv$proteins)
        prot_t <- as.data.frame(prot_t, stringsAsFactors = FALSE)
        prot_t <- tibble::rownames_to_column(prot_t, "sample_id")
        df_merged <- merge(prot_t, rv$annot, by = "sample_id")
        new_names <- gsub(";.*", "", names(df_merged))
        colnames(df_merged) <- new_names

        # apply subset if provided (from file or text)
        userSubsetLines <- character(0)
        if (length(rv$subsetFileIDs) > 0) {
            userSubsetLines <- rv$subsetFileIDs
            appendLog(paste("Applying subset from file with", length(userSubsetLines), "IDs."))
        } else {
            subsetText <- input$proteinSubset
            if (nchar(subsetText) > 0) {
                lines <- unlist(strsplit(subsetText, "\n"))
                lines <- trimws(lines)
                lines <- lines[lines != ""]
                userSubsetLines <- lines
                if (length(userSubsetLines) > 0) {
                    appendLog(paste("Auto-applying subset with:", length(userSubsetLines), "proteins from text area."))
                } else {
                    appendLog("No subset typed in text area (lines are empty?).")
                }
            }
        }

        if (length(userSubsetLines) > 0) {
            keepCols <- unique(c("sample_id", userSubsetLines, input$annotationColumn))
            # only keep columns present in the merged data
            df_merged <- df_merged[, intersect(keepCols, colnames(df_merged)), drop = FALSE]
        } else {
            appendLog("No subset typed or file, using all proteins.")
        }

        # force numeric conversion for columns except 'sample_id' and annotationColumn
        keepColsAlways <- c("sample_id", input$annotationColumn)
        for (cn in colnames(df_merged)) {
            if (!(cn %in% keepColsAlways)) {
                tryNum <- suppressWarnings(as.numeric(df_merged[[cn]]))
                if (all(is.na(tryNum))) {
                    df_merged[[cn]] <- NULL
                    appendLog(paste("Dropping column:", cn, "not numeric or all NA."))
                } else {
                    df_merged[[cn]] <- tryNum
                }
            }
        }
        appendLog(paste("Final columns after numeric conversion:", paste(colnames(df_merged), collapse = ", ")))
        rv$mergedData <- df_merged
        appendLog("Merging complete.")
    }

    # merge annotation and protein if both exist
    attemptMerge <- function() {
        if (!is.null(rv$annot) && !is.null(rv$proteins)) {
            doMerge()
        }
    }

    # load annotation
    observeEvent(input$annotationFile, {
        req(input$annotationFile)
        if (!grepl("\\.tsv$", input$annotationFile$name, ignore.case = TRUE)) {
            appendLog("Error: Annotation file must be a .tsv!")
            return(NULL)
        }
        df <- tryCatch({
            raw <- read.delim(input$annotationFile$datapath, header = TRUE, stringsAsFactors = FALSE)
            colnames(raw)[1] <- "sample_id"
            raw
        }, error = function(e) {
            appendLog(paste("Error reading annotation file:", e$message))
            NULL
        })
        if (!is.null(df)) {
            rv$annot <- df
            appendLog(paste("Annotation file loaded. Rows:", nrow(df), "Cols:", ncol(df)))
            attemptMerge()
        } else {
            appendLog("Annotation file read failed or is empty.")
        }
    })

    # trigger merge if both annotation and protein data are loaded
    observeEvent(list(input$inputContinueButton, input$trainToTestButton), {
        if (is.null(rv$proteins)) {
            req(input$proteinFile)
            req(input$annotationFile)
            if (!grepl("\\.tsv$", input$proteinFile$name, ignore.case = TRUE)) {
                appendLog("Error: Protein file must be a .tsv!")
                return(NULL)
            }
            df <- tryCatch({
                read.delim(input$proteinFile$datapath, header = TRUE, row.names = 1, stringsAsFactors = FALSE)
            }, error = function(e) {
                appendLog(paste("Error reading protein matrix:", e$message))
                NULL
            })
            if (!is.null(df)) {
                rv$proteins <- df
                appendLog(paste("Protein matrix loaded. Dimensions:", paste(dim(df), collapse = " x ")))
                attemptMerge()
            } else {
                appendLog("Protein matrix read failed or is empty.")
            }
        }
    })

    # render Annotation Column UI (exclude sample_id)
    output$annotationColUI <- renderUI({
        req(rv$annot)
        cols <- setdiff(colnames(rv$annot), "sample_id")
        selectInput("annotationColumn", "Annotation Column:", choices = cols)
    })

    # render Class Labels UI
    output$classSelectionUI <- renderUI({
        req(rv$annot, input$annotationColumn)
        colName <- input$annotationColumn
        uniqueVals <- sort(unique(rv$annot[[colName]]))
        if (length(uniqueVals) < 2) {
            return(tags$p("Fewer than 2 unique values in this column. Please select a different column."))
        }
        tagList(
            selectInput("negClass", "Negative Class (0)", choices = uniqueVals),
            selectInput("posClass", "Positive Class (1)", choices = uniqueVals, selected = uniqueVals[length(uniqueVals)])
        )
    })

    # input tab: Continue button – partition data and write subsetted train/test matrices
    observeEvent(input$inputContinueButton, {
        set.seed(1234)
        req(rv$mergedData, input$annotationColumn, input$negClass, input$posClass)
        appendLog("User clicked Continue on Input Tab.")

        # process data, adding the annotationColumn column
        df <- preprocessData(rv$mergedData, input$annotationColumn, input$negClass, input$posClass)

        # split data into training and testing sets
        # if trainSplit is 0, use the entire dataset for both training and testing
        if (input$trainSplit == 0) {
            # empty trainDF
            trainDF <- df[0, ]
            testDF <- df
        } else {
            inTrain <- createDataPartition(df[[input$annotationColumn]], p = input$trainSplit, list = FALSE)
            trainDF <- df[inTrain, ]
            testDF  <- df[-inTrain, ]
        }

        rv$trainData <- trainDF
        rv$testData  <- testDF

        # if a subset is provided, filter trainDF and testDF
        subsetIDs <- character(0)
        if (length(rv$subsetFileIDs) > 0) {
            subsetIDs <- rv$subsetFileIDs
        } else if (nchar(input$proteinSubset) > 0) {
            subsetIDs <- trimws(unlist(strsplit(input$proteinSubset, "\n")))
            subsetIDs <- subsetIDs[subsetIDs != ""]
        }
        if (length(subsetIDs) > 0) {
            trainDF <- trainDF[, c("sample_id", intersect(colnames(trainDF), subsetIDs), input$annotationColumn), drop = FALSE]
            testDF <- testDF[, c("sample_id", intersect(colnames(testDF), subsetIDs), input$annotationColumn), drop = FALSE]
            rv$trainData <- trainDF
            rv$testData <- testDF
            appendLog("Data subset applied to training and testing matrices.")
        }
        # write the training and testing matrices (transposed)
        resolvedOutDir <- resolveOutputDir(input$outputDir)
        if (is.na(resolvedOutDir)) return()  # Stop if output path is invalid
        if (nchar(input$outputDir) > 0 && dir.exists(resolvedOutDir)) {
            fn_ts <- rv$session_ts
            trainMatPath <- file.path(resolvedOutDir, paste0("train_matrix_", fn_ts, ".tsv"))
            trainDF <- trainDF %>%
                        remove_rownames() %>%
                        column_to_rownames("sample_id")
            write.table(trainDF, file = trainMatPath, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            testDF <- testDF %>%
                        remove_rownames() %>%
                        column_to_rownames("sample_id")
            testMatPath <- file.path(resolvedOutDir, paste0("test_matrix_", fn_ts, ".tsv"))
            write.table(testDF, file = testMatPath, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            appendLog("Training and testing matrices written (transposed).")
        }
        updateTabItems(session, "tabs", "train_tab")
    })

    # model training: Continue to train button
    observeEvent(input$trainContinueButton, {
        set.seed(1234)
        appendLog("User clicked Continue on Model Training. Starting Bayesian Optimization and training model...")

        appendLog(sprintf(
              paste(
                "Bayesian-opt bounds:",
                "eta=[%.3f, %.3f], max_depth=[%d, %d], subsample=[%.2f, %.2f],",
                "colsample_bytree=[%.2f, %.2f], min_child_weight=[%d, %d],",
                "gamma=[%.2f, %.2f], alpha=[%.2f, %.2f], lambda=[%d, %d]"
              ),
              input$eta_range[1], input$eta_range[2],
              input$md_range[1], input$md_range[2],
              input$subsample_range[1], input$subsample_range[2],
              input$colsample_range[1], input$colsample_range[2],
              input$child_weight_range[1], input$child_weight_range[2],
              input$gamma_range[1], input$gamma_range[2],
              input$alpha_min, input$alpha_max,
              input$lambda_min, input$lambda_max
            ))

        rv$trainingInProgress <- TRUE
        nrounds <- 1
        withProgress(message = "Training XGBoost model...", value = 0, {
            incProgress(0.1, detail = "Using training data...")
            req(rv$trainData)

            trainDF <- rv$trainData
            md_min <- input$md_range[1]
            md_max <- input$md_range[2]
            cw_min <- input$child_weight_range[1]
            cw_max <- input$child_weight_range[2]
            lam_min <- input$lambda_min
            lam_max <- input$lambda_max

            xgb_cv_bayes <- function(eta, max_depth, subsample, colsample_bytree,
                                       min_child_weight, gamma, alpha, lambda) {
                features <- trainDF[, setdiff(colnames(trainDF), c("sample_id", input$annotationColumn)), drop = FALSE]
                label <- as.numeric(trainDF[[input$annotationColumn]]) - 1
                dtrain <- xgb.DMatrix(as.matrix(features), label = label)

                    cv_result <- tryCatch({
                        xgb.cv(
                            params = list(
                                booster = "gbtree",
                                objective = "binary:logistic",
                                eval_metric = "auc",
                                eta = eta,
                                max_depth = max_depth,
                                subsample = subsample,
                                colsample_bytree = colsample_bytree,
                                min_child_weight = min_child_weight,
                                gamma = gamma,
                                alpha = alpha,
                                lambda = lambda
                            ),
                            data = dtrain,
                            nrounds = 1000,
                            nfold = 5,
                            early_stopping_rounds = 50,
                            verbose = 0,
                            stratified = TRUE
                        )
                    }, error = function(e) {
                      appendLog(paste("ERROR inside xgb.cv:", e$message))
                    })

                if (is.null(cv_result)) {
                    appendLog("ERROR: cv_result is NULL.")
                    return(list(Score = 0.5, Pred = 0))
                } else {
                    bestAUC <- max(cv_result$evaluation_log$test_auc_mean, na.rm = TRUE)
                    nrounds <- nrounds + 1
                    appendLog(paste("xgb_cv_bayes => rounds:", nrounds,
                                    "eta:", eta,
                                    "max_depth:", max_depth,
                                    "subsample:", subsample,
                                    "colsample_bytree:", colsample_bytree,
                                    "min_child_weight:", min_child_weight,
                                    "gamma:", gamma,
                                    "alpha:", alpha,
                                    "lambda:", lambda,
                                    "AUC:", bestAUC))
                    return(list(Score = bestAUC, Pred = 0))
                }
            }

            bo_result <- BayesianOptimization(
                FUN = xgb_cv_bayes,
                bounds = list(
                    eta = c(as.double(input$eta_range[1]), as.double(input$eta_range[2])),
                    max_depth = c(as.integer(md_min), as.integer(md_max)),
                    subsample = c(as.double(input$subsample_range[1]), as.double(input$subsample_range[2])),
                    colsample_bytree = c(as.double(input$colsample_range[1]), as.double(input$colsample_range[2])),
                    min_child_weight = c(as.integer(cw_min), as.integer(cw_max)),
                    gamma = c(as.double(input$gamma_range[1]), as.double(input$gamma_range[2])),
                    alpha = c(as.double(input$alpha_min), as.double(input$alpha_max)),
                    lambda = c(as.double(lam_min), as.double(lam_max))
                ),
                init_points = 5,
                n_iter = 20,
                acq = "ucb",
                kappa = 2.576,
                verbose = TRUE
            )
            best_par <- bo_result$Best_Par
            best_val <- bo_result$Best_Value
            appendLog(paste("Bayesian Optimization done. Best AUC =", round(best_val, 3)))
            incProgress(0.7, detail = "Training final model...")
            rv$bestParams <- list(
                booster = "gbtree",
                objective = "binary:logistic",
                eval_metric = "auc",
                eta = best_par["eta"],
                max_depth = as.integer(best_par["max_depth"]),
                subsample = best_par["subsample"],
                colsample_bytree = best_par["colsample_bytree"],
                min_child_weight = as.integer(best_par["min_child_weight"]),
                gamma = best_par["gamma"],
                alpha = best_par["alpha"],
                lambda_val = best_par["lambda"]
            )
            featCols <- setdiff(colnames(trainDF), c("sample_id", input$annotationColumn))
            features <- trainDF[, featCols, drop = FALSE]
            label <- as.numeric(trainDF[[input$annotationColumn]]) - 1
            appendLog(paste("Final training columns used for XGB:", paste(featCols, collapse = ", ")))
            dtrain <- xgb.DMatrix(as.matrix(features), label = label)
            model <- xgb.train(params = rv$bestParams, data = dtrain, nrounds = 1000, verbose = 0)
            rv$xgbModel <- model
            resolvedOutDir <- resolveOutputDir(input$outputDir)
            if (is.na(resolvedOutDir)) return()  # Stop if output path is invalid
            if (nchar(input$outputDir) > 0 && dir.exists(resolvedOutDir)) {
                fn_ts <- rv$session_ts
                savePath <- file.path(resolvedOutDir, paste0("xgb_model_", fn_ts, ".rds"))
                saveRDS(model, savePath)
                appendLog(paste("Model saved to:", savePath))
                bestParamsDF <- data.frame(Param = names(rv$bestParams), Value = unlist(rv$bestParams))
                bestParamPath <- file.path(resolvedOutDir, paste0("best_params_", fn_ts, ".tsv"))
                write.table(bestParamsDF, file = bestParamPath, sep = "\t", row.names = FALSE, quote = FALSE)
                appendLog("Best hyperparameters written.")
            }
            incProgress(1.0, detail = "Done")
            appendLog("Final XGBoost model trained.")
        })
        rv$trainingInProgress <- FALSE
    })

    observeEvent(input$testContinueButton, {
        updateTabItems(session, "tabs", "test_tab")
    })

    # continue to model testing tab
    observeEvent(input$trainToTestButton, {
        set.seed(1234)
        req(rv$mergedData, input$annotationColumn, input$negClass, input$posClass)
        appendLog("User clicked 'Jump to Model Testing' on Input Tab.")

        # process data, adding the annotationColumn column
        df <- preprocessData(rv$mergedData, input$annotationColumn, input$negClass, input$posClass)

        # split data into training and testing sets
        # if trainSplit is 0, use the entire dataset for both training and testing
        if (input$trainSplit == 0) {
            # empty trainDF
            trainDF <- df[0, ]
            testDF <- df
        } else {
            inTrain <- createDataPartition(df[[input$annotationColumn]], p = input$trainSplit, list = FALSE)
            trainDF <- df[inTrain, ]
            testDF  <- df[-inTrain, ]
        }

        rv$trainData <- trainDF
        rv$testData  <- testDF

        # if a subset is provided, filter trainDF and testDF
        subsetIDs <- character(0)
        if (length(rv$subsetFileIDs) > 0) {
            subsetIDs <- rv$subsetFileIDs
        } else if (nchar(input$proteinSubset) > 0) {
            subsetIDs <- trimws(unlist(strsplit(input$proteinSubset, "\n")))
            subsetIDs <- subsetIDs[subsetIDs != ""]
        }
        if (length(subsetIDs) > 0) {
            trainDF <- trainDF[, c("sample_id", intersect(colnames(trainDF), subsetIDs), input$annotationColumn), drop = FALSE]
            testDF <- testDF[, c("sample_id", intersect(colnames(testDF), subsetIDs), input$annotationColumn), drop = FALSE]
            rv$trainData <- trainDF
            rv$testData <- testDF
            appendLog("Data subset applied to training and testing matrices.")
        }
        # write the training and testing matrices (transposed)
        resolvedOutDir <- resolveOutputDir(input$outputDir)
        if (is.na(resolvedOutDir)) return()  # Stop if output path is invalid
        if (nchar(input$outputDir) > 0 && dir.exists(resolvedOutDir)) {
            fn_ts <- rv$session_ts
            trainMatPath <- file.path(resolvedOutDir, paste0("train_matrix_", fn_ts, ".tsv"))
            trainDF <- trainDF %>%
                        remove_rownames() %>%
                        column_to_rownames("sample_id")
            write.table(trainDF, file = trainMatPath, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            testDF <- testDF %>%
                        remove_rownames() %>%
                        column_to_rownames("sample_id")
            testMatPath <- file.path(resolvedOutDir, paste0("test_matrix_", fn_ts, ".tsv"))
            write.table(testDF, file = testMatPath, sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
            appendLog("Training and testing matrices written (transposed).")
        }
        updateTabItems(session, "tabs", "test_tab")
    })

    # model testing: Continue to Evaluate button
    observeEvent(input$evaluateContinueButton, {
        set.seed(1234)
        if (is.null(rv$xgbModel) && !is.null(input$pretrainedModel)) {
            appendLog("No in-session model found. Attempting to load pretrained model from file.")
            pretrainedPath <- input$pretrainedModel$datapath
            loadedModel <- tryCatch({
                readRDS(pretrainedPath)
            }, error = function(e) {
                appendLog(paste("ERROR loading pretrained model:", e$message))
                shinyjs::alert("Error loading pretrained model, see log.")
                return(NULL)
            })
            if (!is.null(loadedModel)) {
                rv$xgbModel <- loadedModel
                feature_names <- loadedModel$feature_names

                if (!is.null(rv$testData) && nrow(rv$testData) > 0) {
                    # select all feature names of model and if not available in testdata fill with NAs
                    missing_features <- setdiff(feature_names, colnames(rv$testData))
                    for (mf in missing_features) {
                      rv$testData[[mf]] <- NA
                    }

                    # sort rv$testData by feature_names
                    rv$testData <- rv$testData[, c("sample_id", feature_names, input$annotationColumn), drop = FALSE]
                } else {
                    appendLog("No test data available to match feature names.")
                }

              appendLog(paste("Loaded pretrained model from:", pretrainedPath))
            }
        }

        if (is.null(rv$testData)) {
            if (!is.null(rv$mergedData)) {
                rv$testData <- preprocessData(rv$mergedData, input$annotationColumn, input$negClass, input$posClass)
                appendLog(paste("No test data partition found; using entire merged data (", nrow(rv$testData), " rows) as test data."))
            } else {
                appendLog("No test data available.")
                shinyjs::alert("No test data available. Check your input files.")
                return(NULL)
            }
        }
        req(rv$xgbModel, rv$testData)
        appendLog("Evaluating model on test data...")
        testDF <- rv$testData
        featCols <- setdiff(colnames(testDF), c("sample_id", input$annotationColumn))
        features <- testDF[, featCols, drop = FALSE]
        label <- as.numeric(testDF[[input$annotationColumn]]) - 1
        appendLog(paste("Test model columns used for XGB:", paste(featCols, collapse = ", ")))
        dtest <- xgb.DMatrix(as.matrix(features), label = label)
        pred_probs <- predict(rv$xgbModel, dtest)
        rv$pred_probs <- pred_probs
        if (length(pred_probs) != nrow(testDF)) {
            appendLog(paste("ERROR: Mismatch => #predictions =", length(pred_probs), "#rows in testDF =", nrow(testDF)))
            shinyjs::alert("Mismatch in predictions vs. test data rows. Check the log.")
            return(NULL)
        }
        roc_obj <- pROC::roc(response = label, predictor = pred_probs)
        best_thresh <- as.numeric(pROC::coords(roc_obj, "best", ret = "threshold", best.method = "youden"))
        pos_label <- levels(testDF[[input$annotationColumn]])[2]
        neg_label <- levels(testDF[[input$annotationColumn]])[1]
        pred_class <- ifelse(pred_probs > best_thresh, pos_label, neg_label)
        pred_class <- factor(pred_class, levels = c(neg_label, pos_label))
        if (length(pred_class) != length(testDF[[input$annotationColumn]])) {
            appendLog("ERROR: predicted classes length != testDF rows. Could not compute confusionMatrix.")
            shinyjs::alert("Prediction mismatch with test labels, see log.")
            return(NULL)
        }

        cm <- caret::confusionMatrix(pred_class, testDF[[input$annotationColumn]], positive = pos_label)
        rv$confMat <- cm
        rv$rocObj <- roc_obj
        rv$bestThresh <- best_thresh
        appendLog("Test evaluation complete.")
        resolvedOutDir <- resolveOutputDir(input$outputDir)
        if (is.na(resolvedOutDir)) return()  # Stop if output path is invalid
        if (nchar(input$outputDir) > 0 && dir.exists(resolvedOutDir)) {
            fn_ts <- rv$session_ts
            sample_ids <- if ("sample_id" %in% colnames(testDF)) testDF$sample_id else rownames(testDF)
            pred_df <- data.frame(sample_id = sample_ids, predicted_probability = rv$pred_probs)
            write.table(pred_df, file = file.path(resolvedOutDir, paste0("predicted_probabilities_", fn_ts, ".tsv")),
                        sep = "\t", row.names = FALSE, quote = FALSE)
            appendLog("Predicted probabilities written.")
            results_df <- data.frame(
                Accuracy = rv$confMat$overall["Accuracy"],
                Sensitivity = rv$confMat$byClass["Sensitivity"],
                Specificity = rv$confMat$byClass["Specificity"],
                AUC = rv$rocObj$auc,
                Best_Threshold = rv$bestThresh
            )
            write.table(results_df, file = file.path(resolvedOutDir, paste0("evaluation_results_", fn_ts, ".tsv")),
                        sep = "\t", row.names = FALSE, quote = FALSE)
            appendLog("Evaluation results written.")
            cm_df <- as.data.frame(rv$confMat$table)
            write.table(cm_df, file = file.path(resolvedOutDir, paste0("confusion_matrix_", fn_ts, ".tsv")),
                        sep = "\t", row.names = FALSE, quote = FALSE)
            appendLog("Confusion matrix written.")
            roc_data <- data.frame(tpr = rv$rocObj$sensitivities, fpr = 1 - rv$rocObj$specificities)
            p <- ggplot(roc_data[order(roc_data$tpr, decreasing = FALSE), ], aes(x = fpr, y = tpr)) +
                geom_step(color = "#505d44", linewidth = 1) +
                geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
                xlim(0, 1) + ylim(0, 1) +
                labs(title = paste0("ROC Curve (AUC = ", round(rv$rocObj$auc, 3), ")"),
                     x = "1 - Specificity", y = "Sensitivity") +
                theme_minimal()
          ggsave(file = file.path(resolvedOutDir, paste0("roc_curve_", fn_ts, ".png")),
                 plot = p, device = "png", width = 6, height = 6)
          appendLog("ROC curve written as PNG.")
          appendLog("Test outputs written.")
        }
    })

    # continue from testing to performance visualization
    observeEvent(input$testToPerfButton, {
        updateTabItems(session, "tabs", "perf_tab")
    })

    # continue from performance visualization to log
    observeEvent(input$perfToLogButton, {
        updateTabItems(session, "tabs", "log_tab")
    })

    # render UI outputs
    output$bestParamsTable <- renderDataTable({
        req(rv$bestParams)
        df <- data.frame(Param = names(rv$bestParams), Value = unlist(rv$bestParams))
        datatable(df, options = list(dom = 't'))
    })

    output$confusionMatrix <- renderPrint({
        req(rv$confMat)
        rv$confMat
    })

    output$testResults <- renderPrint({
        req(rv$confMat, rv$rocObj)
        cm <- rv$confMat
        cat("Accuracy:", cm$overall["Accuracy"], "\n")
        cat("Sensitivity:", cm$byClass["Sensitivity"], "\n")
        cat("Specificity:", cm$byClass["Specificity"], "\n")
        cat("Best Threshold:", rv$bestThresh, "\n")
    })

    output$rocPlot <- renderPlot({
        req(rv$rocObj)
        r <- rv$rocObj
        roc_data <- data.frame(tpr = r$sensitivities, fpr = 1 - r$specificities)
        ggplot(roc_data[order(roc_data$tpr, decreasing = FALSE), ], aes(x = fpr, y = tpr)) +
            geom_step(color = "#505d44", linewidth = 1) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            xlim(0, 1) + ylim(0, 1) +
            labs(title = paste0("ROC Curve (AUC = ", round(r$auc, 3), ")"),
                 x = "1 - Specificity", y = "Sensitivity") +
            theme_minimal()
    })

    output$trainOutputs <- renderPrint({
        cat("Training outputs (training and testing matrices, best parameters, model) saved to:\n", input$outputDir)
    })

    output$testOutputs <- renderPrint({
        cat("Testing outputs (predictions, evaluation results, confusion matrix, ROC curve) saved to:\n", input$outputDir)
    })

    output$logOutput <- renderText({
        rv$logs
    })

    output$sessionInfo <- renderPrint({
        sessionInfo()
    })
}

shinyApp(ui, server)
