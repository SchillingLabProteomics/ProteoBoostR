suppressWarnings(suppressMessages(library(shiny)))
suppressWarnings(suppressMessages(library(shinydashboard)))
suppressWarnings(suppressMessages(library(shinyjs)))
suppressWarnings(suppressMessages(library(DT)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(tibble)))
suppressWarnings(suppressMessages(library(caret)))
suppressWarnings(suppressMessages(library(xgboost)))
suppressWarnings(suppressMessages(library(pROC)))
suppressWarnings(suppressMessages(library(rBayesianOptimization)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(fontawesome)))

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

get_model_features <- function(model) {
    fn <- tryCatch(model[["feature_names"]], error = function(e) NULL)
    if (!is.null(fn) && length(fn) > 0) return(fn)
    imp <- tryCatch(xgb.importance(model = model), error = function(e) NULL)
    if (!is.null(imp) && "Feature" %in% names(imp)) return(unique(imp$Feature))
    character(0)
}

read_best_threshold <- function(path) {
    df <- tryCatch(read.delim(path, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(df) && "Best_Threshold" %in% names(df)) {
        as.numeric(df$Best_Threshold[1])
    } else {
        NA_real_
    }
}

doMergeStatic <- function(annot_df, protein_df, annotationColumn, subsetIDs = character(0)) {
    prot_t <- as.data.frame(t(protein_df), stringsAsFactors = FALSE)
    prot_t <- tibble::rownames_to_column(prot_t, "sample_id")
    df_merged <- if (!is.null(annot_df)) dplyr::left_join(prot_t, annot_df, by = "sample_id") else prot_t
    colnames(df_merged) <- gsub(";.*", "", names(df_merged))
    # optional subset
    if (length(subsetIDs) > 0) {
        keepCols <- unique(c("sample_id", subsetIDs, annotationColumn))
        df_merged <- df_merged[, intersect(keepCols, colnames(df_merged)), drop = FALSE]
    }
    # numeric conversion for features
    keepColsAlways <- c("sample_id", if (!is.null(annotationColumn)) annotationColumn else character(0))
    for (cn in colnames(df_merged)) {
        if (!(cn %in% keepColsAlways)) {
            tryNum <- suppressWarnings(as.numeric(df_merged[[cn]]))
             if (all(is.na(tryNum))) {
                df_merged[[cn]] <- NULL
            } else {
                df_merged[[cn]] <- tryNum
            }
        }
    }
    df_merged
}

makeAdhocTabUI <- function(id) {
    tabPanel(
        title = tagList(
             span(paste("Test", id)),
             actionButton(paste0("adhoc_remove_", id), NULL, icon = icon("times"),
                          class = "btn-link", style = "margin-left:6px;padding:0;")
        ),
        value = paste0("adhoc_", id),
        fluidRow(
            box(title = "Dataset", status = "primary", width = 6,
                fileInput(paste0("adhoc_annot_", id), "Annotation (.tsv, optional)"),
                fileInput(paste0("adhoc_prot_", id), "Protein Matrix (.tsv)"),
                textInput(paste0("adhoc_outdir_", id), "Output Directory (optional):", value = ""),
                uiOutput(paste0("adhoc_outdir_status_", id)),
                uiOutput(paste0("adhoc_annotColUI_", id)),
                uiOutput(paste0("adhoc_classUI_", id))
            ),
            box(title = "Model & Threshold", status = "primary", width = 6,
                helpText("Uploaded model/eval take precedence over in-session values."),
                uiOutput(paste0("adhoc_model_info_", id)),
                fileInput(paste0("adhoc_model_", id), "Saved Model (.rds)"),
                tags$hr(),
                fileInput(paste0("adhoc_evaltsv_", id), "Evaluation results (.tsv with Best_Threshold)"),
                uiOutput(paste0("adhoc_threshold_info_", id)),
                sliderInput(paste0("adhoc_thresh_band_", id),
                            "\u00B1 range around base threshold",
                            min = 0, max = 0.5, value = 0.1, step = 0.01)
            )
        ),
        fluidRow(
            box(title = "", status = "primary", width = 12,
                actionButton(paste0("adhoc_eval_", id), "Evaluate", icon = icon("play"))
            )
        ),
        fluidRow(
            box(title = "Predicted Samples (ranked)", status = "primary", width = 6,
                plotOutput(paste0("adhoc_predplot_", id), height = "360px")
            ),
            box(title = "Scores", status = "primary", width = 6,
                DT::DTOutput(paste0("adhoc_scoretbl_", id))
            )
        ),
        fluidRow(
            box(title = "Confusion Matrix", status = "info", width = 6,
               verbatimTextOutput(paste0("adhoc_cm_", id))
            ),
            box(title = "ROC Curve", status = "primary", width = 6,
                plotOutput(paste0("adhoc_roc_", id), height = "360px")
            )
        ),
        fluidRow(
            box(title = "Metrics", status = "info", width = 6,
               verbatimTextOutput(paste0("adhoc_metrics_", id))
            )
        )
    )
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

            /* Hover text/icons white */
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li > a:hover,
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li:hover > a,
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li > a:focus {
              color: #fff !important;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li > a:hover .fa,
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li:hover > a .fa,
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li > a:focus .fa {
              color: #fff !important;
            }

            /* Active text/icons white */
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a {
              color: #fff !important;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a .fa {
              color: #fff !important;
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
            menuItem("Input for Training", tabName = "input_tab", icon = icon("upload")),
            menuItem("Model Training", tabName = "train_tab", icon = icon("cogs")),
            menuItem("Model Testing", tabName = "test_tab", icon = icon("play")),
            menuItem("Model Application", tabName = "adhoc_test_tab", icon = icon("flask")),
            menuItem("Log", tabName = "log_tab", icon = icon("clipboard-list"))
        )
    ),
    dashboardBody(
        tabItems(
            # landing page
            tabItem(tabName = "landing_tab",
                fluidRow(
                    box(title = "Welcome to ProteoBoostR (v1.0.0)", status = "primary", width = 12,
                        br(),
                        p("ProteoBoostR is a Shiny-based tool for supervised classification in proteomics data. It leverages XGBoost with Bayesian optimization to train, test, and apply predictive models. The UI supports applying trained models to independent datasets, displays ranked prediction tables, and offers a configurable random seed for reproducibility."),
                        br(),
                        p("Quick Start (Defaults)"),
                        tags$ol(
                          tags$li("Upload your annotation and protein matrix files."),
                          tags$li("Set the annotation column and class labels."),
                          tags$li("(Optionally) Provide a protein subset."),
                          tags$li("Specify your output directory (folder name only in Docker)."),
                          tags$li("Set the random seed (optional) and the train/test split."),
                          tags$li("Click ", tags$b("Start Training"), " to process data and train the model."),
                          tags$li("Go to ", tags$b("Model Testing"), " to evaluate the model. Review the ranked prediction table, metrics, confusion matrix, and ROC curve."),
                          tags$li("Use ", tags$b("Model Application"), " to apply your model to independent datasets (with or without labels).")
                        ),
                        br(),
                        br(),
                        actionButton("goLandingContinue", "Start Training", icon = icon("play")),
                        actionButton("goLandingToAdhoc", "Jump to Model Application", icon = icon("flask"))
                    )
                )
            ),
            # input tab
            tabItem(tabName = "input_tab",
                fluidRow(
                    box(title = "Upload Files", status = "primary", width = 6,
                        fileInput("annotationFile", "Annotation File (.tsv only)"),
                        helpText("Annotation must be a TSV file with first column named sample_id).", style = "margin-top: -20px"),
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
                        helpText("Pick the column for class label (sample_id excluded)."),
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
                    ),
                    box(title = "Random Seed", status = "info", width = 6,
                        numericInput("randSeed", "Seed:", value = 1234, min = 1, step = 1),
                        helpText("Set a random seed for reproducibility.")
                    )
                ),
                fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("inputContinueButton", "Continue to Model Training", icon = icon("arrow-right"))
                    )
                )
            ),
            # model training tab
            tabItem(tabName = "train_tab",
                fluidRow(
					# collapsible box for basic optimization options, default = open
                    box(title = "Bayesian Optimization Settings (General)", status = "primary", width = 12,
						solidHeader = TRUE, collapsible = TRUE, collapse = FALSE,
                        helpText("General parameters control learning rate, tree depth, and sampling."),
						box(status = "primary", width = 6,
							sliderInput("eta_range", "eta range:", min = 0, max = 1, value = c(0.1, 0.3), step = 0.001),
							helpText("Controls how much the model updates at each boosting step. A lower eta slows learning but allows for better convergence."),
							helpText("Use a higher eta (0.1 - 0.3) for small datasets to avoid underfitting. For large datasets, use a lower eta (0.01 - 0.1) with more boosting rounds to improve performance."),
							sliderInput("subsample_range", "subsample range:", min = 0, max = 1, value = c(0.8, 1.0),  step = 0.1),
							helpText("Specifies the fraction of training samples used to grow each tree, introducing randomness and reducing overfitting."),
							helpText("Use lower values (0.6 - 0.9) for small datasets to prevent overfitting. For very large datasets, subsample ~0.8 can speed up training without losing performance.")
						),
						box(status = "primary", width = 6,
							sliderInput("md_range", "max_depth range:", min = 0, max = 20, value = c(3,6), step = 1),
							helpText("Determines the depth of each decision tree. A larger depth captures more complexity but increases overfitting risk."),
							helpText("For small datasets, keep it low (2-6) to avoid overfitting. For large datasets, experiment with higher values (4-10) if needed."),
							sliderInput("colsample_range", "colsample_bytree range:", min = 0, max = 1, value = c(0.8, 1.0),  step = 0.1),
							helpText("Controls how many features are randomly sampled when growing a tree. Helps prevent reliance on specific features."),
							helpText("If you have few features (<10), keep it high (0.8 - 1.0). If you have many features (>50), tune it lower (0.5 - 0.8) to avoid overfitting.")
						),
					),
				),
                fluidRow(
                    # collapsible box for regularization options, default = closed
                    box(title = "Advanced Bayesian Optimization Settings (Regularization)", status = "primary", width = 12,
                        solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
						helpText("Regularization parameters regulate model complexity and reduce overfitting."),
						box(status = "primary", width = 6,
							sliderInput("child_weight_range", "min_child_weight range:", min = 0, max = 20, value = c(0,10), step = 1),
							helpText("Sets the minimum sum of instance weights needed to make a split. Higher values make the model conservative."),
							helpText("Use higher values (3-10) for small datasets to avoid overfitting. For large datasets, lower values (1-3) allow more splits for better learning."),
							numericInput("alpha_min", "alpha (min)", 0),
							numericInput("alpha_max", "alpha (max)", 5),
							helpText("Encourages sparsity in feature selection by penalizing the absolute values of leaf weights, helping with feature selection."),
							helpText("If you have many features, use higher values (1-5) to eliminate irrelevant ones. If you have few features, keep it low (0-1).")
						),
						box(status = "primary", width = 6,
							sliderInput("gamma_range", "gamma range:", min = 0, max = 20, value = c(0, 5), step = 0.1),
							helpText("Forces tree splits to have a significant loss reduction, controlling overfitting."),
							helpText("Use higher values (≥1) for small datasets to avoid unnecessary splits. In large datasets, tune it from 0 to 5 based on performance."),
							numericInput("lambda_min", "lambda (min)", 0),
							numericInput("lambda_max", "lambda (max)", 5),
							helpText("Penalizes large weight values, helping to smooth the model and reduce overfitting."),
							helpText("Use higher values (1-10) for small datasets to control overfitting. For large datasets, set it moderate (1-5) unless overfitting is observed.")
						),
                    ),
                ),
				fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("trainContinueButton", "Train the Model", icon = icon("gears"))
                    )
                ),
                fluidRow(
                    box(title = "Best Hyperparameters", status = "primary", width = 12,
                         dataTableOutput("bestParamsTable")
                    )
                ),
                fluidRow(
                    box(status = "primary", width = 12,
                        actionButton("testContinueButton", "Continue to Model Testing", icon = icon("arrow-right"))
                    )
                )
            ),
            # model testing tab
            tabItem(tabName = "test_tab",
                fluidRow(
                    box(title = "", status = "primary", width = 12,
                        actionButton("evaluateContinueButton", "Test Model", icon = icon("play"))
                    )
                ),
				fluidRow(
					box(title = "Predicted Samples (ranked)", status = "primary", width = 6,
						plotOutput("test_predplot", height = "360px")
					),
					box(title = "Scores", status = "primary", width = 6,
						DT::DTOutput("test_scoretbl")
					)
				),
				fluidRow(
					box(title = "Confusion Matrix", status = "info", width = 6,
						verbatimTextOutput("confusionMatrix")
					),
					box(title = "ROC Curve", status = "primary", width = 6,
						plotOutput("rocPlot", height = "360px")
					)
				),
				fluidRow(
					box(title = "Metrics", status = "info", width = 6,
						verbatimTextOutput("testResults")
					)
				)
            ),
            # log tab
            tabItem(tabName = "log_tab",
                fluidRow(
                    box(title = "Session Info", status = "primary", width = 12,
                        verbatimTextOutput("sessionInfo")
                    )
                )
            ),
            tabItem(tabName = "adhoc_test_tab",
                fluidRow(
                    box(width = 12, status = "primary",
                        tabsetPanel(id = "adhocTabs", type = "tabs",
                            tabPanel(title = HTML("<i class='fa fa-plus'></i> Add test"), value = "adder_tab")
                        )
                    )
                )
            )
        )
    ),
    title = "ProteoBoostR"
)

server <- function(input, output, session) {

    options(shiny.maxRequestSize = 100 * 1024^2) # 100MB limit

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
        logSet = FALSE,
        adhocCounter = 0,
        adhocTabs = list(),
        adhocRegistered = character(0),
        applyBusy = FALSE
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

    # landing page -> go to adhoc tab
    observeEvent(input$goLandingToAdhoc, {
      updateTabItems(session, "tabs", "adhoc_test_tab")
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

    # tab enabling logic
    disableAllExceptLog <- function() {
        shinyjs::addClass(selector = "a[data-value='input_tab']", class = "tab-disabled")
        shinyjs::addClass(selector = "a[data-value='train_tab']", class = "tab-disabled")
        shinyjs::addClass(selector = "a[data-value='test_tab']", class = "tab-disabled")
    }

    enableTabsBasedOnInputs <- function() {
        validInput <- !is.null(rv$annot) && !is.null(rv$proteins) && nchar(input$outputDir) > 0 && dir.exists(input$outputDir)
        if (validInput) {
            shinyjs::removeClass(selector = "a[data-value='train_tab']", class = "tab-disabled")
            shinyjs::removeClass(selector = "a[data-value='test_tab']", class = "tab-disabled")
        } else {
            shinyjs::addClass(selector = "a[data-value='train_tab']", class = "tab-disabled")
            shinyjs::addClass(selector = "a[data-value='test_tab']", class = "tab-disabled")
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
        df_merged <- dplyr::left_join(prot_t, rv$annot, by = "sample_id")
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
        } else {
            appendLog("Annotation file read failed or is empty.")
        }
    })

    # load protein matrix
    observeEvent(input$proteinFile, {
        req(input$proteinFile)
        if (!grepl("\\.tsv$", input$proteinFile$name, ignore.case = TRUE)) {
            appendLog("Error: Protein file must be a .tsv!")
            return(NULL)
        }
        df <- tryCatch({
            read.delim(input$proteinFile$datapath, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
        }, error = function(e) {
            appendLog(paste("Error reading protein matrix:", e$message))
            NULL
        })
        if (!is.null(df)) {
            rv$proteins <- df
            appendLog(paste("Protein matrix loaded. Dimensions:", paste(dim(df), collapse = " x ")))
        } else {
            appendLog("Protein matrix read failed or is empty.")
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
          selectInput("posClass", "Positive Class (1)", choices = uniqueVals, selected = uniqueVals[1]),
          selectInput("negClass", "Negative Class (0)", choices = uniqueVals, selected = uniqueVals[2])
        )
    })

    # input tab: Continue button – partition data and write subsetted train/test matrices
    observeEvent(input$inputContinueButton, {
        set.seed(if (is.null(input$randSeed) || is.na(input$randSeed)) 1234 else input$randSeed)

        # ensure data present and merged BEFORE requiring mergedData
        req(rv$annot, rv$proteins)
        if (is.null(rv$mergedData)) doMerge()

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
            trainDF <- trainDF[, c("sample_id", intersect(subsetIDs, colnames(trainDF)), input$annotationColumn), drop = FALSE]
            testDF <- testDF[, c("sample_id", intersect(subsetIDs, colnames(testDF)), input$annotationColumn), drop = FALSE]
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
        set.seed(if (is.null(input$randSeed) || is.na(input$randSeed)) 1234 else input$randSeed)
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
                lambda = best_par["lambda"]
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

    # model testing: Continue to Evaluate button
    observeEvent(input$evaluateContinueButton, {
        if (isTRUE(rv$applyBusy)) return(invisible(NULL))  # guard double-clicks
        rv$applyBusy <- TRUE
        shinyjs::disable("evaluateContinueButton")
        on.exit({ rv$applyBusy <- FALSE; shinyjs::enable("evaluateContinueButton") }, add = TRUE)

        withProgress(message = "Applying model...", value = 0, {
            incProgress(0.05, detail = "Preparing test set")
            set.seed(if (is.null(input$randSeed) || is.na(input$randSeed)) 1234 else input$randSeed)

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

            incProgress(0.20, detail = "Building DMatrix")
            appendLog("Evaluating model on test data...")
            testDF <- rv$testData
            featCols <- setdiff(colnames(testDF), c("sample_id", input$annotationColumn))
            features <- testDF[, featCols, drop = FALSE]
            label <- as.numeric(testDF[[input$annotationColumn]]) - 1
            appendLog(paste("Test model columns used for XGB:", paste(featCols, collapse = ", ")))
            dtest <- xgb.DMatrix(as.matrix(features), label = label)

            incProgress(0.40, detail = "Predicting")
            pred_probs <- predict(rv$xgbModel, dtest)
            rv$pred_probs <- pred_probs
            if (length(pred_probs) != nrow(testDF)) {
              appendLog(paste("ERROR: Mismatch => #predictions =", length(pred_probs), "#rows in testDF =", nrow(testDF)))
              shinyjs::alert("Mismatch in predictions vs. test data rows. Check the log.")
              return(NULL)
            }

            incProgress(0.60, detail = "Computing ROC / threshold")
            roc_obj <- pROC::roc(response = label, predictor = pred_probs)
            best_thresh <- pROC::coords(roc_obj, "best", ret = "threshold", best.method = "youden")
            # if best_thres is a list take the threshold closest to 0.5
            # common approach when predicted probabilities are expected to behave like calibrated risk scores
            if (nrow(best_thresh) > 1) {
                best_thresh <- as.numeric(best_thresh[which.min(abs(unlist(best_thresh) - 0.5)),])
            } else {
                best_thresh <- as.numeric(best_thresh)
            }
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

            # build ranked table for test_tab
            sample_ids <- if ("sample_id" %in% colnames(testDF)) testDF$sample_id else rownames(testDF)

            pred_tbl_test <- data.frame(
                sample_id = sample_ids,
                ann = testDF[[input$annotationColumn]],
                preds = pred_probs,
                stringsAsFactors = FALSE
            ) %>%
                dplyr::arrange(dplyr::desc(preds)) %>%
                dplyr::mutate(rank = dplyr::row_number())

            output$test_predplot <- renderPlot({
                ggplot(pred_tbl_test, aes(x = rank, y = preds, color = ann)) +
                    geom_point() +
                    geom_hline(yintercept = best_thresh) +
                    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
                    labs(x = "Sample rank (desc prob)", y = "Predicted probability", color = input$annotationColumn) +
                    theme_minimal()
            })

            output$test_scoretbl <- DT::renderDT({
                DT::datatable(
                    pred_tbl_test %>% dplyr::select(rank, sample_id, ann, preds),
                    rownames = FALSE,
                    options = list(pageLength = 10, order = list(list(0, "asc")))
                )
            })

            appendLog("Test evaluation complete.")

            incProgress(0.80, detail = "Saving outputs")
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
            incProgress(1.0, detail = "Done")
        })
    })

    observeEvent(input$adhocTabs, {
        if (identical(input$adhocTabs, "adder_tab")) {
            rv$adhocCounter <- rv$adhocCounter + 1
            id <- rv$adhocCounter
            insertTab(inputId = "adhocTabs",
                      tab = makeAdhocTabUI(id),
                      target = "adder_tab", position = "before", select = TRUE)
            rv$adhocTabs[[as.character(id)]] <- list()
        }
    })

    observe({
        ids <- names(rv$adhocTabs)
            for (id in ids) {
                if (id %in% rv$adhocRegistered) next
                rv$adhocRegistered <- c(rv$adhocRegistered, id)

                local({
                    myid <- id

                    # load annotation
                    observeEvent(input[[paste0("adhoc_annot_", id)]], ignoreInit = TRUE, {
                        f <- input[[paste0("adhoc_annot_", id)]]
                        if (is.null(f)) return(NULL)
                        df <- tryCatch({
                            raw <- read.delim(f$datapath, header = TRUE, stringsAsFactors = FALSE)
                            colnames(raw)[1] <- "sample_id"
                            raw
                        }, error = function(e) {
                            appendLog(paste("Ad-hoc", id, "annotation error:", e$message))
                            NULL
                        })
                        rv$adhocTabs[[id]]$annot <- df
                    })

                    # load protein matrix
                    observeEvent(input[[paste0("adhoc_prot_", id)]], ignoreInit = TRUE, {
                        f <- input[[paste0("adhoc_prot_", id)]]
                        if (is.null(f)) return(NULL)
                        df <- tryCatch({
                            read.delim(f$datapath, header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
                        }, error = function(e) {
                            appendLog(paste("Ad-hoc", id, "protein error:", e$message))
                            NULL
                        })
                        rv$adhocTabs[[id]]$prot <- df
                    })

                    # annotation column selector
                    output[[paste0("adhoc_annotColUI_", id)]] <- renderUI({
                        req(rv$adhocTabs[[id]]$annot)
                        cols <- setdiff(colnames(rv$adhocTabs[[id]]$annot), "sample_id")
                        selectInput(paste0("adhoc_annotCol_", id), "Annotation Column:", choices = cols)
                    })

                    # classes selector
                    output[[paste0("adhoc_classUI_", id)]] <- renderUI({
                        req(rv$adhocTabs[[id]]$annot, input[[paste0("adhoc_annotCol_", id)]])
                        u <- sort(unique(rv$adhocTabs[[id]]$annot[[ input[[paste0("adhoc_annotCol_", id)]] ]]))
                        if (length(u) < 2) return(tags$p("Need ≥2 unique values in this column."))
                        tagList(
                          selectInput(paste0("adhoc_pos_", id), "Positive Class (1)", choices = u, selected = u[1]),
                          selectInput(paste0("adhoc_neg_", id), "Negative Class (0)", choices = u, selected = u[2])
                        )
                    })

                    # check output directory status
                    output[[paste0("adhoc_outdir_status_", id)]] <- renderUI({
                        userPath <- input[[paste0("adhoc_outdir_", id)]]
                        if (is.null(userPath) || nchar(userPath) == 0) {
                            # fallback to main outputDir status
                            mainPath <- if (!is.null(input$outputDir)) input$outputDir else ""
                            resolved <- if (nchar(mainPath) > 0) resolveOutputDir(mainPath) else NA
                            if (is.na(resolved)) {
                                span(icon("times"), " Using main Output Path: invalid / not set", class = "invalid-path")
                            } else if (dir.exists(resolved)) {
                                span(icon("check"), paste(" Using main Output Path:", resolved), class = "valid-path")
                            } else {
                                span(icon("times"), paste(" Using main Output Path (invalid):", resolved), class = "invalid-path")
                            }
                        } else {
                            resolved <- resolveOutputDir(userPath)
                            if (is.na(resolved)) {
                                span(icon("times"), " Invalid Output Path: enter folder name only (Docker)", class = "invalid-path")
                            } else if (dir.exists(resolved)) {
                                span(icon("check"), paste(" Valid Output Path:", resolved), class = "valid-path")
                            } else {
                                span(icon("times"), paste(" Invalid Output Path:", resolved), class = "invalid-path")
                            }
                        }
                    })

                    # check model status
                    output[[paste0("adhoc_model_info_", id)]] <- renderUI({
                        upl <- input[[paste0("adhoc_model_", id)]]
                        hasUpload  <- !is.null(upl)
                        hasSession <- !is.null(rv$xgbModel)
                        sessThr    <- !is.null(rv$bestThresh) && !is.na(rv$bestThresh)

                        tags$div(
                            if (hasUpload) {
                                span(icon("check"), paste(" Uploaded model selected:", upl$name), class = "valid-path")
                            } else if (hasSession) {
                                span(icon("check"), " In-session model used", class = "valid-path")
                            } else {
                                span(icon("times"), " No model available yet", class = "invalid-path")
                            },
                            tags$br(),
                            if (hasUpload && hasSession) {
                                tags$small("Note: Uploaded model overrides in-session model")
                            }
                        )
                    })

                    # check threshold info
                    observeEvent(
                        list(input[[paste0("adhoc_evaltsv_", id)]],
                           input[[paste0("adhoc_thresh_band_", id)]],
                           rv$bestThresh),
                        {
                            base_thr <- NA_real_
                            evalTsv <- input[[paste0("adhoc_evaltsv_", id)]]
                            if (!is.null(evalTsv)) base_thr <- read_best_threshold(evalTsv$datapath)
                            if (is.na(base_thr) && !is.null(rv$bestThresh)) base_thr <- rv$bestThresh

                            band <- input[[paste0("adhoc_thresh_band_", id)]]
                            if (is.null(band)) band <- 0.15
                            tmin <- if (!is.na(base_thr)) max(0, base_thr - band) else NA_real_
                            tmax <- if (!is.na(base_thr)) min(1, base_thr + band) else NA_real_

                            output[[paste0("adhoc_threshold_info_", id)]] <- renderUI({
                                if (is.na(base_thr)) {
                                    tagList(
                                        span(icon("times"),
                                             " No base threshold available yet",
                                             class = "invalid-path")
                                    )
                                } else {
                                    tagList(
                                        span(icon("check"), sprintf(" Base threshold: %f", base_thr), class = "valid-path"),
                                        tags$br(),
                                        tags$span(sprintf("Band: \u00B1%.2f [%f, %f]", band, tmin, tmax))
                                    )
                                }
                            })
                        },
                        ignoreInit = FALSE
                    )

                    # evaluate
                    observeEvent(input[[paste0("adhoc_eval_", id)]], {
                        btn_id <- paste0("adhoc_eval_", myid)
                        shinyjs::disable(btn_id)
                        on.exit(shinyjs::enable(btn_id), add = TRUE)

                        withProgress(message = "Evaluating on new dataset...", value = 0, {
                            incProgress(0.05, detail = "Checking inputs")
                            req(rv$adhocTabs[[id]]$prot)
                            has_ann_file <- !is.null(rv$adhocTabs[[id]]$annot)
                            annCol  <- if (has_ann_file) input[[paste0("adhoc_annotCol_", id)]] else NULL

                            incProgress(0.15, detail = "Merging & preprocessing")
                            merged <- doMergeStatic(
                                if (has_ann_file) rv$adhocTabs[[id]]$annot else NULL,
                                rv$adhocTabs[[id]]$prot,
                                annCol
                            )

                            # is single sample or single class
                            is_single <- nrow(merged) == 1 || if(!is.null(annCol)) length(unique(merged[[annCol]])) == 1 else TRUE

                            has_labels <- (has_ann_file && !is_single) && !is.null(annCol) && annCol %in% colnames(merged)

                            df <- if (has_labels && !is_single) {
                                preprocessData(
                                    merged,
                                    annCol,
                                    input[[paste0("adhoc_neg_", id)]],
                                    input[[paste0("adhoc_pos_", id)]]
                                )
                            } else {
                                merged
                            }

                            incProgress(0.30, detail = "Loading model")
                            model <- NULL
                            # prefer uploaded model if provided
                            upl <- input[[paste0("adhoc_model_", id)]]
                            if (!is.null(upl)) {
                                model <- tryCatch(readRDS(upl$datapath),
                                                  error = function(e) {
                                                      shinyjs::alert("Error loading uploaded .rds model.")
                                                      appendLog(e$message)
                                                      NULL
                                                  })
                            } else if (!is.null(rv$xgbModel)) {
                                model <- rv$xgbModel
                            }
                            req(model)

                            incProgress(0.45, detail = "Aligning features")
                            model_feats <- get_model_features(model)
                            if (length(model_feats) == 0) {
                                shinyjs::alert("Model has no feature_names; cannot align.")
                                return(NULL)
                            }
                            featColsDf <- setdiff(colnames(df), c("sample_id", annCol))
                            present <- intersect(model_feats, featColsDf)
                            missing <- setdiff(model_feats, featColsDf)
                            X <- df[, present, drop = FALSE]
                            for (m in missing) X[[m]] <- NA_real_
                            X <- X[, model_feats, drop = FALSE]
                            X[] <- lapply(X, function(v) as.numeric(v))

                            incProgress(0.60, detail = "Predicting")
                            dtest <- if (has_ann_file && !is_single) {
                                        xgb.DMatrix(as.matrix(X),
                                                    label = as.numeric(df[[annCol]]) - 1,
                                                    missing = NA)
                                    } else {
                                        xgb.DMatrix(as.matrix(X),
                                                    missing = NA)
                                    }
                            pred_probs <- predict(model, dtest)

                            incProgress(0.70, detail = "Scoring & thresholds")
                            pred_tbl <- tibble::tibble(
                                sample_id = df$sample_id,
                                preds     = as.numeric(pred_probs)
                            )
                            if (has_labels) pred_tbl$ann <- df[[annCol]]

                            pred_tbl <- pred_tbl %>%
                                dplyr::arrange(dplyr::desc(preds)) %>%
                                dplyr::mutate(rank = dplyr::row_number())

                            # set base threshold
                            base_thr <- NA_real_
                            evalTsv <- input[[paste0("adhoc_evaltsv_", id)]]
                            if (!is.null(evalTsv)) base_thr <- read_best_threshold(evalTsv$datapath)
                            if (is.na(base_thr) && !is.null(rv$bestThresh)) base_thr <- rv$bestThresh

                            band <- input[[paste0("adhoc_thresh_band_", id)]]
                            tmin <- if (!is.na(base_thr)) max(0, base_thr - band) else NA_real_
                            tmax <- if (!is.na(base_thr)) min(1, base_thr + band) else NA_real_

                            cm <- NULL
                            roc_obj <- NULL
                            n_kept <- 0
                            n_gray <- 0

                            # show threshold info
                            output[[paste0("adhoc_threshold_info_", id)]] <- renderUI({
                                tagList(
                                    p(sprintf("Base threshold: %f", base_thr)),
                                    p(sprintf("Band: \u00B1%.2f [%f, %f]", band, tmin, tmax))
                                )
                            })

                            # zone classification
                            if (!is.na(base_thr)) {
                                levs <- if (has_labels) levels(df[[annCol]]) else c("neg","pos")
                                pos_label <- levs[length(levs)]
                                neg_label <- levs[1]
                                pred_tbl <- pred_tbl %>%
                                                mutate(
                                                    zone = dplyr::case_when(
                                                        preds < tmin ~ "below",
                                                        preds > tmax ~ "above",
                                                        TRUE         ~ "not classified"
                                                    ),
                                                    pred_label = dplyr::case_when(
                                                        zone == "below" ~ neg_label,
                                                        zone == "above" ~ pos_label,
                                                        TRUE            ~ NA_character_
                                                    )
                                                )
                            }

                            if (has_labels) {
                                keep_idx <- pred_tbl$zone != "not classified"
                                n_kept <- sum(keep_idx)
                                n_gray <- nrow(pred_tbl) - n_kept

                                # filtered y and preds
                                y_kept      <- pred_tbl$ann[keep_idx]
                                preds_kept  <- pred_tbl$preds[keep_idx]

                                # confusion matrix
                                if (n_kept > 0) {
                                    pred_class_kept <- ifelse(preds_kept > base_thr, pos_label, neg_label)
                                    pred_class_kept <- factor(pred_class_kept, levels = c(neg_label, pos_label))
                                    y_kept_fac <- factor(y_kept, levels = c(neg_label, pos_label))
                                    if (length(pred_class_kept) == length(y_kept_fac)) {
                                        cm <- caret::confusionMatrix(pred_class_kept, y_kept_fac, positive = pos_label)
                                    }
                                }

                                # ROC
                                if (n_kept >= 2 && length(unique(y_kept)) == 2) {
                                    roc_obj <- pROC::roc(response = as.numeric(y_kept) - 1, predictor = preds_kept)
                                }
                            }

                            incProgress(0.85, detail = "Rendering outputs")

                            # ranked with plot + score table
                            output[[paste0("adhoc_predplot_", id)]] <- renderPlot({
                                p <- ggplot(pred_tbl, aes(x = rank, y = preds)) +
                                        geom_point() +
                                        scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
                                        labs(x = "Sample rank (desc prob)", y = "Predicted probability") +
                                        theme_minimal()

                                if (has_labels) {
                                     p <- p + aes(color = ann) + labs(color = annCol)
                                }
                                if (!is.na(base_thr)) {
                                    p <- p +
                                        annotate("rect", xmin = -Inf, xmax = Inf, ymin = tmin, ymax = tmax,
                                                 alpha = 0.15, fill = "grey70") +
                                        geom_hline(yintercept = base_thr)
                                }
                                p
                            })

                            # score table
                            output[[paste0("adhoc_scoretbl_", id)]] <- DT::renderDT({
                                cols <- c("rank", "sample_id", "preds")
                                if (has_labels) cols <- c(cols, "ann")
                                if (!is.null(pred_tbl$pred_label)) cols <- c(cols, "pred_label", "zone")
                                DT::datatable(pred_tbl[, cols, drop = FALSE],
                                    rownames = FALSE,
                                    options = list(pageLength = 10, order = list(list(0, "asc")))
                                )
                            })

                            if (has_labels && !is.null(cm)) {
                                output[[paste0("adhoc_cm_", id)]] <- renderPrint(cm)
                            } else {
                                output[[paste0("adhoc_cm_", id)]] <- renderPrint({
                                                                        cat("No annotation provided or no evaluable samples.")
                                })
                            }

                            if (has_labels && !is.null(roc_obj)) {
                                # ROC (filtered)
                                output[[paste0("adhoc_roc_", id)]] <- renderPlot({
                                    roc_data <- data.frame(tpr = roc_obj$sensitivities, fpr = 1 - roc_obj$specificities)
                                    ggplot(roc_data[order(roc_data$tpr, decreasing = FALSE), ], aes(x = fpr, y = tpr)) +
                                        geom_step() +
                                        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
                                        xlim(0, 1) + ylim(0, 1) +
                                        labs(title = paste0("ROC (AUC = ", round(as.numeric(roc_obj$auc), 3), ")"),
                                             x = "1 - Specificity", y = "Sensitivity") +
                                        theme_minimal()
                                })
                            } else {
                                output[[paste0("adhoc_roc_", id)]] <- renderPlot({
                                  plot.new(); text(0.5, 0.5, "No ROC available (no annotation / insufficient samples).")
                                })
                            }
                            # metrics (filtered)
                            output[[paste0("adhoc_metrics_", id)]] <- renderPrint({
                                if (has_labels && !is.null(cm)) {
                                    acc <- if (!is.null(cm)) cm$overall["Accuracy"] else NA
                                    sen <- if (!is.null(cm)) cm$byClass["Sensitivity"] else NA
                                    spe <- if (!is.null(cm)) cm$byClass["Specificity"] else NA
                                    auc <- if (!is.null(roc_obj)) as.numeric(roc_obj$auc) else NA
                                    cat("Kept:", n_kept, "\n")
                                    cat("Ignored (treshold range):", n_gray, "\n")
                                    cat("Applied Threshold:", base_thr, "\n")
                                    cat("Accuracy:", acc, "\n")
                                    cat("Sensitivity:", sen, "\n")
                                    cat("Specificity:", spe, "\n")
                                    cat("AUC:", auc, "\n")
                                } else {
                                    cat("No annotation provided. Cannot compute metrics.")
                                }
                            })

                            incProgress(0.95, detail = "Saving files")
                            outdir_user <- input[[paste0("adhoc_outdir_", id)]]
                            resolvedOutDir <- if (nchar(outdir_user) > 0) resolveOutputDir(outdir_user) else resolveOutputDir(input$outputDir)
                            if (!is.na(resolvedOutDir) && dir.exists(resolvedOutDir)) {
                                fn_ts <- paste0(rv$session_ts, "_adhoc", id)

                                write.table(data.frame(sample_id = pred_tbl$sample_id,
                                                       predicted_probability = pred_tbl$preds),
                                            file = file.path(resolvedOutDir, paste0("predicted_probabilities_", fn_ts, ".tsv")),
                                            sep = "\t", row.names = FALSE, quote = FALSE)
                                if (has_labels && !is.null(cm)) {
                                    write.table(data.frame(
                                        Accuracy = cm$overall["Accuracy"],
                                        Sensitivity = cm$byClass["Sensitivity"],
                                        Specificity = cm$byClass["Specificity"],
                                        AUC = if (!is.null(roc_obj)) as.numeric(roc_obj$auc) else NA,
                                        Applied_Threshold = base_thr
                                    ),
                                    file = file.path(resolvedOutDir, paste0("evaluation_results_", fn_ts, ".tsv")),
                                    sep = "\t", row.names = FALSE, quote = FALSE)

                                    write.table(as.data.frame(cm$table),
                                                file = file.path(resolvedOutDir, paste0("confusion_matrix_", fn_ts, ".tsv")),
                                                sep = "\t", row.names = FALSE, quote = FALSE)
                                }

                                p_save <- ggplot(pred_tbl, aes(x = rank, y = preds)) +
                                    geom_point() +
                                    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
                                    labs(x = "Sample rank", y = "Predicted class probability") +
                                    theme_minimal()
                                if (has_labels) p_save <- p_save + aes(color = ann) + labs(color = annCol)
                                if (!is.na(base_thr)) {
                                    p_save <- p_save +
                                        annotate("rect", xmin = -Inf, xmax = Inf, ymin = tmin, ymax = tmax,
                                                 alpha = 0.15, fill = "grey70") +
                                        geom_hline(yintercept = base_thr)
                                }
                                ggsave(file = file.path(resolvedOutDir, paste0("predicted_samples_", fn_ts, ".png")),
                                       plot = p_save, device = "png", width = 7, height = 5)
                                appendLog(paste("Ad-hoc", id, "outputs written to:", resolvedOutDir))
                            } else {
                                appendLog(paste("Ad-hoc", id, "no valid output directory. Skipping writes."))
                            }

                            incProgress(1.0, detail = "Done")
                        })
                    })

                    # remove tab
                    observeEvent(input[[paste0("adhoc_remove_", id)]], {
                        removeTab("adhocTabs", target = paste0("adhoc_", id))
                        rv$adhocTabs[[id]] <- NULL
                    })
                })
            }
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
