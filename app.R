# app.R
# Shiny dashboard for inspecting GLM's fitted on freqMPTLfreq2 (or any data with the same columns)
# To use the model, simply run the app.R (through your preferred means) and put your models in a models/ folder in the same directory
# Assumptions and Notes:
# - Models should be poisson, quasipoisson, negative binomial (such as from MASS), or inherit from zeroinfl (from library pscl)
# - Place .rds model files in a "models/" folder next to this app.R (use saveRDS(myModel, "filename.rds") to save the model)
# - Models should be saved with their model/data embedded. if myModel <- glm(...), then defaults should work
# - If myModel <- glm.nb(...), or myModel <- zeroinfl(...) you may need to manually attach data with myModel$data <- my_data
# - Using saveRDS(myModel, "filename.rds") should work for saving the objects properly
# - Make sure to include all columns and not just the ones you are fitting on, so that you can inspect residuals and fits across those predictors
# - Expected columns in the model data (case sensitive, in any order):
#   IDpol, ClaimNb, Exposure, Area, VehPower, VehAge, DrivAge, BonusMalus, VehBrand, VehGas, Density, Region
# - Categorical predictors: Area, VehBrand, VehGas, Region
# - Continuous predictors : VehPower, VehAge, DrivAge, BonusMalus, Density
# - You should be able to easily generalize this project to any set of predictors by changing the list of expected predictors in the R script

library(shiny)
library(ggplot2)
library(dplyr)
library(broom)
library(DT)
library(plotly)
library(readr)
library(MASS) # for glm.nb detection
library(tools)
library(pscl) #for zero inflated models

# Helper: scan models folder for .rds
list_model_files <- function(models_dir = "models") {
  if (!dir.exists(models_dir)) return(character(0))
  files <- list.files(models_dir, pattern = "\\.rds$", full.names = TRUE)
  names(files) <- basename(files)
  files
}


# Build grouped summary
group_summary <- function(df, predictor, nbins = 25,
                          cont_vars = c("VehPower","VehAge","DrivAge","BonusMalus","Density"),
                          resid_type = "pearson", mod, conf = 0.95, dispersion = 1) {
  
  bin_col <- ".__bin"
  df[[bin_col]] <- NA_character_
  
  # Make predictions for frequencies
  pred_count <- tryCatch(predict(mod, type = "response", newdata = df),
                         error = function(e) {
                           if (!is.null(mod$fitted.values) && length(mod$fitted.values) == nrow(df)) return(mod$fitted.values)
                           stop("Cannot get predictions from model: ", e$message)
                         })
  df$.pred_count <- pred_count
  df$.pred_freq  <- df$.pred_count / df$Exposure
  df$.obs_freq   <- df$ClaimNb / df$Exposure
  
  
  # select residual vector requested
  if (resid_type == "response") {
    res_vec = residuals(mod, type = "response");
  } else if (resid_type == "pearson") {
    if (inherits(mod, "zeroinfl")) {
      res_vec = residuals(mod, type = "pearson")
    } else {
      res_vec = residuals(mod, type = "pearson") / summary(mod)$dispersion #for quasipoisson. negbin and poisson have dispersion = 1
    }
  } else if (resid_type == "deviance") {
    #deviance residuals unavailable for zero inflated models
    res_vec = tryCatch(residuals(mod, type = "deviance"), error = function(e) rep(NA, nrow(df))) 
  }
  df$.resid <- res_vec
  
  # aggregate by the requested predictor
  z <- qnorm(1 - (1 - conf) / 2)
  if (predictor %in% cont_vars) {
    #cut the continuous predictor into bins
    grouped <- df %>%
      group_by(level = cut(df[[predictor]], nbins, include.lowest = TRUE)) %>%
      summarise(
        n_policies = n(),
        sum_exposure = sum(Exposure, na.rm = TRUE),
        sum_claims = sum(ClaimNb, na.rm = TRUE),
        observed_freq = ifelse(sum_exposure > 0, sum_claims / sum_exposure, NA_real_),
        predicted_count = sum(.pred_count, na.rm = TRUE),
        predicted_freq = ifelse(sum_exposure > 0, predicted_count / sum_exposure, NA_real_),
        mean_resid = mean(.resid, na.rm = TRUE),
        sd_resid = sd(.resid, na.rm = TRUE), 
        .groups = "drop"
      )
  } else {
    #use the categorical variable as bins
    grouped <- df %>%
      group_by(level = df[[predictor]]) %>%
      summarise(
        n_policies = n(),
        sum_exposure = sum(Exposure, na.rm = TRUE),
        sum_claims = sum(ClaimNb, na.rm = TRUE),
        observed_freq = ifelse(sum_exposure > 0, sum_claims / sum_exposure, NA_real_),
        predicted_count = sum(.pred_count, na.rm = TRUE),
        predicted_freq = ifelse(sum_exposure > 0, predicted_count / sum_exposure, NA_real_),
        mean_resid = mean(.resid, na.rm = TRUE),
        sd_resid = sd(.resid, na.rm = TRUE), 
        .groups = "drop"
      )
  }
  
  
  # CIs (normal approx). Use dispersion factor for observed/predicted count SE.
  grouped <- grouped %>%
    mutate(
      dispersion = dispersion,
      se_obs = ifelse(sum_exposure > 0, sqrt(pmax(sum_claims, 0)) * sqrt(dispersion) / sum_exposure, NA_real_),
      obs_lo = observed_freq - z * se_obs,
      obs_hi = observed_freq + z * se_obs,
      se_pred = ifelse(sum_exposure > 0, sqrt(pmax(predicted_count, 0)) * sqrt(dispersion) / sum_exposure, NA_real_),
      pred_lo = predicted_freq - z * se_pred,
      pred_hi = predicted_freq + z * se_pred,
      se_mean_resid = ifelse(n_policies > 1, sd_resid / sqrt(n_policies), NA_real_),
      resid_lo = mean_resid - z * se_mean_resid,
      resid_hi = mean_resid + z * se_mean_resid
    )
  grouped
}

ui <- fluidPage(
  titlePanel("GLM Frequency Model Explorer"),
  
  fluidRow(
    column(
      width = 3,
      h4("Model selection"),
      actionButton("rescan", "Rescan models/ folder"),
      br(), br(),
      uiOutput("model_picker_ui"),
      verbatimTextOutput("model_info"),
      hr(),
      h4("Graph Controls"),
      selectInput("predictor", "Predictor",
                  choices = c("Area","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region"),
                  selected = "Area"),
      numericInput("nbins", "Number of bins (equal-width for continuous)", value = 15, min = 2, max = 200, step = 1),
      radioButtons("resid_type", "Residual type", choices = c("pearson","response","deviance"), selected = "pearson"),
      numericInput(
        inputId = "n_points",
        label = "Number of individual residuals to plot",
        value = 1000,   # default
        min = 100,      # minimum
        max = 100000,    # Your computer will probably crash
        step = 100
      ),
      selectInput(
        inputId = "resid_xaxis",
        label = "Plot individual residuals against:",
        choices = c("Fitted frequency" = "fitted",
                    "Selected predictor" = "predictor"),
        selected = "predictor"
      ),
      sliderInput("conf_level", "CI level", min = 0.8, max = 0.99, value = 0.95, step = 0.01),
      hr(),
      downloadButton("download_grouped_csv", "Download grouped CSV"),
      br(), br(),
      p("Place one or more .rds model files in a folder named 'models/' next to this app.R. Models must include their training data in myModel$data")
    ),
    column(
      width = 9,
      fluidRow(
        column(8,
               uiOutput("metrics_ui")
        )
      ),
      fluidRow(
        column(12,
               h4("Grouped Observed vs Predicted"),
               plotlyOutput("grp_obs_pred_plot", height = "250px"),
               downloadButton("download_grp_plot", "Download grouped plot PNG")
        )
      ),
      fluidRow(
        column(12,
               h4("Grouped Residuals (mean ± CI)"),
               plotlyOutput("grp_resid_plot", height = "250px"),
               downloadButton("download_grp_resid_plot", "Download grouped residuals PNG")
        )
      ),
      fluidRow(
        column(12,
               h4("Individual Residuals"),
               plotlyOutput("indiv_resid_plot", height = "250px"),
               downloadButton("download_indiv_plot", "Download residuals scatter PNG")
        )
      )
    )
  ),
  
  fluidRow(
    column(
      width = 12,
      h4("Model Coefficients & Significance"),
      DTOutput("coef_table"),
      br(),
      verbatimTextOutput("notes")
    )
  )
)


server <- function(input, output, session) {
  models_dir <- reactiveVal("models")
  # reactive list of files
  model_files <- reactive({
    input$rescan
    Sys.sleep(0.1)
    list_model_files(models_dir())
  })
  
  # model picker UI
  output$model_picker_ui <- renderUI({
    files <- model_files()
    if (length(files) == 0) {
      tags$div(
        tags$b("No models found in 'models/' folder."),
        tags$div("Put .rds files in models/ and click 'Rescan models/ folder'.")
      )
    } else {
      selectInput("model_file", "Choose model (.rds)", choices = files, selected = files[1])
    }
  })
  
  # load selected model object and its data
  loaded <- reactive({
    req(input$model_file)
    path <- input$model_file
    mod <- tryCatch(readRDS(path), error = function(e) {
      showNotification(paste("Error reading", path, e$message), type = "error")
      return(NULL)
    })
    if (is.null(mod)) return(NULL)
    mf <- mod
    if (is.null(mf)) {
      # Notify user - model lacks embedded data
      showModal(modalDialog(
        title = "Model data not found in object",
        paste0("The model '", basename(path), "' does not contain an accessible model frame/data. ",
               "Please save models with their data (e.g., `saveRDS(mod, 'model.rds')` where mod was fit with model=TRUE)"),
        easyClose = TRUE
      ))
      return(list(mod = mod, data = NULL, error = TRUE))
    }
    # ensure correct columns present
    expected_cols <- c("IDpol","ClaimNb","Exposure","Area","VehPower","VehAge","DrivAge","BonusMalus","VehBrand","VehGas","Density","Region")
    missing <- setdiff(expected_cols, names(as.data.frame(mf$data)))
    if (length(missing) > 0) {
      showModal(modalDialog(
        title = "Model data missing expected columns",
        paste("The model data is missing these required columns:", paste(missing, collapse = ", ")),
        easyClose = TRUE
      ))
      return(list(mod = mod, data = as.data.frame(mf$data), error = TRUE))
    }
    list(mod = mod, data = as.data.frame(mf$data), error = FALSE)
  })
  
  output$model_info <- renderText({
    ld <- loaded()
    if (is.null(ld)) return("No model loaded.")
    if (isTRUE(ld$error)) return("Model loaded but missing data or columns; see modal.")
    mod <- ld$mod
    cls <- paste(class(mod), collapse = ", ")
    fam <- NA
    if (inherits(mod, "zeroinfl")) {
      fam <- "Zero Inflated Poisson"
    } else {
      try({ fam <- paste0(family(mod)$family, "(", family(mod)$link, ")") }, silent = TRUE)
    }
    paste0("Loaded: ", basename(input$model_file), "\nClass: ", cls, "\nFamily/link: ", fam)
  })
  
  
  # reactive: compute grouped summary
  grouped_df <- reactive({
    ld <- loaded()
    req(ld)
    if (isTRUE(ld$error)) return(NULL)
    df <- ld$data
    mod <- ld$mod
    # compute dispersion estimate (Pearson)
    if (inherits(mod, "zeroinfl")) {
      disp <- 1
    } else if (family(mod)$family == "poisson" || family(mod)$family == "quasipoisson") { #poisson or quasipoisson
      disp <- summary(mod)$dispersion
    } else { #negative binomial theta
      disp <- mod$theta
    }
    # group
    grp <- group_summary(df = df,
                         predictor = input$predictor,
                         nbins = input$nbins,
                         resid_type = input$resid_type,
                         mod = mod,
                         conf = input$conf_level,
                         dispersion = ifelse(is.na(disp) || disp <= 0, 1, disp))
    grp
  })
  
  # model metrics UI
  output$metrics_ui <- renderUI({
    ld <- loaded()
    req(ld)
    if (isTRUE(ld$error)) return(NULL)
    mod <- ld$mod
    out <- tagList()
    # BIC, deviance, logLik, pearson sum of squares as available
    bic <- tryCatch(BIC(mod), error = function(e) NA)
    dev <- tryCatch(deviance(mod), error = function(e) NA)
    logl <- tryCatch(logLik(mod), error = function(e) NA)
    disp <- tryCatch({prs <- residuals(mod, type="pearson"); sum(prs^2, na.rm=TRUE)/df.residual(mod)}, error = function(e) NA)
    out <- tagAppendChildren(out, tags$p(tags$b("BIC:"), ifelse(is.na(bic), "N/A", round(bic,3))))
    out <- tagAppendChildren(out, tags$p(tags$b("Deviance:"), ifelse(is.na(dev), "N/A", round(dev,3))))
    out <- tagAppendChildren(out, tags$p(tags$b("Log Likelihood:"), ifelse(is.na(logl), "N/A", as.numeric(round(logl,3)))))
    out <- tagAppendChildren(out, tags$p(tags$b("sum(Pearson^2) / df:"), ifelse(is.na(disp), "N/A", round(disp,3))))
    out
  })
  
  # coefficient table
  output$coef_table <- renderDT({
    ld <- loaded()
    req(ld)
    if (isTRUE(ld$error)) return(datatable(data.frame(Message="Model missing data or columns."), options = list(dom='t')))
    mod <- ld$mod
    tidy_mod <- tryCatch(broom::tidy(mod, conf.int = TRUE), error = function(e) {
      # fallback to summary table
      coefs <- coef(summary(mod))
      if (inherits(mod, "zeroinfl")) {
        df = as.data.frame(coefs)
        return(df)
      }
      df <- as.data.frame(coefs)
      df$term <- rownames(df)
      rownames(df) <- NULL
      cis <- tryCatch(confint(mod), error = function(e) NA)
      if (is.matrix(cis)) {
        df$conf.low <- cis[,1]
        df$conf.high <- cis[,2]
      }
      df
    })
    datatable(as.data.frame(tidy_mod), options = list(pageLength = 20, scrollX = TRUE))
  })
  
  
  # notes
  output$notes <- renderText({
    "Notes:
- Grouped observed and predicted frequencies are plotted as frequency = sum_claims / sum_exposure and predicted_count / sum_exposure.
- CIs use normal approximation with Pearson dispersion scaling (for quasipoisson).
- Models must include training/test data when saved with appropriate columns"
  })
  
  
  
  # Plot 1: grouped observed vs predicted
  grp_plot_gg <- reactive({
    grp <- grouped_df() 
    req(grp)
    p <- ggplot(grp, aes(x = level)) +
      geom_point(aes(y = observed_freq, color = "Observed"), size = 2, position = position_nudge(x = -0.15)) +
      geom_errorbar(aes(ymin = obs_lo, ymax = obs_hi, color = "Observed"), width = 0.2, position = position_nudge(x = -0.15)) +
      geom_point(aes(y = predicted_freq, color = "Predicted"), size = 2, position = position_nudge(x = 0.15)) +
      geom_errorbar(aes(ymin = pred_lo, ymax = pred_hi, color = "Predicted"), width = 0.2, position = position_nudge(x = 0.15)) +
      scale_color_manual(values = c("Observed" = "black", "Predicted" = "blue")) +
      labs(x = input$predictor, y = "Frequency", color = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p
  })
  
  output$grp_obs_pred_plot <- renderPlotly({
    p <- grp_plot_gg()
    ggplotly(p, tooltip = c("x","y","color"))
  })
  
  # download grouped plot
  output$download_grp_plot <- downloadHandler(
    filename = function() {
      paste0("grouped_obs_pred_", input$predictor, ".png")
    },
    content = function(file) {
      ggsave(file, plot = grp_plot_gg(), width = 10, height = 5, dpi = 150)
    }
  )
  
  # Plot 2: grouped residuals
  grp_resid_gg <- reactive({
    grp <- grouped_df()
    req(grp)
    grp$level <- factor(as.character(grp$level), levels = unique(grp$level))
    p <- ggplot(grp, aes(x = level, y = mean_resid)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = resid_lo, ymax = resid_hi), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = input$predictor, y = paste0("Mean ", input$resid_type, " residual")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p
  })
  
  output$grp_resid_plot <- renderPlotly({
    p <- grp_resid_gg()
    ggplotly(p, tooltip = c("x","y"))
  })
  
  output$download_grp_resid_plot <- downloadHandler(
    filename = function() {
      paste0("grouped_resid_", input$predictor, ".png")
    },
    content = function(file) {
      ggsave(file, plot = grp_resid_gg(), width = 10, height = 4, dpi = 150)
    }
  )
  
  # Plot 3: individual residuals vs fitted frequency
  indiv_gg <- reactive({
    ld <- loaded()
    req(ld)
    if (isTRUE(ld$error)) return(NULL)
    df <- ld$data
    mod <- ld$mod
    # predictions and residuals
    pred <- tryCatch(predict(mod, type = "response", newdata = df), error = function(e) {
      if (!is.null(mod$fitted.values) && length(mod$fitted.values) == nrow(df)) return(mod$fitted.values)
      stop("Cannot get predictions from model.")
    })
    df$.pred_freq <- pred / df$Exposure
    df$.obs_freq <- df$ClaimNb / df$Exposure
    df$.resid <- tryCatch(residuals(mod, type = input$resid_type), error = function(e) rep(NA, nrow(df)))
    
    max_points <- input$n_points
    if (nrow(df) > max_points) {
      indices <- sample(seq_len(nrow(df)), max_points)
    } else {
      indices <- seq_len(nrow(df))
    }
    if (input$resid_xaxis == "fitted") {
      p <- ggplot(df[indices, ], aes(x = .pred_freq, y = .resid)) +
        geom_point(aes(text = paste0("IDpol: ", IDpol, "<br>ClaimNb: ", ClaimNb, "<br>Exposure: ", Exposure)), alpha = 0.4) +
        labs(x = "Fitted frequency (predicted count / Exposure)", y = paste0(input$resid_type, " residual")) +
        theme_minimal()
    } else {
      p <- ggplot(df[indices, ], aes_string(x = input$predictor, y = ".resid")) +
        geom_point(aes(text = paste0("IDpol: ", IDpol, "<br>ClaimNb: ", ClaimNb, "<br>Exposure: ", Exposure)), alpha = 0.4) +
        labs(x = input$x_resid, y = paste0(input$resid_type, " residuals")) +
        theme_minimal()
    }
    p
  })
  
  output$indiv_resid_plot <- renderPlotly({
    p <- indiv_gg()
    ggplotly(p, tooltip = "text")
  })
  
  output$download_indiv_plot <- downloadHandler(
    filename = function() {
      paste0("individual_resid_", input$predictor, ".png")
    },
    content = function(file) {
      ggsave(file, plot = indiv_gg(), width = 8, height = 5, dpi = 150)
    }
  )
  
  # grouped CSV download
  output$download_grouped_csv <- downloadHandler(
    filename = function() {
      paste0("grouped_summary_", input$predictor, ".csv")
    },
    content = function(file) {
      grp <- grouped_df()
      req(grp)
      write_csv(grp, file)
    }
  )
  
}

shinyApp(ui, server)