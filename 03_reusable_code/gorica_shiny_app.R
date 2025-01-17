library(shiny)
library(restriktor)
library(Matrix)  # For matrix validation

# Define the Shiny app
ui <- fluidPage(
  
  titlePanel("Hypothesis Evaluation with GORICA"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Input Parameters"),
      
      textInput("num_params", "Number of Parameters:", value = "2", placeholder = "Enter a number"),
      
      uiOutput("parameter_inputs"),
      
      h3("Covariance Matrix"),
      
      uiOutput("covariance_inputs"),
      
      h3("Hypotheses"),
      
      textInput("num_hypotheses", "Number of Hypotheses:", value = "1", placeholder = "Enter a number"),
      
      uiOutput("hypotheses_inputs"),
      
      checkboxInput("add_null", "Include Null Hypothesis (Heq = TRUE):", value = FALSE),
      
      actionButton("run_analysis", "Run Analysis")
    ),
    
    mainPanel(
      h3("Results"),
      
      verbatimTextOutput("gorica_results")
    )
  )
)

server <- function(input, output, session) {
  
  # Reactive values to store hypotheses
  stored_hypotheses <- reactiveValues(values = list())
  
  # Dynamically generate inputs for parameters
  output$parameter_inputs <- renderUI({
    n_params <- as.numeric(input$num_params)
    if (is.na(n_params) || n_params < 2) return(NULL)
    
    lapply(1:n_params, function(i) {
      textInput(
        inputId = paste0("param_", i),
        label = paste("Value for Parameter", i, "(Name: ", letters[i], "):"),
        value = "0"
      )
    })
  })
  
  # Dynamically generate inputs for covariance matrix
  output$covariance_inputs <- renderUI({
    n_params <- as.numeric(input$num_params)
    if (is.na(n_params) || n_params < 2) return(NULL)
    
    matrixInputs <- lapply(1:n_params, function(i) {
      lapply(1:n_params, function(j) {
        textInput(
          inputId = paste0("cov_", i, "_", j),
          label = paste("Covariance [", letters[i], ",", letters[j], "]: " ),
          value = ifelse(i == j, "1", "0")
        )
      })
    })
    
    do.call(tagList, unlist(matrixInputs, recursive = FALSE))
  })
  
  # Dynamically generate inputs for hypotheses
  output$hypotheses_inputs <- renderUI({
    n_hypotheses <- as.numeric(input$num_hypotheses)
    if (is.na(n_hypotheses) || n_hypotheses < 1) return(NULL)
    
    lapply(1:n_hypotheses, function(i) {
      textInput(
        inputId = paste0("hypothesis_", i),
        label = paste("Hypothesis", i, ":"),
        value = if (!is.null(stored_hypotheses$values[[paste0("hypothesis_", i)]])) stored_hypotheses$values[[paste0("hypothesis_", i)]] else "",
        placeholder = "Example: a > b"
      )
    })
  })
  
  # Observe changes in the number of hypotheses and update stored values
  observeEvent(input$num_hypotheses, {
    n_hypotheses <- as.numeric(input$num_hypotheses)
    if (is.na(n_hypotheses) || n_hypotheses < 1) return()
    
    # Update stored values and initialize new ones as empty strings
    for (i in 1:n_hypotheses) {
      if (is.null(stored_hypotheses$values[[paste0("hypothesis_", i)]])) {
        stored_hypotheses$values[[paste0("hypothesis_", i)]] <- ""
      }
    }
    
    # Remove excess stored values if the number decreases
    stored_hypotheses$values <- stored_hypotheses$values[1:n_hypotheses]
  })
  
  # Synchronize text inputs with stored values without re-rendering UI
  observe({
    lapply(names(stored_hypotheses$values), function(id) {
      updateTextInput(session, inputId = id, value = stored_hypotheses$values[[id]])
    })
  })
  
  # Extract parameters and covariance matrix
  get_parameters <- reactive({
    n_params <- as.numeric(input$num_params)
    if (is.na(n_params) || n_params < 2) stop("Number of parameters must be at least 2 and numeric.")
    
    params <- sapply(1:n_params, function(i) {
      val <- as.numeric(input[[paste0("param_", i)]])
      if (is.na(val)) stop("All parameter values must be numeric.")
      val
    })
    names(params) <- letters[1:n_params]
    params
  })
  
  get_covariance_matrix <- reactive({
    n_params <- as.numeric(input$num_params)
    if (is.na(n_params) || n_params < 2) stop("Number of parameters must be at least 2 and numeric.")
    
    cov_matrix <- matrix(
      sapply(1:n_params, function(i) {
        sapply(1:n_params, function(j) {
          val <- as.numeric(input[[paste0("cov_", i, "_", j)]])
          if (is.na(val)) stop("All covariance values must be numeric.")
          val
        })
      }),
      nrow = n_params
    )
    colnames(cov_matrix) <- letters[1:n_params]
    rownames(cov_matrix) <- letters[1:n_params]
    cov_matrix
  })
  
  # Extract hypotheses
  get_hypotheses <- reactive({
    n_hypotheses <- as.numeric(input$num_hypotheses)
    if (is.na(n_hypotheses) || n_hypotheses < 1) stop("Number of hypotheses must be at least 1 and numeric.")
    
    hypotheses <- lapply(1:n_hypotheses, function(i) {
      hypothesis <- input[[paste0("hypothesis_", i)]]
      if (is.null(hypothesis) || hypothesis == "") {
        stop(paste("Hypothesis", i, "is empty. Please provide a valid hypothesis."))
      }
      hypothesis
    })
    names(hypotheses) <- paste0("h", 1:n_hypotheses)
    
    hypotheses
  })
  
  # Validate covariance matrix
  validate_covariance <- reactive({
    cov_matrix <- get_covariance_matrix()
    isSymmetric(cov_matrix) && all(eigen(cov_matrix)$values >= 0)
  })
  
  # Run GORICA analysis
  observeEvent(input$run_analysis, {
    if (!validate_covariance()) {
      showModal(modalDialog(
        title = "Invalid Covariance Matrix",
        "The covariance matrix must be symmetric and positive semi-definite.",
        easyClose = TRUE
      ))
      return()
    }
    
    params <- get_parameters()
    cov_matrix <- get_covariance_matrix()
    hypotheses <- get_hypotheses()
    
    # Run GORICA analysis
    gorica_res <- goric(
      object = params,
      VCOV = cov_matrix,
      hypotheses = hypotheses,
      Heq = input$add_null
    )
    
    output$gorica_results <- renderPrint({
      gorica_res
    })
  })
}

shinyApp(ui = ui, server = server)
