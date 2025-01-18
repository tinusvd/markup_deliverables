# app.R
library(shiny)
library(restriktor)
library(Matrix)
library(ggplot2)

# Module UI for the analysis page
analysisUI <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      h3("Input Parameters"),
      numericInput(ns("num_params"), "Number of Parameters:", value = 2, min = 2),
      uiOutput(ns("parameter_inputs")),
      
      h3("Covariance Matrix"),
      uiOutput(ns("covariance_inputs")),
      
      h3("Hypotheses"),
      numericInput(ns("num_hypotheses"), "Number of Hypotheses:", value = 1, min = 1),
      uiOutput(ns("hypotheses_inputs")),
      
      checkboxInput(ns("add_null"), "Include Null Hypothesis (Heq = TRUE):", value = FALSE),
      actionButton(ns("run_analysis"), "Run Analysis")
    ),
    
    mainPanel(
      h3("Results"),
      verbatimTextOutput(ns("gorica_results")),
      h3("Visualization of GORICA Weights"),
      plotOutput(ns("gorica_plot"))
    )
  )
}

# Module server for the analysis page
analysisServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Parameter inputs
    output$parameter_inputs <- renderUI({
      req(input$num_params)
      lapply(1:input$num_params, function(i) {
        numericInput(
          session$ns(paste0("param_", i)),
          label = paste("Value for Parameter", i, "(Name: ", letters[i], "):"),
          value = 0
        )
      })
    })
    
    # Covariance matrix inputs
    output$covariance_inputs <- renderUI({
      req(input$num_params)
      div(
        style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px;",
        lapply(1:input$num_params, function(i) {
          lapply(1:input$num_params, function(j) {
            numericInput(
              session$ns(paste0("cov_", i, "_", j)),
              label = paste("Cov [", letters[i], ",", letters[j], "]"),
              value = ifelse(i == j, 1, 0)
            )
          })
        })
      )
    })
    
    # Hypotheses inputs
    output$hypotheses_inputs <- renderUI({
      req(input$num_hypotheses)
      lapply(1:input$num_hypotheses, function(i) {
        textInput(
          session$ns(paste0("hypothesis_", i)),
          label = paste("Hypothesis", i, ":"),
          placeholder = "Example: a > b"
        )
      })
    })
    
    # Get parameters
    get_parameters <- reactive({
      req(input$num_params)
      params <- sapply(1:input$num_params, function(i) {
        input[[paste0("param_", i)]]
      })
      names(params) <- letters[1:input$num_params]
      params
    })
    
    # Get covariance matrix
    get_covariance_matrix <- reactive({
      req(input$num_params)
      cov_matrix <- matrix(
        sapply(1:input$num_params, function(i) {
          sapply(1:input$num_params, function(j) {
            input[[paste0("cov_", i, "_", j)]]
          })
        }),
        nrow = input$num_params
      )
      colnames(cov_matrix) <- letters[1:input$num_params]
      rownames(cov_matrix) <- letters[1:input$num_params]
      cov_matrix
    })
    
    # Get hypotheses
    get_hypotheses <- reactive({
      req(input$num_hypotheses)
      hypotheses <- lapply(1:input$num_hypotheses, function(i) {
        input[[paste0("hypothesis_", i)]]
      })
      names(hypotheses) <- paste0("h", 1:input$num_hypotheses)
      hypotheses[sapply(hypotheses, function(x) !is.null(x) && x != "")]
    })
    
    # Validate matrix
    validate_covariance <- reactive({
      req(get_covariance_matrix())
      tryCatch({
        cov_matrix <- get_covariance_matrix()
        isSymmetric(cov_matrix) && all(eigen(cov_matrix)$values >= -1e-10)
      }, error = function(e) FALSE)
    })
    
    # Run analysis
    observeEvent(input$run_analysis, {
      req(get_parameters(), get_covariance_matrix(), get_hypotheses())
      
      if (!validate_covariance()) {
        showNotification(
          "The covariance matrix must be symmetric and positive semi-definite.",
          type = "error"
        )
        return()
      }
      
      tryCatch({
        gorica_res <- goric(
          object = get_parameters(),
          VCOV = get_covariance_matrix(),
          hypotheses = get_hypotheses(),
          Heq = input$add_null
        )
        
        # Results
        output$gorica_results <- renderPrint({
          print(gorica_res)
        })
        
        # Plot
        output$gorica_plot <- renderPlot({
          weights <- data.frame(
            Hypothesis = rownames(gorica_res$result),
            Weight = gorica_res$result[, 7]
          )
          
          ggplot(weights, aes(x = Hypothesis, y = Weight, fill = Hypothesis)) +
            geom_bar(stat = "identity") +
            geom_text(aes(label = sprintf("%.3f", Weight)),
                      position = position_dodge(width = 0.9),
                      vjust = -0.5) +
            labs(title = "GORICA Weights",
                 x = "Hypothesis",
                 y = "Weight") +
            scale_fill_brewer(palette = "Set3") +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none"
            )
        })
      }, error = function(e) {
        showNotification(paste("Error in analysis:", e$message), type = "error")
      })
    })
  })
}

# Main UI
ui <- navbarPage(
  title = "GORICA Analysis Tool",
  
  # Landing page
  tabPanel(
    "Home",
    div(
      class = "container",
      style = "max-width: 800px; margin: 40px auto; padding: 0 20px;",
      
      h1("GORICA Analysis Tool", align = "center"),
      p("Welcome to the GORICA (Generalized Order-Restricted Information Criterion Approximation) 
        analysis tool. This application helps you evaluate informative hypotheses without writing code. The goal of this Shiny app is to help you understand the influence of parameters and covariances on the support for your hypothesis. This app uses the goric function (Vanbrabant & Kuiper, 2024) of the restriktor (Vanbrabant & Rosseel, 2020) package in R. The author of this Shiny app is Martijn L.G. van Dam.", 
        style = "font-size: 1.2em;"),
      
      hr(),
      
      h3("What is the GORICA?"),
      p("GORICA is an abbreviation of 'Generalized Order Restricted Information Criterion Approximation'. It is an extension of the well-known Akaike's Information Criterion (AIC) and allows researchers to evaluate informative hypotheses in a broad range of statistical models.
        The GORICA quantifies the support in the data for the hypotheses at hand and indicates which hypothesis is to prefer. If weights are applied to the GORICA, one can determine how much more support there is in the data for one hypothesis over another."),

      h3("How does this Shiny app work?"),
      p("On the Information page, a description of how this app can be used is given. You will also find a description with interpretations of the output.")
    )
  ),
  
  # Information page
  tabPanel(
    "Information",
    titlePanel("Information on the App and the Output"),
    div(
      class = "container",
      style = "max-width: 1000px; margin: 0 auto; padding: 20px;",
      
      # Input Parameters Section
      h2("Input Parameters", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        h3("Number of Parameters"),
        p("Enter the number of parameters you want to analyze. These parameters represent the estimates from your analysis, such as means, regression coefficients, or other statistical parameters."),
        h3("Parameter Values"),
        p("For each parameter, you need to input:"),
        tags$ul(
          tags$li("The actual value of the parameter (e.g., mean, regression coefficient)"),
          tags$li("Parameters are automatically named using letters (a, b, c, etc.)")
        ),
        p("Please note that if you are comparing regression coefficients, you need the standardized estimates!")
      ),
      
      # Covariance Matrix Section
      h2("Covariance Matrix", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        p("The covariance matrix represents the uncertainty in your parameter estimates. It should be:"),
        tags$ul(
          tags$li("Symmetric - the value in position [i,j] equals the value in position [j,i]"),
          tags$li("Positive semi-definite - a mathematical property ensuring the matrix represents valid variances and covariances. Please note that the Shiny App also has a built-in check for positive semi-definite matrices."),
          tags$li("Diagonal elements (variances) should be positive")
        ),
        p("If you're unsure about your covariance matrix, consult the output of your statistical analysis or statistical software.")
      ),
      
      # Hypotheses Section
      h2("Defining Hypotheses", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        p("You can specify multiple hypotheses to compare. Each hypothesis should be written using:"),
        tags$ul(
          tags$li("Parameter names (a, b, c, etc.), where a corresponds to parameter 1, b to parameter 2, etc."),
          tags$li("(In)equality constraints (>, <, =)"),
          tags$li("A plus (+) or minus (-) sign accompanied by an (in)equality constraint"),
          tags$li("Numbers if needed")
        ),
        h4("Example Hypotheses:"),
        tags$ul(
          tags$li("a > b - Parameter a is greater than parameter b"),
          tags$li("a > b > c - Ordered parameters"),
          tags$li("a = b - Parameters are equal"),
          tags$li("a > 0 - Parameter is positive"),
          tags$li("-0.05 < a - b > 0.05 - The difference between two parameters is about 0")
        )
      ),
      
      # Understanding Results Section
      h2("Understanding the Results", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        h3("GORICA Values"),
        p("The output provides several key pieces of information:"),
        tags$ul(
          tags$li(strong("Log likelihood:"), " The approximated log-likelihood value of the model"),
          tags$li(strong("penalty:"), " The penalty term applied to the log likelihood. It represents the expected number of distinct parameters"),
          tags$li(strong("GORICA values:"), " Lower values indicate better support for the hypothesis"),
          tags$li(strong("Log likelihood weights:"), " The weights of the log likelihood"),
          tags$li(strong("Penalty weights:"), " The weights of the penalty (complexity)"),
          tags$li(strong("GORICA weights:"), " Represent the relative support for each hypothesis, summing to 1")
        ),
        
        h3("Interpreting the Results"),
        p("To interpret the results, consider:"),
        tags$ul(
          tags$li("The hypothesis with the lowest GORICA value is best supported by the data"),
          tags$li("GORICA weights can be divided by each other to obtain a ratio of the weights. E.g., if hypothesis 1 and 2 have weights of 0.8 and 0.2 respectively, then hypothesis 1 has 0.8/0.2 = 4 times more support in the data than hypothesis 2."),
          tags$li("If the null hypothesis (Heq) has the highest weight, your informative hypotheses might not be well-supported by the data")
        ),
      ),
      
      # Tips and Recommendations Section
      h2("Tips and Recommendations", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        tags$ul(
          tags$li("Always verify your input data carefully"),
          tags$li("Start with simple hypotheses before testing more complex ones"),
          tags$li("Consider including the null hypothesis as a reference if it is reasonable"),
          tags$li("Pay attention to error messages about the covariance matrix"),
          tags$li("Save your results for later reference")
        )
      ),
      
      h2("Questions or Remarks", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        p("If you have any questions or remarks about this Shiny app, you can reach out to the author via email: m.l.g.vandam@uu.nl")
      ),
      
      # References Section
      h2("References", style = "color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px;"),
      div(
        style = "margin-bottom: 30px;",
        tags$ul(
          tags$li(HTML("Kuiper, R. M., & Hoijtink, H. (2013). A Fortran 90 program for the generalization of the order-restricted information criterion. <i>Journal of Statistical Software, 54</i>(8), 1-19. https://doi.org/10.18637/jss.v054.i08")),
          tags$li(HTML("Vanbrabant, L., & Rosseel, Y. (2020). An introduction to restriktor: Evaluating informative hypotheses for linear models. <i>R package</i> version 0.2-800. https://doi.org/10.4324/9780429273872-14")),
          tags$li(HTML("Vanbrabant, L., & Kuiper, R. (2024). <i>restriktor: Restricted statistical estimation and inference for linear models</i>. R package version 0.6-10. https://www.restriktor.org"))
        )
      )
    )
  ),
  
  # Analysis page
  tabPanel(
    "Analysis",
    titlePanel("GORICA Analysis"),
    analysisUI("analysis")
  )
)

# Main server
server <- function(input, output, session) {
  observeEvent(input$goto_analysis, {
    updateNavbarPage(session, "GORICA Analysis Tool", selected = "Analysis")
  })
  
  analysisServer("analysis")
}

shinyApp(ui = ui, server = server)