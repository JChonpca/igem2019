# Shiny App for burden-model

# Authors: Ginny Mortensen / Jeff Barrick

library(tidyverse)
library(shiny)
library(deSolve)
library(adaptivetau)

## ODE Model

burdenode <- function(t, y, pars) {
  with(as.list(c(y, pars)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells (Et) grow/mutate
    dFt <- Ft + (u * (1-b) * Et)            # rate at which failed cells (Ft) accumulate/grow
    return(list(c(dEt, dFt)))
  })
}

ode_model_results <- function(burden, mutation_rate, simulation.time) {
  yini  <- c(Et = 1, Ft = 0)  # initial values of Et and Ft. Begin with one engineered cell.
  times <- seq(0,simulation.time, 0.1)    # time sequence for ODE solver
  pars <- c(b = burden, u = mutation_rate)
  out <- ode(yini, times, burdenode, pars) # run solver
  out <- data.frame(out)                     # save output into dataframe
  out$fraction <- (out$Et/(out$Et+out$Ft))   # fraction of engineered cells remaining in population
  out$generation <- log2(out$Et+out$Ft)      # number of cell doublings
  return(out)
}

##########################################################################################################

transition.list = list(
  c(e = +1),           # trans 1: e cells grow
  c(f = +1),           # trans 3: f cells grow
  c(e = -1, f = +1)    # trans 2: e cells mutate into f cells
)

rates.function <- function(x, pars, t) {
  return(c(
    (1 - pars$b) * x["e"],             # rate of e cells growing
    x["f"],                              # rate of f cells growing
    (1 - pars$b) * x["e"] * pars$u   # rate of e cells mutating into f cells
  )) 
}

sim.results.df = data.frame() # save as a dataframe

stochastic_model_results <- function(burden, mutation_rate, n, simulation.time){
  for (this.seed in 1:n) {
    pars = list(b = burden, u = mutation_rate);
    sim.results = ssa.adaptivetau(
      init.values = c(e = 1, f = 0),
      transition.list,
      rates.function,
      pars,
      tf = simulation.time
    )
    
    # reformat and calculate generations
    this.sim.results.df =   data.frame (
      time = sim.results[,c("time")],
      generation = log2(sim.results[,c("e")] + sim.results[,c("f")]),
      e = sim.results[,c("e")],
      f = sim.results[,c("f")],
      fr.e = sim.results[,c("e")] / (sim.results[,c("e")] + sim.results[,c("f")]),
      seed = this.seed,
      b = burden,
      u = mutation_rate,
      method = "stochastic"
    )
    sim.results.df = rbind(sim.results.df, this.sim.results.df)
  }
  sim.results.df$seed = as.factor(sim.results.df$seed)
  return(sim.results.df)
}

########################################################################################################

## Plot model

# Define UI for application
ui <- fluidPage(
  titlePanel("Evolutionary Failure of an Engineered Cell Population"),
  tags$p(
    "For  information about interpreting these simulations see the ",
    tags$a(href="https://2019.igem.org/Team:Austin_UTexas", 
           "2019 UT Austin iGEM Team wiki"),
    " and the ",
    tags$a(href="http://parts.igem.org/Part:BBa_K3174002",
      "guide on the iGEM Registry")
  ),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = 30),
    sliderInput("u", label = "Log10 mutation rate", min = -12, max = 0, value = -5, step=0.5),
    sliderInput("n", label = "Number of stochastic simulations (red)", min = 1, max = 100, value = 10),
    checkboxGroupInput(
      "scale_guides", 
      label = "Culture scale guides (blue)",
      selected = list("colony", "test tube", "flask", "bioreactor"),
      choices = list("colony", "test tube", "flask", "bioreactor")
    ),
    sliderInput("time", label = "Max cell doublings to graph", min = 50, max = 500, value = 100, step=50),
    checkboxInput("ode", label = "Show ODE model for comparison (black)", value=TRUE),
    textInput("seed", "Random seed", value = "2716", width = NULL, placeholder = NULL)
  ),
  mainPanel(
    plotOutput(outputId = "produce")
  )
)

server <- function(input, output) {
  
  output$produce <- renderPlot(
    expr = {
      
      set.seed(input$seed)
      b <- as.numeric(input$b)/100
      u <- as.numeric(10^input$u)
      volumes <-input$volumes
      stochastic.model.results <- stochastic_model_results(burden = b, mutation_rate = u, n = input$n, simulation.time = input$time)
      ode.model.results <- ode_model_results(burden = b, mutation_rate = u, simulation.time = input$time)
      p = ggplot(data = stochastic.model.results, 
                 aes(x = generation, 
                     y = fraction), clip="off") +
        geom_line(data = stochastic.model.results, aes(x = generation, y = fr.e, group = seed), color="red", alpha=0.5, size = 1) +
        
        labs(x = "Cell Doublings", 
             y = "Fraction of Engineered Cells Remaining") +
        scale_x_continuous(limits = c(0,input$time), expand = expansion(mult=0.01) )+
        scale_y_continuous(limits = c(0,1), expand = expansion(add=0.02), breaks=seq(0, 1, by=0.1) ) + 
        theme_classic() +
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
        theme(text = element_text(size = 18)) +
        NULL
      
      if (input$ode) {
        p = p + geom_line(data = ode.model.results, color = "black", size = 1.5)
      }
      
      if ("colony" %in% input$scale_guides) {
        p = p + geom_vline(xintercept = 22.93156857, 
                           color = "#6BAED5", 
                           size = 1)
      }
      
      if ("test tube" %in% input$scale_guides) {
        p = p + geom_vline(xintercept = 34.21928095, 
                           color = "#4392C5", 
                           size = 1)
      }
      
      if ("flask" %in% input$scale_guides) {
        p = p + geom_vline(xintercept = 41.12617154, 
                           color = "#2572B5", 
                           size = 1)
      }
      
      if ("bioreactor" %in% input$scale_guides) {
        p = p + geom_vline(xintercept = 56.15084952, 
                           color = "#1A4990", 
                           size = 1)
      }

      p
    }, 
    height = 720)
}

# Run the application 
shinyApp(ui = ui, server = server)
