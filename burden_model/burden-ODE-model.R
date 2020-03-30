
#########################
# SHINY * IN PROGRESS * #
#   -GAM                #
#########################


library(tidyverse)
library(shiny)
library(deSolve)
library(adaptivetau)

burdenode <- function(t, y, pars) {
  with(as.list(c(y, pars)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells (Et) grow/mutate
    dFt <- Ft + (u * (1-b) * Et)            # rate at which failed cells (Ft) accumulate/grow
    return(list(c(dEt, dFt)))
  })
}

all_the_stuff <- function(burden, mutation_rate) {
  yini  <- c(Et = 1, Ft = 0)  # initial values of Et and Ft. Begin with one engineered cell.
  times <- seq(0,200, 0.1)    # time sequence for ODE solver
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

all_the_things <- function(burden, mutation_rate){
 for (this.seed in 1:20) {
   set.seed(this.seed) # set random number generator seed to be reproducible
   pars = list(b = burden, u = mutation_rate);
   sim.results = ssa.adaptivetau(
     init.values = c(e = 1, f = 0),
     transition.list,
     rates.function,
     pars,
     tf = 200
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
# create vector of cell doublings to saturation
vol = {log2((5 * 10^9) * (.005 * 1000))}
vol2 = {log2((5 * 10^9) * (2 * 1000))}
vol200 = {log2((5 * 10^9) * (200 * 1000))}
vol200k = {log2((5 * 10^9) * (200000 * 1000))}
volumes <- c(vol, vol2, vol200, vol200k)

# Define UI for application
ui <- fluidPage(
  titlePanel("Modeling Burden"),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = 10),
    radioButtons("u", label = "Escape Rate", 
                 choiceNames = list("1E-5", "1E-6", "1E-7", "1E-8"), 
                 choiceValues = list(1e-5, 1e-6, 1e-7, 1e-8)),
    checkboxGroupInput("volumes", label = "Generations to Saturate Volumes:",
                       choices = c("5mL"=vol, "2L"= vol2, "200L"=vol200,"200,000L"=vol200k))
  ),
  mainPanel(
    plotOutput(outputId = "produce")
  )
)

server <- function(input, output) {
  output$produce <- renderPlot({
    b <- as.numeric(input$b)/100
    u <- as.numeric(input$u)
    volumes <-input$volumes
    out.here <- all_the_stuff(burden = b, mutation_rate = u)
    out.there <- all_the_things(burden = b, mutation_rate = u)
    ggplot(data = out.here, 
           aes(x = generation, 
               y = fraction)) +
      geom_line(color = "red", size = 1.5) +
      geom_line(data = out.there, aes(x = generation, y = fr.e, group = seed)) +
      geom_line()+
      geom_vline(xintercept = as.numeric(volumes), 
                 color = "#696969", 
                 size = 0.8) +
      labs(title = "Production Curve", 
           x = "Generations", 
           y = "Fraction of Producing Cells") +
      NULL
  }, 
  height = 600)
}

# Run the application 
shinyApp(ui = ui, server = server)
