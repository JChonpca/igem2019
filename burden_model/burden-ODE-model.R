library(tidyverse)
library(shiny)
library(deSolve)

#These lines should be deleted, they are global variables but we don't want them as global
# pars<- c(
#   b = 0.4,   # burden value
#   u = 1e-5   # mutation rate
# )
#I added some inputs to the 

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
  ## Run model
  
  out   <- ode(yini, times, burdenode, pars)
  out <- data.frame(out)                    # save output into dataframe
  
  out$fraction <- (out$Et/(out$Et+out$Ft))  # fraction of engineered cells remaining in population
  
  out$generation <- log2(out$Et+out$Ft)     # number of cell doublings
  return(out)
}

## Plot model
# create vector of generations to saturate different volumes
vol={log2((5*10^9)*(.005*1000))}
vol2={log2((5*10^9)*(2*1000))}
vol200={log2((5*10^9)*(200*1000))}
vol200k={log2((5*10^9)*(200000*1000))}
volumes<-c(vol, vol2, vol200, vol200k)

## Plot output of model
# plot.mod <- ggplot(out, aes(x=generation, y=fraction)) +
#   geom_line() +
#   ylim(c(0, 1)) +
#   xlim(c(0, 200)) +
#   geom_vline(xintercept = as.numeric(volumes), color = "#696969")
# plot.mod



# Define UI for application
ui <- fluidPage(
  titlePanel("Modeling Burden"),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = 10),
    radioButtons("k", label = "Escape Rate", 
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
    k <- as.numeric(input$k)
    volumes <-input$volumes
    out.here <- all_the_stuff(burden = b, mutation_rate = k)
    ggplot(data= out.here, 
           aes(x=generation, 
               y = fraction))+
    geom_line(color = "#ff3232", size = 1.5)+
    geom_vline(xintercept = as.numeric(volumes), 
                color = "#696969", 
                size = 0.8)+
    labs(title = "Production Curve", 
         x = "Generations", 
         y = "Fraction of Producing Cells") +
    NULL
  }, 
  height = 600)
}


# Run the application 
shinyApp(ui = ui, server = server)
