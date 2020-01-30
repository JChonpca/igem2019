library(tidyverse)
library(shiny)
library(deSolve)

pars<- c(
  b = 0.4,   # burden value
  u = 1e-5   # mutation rate
)

burdenode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells (Et) grow/mutate
    dFt <- Ft + (u * (1-b) * Et)            # rate at which failed cells (Ft) accumulate/grow
    return(list(c(dEt, dFt)))
  })
  
  
  yini  <- c(Et = 1, Ft = 0)  # initial values of Et and Ft. Begin with one engineered cell.
  times <- seq(0,200, 0.1)    # time sequence for ODE solver
  
  ## Run model
  out   <- ode(yini, times, burdenode, pars)
  out <- data.frame(out)                    # save output into dataframe
  out$fraction <- (out$Et/(out$Et+out$Ft))  # fraction of engineered cells remaining in population
  out$generation <- log2(out$Et+out$Ft)     # number of cell doublings
}


## Plot model
# create vector of generations to saturate different volumes
vol={log2((5*10^9)*(.005*1000))}
vol2={log2((5*10^9)*(2*1000))}
vol200={log2((5*10^9)*(200*1000))}
vol200k={log2((5*10^9)*(200000*1000))}
volumes<-c(vol, vol2, vol200, vol200k)

## Plot output of model
plot.mod <- ggplot(out, aes(x=generation, y=fraction)) +
  geom_line() +
  ylim(c(0, 1)) +
  xlim(c(0, 200)) +
  geom_vline(xintercept = as.numeric(volumes), color = "#696969")
plot.mod



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
    b <- input$b/100
    k <- input$k
    volumes<-input$volumes
    out.here<- burdenode(t, y, p) 
    ggplot(data= out.here, aes(x=generation, y = fraction))+geom_line(color = "#ff3232", size = 1.5)+
      geom_vline(xintercept = as.numeric(volumes), color = "#696969", size = 0.8)+
      labs(title = "Production Curve", x = "Generations", y = "Fraction of Producing Cells") +
      NULL
  }, 
  height = 600)
}


# Run the application 
shinyApp(ui = ui, server = server)
