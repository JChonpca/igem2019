
######################################## BURDEN MODEL ###############################################
###### Fraction of engineered cells in time as predicted by burden (b) and mutation rate (u) ########
#####################################################################################################

library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)
library(deSolve)

## Build ODE
burdenmod <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells decrease while growing
    dFt <- Ft + (u * (1-b) * Et)  # rate at which failed cells accumulate
    return(list(c(dEt, dFt)))
  })
}

plot.ode.mod<- function(pars) {

  print(pars)
  #initial values for Et and Ft
  yini  <- c(Et = 1, Ft = 0)
  
  #build a sequence of times for the solver to use.
  times <- seq(0, pars$this.generation[1], by = 1)
  
## Run model
  #Use the ode solver with initial values, time sequence, function, and parameter dataframe
  out   <- ode(yini, times, burdenmod, pars)
  
  #convert ODE output into data frame for further processing
  out <- data.frame(out)
  
  #calculate fraction of engineered cells 
  out$fraction <- (out$Et/(out$Et+out$Ft))
  
  #calculate generations from number of cell doublings
  out$generation <- log2(out$Et+out$Ft) 

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

}

## Build Shiny app

## Define UI for application
ui <- fluidPage(
  titlePanel("Modeling Burden"),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = b),
    radioButtons("u", label = "Mutation Rate", 
                 choiceNames = list("1E-5", "1E-6", "1E-7", "1E-8"), 
                 choiceValues = list(u, 1e-6, 1e-7, 1e-8)),
    checkboxGroupInput("volumes", label = "Saturation Volumes:",
                     choices = c("5mL"=vol, "2L"= vol2, "200L"=vol200,"200,000L"=vol200k))
  ),
  mainPanel(
    plotOutput("plotmod")
  )
)

server <- function(input, output) {
  output$plotmod <- renderPlot({
    pars <- data.frame(
      b = input$b / 100,
      u = as.numeric(input$u),
      volumes = input$volumes
    )
    plot.ode.mod(pars)
  }, 
  height = 600)
}


## Run the application 
shinyApp(ui = ui, server = server)
