
library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)
library(deSolve)
library(adaptivetau)

  ######################################## BURDEN MODEL ###############################################
  ##### Plot fraction of engineered cells in time as dictated by burden (b) and mutation rate (u) #####
  #####################################################################################################

## this syntax is weird, do you want all of this commented out...?
#pars<- c(
#  b  = 0.1  # burden as reduction in growth rate
#  u  = 1e-5  # mutation rate
#)
  
#I need more comments to help me understand what the code is doing -SPL
burdenmod <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells decrease while growing
    dFt <- Ft + (u * (1-b) * Et)  # rate at which failed cells accumulate
    return(list(c(dEt, dFt)))
  })
}

plot.ode.mod<- function(p) {

pars <- c(
  b=p$b,
  u=p$u,
  fraction=p$fraction,
  this.generation=p$this.generation
)
yini  <- c(Et = 1, Ft = 0)

times <- seq(0, p$this.generation, 0, 01)
#times <- seq(0,200)

## Run model
out   <- ode(yini, times, burdenmod, pars)
summary(out)

out <- data.frame(out)
out$fraction <- (out$Et/(out$Et+out$Ft))  # fraction of engineered cells
out$generation <- log2(out$Et+out$Ft)  # number of cell doublings
outout = outout %>% filter(generation<=p$this.generation)

## Plot model
plot.mod = ggplot(outout, aes(x=generation, y=fraction)) + geom_line() + ylim(c(0, 1)) + xlim(c(0, 150))
plot.mod
}

## Build Shiny app
v= 0.005
v2=2
v200=200
v200k=200000
vol={log2((5*10^9)*(v*1000))}
vol2={log2((5*10^9)*(v2*1000))}
vol200={log2((5*10^9)*(v200*1000))}
vol200k={log2((5*10^9)*(v200k*1000))}
volumes<-c(vol, vol2, vol200, vol200k)  # generation of saturation for volumes

## Define UI for application
ui <- fluidPage(
  titlePanel("Modeling Burden"),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = b),
    radioButtons("u", label = "Mutation Rate", 
                 choiceNames = list("1E-5", "1E-6", "1E-7", "1E-8"), 
                 choiceValues = list(u, 1e-6, 1e-7, 1e-8)),
    checkboxGroupInput("volumes", label = "Saturation Volumes:",
                       choices = c("5mL"=vol, "2L"= vol2, "200L"=vol200,"200,000L"=vol200k)),
    sliderInput("this.generation", label = "Generations", min = 50, max = 200, value = 100)
  ),
  mainPanel(
    plotOutput("plotmod")
  )
)

server <- function(input, output) {
  output$plotmod <- renderPlot({
    pars = data.frame(
      b = input$b,
      u = input$u,
      this.generation = input$this.generation,
      volumes = input$volumes
      )
    plot.ode.mod(pars)
  }, 
  height = 600)
}


## Run the application 
shinyApp(ui = ui, server = server)


############################# IN PROGRESS - STOCHASTIC SIMULATION ####################################
# commenting out unused code -SPL
#parms<-list(b=0.1, u=1e5);
#stochastic_burden=ssa.adaptivetau(c(1,0),
#                  matrix(c(,0, -1,1), nrow=2), plot.ode.mod, parms, tf=200)
#stochastic_burden_plot

