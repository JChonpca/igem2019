
library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)
library(deSolve)
library(adaptivetau)

  ######################################## BURDEN MODEL ###############################################
  ##### Plot fraction of engineered cells in time as dictated by burden (b) and mutation rate (u) #####
  #####################################################################################################

#example parameters to test ODE function comment out and 
#run with burden mod function to test 
# test_pars<- tibble(
#   b  = 0.1,  # burden as reduction in growth rate
#   u  = 1e-5,  # mutation rate
#   fraction = .5,
#   this.generation = 200
# )
  
#I need more comments to help me understand what the code is doing -SPL
#Okay I get it now. But I had to change the structure of vectors / data frames
#for this different things to get it to work
burdenmod <- function(t, y, p) {
  ## This function contains the ODE to be passed to the ODE solver
  with(as.list(c(y, p)), {
    dEt <- ((1-b) * Et) - (u * (1-b) * Et)  # rate at which engineered cells decrease while growing
    dFt <- Ft + (u * (1-b) * Et)  # rate at which failed cells accumulate
    return(list(c(dEt, dFt)))
  })
}

plot.ode.mod<- function(pars) {
#this function takes a DATA FRAME of parameters
#this function reparses that data frame into a new one. Maybe not necessary now
# pars <- tibble(
#   b=p$b,
#   u=p$u,
#   fraction=p$fraction,
#   this.generation=p$this.generation
# )
print(pars)
#initial values for Et and Ft
yini  <- c(Et = 1, Ft = 0)

#build a sequence of times for the solver to use. I removed and extra element in this function
times <- seq(0, pars$this.generation[1], by = 1)

## Run model
#Use the ode solver with yini, times, burdenmod function, and parameter df
out   <- ode(yini, times, burdenmod, pars)

#convert ODE output into data frame for further processing
out <- data.frame(out)

#calculate fraction of engineered cells 
out$fraction <- (out$Et/(out$Et+out$Ft))
#calculate generations from number of cell doublings
out$generation <- log2(out$Et+out$Ft) 
#filter output into new dataframe outout below certain # of generations
#why is this necessary? -SPL
outout <- out %>% filter(generation<=pars$this.generation[1])

## Plot output of model
plot.mod <- ggplot(outout, aes(x=generation, y=fraction)) +
  geom_line() +
  ylim(c(0, 1)) +
  xlim(c(0, 150))
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
    #checkboxGroupInput("volumes", label = "Saturation Volumes:",
    #                   choices = c("5mL"=vol, "2L"= vol2, "200L"=vol200,"200,000L"=vol200k)),
    sliderInput("this.generation", label = "Generations", min = 50, max = 200, value = 100)
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
      this.generation = input$this.generation
      #volumes = input$volumes
      )
    plot.ode.mod(pars)
    #ggplot(data = mtcars, aes(x = mpg, y = disp))
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

