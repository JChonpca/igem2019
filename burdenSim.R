#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(deSolve)
library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

################################# ODE CONSTRUCTION ######################################
## Capacity monitor parameters

## X2 = mRNA size (bp)
## X3 = RBS strength **Adjustable
## X4 = "lumped parameter" of translational resources consumed (s^-1) **Adjustable

## G = average concentration of free ribosomes
## M = average concentration of mRNA transcripts in cytoplasm
## n = number of positions as mRNA size/30

X2 = 90
X3 = 2
X4 = 2

G = 312.5
M = 112.5
n = X2/30

## Ribosome parameters

## X5 = binding rate
## X6 = unbinding rate
## X7 = first synthesis rate

a1 = .00001
a2 = 200
b = 1

X5 = a1*X3
X6 = a2/X3
X7 = b*X3

## Y0 = average concentration for ribosome bound to RBS of all synthetic gene mRNAs
## Yn = average concentration of ribosomes bound to RBS across all mRNAs in position n

##Model kinetics: n+3 number of states.

translational_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dG <- -(M*X5*G)*(1-(Y0/M)) + X6*Y0 + X4*Yn 
    dY0 <- (M*X5*G)*(1-(Y0/M)) - X6*Y0 - X7*Y0*(1-(Y0/M))
    dY1 <- X7*Y0*(1-(Y1/M)) - X4*Y1*(1-(Y2/M))
    dY2 <- X4*Y1*(1-(Y2/M)) - X4*Y2*(1-(Yn/M))
    dYn <- X4*Y2*(1-(Yn/M)) - X4*Yn
    dProt <- X4*Yn
    return(list(c(dG, dY0, dY1, dY2, dYn, dProt)))
  })
}

plot_ode <- function(p) 

parms <- c(
          X2=X2,
          X3=X3,
          X4=X4,
          X5=X5,
          X6=X6,
          X7=X7,
          M=M
          )
  
print(parms)
  
##Run model
  
  out <- ode(y = c(G=312.5, Y0=0, Y1=0, Y2=0, Yn=0, Prot=0), times = seq(0,1000), translational_ode, parms)
  out = data.frame(out)
  out
  write.csv(out, paste0("translational.csv"))

  ########################  CREATE SHINY  ################################################################
 
 
#Create levels  
  outout = gather(out, key = "type", value = "concentration", G, Y0, Y1, Y2, Yn, Prot)
  outout$type = factor(outout$type)
  levels(outout$type) = c("free ribosome concentration", 
                          "initial state", 
                          "state 1", 
                          "state 2", 
                          "final state",
                          "protein output"
  )
  
#Create plot  
  linear_plot = ggplot(data = outout, aes(x = times, y = concentration, color = type)) + geom_line() + 
    labs(title = "Translational Burden Simulator", x = "RBS Strength", y = "Translational Resources Consumed")

# Define UI for application
ui <- fluidPage(
  sidebarPanel(
    sliderInput("X2", label = "mRNA size (bp)", min = -10, max = 100, round = FALSE, value = X2),
    sliderInput("X3", label = "RBS strength", min = -10, max = 100,  round = FALSE, value = X3),
    sliderInput("X4", label = "translational resources consumed", min = -10, max = 100, round = FALSE, value = X4)
  ),
  mainPanel(
    plotOutput("plot")
  )
)

# Define server logic 
server <- function(input, output) {
  output$plot <- renderPlot({
    parms <- data.frame( 
      X2 = input$X2,
      X3 = input$X3,
      X4 = input$X4,
      X5=input$X5,
      X6=input$X6,
      X7=input$X7,
      M=input$M
    )
    plot_ode(parms) 
  }, 
  height = 600)
}

# Run the application 
shinyApp(ui = ui, server = server)



