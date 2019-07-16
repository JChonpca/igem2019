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


##  Trying to model a situation where a capacity monitor mRNA and a synthetic construct mRNA compete for
##  a fixed concentration of ribosomes within a host cell. 

##  The app should allow users to adjust translational strength (RBS strength)
##  and the mRNA size of the construct to output the predicted production of the contruct and
##  burden (translational resources consumed) as dictated by the model.

################################# ODE Construction ######################################
## Capacity monitor parameters

## X2 = mRNA size (bp) **Adjustable
## X3 = RBS strength **Adjustable
## X4 = "lumped parameter" gamma, translational resources consumed (s^-1) **Adjustable
## G = average concentration of free ribosomes **Adjustable
## M = average concentration of mRNA transcripts in cytoplasm **Adjustable
## n = number of positions as mRNA size/30

G = 2500
M = 900
n = X2/30

X2 = 720 
X3 = 1
X4 = 1

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
  
  out <- ode(y = c(G=2500, Y0=0, Y1=0, Y2=0, Yn=0, Prot=0), times = seq(0,1000), translational_ode, parms)
  out = data.frame(out)
  out
  write.csv(out, paste0("translational.csv"))

  ########################  CREATE SHINY  ################################################################
  
  
  outout = gather(out, key = "type", value = "n" X2, X3, X4, X5, X6, X7, X8) %>% 
    select(-total, -time)
  outout$type = factor(outout$type)
  levels(outout$type) = c("mRNA size", 
                          "", 
                          "deletion-plasmid", 
                          "wt-plasmid-integration", 
                          "satellite-plasmid-integration", 
                          "deletion-plasmid-integration", 
                          "no-plasmid-integration",
                          "no-plasmid"
  )
  outout = outout %>% filter(generation<=p$end_generation)
  
  log_plot = ggplot(outout, aes(x=generation,y=n, color=type)) + geom_line() + scale_y_log10(limits=c(1E-10, 1))
  
  linear_plot = ggplot(outout, aes(x=generation,y=n, color=type)) + geom_line() + ylim(c(0, 1))
  
  # Can only return one plot object, so combine them  
  grid.arrange(linear_plot, log_plot, nrow = 2)
}


# Define UI for application
ui <- fluidPage(
  sidebarPanel(
    sliderInput("w1", label = "wt-plasmid-fitness", value = w1),
    sliderInput("w2", label = "satellite-plasmid-fitness", value = w2),
    sliderInput("w3", label = "deletion-plasmid-fitness", value = w3),
    sliderInput("w4", label = "wt-plasmid-integration-fitness", value = w4),
    sliderInput("w5", label = "satellite-plasmid-integration-fitness", value = w5),
    sliderInput("w6", label = "deletion-plasmid-integration-fitness", value = w6),
    sliderInput("w7", label = "no-plasmid-integration-fitness", value = w7),
    sliderInput("w8", label = "no-plasmid-fitness", value = w8),
    numericInput("u_del", label = "deletion-mutation-rate", value = u_del),
    numericInput("u_sat", label = "satellite-mutation rate", value = u_sat),
    numericInput("u_int", label = "integration-mutation-rate", value = u_int),
    numericInput("u_wt_del_loss", label = "wt-or-deletion-plasmid-loss-rate", value = u_wt_del_loss),
    numericInput("u_sat_loss", label = "satellite-plasmid-loss-rate", value = u_sat_loss),
    numericInput("end_generation", label = "end generation", value = 100)
  ),
  mainPanel(
    plotOutput("plot")
  )
)

# Define server logic 
server <- function(input, output) {
  
  #  output$satellite_plasmid <- renderPlot({
  #    parms <- c(w1=input$w1, w2=input$w2)
  #    out <- ode(y = c(X1=1, X2=1), times=seq(0, 10, .1), satellite_plasmid, parms)
  #    matplot.0D(out)
  #  })
  
  output$plot <- renderPlot({
    if (input$w8 != 0) {
      min_w = 0.5*min(input$w1, input$w2, input$w3, input$w4, input$w5, input$w6, input$w7, input$w8)
      fixed_w8 = input$w8 / min_w
    } else {
      min_w = 0.5*min(input$w1, input$w2, input$w3, input$w4, input$w5, input$w6, input$w7)
      fixed_w8 = 0
    }
    
    parms <- data.frame( 
      w1=input$w1/min_w, 
      w2=input$w2/min_w, 
      w3=input$w3/min_w,
      w4=input$w4/min_w,
      w5=input$w5/min_w,
      w6=input$w6/min_w,
      w7=input$w7/min_w,
      w8=fixed_w8,
      u_del=input$u_del,
      u_sat=input$u_sat,
      u_int=input$u_int,
      u_wt_del_loss=input$u_wt_del_loss,
      u_sat_loss=input$u_sat_loss,
      end_generation=input$end_generation)
    plot_ode(parms) 
  }, 
  height = 600)
}

# Code for testing outside of Shiny

if (w8 != 0) {
  min_w = 0.5*min(w1, w2, w3, w4, w5, w6, w7, w8)
  fixed_w8 = w8 / min_w
} else {
  min_w = 0.5*min(w1, w2, w3, w4, w5, w6, w7)
  fixed_w8 = 0
}
parms <- data.frame( 
  w1=w1/min_w, 
  w2=w2/min_w, 
  w3=w3/min_w,
  w4=w4/min_w,
  w5=w5/min_w,
  w6=w6/min_w,
  w7=w7/min_w,
  w8=fixed_w8,
  u_del=u_del,
  u_sat=u_sat,
  u_int=u_int,
  u_wt_del_loss=u_wt_del_loss,
  u_sat_loss=u_sat_loss,
  end_generation=100
)
plot_ode(parms) 

# Run the application 
shinyApp(ui = ui, server = server)
