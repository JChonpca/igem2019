#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)
library(RColorBrewer)

################################### GRAPHIC MODEL ###########################################
######Plot fraction of producers in time as dictated by burden (b) and escape rate (k)#######

##Create parameters
u = 1 #specific growth rate of producers
k = 1*10^-5
b = .1
e = exp(1)
t = seq(0,100)

##plug into curves
eq = function (t){(k+0*u)/(k*e^((k+u*0)*t)+0*u)}
eq1 = function (t){(k+.1*u)/(k*e^((k+u*.1)*t)+.1*u)}
eq2 = function (t){(k+.2*u)/(k*e^((k+u*.2)*t)+.2*u)}
eq3 = function (t){(k+.3*u)/(k*e^((k+u*.3)*t)+.3*u)}
eq4 = function (t){(k+.4*u)/(k*e^((k+u*.4)*t)+.4*u)}
eq4 = function (t){(k+.5*u)/(k*e^((k+u*.5)*t)+.5*u)}

##create generation @ volume reference lines (liters)
v= 0.005
v2=2
v200=200
v200k=200000

vol={log2((5*10^9)*(v*1000))}
vol2={log2((5*10^9)*(v2*1000))}
vol200={log2((5*10^9)*(v200*1000))}
vol200k={log2((5*10^9)*(v200k*1000))}

##plot model
producers= ggplot(data.frame(x=t), aes(x=x)) +
  stat_function(fun=eq, geom="line", color= "#ffcccc", size = 1.5)+
  stat_function(fun=eq1, geom="line", color= "#ff9999", size = 1.5)+
  stat_function(fun=eq2, geom="line", color= "#ff3232", size = 1.5)+
  stat_function(fun=eq3, geom="line", color= "#cc0000", size = 1.5)+
  stat_function(fun=eq4, geom="line", color= "#7f0000", size = 1.5)+
  geom_vline(xintercept = vol, "#3f3f3f", size = 0.8)+
  geom_vline(xintercept = vol2, color = "#696969", size = 0.8)+
  geom_vline(xintercept = vol200, color = "#a8a8a8", size = 0.8)+
  geom_vline(xintercept = vol200k, color = "#d3d3d3", size = 0.8)+
  labs(title = "Production Curves for Differing Burden Values (Escape Rate = 10^-5)")+
  xlab("Generations") + ylab("Fraction of Producing Cells") +
  NULL
producers

################################## SIMULATION ################################################

#initialize fixed parameters and dataframe for function to leave values
u = 1
e = exp(1)
t = seq(0, 100)

produce.this=function(b, k){
  burden <- data.frame("Generation" = t)
  burden %>% 
    mutate(Fraction = (k+b*u)/(k*e^((k+u*b)*Generation)+b*u))
}

# Define UI for application
ui <- fluidPage(
  titlePanel("Modeling Burden"),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = 10),
    sliderInput("k", label = "Escape Rate", min = 10^-8, max = 10^-5, value = 10^-5)
),
  mainPanel(
    plotOutput(outputId = "produce")
  )
)

server <- function(input, output) {
  output$produce <- renderPlot({
      b <- input$b/100
      k <- input$k
    plot_this<- produce.this(b, k) 
    ggplot(data= plot_this, aes(x=Generation, y = Fraction))+geom_line(color = "#ff3232", size = 1.5)+
      geom_vline(xintercept = vol, "#3f3f3f", size = 0.8)+
      geom_vline(xintercept = vol2, color = "#696969", size = 0.8)+
      geom_vline(xintercept = vol200, color = "#a6a6a6", size = 0.8)+
      geom_vline(xintercept = vol200k, color = "#a8a8a8", size = 0.8)+
      labs(title = "Production Curves for Differing Burden Values",x = "Generations", y = "Fraction of Producing Cells") +
      NULL
  }, 
  height = 600)
}


# Run the application 
shinyApp(ui = ui, server = server)
