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
k = 1*10^-6
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
  geom_vline(xintercept = vol, color = "#7fbff5", size = 0.8)+
  geom_vline(xintercept = vol2, color = "#619ed2", size = 0.8)+
  geom_vline(xintercept = vol200, color = "#3c66da", size = 0.8)+
  geom_vline(xintercept = vol200k, color = "#0c34a1", size = 0.8)+
  labs(title = expression(paste("Failure Curves (", mu, " = 1E-6)")))+
  xlab("Generations") + ylab("Fraction of Engineered Cells") +
  NULL
producers

################################## SIMULATION ################################################

#initialize fixed parameters and dataframe for function to leave values
t = seq(0, 100)
volumes<-c(vol, vol2, vol200, vol200k)

produce.this=function(b, k, t){
  burden <- data.frame("Generation" = seq(0, t))
  burden %>% 
    mutate(Fraction = (k+b*1)/(k*(exp(1))^((k+1*b)*Generation)+b*u))
}

# Define UI for application
ui <- fluidPage(
  titlePanel("Modeling Burden"),
  sidebarPanel(
    sliderInput("b", label = "Burden (%)", min = 0, max = 100, value = 10),
    radioButtons("k", label = "Escape Rate", 
                 choiceNames = list("1E-5", "1E-6", "1E-7", "1E-8"), 
                 choiceValues = list(1e-5, 1e-6, 1e-7, 1e-8)),
    checkboxGroupInput("volumes", label = "Generations to Saturate Volumes:",
                       choices = c("5mL"=vol, "2L"= vol2, "200L"=vol200,"200,000L"=vol200k)),
    sliderInput("t", label = "Number of Generations Elapsed", min = 50, max = 150, value = 100)
  ),
  mainPanel(
    plotOutput(outputId = "produce")
  )
)

server <- function(input, output) {
  output$produce <- renderPlot({
    b <- input$b/100
    k <- input$k
    t<- input$t
    volumes<-input$volumes
    plot_this<- produce.this(b, as.numeric(k), t) 
    ggplot(data= plot_this, aes(x=Generation, y = Fraction))+geom_line(color = "#ff3232", size = 1.5)+
      geom_vline(xintercept = as.numeric(volumes), color = "#696969", size = 0.8)+
      labs(title = "Production Curve", x = "Generations", y = "Fraction of Producing Cells") +
      NULL
  }, 
  height = 600)
}


# Run the application 
shinyApp(ui = ui, server = server)