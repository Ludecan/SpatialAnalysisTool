#setwd('D:/Workspace/TutorialShiny/spatialAnalysisTool')
# Load packages ----
library(shiny)
library(gstat)
library(sp)
source('motionLib/mapping/mapearEx.r')
data(meuse)
data(meuse.grid)
coordinates(meuse) = ~x+y
proj4string(meuse) <- CRS("+init=epsg:28992")
coordinates(meuse.grid) = ~x+y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
gridded(meuse.grid) <- TRUE
meuse$logZinc <- log(meuse$zinc)

obsPlot <- mapearPuntosGGPlot(meuse, continuo = TRUE, zcol = 'logZinc', dibujar = FALSE, 
                              titulo = 'Logarithm of Zinc Concentration[ppm] for the Meuse Dataset')

lzn.vgm = variogram(object = logZinc~1, data = meuse)
maxGamma <- round(max(lzn.vgm$gamma), digits = 1)
maxDist <- round(max(lzn.vgm$dist))
lzn.vgm.map = variogram(object = logZinc~1, data = meuse, map=T, cutoff=maxDist, width=maxDist/15)
# This is an abuse, it's acutally assigning any cartesian CRS that we need. We use Uruguay's UTM projection
proj4string(lzn.vgm.map$map) <- CRS('+proj=utm +zone=21 +south +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0 +units=m +no_defs')
# We take away points in the map calculated with less than 5 observations to avoid excessive noise
lzn.vgm.map$map$var1[lzn.vgm.map$map$np.var1 < 5] <- NA
# I don't like sp's base plotting methods so I'm using my own based on ggplot2
lzn.vgm.map <- mapearGrillaGGPlot(grilla = lzn.vgm.map$map, continuo = T, titulo = 'Variogram map')

#For debugging
#lzn.fit = fit.variogram(object = lzn.vgm, model = vgm(psill = 0.7, model = 'Exp', range = 500, nugget = 0))
#lzn.kriged = krige(logZinc~1, meuse, meuse.grid, model = lzn.fit)
#mapearGrillaGGPlot(lzn.kriged, continuo = T)

availableModels <- list('Exp', 'Sph', 'Gau', 'Cir', 'Pen')
names(availableModels) <- c('Exponential', 'Spherical', 'Gaussian', 'Circular', 'Pentaspherical')

# User interface ----
ui <- fluidPage(
  titlePanel(h1("Spatial Analysis Tool")),
  sidebarLayout(
    sidebarPanel(
      helpText("This tool allows modelling and interpolating Spatial Data through an Ordinary Kriging approach.
               Select a Variogram Model and Parameters to fit the Empirical Variogram in the Top-Left Panel. 
               The Bottom-Right panel will update the interpolation results to reflect your chosen model.
               You can alter the variogram model/parameters and see how they affect the interpolated values.
               Currently WIP. Next steps will allow you to import your own data."),
      selectInput('model', label = 'Variogram Model', 
                  choices = availableModels, selected = availableModels[1]),
      sliderInput('Nugget', label = 'Nugget', min = 0, max = maxGamma, value = 0),
      sliderInput('Range', label = 'Range', min = 1, max = maxDist, value = maxDist / 3),
      sliderInput('Psill', label = 'Psill', min = 0.1, max = maxGamma, value = maxGamma * 0.7),
      actionButton('AutoFit', label = 'Auto Fit'),
      helpText(" Gstat does not fit anisotropy parameters automatically so these must be specified manually."),
      sliderInput('AnisDir', label = 'Anisotropy Direction', min = 0, max = 179, value = 55),
      sliderInput('AnisRatio', label = 'Anisotropy Ratio', min = 0.01, max = 1, value = 0.5)
      ),
    mainPanel(
      fluidRow(
        column(width = 6, 
               plotOutput("vgmplot")),
        column(width = 6, 
               plotOutput("vgmmapplot"))
      ),
      fluidRow(
        column(width = 6, 
               plotOutput("obsplot")),
        column(width = 6, 
               plotOutput("krigingplot"))
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  formState <- reactive({
    if (input$model == 'Pow') {
      minRange <- 0.01
      maxRange <- 2
      rangeValue <- 1
      rangeStep <- 0.1
    } else {
      minRange <- 1
      maxRange <- maxDist
      rangeValue <- maxDist / 3
      rangeStep <- NULL
    }
    
    return(list(minRange=minRange, maxRange=maxRange, rangeValue=rangeValue, rangeStep=rangeStep))
  })
  
  updateVGM <- reactive({
    vgm(psill = input$Psill, model = input$model, range = input$Range, nugget = input$Nugget)
  })
  
  updateVGMWithAnis <- reactive({
    vgm(psill = input$Psill, model = input$model, range = input$Range, nugget = input$Nugget, anis = c(input$AnisDir, input$AnisRatio))
  })
  
  observeEvent(input$model, {
    state <- formState()
    updateSliderInput(session, inputId = 'Range', label = 'Range', min = state$minRange, max = state$maxRange, 
                      value = state$rangeValue, step = state$rangeStep)  
  })
  
  observeEvent(input$AutoFit, {
    lzn.fit = fit.variogram(object = lzn.vgm, model = updateVGM())
    nugget <- lzn.fit$psill[1]
    range <- lzn.fit$range[2]
    psill <- lzn.fit$psill[2]
    
    state <- formState()
    updateSliderInput(session, inputId = 'Range', label = 'Range', min = state$minRange, max = state$maxRange, 
                      value = state$rangeValue, step = state$rangeStep)  
    
    updateSliderInput(session, inputId = 'Nugget', label = 'Nugget', min = 0, max = maxGamma, value = nugget)
    updateSliderInput(session, inputId = 'Psill', label = 'Psill', min = 0.1, max = maxGamma, value = psill)
  })
  
  output$dbg <- renderText(paste('Model=', input$model, ', Nugget=' ,input$Nugget, 
                                 ', Range=', input$Range, ', Psill=', input$Psill, '.'))
  output$vgmplot <- renderPlot({
    # lzn.fit = fit.variogram(object = lzn.vgm, model = vgm(psill = input$Psill, model = model, range = input$Range, nugget = input$Nugget))
    lzn.fit = updateVGM()
    plot(lzn.vgm, lzn.fit, main=paste('Empirical and Model Variogram. ', input$model, ', Nugget=', input$Nugget, 
                                      ', Range=', input$Range, ', Psill=', input$Psill, '.', sep=''))
  })
  
  output$vgmmapplot <- renderPlot({
    print(lzn.vgm.map)
  })
  
  output$obsplot <- renderPlot({
    print(obsPlot)
  })
  
  output$krigingplot <- renderPlot({
    lzn.kriged = krige(logZinc~1, meuse, meuse.grid, model = updateVGMWithAnis())
    print(mapearGrillaGGPlot(lzn.kriged, continuo = T, dibujar= FALSE, titulo = 'Kriging Interpolation'))
  })
}

# Run the app
shinyApp(ui, server)
