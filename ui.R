#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)


# Define UI for application that draws a histogram
shinyUI(
    fluidPage(
        theme = shinytheme("flatly"),
        tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")),   
        navbarPage(
            
            title = strong("translational PKPD"),
            
            tabPanel(
                strong("human PK prediction"),
                fluidRow(
                    column(4, 
                           # h3("Dose regimen"),
                           fluidRow(
                               column(6),
                               column(6,
                                      sliderInput("cycle_pk", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("amt_pk1", label = "dose (mg)", min = 0, max = 15, step = 0.5, value = c(5))
                               ),
                               column(6,
                                      sliderInput("amt_pk2", label = "dose interruption (mg)", min = 0, max = 15, step = 0.5, value = c(0))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("interval_pk1", label = "dose interval (h)", min = 0, max = 48, step = 2, value = c(24))
                               ),
                               column(6,
                                      sliderInput("interval_pk2", label = "interrupted interval (h)", min = 0, max = 48, step = 2, value = c(24))
                               )
                           ),
                           fluidRow(
                               column(6,
                                      sliderInput("n_pk1", label = "dose times", min = 1, max = 28, step = 1, value = c(1))
                               ),
                               column(6,
                                      sliderInput("n_pk2", label = "interruption times", min = 1, max = 28, step = 1, value = c(7))
                               )
                           )
                    ),
                column(8,
                       plotOutput("distPlotPK"))
                ),
                
                fluidRow(
                    column(3,
                           sliderInput("cl_pk", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(5))),
                    column(3,
                           sliderInput("v_pk", label = "V (L)", min = 10, max = 500, step = 10, value = c(100))),
                    column(3,
                           sliderInput("ka_pk", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.5))),
                    column(3,
                           sliderInput("f_pk", label = "F", min = 0.1, max = 1, step = 0.05, value = c(0.5)))
                ),
                # absolutePanel(
                #     top = 80, right = 250, width = 150, height = 10, draggable = TRUE,
                #     HTML(
                #         paste0("<strong>t<sub>", "1/2", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                #                 textOutput(outputId = "halflife", inline = T), "</code> hr</strong>")
                #          )
                # ),
                absolutePanel(
                    top = 150, right = 60, width = 130, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>C<sub>", "max", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkcmax1", inline = T), "</code> ng/mL</strong>")
                    )
                ),
                absolutePanel(
                    top = 190, right = 50, width = 180, height = 10, draggable = TRUE,
                    HTML(
                        paste0("<strong>AUC<sub>", "tau", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                               textOutput(outputId = "pkauc1", inline = T), "</code> hr*ng/mL</strong>")
                    )
                ),
                
                
                # absolutePanel(
                #     top = 123, right = 8, width = 220, height = 10, draggable = FALSE,
                #     HTML(paste0("<strong>AUC<sub>", textOutput(outputId = "auclower", inline = T), "-",
                #                 textOutput(outputId = "aucupper", inline = T), "hr</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                #                 textOutput(outputId = "pkauc1", inline = T), "</code> hr*mg/L</strong>"))
                # )
                
                # absolutePanel(
                #     top = 35, right = 150, width = 200, height = 10, draggable = TRUE,
                #     sliderInput("cycle_pk", label = "Cycle number", min = 1, max = 10, step = 1, value = c(2))
                # ),
                # absolutePanel(
                #     top = 65, right = 30, width = 100, height = 10, draggable = TRUE,
                #     img(src="LOGOdMed.png", height = 50)
                # )
                # absolutePanel(
                #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                #     h3("ing")
                # ),
                absolutePanel(
                    top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                    img(src="LOGOdMed.png", height = 50)
                )
            ),
            
            tabPanel(
                strong("tumor growth"),
                fluidRow(
                    fluidRow(
                        column(2, 
                               # h3("Dose regimen"),
                               fluidRow(
                                   column(12),
                                   br(),
                                   br(),
                                   br(),
                                   br(),
                                   br()
                                   # column(6,
                                   #        sliderInput("cycle_tumor", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                                   # )
                               ),
                               fluidRow(
                                   column(12,
                                          sliderInput("amt_tumor1", label = "dose (mg)", min = 0, max = 15, step = 0.5, value = c(5))
                                   )
                                   # column(6,
                                   #        sliderInput("amt_tumor2", label = "dose interruption (mg)", min = 0, max = 15, step = 0.5, value = c(5))
                                   # )
                               ),
                               fluidRow(
                                   column(12,
                                          sliderInput("interval_tumor1", label = "dose interval (h)", min = 0, max = 48, step = 2, value = c(24))
                                   )
                                   # column(6,
                                   #        sliderInput("interval_tumor2", label = "interrupted interval (h)", min = 0, max = 48, step = 2, value = c(24))
                                   # )
                               ),
                               fluidRow(
                                   column(12,
                                          sliderInput("n_tumor1", label = "dose times", min = 1, max = 28, step = 1, value = c(28))
                                   )
                                   # column(6,
                                   #        sliderInput("n_tumor2", label = "interruption times", min = 1, max = 28, step = 1, value = c(7))
                                   # )
                               )
                        ),
                        column(5, plotOutput("distPlotTumorPK")),
                        column(5, plotOutput("distPlotTumorPD"))
                    ),    
                    
                    fluidRow(
                        column(3,
                               sliderInput("cl_tumor", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(5))),
                        column(3,
                               sliderInput("v_tumor", label = "V (L)", min = 10, max = 500, step = 10, value = c(100))),
                        column(3,
                               sliderInput("ka_tumor", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.5))),
                        column(3,
                               sliderInput("f_tumor", label = "F", min = 0.1, max = 1, step = 0.05, value = c(0.5)))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("emax_tumor", label = "EMAX", min = 0.0001, max = 0.01, step = 0.0001, value = c(0.005))),
                        column(3,
                               sliderInput("ec50_tumor", label = "EC50 (ng/mL)", min = 1, max = 200, step = 5, value = c(50))),
                        column(2,
                               sliderInput("lamda0_tumor", label = "Lambda0", min = 0.0001, max = 0.01, step = 0.0001, value = c(0.005)))
                        # column(2,
                        #        sliderInput("w0_tumor", label = "W0 (mm^3)", min = 100, max = 100, step = 10, value = c(100))),
                        # column(2,
                        #        sliderInput("frac_tumor", label = "frac", min = 0.05, max = 0.05, step = 0.05, value = c(0.05)))
                    ),
                    absolutePanel(
                        top = 150, right = 120, width = 80, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>TGI<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "tgi", inline = T), "</code> %</strong>")
                        )
                    ),
                    # absolutePanel(
                    #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                    #     h3("ing")
                    # ),
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            
            
            tabPanel(
                strong("efficacy biomarker"),
                fluidRow(
                    fluidRow(
                        column(4, 
                               # h3("Dose regimen"),
                               fluidRow(
                                   column(6),
                                   column(6,
                                          sliderInput("cycle_biomarker", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("amt_biomarker1", label = "dose (mg)", min = 0, max = 15, step = 0.5, value = c(3))
                                   ),
                                   column(6,
                                          sliderInput("amt_biomarker2", label = "dose interruption (mg)", min = 0, max = 15, step = 0.5, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("interval_biomarker1", label = "dose interval (h)", min = 0, max = 48, step = 2, value = c(24))
                                   ),
                                   column(6,
                                          sliderInput("interval_biomarker2", label = "interrupted interval (h)", min = 0, max = 48, step = 2, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("n_biomarker1", label = "dose times", min = 1, max = 28, step = 1, value = c(28))
                                   ),
                                   column(6,
                                          sliderInput("n_biomarker2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                                   )
                               )
                        ),
                        column(4, plotOutput("distPlotEfficacyPKPDtime")),
                        column(4, plotOutput("distPlotEfficacyPKPD"))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("cl_biomarker", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(5))),
                        column(3,
                               sliderInput("v_biomarker", label = "V (L)", min = 10, max = 500, step = 10, value = c(100))),
                        column(3,
                               sliderInput("ka_biomarker", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.5))),
                        column(3,
                               sliderInput("f_biomarker", label = "F", min = 0.1, max = 1, step = 0.05, value = c(0.5)))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("gamma_biomarker", label = "Gamma", min = 1, max = 5, step = 0.1, value = c(1.2))),
                        column(3,
                               sliderInput("emax_biomarker", label = "EMAX", min = 0.1, max = 1, step = 0.05, value = c(1))),
                        column(3,
                               sliderInput("ec50_biomarker", label = "EC50 (ng/mL)", min = 1, max = 50, step = 1, value = c(15)))
                    ),
                    absolutePanel(
                        top = 200, right = 180, width = 130, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>Inhibition<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "inhi", inline = T), "</code> %</strong>")
                        )
                    ),
                    # absolutePanel(
                    #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                    #     h3("ing")
                    # ),
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            
            
            
            tabPanel(
                strong("safety biomarker"),
                fluidRow(
                    fluidRow(
                        column(4, 
                               # h3("Dose regimen"),
                               fluidRow(
                                   column(6),
                                   column(6,
                                          sliderInput("cycle_anc", label = "cycles", min = 1, max = 10, step = 1, value = c(1))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("amt_anc1", label = "dose (mg)", min = 0, max = 15, step = 0.5, value = c(12))
                                   ),
                                   column(6,
                                          sliderInput("amt_anc2", label = "dose interruption (mg)", min = 0, max = 15, step = 0.5, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("interval_anc1", label = "dose interval (h)", min = 0, max = 48, step = 2, value = c(24))
                                   ),
                                   column(6,
                                          sliderInput("interval_anc2", label = "interrupted interval (h)", min = 0, max = 48, step = 2, value = c(0))
                                   )
                               ),
                               fluidRow(
                                   column(6,
                                          sliderInput("n_anc1", label = "dose times", min = 1, max = 28, step = 1, value = c(28))
                                   ),
                                   column(6,
                                          sliderInput("n_anc2", label = "interruption times", min = 1, max = 28, step = 1, value = c(1))
                                   )
                               )
                        ),
                        
                        column(8, 
                               plotOutput("distPlotSafetyPKPD"))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("cl_anc", label = "CL (L/h)", min = 1, max = 10, step = 0.5, value = c(5))),
                        column(3,
                               sliderInput("v_anc", label = "V (L)", min = 10, max = 500, step = 10, value = c(100))),
                        column(3,
                               sliderInput("ka_anc", label = "Ka (1/h)", min = 0.1, max = 1, step = 0.05, value = c(0.5))),
                        column(3,
                               sliderInput("f_anc", label = "F", min = 0.1, max = 1, step = 0.05, value = c(0.5)))
                    ),
                    fluidRow(
                        column(3,
                               sliderInput("slope_anc", label = "Slope", min = 0.0001, max = 0.01, step = 0.00005, value = c(0.0039))),
                        # column(3,
                        #        sliderInput("circ0_anc", label = "Circ0 (*10^9)", min = 5, max = 5, step = 0.05, value = c(5))),
                        column(3,
                               sliderInput("gamma_anc", label = "Gamma", min = 0.001, max = 1, step = 0.001, value = c(0.122))),
                        column(3,
                               sliderInput("mtt_anc", label = "MTT (h)", min = 50, max = 200, step = 5, value = c(125)))
                    ),
                    absolutePanel(
                        top = 290, right = 30, width = 170, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong>ANC nadir<sub>", "", "</sub> = <code style='color:#ec4c3c;background-color:#F8F9F9'>",
                                   textOutput(outputId = "nadir", inline = T), "</code> *10^9</strong>")
                        )
                    ),
                    
                    absolutePanel(
                        top = 180, right = 10, width = 210, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong> Gr3 AE: <sub>", "</sub>  <code style='color:#ec4c3c;background-color:#F8F9F9'>","0.5~<1.0", 
                                   "</code> *10^9</strong>")
                        )
                    ),
                    absolutePanel(
                        top = 210, right = 0, width = 190, height = 10, draggable = TRUE,
                        HTML(
                            paste0("<strong> Gr4 AE: <sub>", "</sub>  <code style='color:#ec4c3c;background-color:#F8F9F9'>","<0.5", 
                                   "</code> *10^9</strong>")
                        )
                    ),
                    
                    # absolutePanel(
                    #     top = 100, left = 125, width = 100, height = 10, draggable = TRUE,
                    #     h3("ing")
                    # ),
                    absolutePanel(
                        top = 110, left = 55, width = 100, height = 10, draggable = TRUE,
                        img(src="LOGOdMed.png", height = 50)
                    )
                )
            ),
            footer = h5(HTML("dMed Copyright 2019 : 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> E </strong>
                       arly <strong style='color:#ec4c3c;background-color:#F8F9F9'> D </strong> evelepment and 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> C </strong> linical 
                       <strong style='color:#ec4c3c;background-color:#F8F9F9'> P </strong>harmacology"), align = "right")
        )
    )
)
