library(shiny)

shinyUI(fluidPage(
  titlePanel("Microbiome plots"),   
      sidebarLayout(
        sidebarPanel(
          radioButtons("select_data", "Use example or import data", choices = c("Example","Import data"), selected = character(0)),
          conditionalPanel(
            condition="input.select_data=='Import data'",
            fileInput("file.bion", "Import bion .tsv File",
                multiple = F,
                accept = c(".tsv","text/csv", "text/comma-separated-values,text/plain"))),
          uiOutput("choose_level"),
          uiOutput("select_allsamples"),   
          uiOutput("select_sample"),
          conditionalPanel(
            condition="input.plot_type!='NMDS'",
            uiOutput("select_top10"),  
            uiOutput("search")),
          selectInput('plot_type', 'Choose plot', c("geom_bar", "geom_tile", "geom_point", "NMDS")),
          conditionalPanel(
            condition="input.plot_type=='geom_tile'",
            checkboxInput("values_in_tiles", "Add values to heatmap", FALSE)),
          checkboxInput("mod_fig", "Figure options", FALSE),
          conditionalPanel(
            condition="input.mod_fig==true",
            selectInput(inputId = "x_label_size",label = "x-axis label size:", choices = c(1:16), selected = 10)),
          conditionalPanel(
            condition="input.mod_fig==true",
            selectInput(inputId = "y_label_size",label = "y-axis label size:", choices = c(1:16), selected = 10)),
          conditionalPanel(
            condition="input.mod_fig==true & input.values_in_tiles==true",
            selectInput(inputId = "value_size",label = "Tile value size:", choices = c(1:5), selected = 3)),
          conditionalPanel(
            condition="input.mod_fig==true",
            selectInput(inputId = "fig_height",label = "Figure height (pixels)", choices = c(500,1000), selected = 500)),
          conditionalPanel(
            condition="input.mod_fig==true",
            selectInput(inputId = "fig_width",label = "Figure width (pixels)", choices = c(500,1000), selected = 500)),
          conditionalPanel(
            condition="input.mod_fig==true",
            selectInput(inputId = "fig_format",label = "Figure format:", choices = c("png","pdf"), selected = "png")),
          downloadButton("downloadFig", "Download")
        ),
        mainPanel(
          plotOutput(outputId = "ggplot", width = 500)
          #DT::dataTableOutput("React_Out")
        )
      )
)
)


