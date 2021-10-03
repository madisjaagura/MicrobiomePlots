
library(shiny)
library(ggplot2)
source(file.path("source","source.R"), local = T)


shinyServer(function(input, output, session){


  reactive_level.select = reactive({
    req(input$c_level)
    input$c_level
    })

  reactive_data = reactive({
    req(reactive_level.select())
    if(input$select_data=="Example" | is.null(input$file.bion$datapath)){
      get.data(level = reactive_level.select(), file.bion = "raw_data/BION.taxmap.species.tsv")} else {
      get.data(level = reactive_level.select(), file.bion = ifelse(input$select_data=="Import data" & !is.null(input$file.bion$datapath), input$file.bion$datapath,"raw_data/BION.taxmap.species.tsv"))}
    })


  reactive_data.long = reactive({
    reactive_data() %>% gather(-tax, key = "sample", value = "proportion") %>% mutate(sample = factor(sample, levels = unique(sample))) #alternatiiv: species.long = melt(species) #pakett - reshape
  })  

  reactive_data.long.select = reactive({
    x = reactive_data.long()
    x <- x %>% filter(tax %in% input$search)
    x <- x %>% filter(sample %in% input$s_sample)
    x
  })  

  reactive_top10data = reactive({
    reactive_data.long() %>% filter(sample %in% input$s_sample) %>%
      group_by(tax) %>%
      mutate(avg.tax = mean(proportion)) %>%
      arrange(desc(avg.tax)) %>% 
      group_by(sample) %>%
      slice(1:10) #%>% mutate(tax= gsub("g__|s__","",tax))
  })  

  observe({
    req(input$select_data) #UI plokk joonistatakse vaid andmeallika selekteerimise korral
    output$choose_level <- renderUI({
      selectizeInput("c_level", "Select taxonomy level:", multiple=F, choices = c("species", "genus", "phylum"))
    })  
  })

  observe({
    req(reactive_data()) #UI plokk joonistatakse vaid andmete olemasolu korral
    
    output$select_top10 <- renderUI({
      checkboxInput("s_top10", "Select TOP10 taxa", FALSE)
    })  

    output$select_allsamples <- renderUI({
      checkboxInput("s_allsamples", "Select all samples", FALSE)
    })  

    output$select_sample <- renderUI({
      selectizeInput("s_sample", "Select samples:", multiple=TRUE, choices = colnames(reactive_data())[-1])
    })  


    output$search <- renderUI({
      selectizeInput("search", "Search taxa: ", multiple=TRUE, choices = unique(reactive_data.long()$tax))
    })
  })

  observe({
    #req(reactive_data.long.select())
    output$React_Out = DT::renderDataTable({
      DT::datatable(reactive_data.long.select(), filter = "top")
    })  
  })


  observeEvent(input$s_allsamples,{
    
    if (input$s_allsamples == T){
      x <- colnames(reactive_data())[-1]
    updateSelectizeInput(session, "s_sample", selected = x, options = list())
    } else {updateSelectizeInput(session, "s_sample", selected = "", options = list())}
  })


  observeEvent(input$s_top10,{
    
    if (input$s_top10 == T){
      x <- unique(reactive_top10data()$tax)
      updateSelectizeInput(session, "search", selected = x, options = list())
    } else {updateSelectizeInput(session, "search", selected = "", options = list())}
  })


  output$out <- renderPrint(paste0(
    length(input$select)
  ))

  observe({
    output$ggplot <- renderPlot({
      input.fig = reactive_data.long.select()
      input.fig$tax = gsub("g__|s__","",input.fig$tax)

         a.ratio = length(unique(input.fig$tax)) / (length(unique(input.fig$sample)))
         p <- ggplot(input.fig, aes(x = tax, y = sample)) + 
             theme_classic() + 
             theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5, size=input$x_label_size), 
               axis.text.y=element_text(size=input$y_label_size), 
               aspect.ratio = a.ratio) + 
             xlab("") + ylab("") 
        
        if(input$plot_type == "geom_tile"){
          p <- p + geom_tile(aes(fill = proportion), colour = "white") + coord_flip() + scale_fill_gradientn(colours = c("white", "steelblue"), limits = c(0.000001,max(input.fig$proportion)), name = "Proportion", trans = "sqrt", na.value="lightgray")
          }     
        if(input$plot_type == "geom_point"){ 
          p <- p + geom_point(aes(size = proportion)) + coord_flip()
          }
        if(input$values_in_tiles){
          p <- p + geom_text(aes(label=signif(proportion, digits = 2)), size = as.numeric(input$value_size))
          }
        if(input$plot_type == "geom_bar"){ 
          library(RColorBrewer)
          cols <- colorRampPalette(brewer.pal(12, "Set3"))
          myPal <- cols(length(unique(input.fig$tax)))
          p <- p + geom_bar(aes(x = sample, y = proportion, fill = tax), stat='identity') + ylab("Proportion") + scale_fill_manual(values = myPal) 
        }




        if(input$plot_type == "NMDS"){ 
          set.seed(1234)
          data = input.fig %>% spread(sample, proportion)
          meta.nmds.data <- metaMDS(t(data[,-1])) 
          NMDS.data =  as.data.frame(t(data[0,]))
          NMDS.data = data.frame(NMDS1 = meta.nmds.data$points[,1], NMDS2 = meta.nmds.data$points[,2],sample= rownames(meta.nmds.data$points))
         #NMDS.data$prevotella = colSums(data[grepl("Prevotella", data$tax),-1])
         #NMDS.data$bacteroides = colSums(data[grepl("Bacteroides", data$tax),-1])

          library(ggrepel)

          p <- ggplot(NMDS.data, aes(x = NMDS1, y = NMDS2, label = sample))  + 
                theme_classic() + 
                geom_point(alpha = 0.5,size = 5) + geom_text_repel() + 
                theme(axis.text.x=element_text(size=input$x_label_size), 
                  axis.text.y=element_text(size=input$y_label_size)) + 
                ggtitle("") + coord_fixed() #+background_grid(major = "xy", minor = "none")  
        }         
      print(p)  
    })
  })
})      

  #     }
  #     if(input$plot_type == "geom_bar_TOP10" | input$plot_type == "geom_bar_TOP10_sample" ){ 
  #       if(length(input$select) > 0){ species.long = species.long[species.long$sample %in% c(input$select,"tax"),]}
  #         species.long$sample = factor(species.long$sample, levels = unique(species.long$sample))
  #         input.fig = species.long[gsub("g__| s__.*$|_..*","",species.long$tax) %in% input$genus,]
  #     }   
  #     if(input$plot_type == "geom_bar_TOP10"){ 
  #       input.fig <- species.long %>% 
  #         group_by(tax) %>%
  #         mutate(avg.tax = mean(proportion)) %>%
  #         arrange(desc(avg.tax)) %>% 
  #         group_by(sample) %>%
  #         slice(1:10) %>%
  #         mutate(tax= gsub("g__|s__","",tax))
        
  #       cols <- colorRampPalette(brewer.pal(12, "Set3"))
  #       myPal <- cols(length(unique(input.fig$tax)))
  #       a.ratio = length(unique(input.fig$tax)) / (length(unique(input.fig$sample)))
  #       p <- ggplot(input.fig, aes(x = tax, y = sample)) +  
  #            theme_classic() + 
  #            geom_bar(aes(x = sample, y = proportion, fill = tax), stat='identity') + 
  #            scale_fill_manual(values = myPal) + 
  #            theme(axis.text.x=element_text(angle=90, hjust=1, vjust = 0.5, size=input$x_label_size), 
  #             axis.text.y=element_text(size=input$y_label_size), 
  #             aspect.ratio = a.ratio)  + ylab("Proportion") + xlab("") #+ theme(aspect.ratio=1/3) 
  #       print(p)
  #     }
  #   }, height = as.numeric(input$fig_height), width = as.numeric(input$fig_width)) #

  

    # output$downloadFig <- downloadHandler(
    #   filename = function() {
    #     paste0(paste0(input$genus, collapse = "_"),".", input$plot_type, ifelse(input$fig_format == "png",".png",".pdf"))
    #   },
    #   content = function(file) {
    #     ggsave(file, width = as.numeric(input$fig_width)/72, height = as.numeric(input$fig_height)/72, units = "in")
    #   }
    # )

