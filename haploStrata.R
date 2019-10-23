library(shiny)
library(ggplot2)
library(pheatmap)


d = read.delim("data/plotWithDecade.txt",
               header = T,
               sep = "\t")
d$decade = factor(d$decade, levels=c('early','middle','late'))
mid <- d$start + ((d$stop - d$start) / 2)
wide <- d$stop - d$start
d = cbind(d, mid, wide)

# varieties = read.delim("data/plotWithDecade.txt",
#                        header = T,
#                        sep = "\t")
# 
haplos = read.delim("data/allelesWithOrder.txt", sep = "\t")


# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("HaploStrata: US Rice"),
  tags$hr(),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      numericInput("chrom", "chromosome to view (uses numeric values)", value = 1),
      #numericInput("winSize", "window size", value = 200),
      # Horizontal line ----
      tags$hr(),
      
      # # Input: Select a file ----
      # fileInput("annotations", "Load your bed file here",
      #           multiple = FALSE,
      #           accept = c("text/csv",
      #                      "text/comma-separated-values,text/plain",
      #                      ".csv")),
      uiOutput("tab2")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      h4("Select a window within any panel and zoom by double-clicking it (all panels will adjust accordingly) | Double click outside box (or with no box) to return to full view | Single click to show individual LD block details below"),
      # Output: Data file ----
      plotOutput(outputId = "plot01", width = "1200px", height = "600px",dblclick = "plot01_dblclick",click="plot01_click",brush = brushOpts(id = "plot01_brush",resetOnNew = TRUE)),
      h5("Heatmap showing similiarity of haplotypes.  The column labels match the haplotypes above (see key at far right)"),
      plotOutput(outputId = "plot02", width = "1200px", height = "300px",dblclick = "plot02_dblclick",click="plot02_click",brush = brushOpts(id = "plot02_brush",resetOnNew = TRUE)),
      tags$head(tags$style(HTML("#info{font-size: 10px;}"))),
      downloadButton("button2",label="export the SNP matrix (in hapmap format) for the selected LD block (shown in heatmap and defined below)",class="butt2"), tags$head(tags$style(".butt2{background-color:black;} .butt2{color: white;} .butt2{font-style: italic;} .butt2{font-size:100%}")),
      downloadButton("button3",label="export sample information and haplotype calls for the selected LD block (shown in heatmap and defined below)",class="butt2"), tags$head(tags$style(".butt2{background-color:black;} .butt2{color: white;} .butt2{font-style: italic;} .butt2{font-size:100%}")),
      #plotOutput(outputId = "plot02", width = "1200px", height = "600px",dblclick = "plot02_dblclick",click="plot02_click",brush = brushOpts(id = "plot02_brush",resetOnNew = TRUE)),tags$head(tags$style(HTML("#info{font-size: 10px;}")))
      verbatimTextOutput(outputId ="contents")
    )
  )
)



# Define server logic to read selected file ----
server <- function(input, output) {
  options(shiny.maxRequestSize=50*1024^2) 
  
  # url1 <- a("QTL-seq file format description and example (test.freq)", href="https://github.com/USDA-ARS-GBRU/QTLsurge/tree/master/test")
  # url2 <- a("Cycle file format description and example (cycle1.freq and cycle2.freq)", href="https://github.com/USDA-ARS-GBRU/QTLsurge/tree/master/test")
  # output$tab <- renderUI({
  #   tagList(url1)
  # })
  # output$tab2 <- renderUI({
  #   tagList(url2)
  # })

  
  dfChrom <- reactive({
    #levels = levels(d$chr)
    d[which(d$chrom == input$chrom),]  
  })
  
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot01 <- renderPlot({
    centr = c(16.7,13.6,19.4,9.7,12.4,15.3,12.1,12.9,2.8,8.2,12,11.9)
    ggplot(data = dfChrom(),aes(y = count, x = mid, width = wide, fill = as.character(order))) + 
      geom_bar(stat="identity",colour="white") + 
      labs(x="position on genome (bp)") +
      coord_cartesian(xlim = ranges$x, expand = FALSE) + #ylim = ranges$y
      facet_grid(decade~.)
  })
  
  
  output$plot02 <- renderPlot({
    validate(
      need(input$plot01_click != "","To see something here, click on an LD block of interest.")
    )
    d <- dfChrom()
    d <- d[which(focalPoint$pos < d$stop & focalPoint$pos >d$start),]
    id = d$id[1] #all ids should be same
    temp = haplos[which(haplos$Locus == as.character(id)),]
    temp = na.omit(temp)
    mat = adist(temp$Allele)
    #par(mar=c(0,0,0,0))
    l = length(temp$Allele[1])
    heatmap(mat, col = paste("gray",1:99,sep=""),symm=T,labRow = temp$Allele, labCol=temp$order, cexRow = 1/l)
  })
  
  
  ################################
  observeEvent(input$plot01_dblclick, {
    brush <- input$plot01_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  focalPoint <- reactiveValues(pos = NULL, order = NULL)
  
  observeEvent(input$plot01_click, {
    brush <- input$plot01_click
    if (!is.null(brush)) {
      focalPoint$pos <- brush$x
      focalPoint$order <- brush$y
    } else {
      focalPoint$pos <- NULL
      focalPoint$order <- NULL
    }
  })
  
  output$button2 <- downloadHandler(
    filename = function() {
      d <- dfChrom()
      d <- d[which(focalPoint$pos < d$stop & focalPoint$pos >d$start),]
      id = d$id[1] #all ids should be same
      paste('snpsInHaplotype_', id, '.txt', sep='')
    },
    content = function(con) {
      d <- dfChrom()
      d <- d[which(focalPoint$pos < d$stop & focalPoint$pos >d$start),]
      id = d$id[1] #all ids should be same
      snps = read.delim("data/common.nohet.thin3k.hmp.txt")
      snps = snps[which(snps$chrom == d$chr[1]),]
      snps = snps[which(snps$pos < d$stop[1] & snps$pos > d$start[1]),]
      write.table(snps,con,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
  )
  
  output$button3 <- downloadHandler(
    filename = function() {
      d <- dfChrom()
      d <- d[which(focalPoint$pos < d$stop & focalPoint$pos >d$start),]
      id = d$id[1] #all ids should be same
      paste('haplotypeAndVarInfo_', id, '.txt', sep='')
    },
    content = function(con) {
      d <- dfChrom()
      d <- d[which(focalPoint$pos < d$stop & focalPoint$pos >d$start),]
      id = d$id[1] #all ids should be same
      var = read.delim("data/varietyLookup.txt")
      selectedLDBlk_order1st_code2nd = var[,id]
      new = cbind(var[,1:7], selectedLDBlk_order1st_code2nd)
      write.table(new,con,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)
    }
  )
  
  
  output$contents <- renderPrint({
    d <- dfChrom()
    d <- d[which(focalPoint$pos < d$stop & focalPoint$pos >d$start),]
    d
    #focalPoint$order
    #head(brushedPoints(dfChrom(), input$plot01_brush, allRows = TRUE))
  })
  
  # output$contents <- renderText({
  #   if(is.null(input$plot01_brush)) return("NULL\n")
  #  temp = input$plot01_brush
  #    #cat(str(temp))
  #    {hoverValue(temp)}
  #   #paste0(str(ranges))
  # })
  
  
  
}

# Create Shiny app ----
shinyApp(ui, server)