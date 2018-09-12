#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(ggplot2)
library(ggrepel)
library(reshape2)

#options(encoding="UTF-8")

expr <- read.table("Hinze_normalized_expression.txt",as.is=TRUE)
info <- read.delim("SampleInfo_hinze.txt",encoding = "UTF-8")
res <- read.delim("DifferentialExpressionResults_TCells.txt")

go <- read.delim("go_to_genes.txt",as.is=TRUE)

go$term <- unlist(lapply(strsplit(go$GO.Term," \\("),"[",1))
go$id <- substr(unlist(lapply(strsplit(go$GO.Term," \\("),"[",2)),1,10)

# limit for T cells
expr <- expr[,info$Cells == "T Cells"]
info <- info[info$Cells == "T Cells",]

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Expression for GO Term"),
  
  # Sidebar with a slider input for number of bins 
  
  #http://127.0.0.1:5428/?goTerm=GO:0022408
  textInput("goTerm","GO Term"),
  
  plotOutput("volcano"),
  plotOutput("boxplots"),
  actionButton("prev","Prev"),
  actionButton("next_","Next"),
  dataTableOutput('table')
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  mytitle = "default"

  rValues <- reactiveValues()
  rValues$page <- 0
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    #cat(paste(query))
    
    if (!is.null(query[['goTerm']])) {
      mytitle = query[['goTerm']]
      updateTextInput(session, "goTerm", value = query[['goTerm']])
    }
    
    #cat("title: ",mytitle)
  })
  observeEvent(input$goTerm,{
    if(input$goTerm != ''){
      i <- which (go$id ==input$goTerm)
      genes <- strsplit(go$Genes[i] ,",")[[1]]
      res.genes <- res[match(genes,row.names(res)),]
      res.genes <- na.omit(res.genes)
      res.genes <- res.genes[!is.na(res.genes$padj) & res.genes$padj <= 0.05,] # results are filtered for significance
      
      rValues$numGenes <- nrow(res.genes)
      rValues$sigGenes <- row.names(res.genes)
    }    
  })
  
  observeEvent(input$prev, {
    if (rValues$page > 0)
      rValues$page <- rValues$page -1
  })
  
  observeEvent(input$next_, {
    nGenes<- rValues$numGenes
    cat("number of genes: ",nGenes)
    if ((1+rValues$page)*24 <= nGenes){
      rValues$page <- rValues$page + 1
    }

  })

  output$volcano <- renderPlot(getVolcano(input$goTerm))
  output$boxplots <- renderPlot({
      #cat("next,",input$next_," prev", input$prev)
      page <- rValues$page
      getBoxPlots(input$goTerm, page)
  })
  output$table <- renderDataTable({
    if(input$goTerm != ''){
      sigTable <- res[row.names(res) %in% rValues$sigGenes,]
      
      geneLinks <- paste("<a target ='_new' href='http://www.informatics.jax.org/searchtool/Search.do?query=",
                         row.names(sigTable),
                         "'>",row.names(sigTable),"</a>")
      data.frame(Genes=geneLinks,sigTable)
    }
  },escape = FALSE)
}


getVolcano <- function(goID) {
  if (goID != ""){
    i <- which (go$id ==goID)
    genes <- strsplit(go$Genes[i] ,",")[[1]]
    res.genes <- res[match(genes,row.names(res)),]
    res.genes <- na.omit(res.genes)
    #get the top 100 significant genes
    top <- row.names(res.genes)[order(res.genes$pvalue)]
    top <- top[1:min(100,length(top))]
    top <- row.names(res.genes) %in% top & res.genes$padj < 0.05 # only show top genes if they are sig
    
    ggplot(mapping=aes(res.genes$log2FoldChange,-log10(res.genes$pvalue))) + geom_point(aes(color=res.genes$padj<0.05)) + 
      labs(x="log 2 fold change", y="-log10 p-value", title=paste("Volcano Plot for",go$GO.Term[i])) +
      geom_text_repel(mapping=aes(res.genes$log2FoldChange[top],
                                  -log10(res.genes$pvalue[top]),
                                  label=row.names(res.genes)[top])) + 
      theme_classic() +labs(color="Significant")
    
  }  
}

getBoxPlots <- function(goID, page){
  if (goID != ""){
    cat("page",page,"\n")
    i <- which (go$id ==goID)
    genes <- strsplit(go$Genes[i] ,",")[[1]]
    res.genes <- res[match(genes,row.names(res)),]
    res.genes <- na.omit(res.genes)
    
    res.genes <- res.genes[!is.na(res.genes$padj) & res.genes$padj <= 0.05,] # results are filtered for significance
    genes.expr <- expr[match(row.names(res.genes),row.names(expr)),] # expression values for thos in our filtered results
    
    #genes.expr <- genes.expr[1:min(24,nrow(genes.expr)),] # reduce to the top 20 most signficant
    
    # now sort by mean
    genes.expr <- genes.expr[order(apply(genes.expr,1,mean),decreasing = TRUE),]
    genes.expr$Gene <- row.names(genes.expr)
    facets <- factor(sort(rep(1:6,4))) # creates 6 blocks of 4 genes per facet
    
    start = 1 + 24*page
    if (start > nrow(genes.expr))  start <- nrow(genes.expr)-24
    if (start < 1) start <- 1
    
    max = 2000
    
    
    #cat("\tBoxplots for",nrow(genes.expr),"\n")
    genes.exprBlock <- genes.expr[start:min((start+23),nrow(genes.expr)),] # reduce to the top 20 most signficant
    #cat("\tBoxplots for",nrow(genes.exprBlock),"\n")
    
    melted <- melt(
      droplevels(cbind(genes.exprBlock,facet=facets[1:nrow(genes.exprBlock)])),
      id.vars = c("Gene","facet"))
    melted$Category <- info$Category[match(melted$variable,info$Sample)]

   # plot the expression levels for the specific genes in the pathway

    ggplot(data=melted, aes(Gene,value,fill=Category)) + geom_boxplot() +# scale_y_log10() +
      labs(y="log2 expression",title=go$GO.Term[i])+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      facet_wrap(~facet, scales = "free")
    #start <- start + 24
    
    
  }  
}
# Run the application 
shinyApp(ui = ui, server = server)


