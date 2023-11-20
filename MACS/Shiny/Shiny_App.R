#
#
#     Methylome and transcriptome shiny database
#
#

library(shiny)





myo_data <- read.csv("./myo_data.csv")
myo_int_data <- read.csv("./myo_int_data.csv")
GO_MYO <- read.csv("./GO_MYO.csv")
GO_MYOINT <- read.csv("./GO_MYOINT.csv")
KEGG_MYO <- read.csv("./KEGG_MYO.csv")
KEGG_MYOINT <- read.csv("./KEGG_MYOINT.csv")
gene_lists <- readRDS("./gene_lists.RData")






###################################################################################

# V5 DMPs + KEGG

#################################################################################



# Shiny UI
ui <- fluidPage(
        titlePanel("Epigenetics Database Search"),
        
        tabsetPanel(
                # DMPs Tab
                tabPanel("DMPs", 
                         sidebarLayout(
                                 sidebarPanel(
                                         selectInput("geneSelection", "Select Gene or Gene List:",
                                                     choices = c("type1", "type2a", "Myonuclei", "Individual Gene"),
                                                     selected = "Individual Gene"),
                                         uiOutput("individualGeneInput"),
                                         actionButton("search", "Search"),
                                         checkboxGroupInput("dataset", "Choose Dataset:",
                                                            choices = list("MYO" = "myo", "MYO+INT" = "myo_int"),
                                                            selected = c("myo", "myo_int")),
                                         sliderInput("deltaMRange", "Delta_M Range:",
                                                     min = min(min(myo_data$delta_M), min(myo_int_data$delta_M)),
                                                     max = max(max(myo_data$delta_M), max(myo_int_data$delta_M)),
                                                     value = c(-1, 1), step = 0.01),
                                         checkboxGroupInput("regulatoryFeature", "Choose Regulatory Feature Group:",
                                                            choices = unique(c(myo_data$Regulatory_Feature_Group, myo_int_data$Regulatory_Feature_Group)),
                                                            selected = unique(c(myo_data$Regulatory_Feature_Group, myo_int_data$Regulatory_Feature_Group))),
                                         checkboxGroupInput("relationIsland", "Choose Relation to Island:",
                                                            choices = unique(c(myo_data$Relation_to_Island, myo_int_data$Relation_to_Island)),
                                                            selected = unique(c(myo_data$Relation_to_Island, myo_int_data$Relation_to_Island))
                                         )
                                 ),
                                 mainPanel(
                                         tableOutput("myoResults"),
                                         tableOutput("myoIntResults")
                                 )
                         )
                ),
                
                # KEGG/Go Pathways Tab
                tabPanel("Pathway Analysis",
                         sidebarLayout(
                                 sidebarPanel(
                                         radioButtons("pathwayType", "Select Pathway Type:",
                                                      choices = list("GO" = "go", "KEGG" = "kegg"),
                                                      selected = "go"),
                                         checkboxGroupInput("pathwayDataset", "Choose Dataset:",
                                                            choices = list("MYO" = "myo", "MYO+INT" = "myo_int"),
                                                            selected = c("myo")),
                                         uiOutput("subsetUI"), # Dynamic UI for GO subset
                                         actionButton("analyzePathways", "Analyze")
                                 ),
                                 mainPanel(
                                         tableOutput("myoPathwayResults"),
                                         tableOutput("myoIntPathwayResults")
                                 )
                         )
                )
        )
)

# Shiny Server
server <- function(input, output) {
        
        # The code you've provided earlier for the DMPs goes here...
        output$individualGeneInput <- renderUI({
                if(input$geneSelection == "Individual Gene") {
                        textInput("geneName", "Enter Gene Name:", value = "")
                }
        })
        
        observeEvent(input$search, {
                if(input$geneSelection == "Individual Gene") {
                        gene <- input$geneName
                        gene_filter_myo <- grepl(gene, myo_data$UCSC_RefGene_Name, ignore.case = TRUE)
                        gene_filter_myo_int <- grepl(gene, myo_int_data$UCSC_RefGene_Name, ignore.case = TRUE)
                } else {
                        genes_in_list <- gene_lists[[input$geneSelection]]
                        gene_filter_myo <- myo_data$UCSC_RefGene_Name %in% genes_in_list
                        gene_filter_myo_int <- myo_int_data$UCSC_RefGene_Name %in% genes_in_list
                }
                
                # Filter Data based on input for myo_data
                myo_filtered_data <- myo_data[
                        gene_filter_myo &
                                myo_data$delta_M >= input$deltaMRange[1] &
                                myo_data$delta_M <= input$deltaMRange[2] &
                                myo_data$Regulatory_Feature_Group %in% input$regulatoryFeature &
                                myo_data$Relation_to_Island %in% input$relationIsland,
                ]
                
                # Filter Data based on input for myo_int_data
                myo_int_filtered_data <- myo_int_data[
                        gene_filter_myo_int &
                                myo_int_data$delta_M >= input$deltaMRange[1] &
                                myo_int_data$delta_M <= input$deltaMRange[2] &
                                myo_int_data$Regulatory_Feature_Group %in% input$regulatoryFeature &
                                myo_int_data$Relation_to_Island %in% input$relationIsland,
                ]
                
                # Render Tables based on selected datasets
                if ("myo" %in% input$dataset) {
                        output$myoResults <- renderTable({ myo_filtered_data }, caption = "MYO Results")
                } else {
                        output$myoResults <- renderTable({ NULL })
                }
                
                if ("myo_int" %in% input$dataset) {
                        output$myoIntResults <- renderTable({ myo_int_filtered_data }, caption = "MYO+INT Results")
                } else {
                        output$myoIntResults <- renderTable({ NULL })
                }
        })
        
        # Additional code for KEGG/Go pathways can be added here...
        output$subsetUI <- renderUI({
                if(input$pathwayType == "go") {
                        selectInput("pathwaySubset", "Select GO Subset:",
                                    choices = unique(GO_MYO$subset),
                                    selected = NULL, multiple = TRUE)
                } else {  # Assuming input$pathwayType == "kegg"
                        selectInput("pathwaySubset", "Select KEGG Subset:",
                                    choices = unique(KEGG_MYO$subset),
                                    selected = NULL, multiple = TRUE)
                }
        })
        
        # KEGG/Go Pathways
        observeEvent(input$analyzePathways, {
                data_myos <- NULL
                data_myoint <- NULL
                
                if (input$pathwayType == "go") {
                        if ("myo" %in% input$pathwayDataset) {
                                data_myos <- GO_MYO
                        }
                        if ("myo_int" %in% input$pathwayDataset) {
                                data_myoint <- GO_MYOINT
                        }
                } else {  # Assuming input$pathwayType == "kegg"
                        if ("myo" %in% input$pathwayDataset) {
                                data_myos <- KEGG_MYO
                        }
                        if ("myo_int" %in% input$pathwayDataset) {
                                data_myoint <- KEGG_MYOINT
                        }
                }
                
                # Filter by subset if selected
                if(!is.null(input$pathwaySubset)) {
                        if(!is.null(data_myos)) {
                                data_myos <- data_myos[data_myos$subset %in% input$pathwaySubset, ]
                        }
                        if(!is.null(data_myoint)) {
                                data_myoint <- data_myoint[data_myoint$subset %in% input$pathwaySubset, ]
                        }
                }
                
                # Render Tables
                output$myoPathwayResults <- renderTable(data_myos, caption = "MYO Pathway Results")
                output$myoIntPathwayResults <- renderTable(data_myoint, caption = "MYO+INT Pathway Results")
        })
}

# Run Shiny App
shinyApp(ui = ui, server = server)




# host shiny app in shinyapp.io

install.packages("rsconnect")

library(rsconnect)

rsconnect::setAccountInfo(name='maxullrich',
                          token='CABAE5B2F33FDB2604DF2CAA21D0C92A',
                          secret='v5OtZDFpumsUSSuPag5N0MYrVZShCfq/FKSjEuif')


