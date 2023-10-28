#
#
#     Methylome and transcriptome shiny database
#
#

library(shiny)
library(ggplot2)

# Sample data for illustration purposes. Replace with your actual datasets.
myo_data <- DMPs_PM_vs_BM %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1))


myo_int_data <- DMPs_PH_vs_BH %>% 
        filter(p.value < 0.05) %>% 
        merge(., anno, by = "cpg") %>% 
        mutate(UCSC_RefGene_Name = sapply(strsplit(UCSC_RefGene_Name, split = ";"), `[`, 1))




###################################################################################

# V5

#################################################################################


type1 <- c("ANKRD2",
           "ATP2A2",
           "CA3",
           "CASQ2",
           "CD36",
           "CYB5R1",
           "FABP3",
           "LDHB",
           "MYH7",
           "MYL12A",
           "MYL2",
           "MYL3",
           "MYL6B",
           "MYOZ2",
           "PDLIM1",
           "PLN",
           "TNNC1",
           "TNNI1",
           "TNNT1",
           "TPM3"
)

type2a <- c("ALDOA",
            "ATP2A1",
            "DDIT4L",
            "ENO3",
            "G0S2",
            "GAPDH",
            "LDHA",
            "MYBPC2",
            "MYH1",
            "MYH2",
            "MYL1",
            "MYLPF",
            "PFKM",
            "PGM1",
            "PKM",
            "SLN",
            "TNNC2",
            "TNNI2",
            "TNNT3",
            "TPM1"
)

Myonuclei = Myonuclei %>% pull(symbol)

gene_lists <- list(
        type1 = type1,
        type2a = type2a,
        Myonuclei = Myonuclei
)



# Shiny UI
ui <- fluidPage(
        titlePanel("Epigenetics Database Search"),
        
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
                                           selected = unique(c(myo_data$Relation_to_Island, myo_int_data$Relation_to_Island)))
                ),
                
                mainPanel(
                        tableOutput("myoResults"),
                        tableOutput("myoIntResults")
                )
        )
)

# Shiny Server
server <- function(input, output) {
        
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
}
# Run Shiny App
shinyApp(ui = ui, server = server)




###################################################################################

# V5 DMPs + KEGG

#################################################################################


GO_MYOINT <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_homogenate.csv")
GO_MYO <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/GO_all_pathways_myonuclei.csv")
KEGG_MYOINT <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_all_pathways_homogenate.csv")
KEGG_MYO <- read.csv("/Users/maxul/Documents/Skole/Master 21-22/Master/DATA/Supplementary files/kegg_all_pathways_myonuclei.csv")

GO_MYOINT <- GO_MYOINT %>% filter(P.DE < 0.05) %>% dplyr::select(2:9)
GO_MYO <- GO_MYO %>% filter(P.DE < 0.05) %>% dplyr::select(2:9)
KEGG_MYO <- KEGG_MYO %>% filter(P.DE < 0.05) %>% dplyr::select(2:8)
KEGG_MYOINT <- KEGG_MYOINT %>% filter(P.DE < 0.05) %>% dplyr::select(2:8)

# add subset to KEGG pathways

names(KEGG_new$kg.sets) %>% as.data.frame() %>% 
        dplyr::select("pathway" = 1) -> subset

subset[KEGG_new$sig.idx,"subset"] <- "Signalling"
subset[KEGG_new$met.idx,"subset"] <- "Metabolic"
subset[KEGG_new$dise.idx,"subset"] <- "Disease"

KEGG_MYO <- KEGG_MYO %>% merge(., subset, by = "pathway")
KEGG_MYOINT <- KEGG_MYOINT %>% merge(., subset, by = "pathway")


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
