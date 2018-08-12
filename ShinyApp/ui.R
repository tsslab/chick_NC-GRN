#ui.R

# load useful packages/libraries ----------------------

library(shiny)
library(shinythemes)
library(d3heatmap)
library(DT)
library(plotly)
library(viridis)
library(knitr)

# User Interface  ----------------------

shinyUI(fluidPage(theme = shinytheme("cosmo"),
                  titlePanel("Chick Neural Crest Gene Regulatory Network, TSS-lab"),
                  sidebarLayout(
                    # Sidebar -----------------------------------------
                    sidebarPanel(position = "left",
                                 radioButtons("main_selection", "Select one:",
                                              choices=c("Survey gene expression"=1,
                                                        "Visualize WGCNA clusters"=2,
                                                        "Search for regulatory elements"=3,
                                                        "Explore the regulatory networks"=4,
                                                        "Explore coexpression genes"=5 
                                                        ))
                    ),
                    # Main  -----------------------------------------
                    mainPanel(
                      tabsetPanel(
                        tabPanel("Data", value=1, 
                                 # RNA-seq  -----------------------------------------
                                 
                                 conditionalPanel(condition="input.main_selection==1", 
                                                  h3("Select a Gene"),
                                                  selectizeInput('Gene_Name', 
                                                                 label = "",
                                                                 choices = NULL),
                                                  
                                                  h3("Gene Information"),
                                                  textOutput("SelectedGene"),
                                                  textOutput("Ensembl_Geneid"),
                                                  textOutput("Full_Gene_Name"),
                                                  br(),
                                                  h3("Expression levels"),
                                                  helpText("Expression levels measured in FPKMs (Fragments Per Kilobase of 
                                                           transcript per Million reads)"),
                                                  plotOutput("plot_fpkm", width = "80%", height = "400px"),
                                                  #verbatimTextOutput("event"),
                                                  br(),
                                                  h3("Differential Expression"),
                                                  helpText("Differential expression was performed using R package DESeq2. Genes were 
                                                            deemed significantly differentially expressed if the adjusted p-value was less than 0.05."),
                                                  tableOutput("table_de")
                                 ),
                                 conditionalPanel(condition="input.main_selection==2", 
                                                  helpText("Select a gene or a cluster to view its expression:"),
                                                  selectizeInput('WGCNA_Choice', 
                                                                 label = NULL,
                                                                 choices = NULL),
                                                  htmlOutput("text_wgcna"),
                                                  #                                                   #uiOutput("wgcnaimage"),
                                                  imageOutput("wgcnaimage",height = "600px"),
                                                  h5("Genes in cluster:"),
                                                  DT::dataTableOutput("WGCNAtable", width="70%")
                                 ),
                                 # Regulatory elements -----------------------------------------
                                 conditionalPanel(condition="input.main_selection==3", 
                                                  h3("Select a Gene"),
                                                  selectizeInput('Gene_Name_reg', 
                                                                 label = "",
                                                                 choices = NULL),
                                                  h3("Associated putative regulatory elements identified in Clusters-3-4-9"),
                                                  tableOutput('table_enh_clust'),
                                                  h3("Transcription factor binding motifs found in Cluster 3 and 9 elements"),
                                                  helpText("Motif positions were identified by screening TF PWMs from the HOCOMOCO v.4 human database against enhancers using the annotatePeaks.pl script from Homer 4.7."),
                                                  d3heatmapOutput("TFBS",  width = "100%", height = "1200px"),     
                                                  h3("Associated putative regulatory elements identified using DiffBind"),
                                                  helpText("ATAC-seq peaks were found differentially bound using DiffBind."),
                                                  tableOutput('table_diffbind'),
                                                  h3("Transcription factor binding motifs found in DiffBind elements"),
                                                  d3heatmapOutput("TFBS_diffbind",  width = "100%", height = "1200px")
                                 ),
                                 # Regulatory networks -----------------------------------------
                                 conditionalPanel(condition="input.main_selection==4", 
                                                  h3("Regulatory networks"),
                                                  selectizeInput('Gene_Name_RN', 
                                                                 label = "",
                                                                 choices = NULL),
                                                  #helpText("More details here"),
                                                  h4("Cluster 3"),
                                                  plotOutput("plot_rn_3", width = "90%", height = "800px"),
                                                  h4("Cluster 9"),
                                                  plotOutput("plot_rn_9", width = "90%", height = "800px")
                                 ),
                                 # Single cell coexpression -----------------------------------------
                                 conditionalPanel(condition="input.main_selection==5", 
                                                  helpText("Select a gene to view the coexpression heatmap or multiple genes to see their coexpression:"),
                                                  selectizeInput('Cor_Gene',
                                                                 label = "",
                                                                 choices = NULL, multiple = TRUE),
                                                  conditionalPanel(condition="output.corgenelen",
                                                    helpText("Select a value for the minimum of the correlation:"),
                                                    numericInput('cor_thres', label=NULL, value=0.7,
                                                               min=0.05, max=1, step=0.05),
                                                    textOutput("text_cor")
                                                  ),
                                                  d3heatmapOutput("heatmap_coexpr",  width = "80%", height = "800px")
                                 )
                        ),
                        # About us -----------------------------------------
                        tabPanel("About us", value=2,
                                 includeMarkdown(knitr::knit("Chicken_paper.Rmd"))),
                        id = "tabselected"))
)
)
)
