#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(XML)
library("methods")
library(DT)
library(stringr)
library(reticulate) # to source python code
use_virtualenv("env")
library(readr)
library(shinythemes)
library(crosstalk)
library(leaflet)
library(patchwork)
library(scales)
library(tidyverse)
library(shinypanels)
library(readxl)

############################################################################################################
# Read foundation xml function

foundationXML_general <- function(inputFile){
    
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    # Extract General Information
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    SubmittedDiagnosis <- gsub('^.*<SubmittedDiagnosis>\\s*|\\s*</SubmittedDiagnosis>.*$', '', xmlText)
    Gender <- gsub('^.*<Gender>\\s*|\\s*</Gender>.*$', '', xmlText)
    DOB <- gsub('^.*<DOB>\\s*|\\s*</DOB>.*$', '', xmlText)
    
    if (str_length(DOB)> 11) {
        DOB <- gsub('^.*<DOB xmlns:xsd="http://www.w3.org/2001/XMLSchema">\\s*|\\s*</DOB>.*$', '', xmlText)
    }
    
    SpecSite <- gsub('^.*<SpecSite>\\s*|\\s*</SpecSite>.*$', '', xmlText)
    if (str_length(SpecSite) > 50) {
        SpecSite <- gsub('^.*<SpecFormat>\\s*|\\s*</SpecFormat>.*$', '', xmlText)
    }
    
    CollDate <- gsub('^.*<CollDate>\\s*|\\s*</CollDate>.*$', '', xmlText)
    if (str_length(CollDate) > 11) {
        CollDate <- gsub('^.*<CollDate xmlns:xsd="http://www.w3.org/2001/XMLSchema">\\s*|\\s*</CollDate>.*$', '', xmlText)
    }
    
    ReceivedDate <- gsub('^.*<ReceivedDate>\\s*|\\s*</ReceivedDate>.*$', '', xmlText)
    if (str_length(ReceivedDate) > 11) {
        ReceivedDate <- gsub('^.*<ReceivedDate xmlns:xsd="http://www.w3.org/2001/XMLSchema">\\s*|\\s*</ReceivedDate>.*$', '', xmlText)
    }
    CountryOfOrigin <- gsub('^.*<CountryOfOrigin>\\s*|\\s*</CountryOfOrigin>.*$', '', xmlText)
    if (str_length(CountryOfOrigin) > 50) {
        CountryOfOrigin <- c("NA")
    }
    # Purity estimates
    Pathological_Purity <- as.numeric(gsub('^.*percent-tumor-nuclei="\\s*|\\s*".*$', '', xmlText))
    Computational_Purity <- as.numeric(gsub('^.*purity-assessment="\\s*|\\s*".*$', '', xmlText))
    # Biomarkers
    TMB_Score <- as.numeric(gsub('^.*tumor-mutation-burden score="\\s*|\\s*".*$', '', xmlText))
    MSS <- gsub('^.*microsatellite-instability status="\\s*|\\s*".*$', '', xmlText)
    if (str_length(MSS)>9) {
        MSS = ""
    }
    
    # Short Variants
    ShortVariants <- gsub('^.*<short-variants>\\s*|\\s*</short-variants>.*$', '', xmlText)
    ShortVariants <- str_split(ShortVariants, " </short-variant> ")
    ShortVariants<- ShortVariants[[1]][] # all short variants with attributes
    Genes <- c()
    Allele_Freqs <- c()
    Depths <- c()
    Protein_Change <- c()
    Positions <- c()
    Functional_Effect <- c()
    Mutations <- c()
    Strand <- c()
    for (i in 1: length(ShortVariants)) {
        line <- ShortVariants[i]
        geneName <- gsub('^.*gene=\"\\s*|\\s*\".*$', '', line)
        Genes[i] <- geneName
        alleleFreq <- gsub('^.*allele-fraction="\\s*|\\s*".*$', '', line)
        Allele_Freqs[i] <- as.numeric(alleleFreq)
        depth <- gsub('^.*depth="\\s*|\\s*".*$', '', line)
        Depths[i] <- as.numeric(depth)
        proteinchange <- gsub('^.*protein-effect="\\s*|\\s*".*$', '', line)
        Protein_Change[i] <- proteinchange
        position <- gsub('^.*position="\\s*|\\s*".*$', '', line)
        Positions[i] <- position
        funct_effect <- gsub('^.*functional-effect="\\s*|\\s*".*$', '', line)
        Functional_Effect[i] <- funct_effect
        mutation <- gsub('^.*cds-effect="\\s*|\\s*".*$', '', line)
        Mutations[i] <- mutation
        strand <- gsub('^.*strand="\\s*|\\s*".*$', '', line)
        Strand[i] <- strand
    }
    
    # find short variants that overlap with copy-number-alterations
    # Copy Number Alterations
    CopyNumberAlterations <- gsub('^.*<copy-number-alterations>\\s*|\\s*</copy-number-alterations>.*$', '', xmlText)
    CopyNumberAlterations <- str_split(CopyNumberAlterations, " </copy-number-alteration> ")
    CopyNumberAlterations<- CopyNumberAlterations[[1]][] # all CNAs with attributes
    CNA_geneName <- c()
    CNA_copyNumber <- c()
    CNA_position <- c()
    for (i in 1 : length(CopyNumberAlterations)) {
        line <- CopyNumberAlterations[i]
        cna_geneName <- gsub('^.*gene=\"\\s*|\\s*\".*$', '', line)
        
        cna_copynumber <- gsub('^.*copy-number=\"\\s*|\\s*".*$', '', line)
        
        cna_pos <- gsub('^.*position=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(cna_copynumber)>5) {
            cna_geneName = ""
            cna_copynumber = ""
            cna_pos = ""
        }
        CNA_geneName[i] <- cna_geneName
        CNA_copyNumber[i] <- cna_copynumber
        CNA_position[i] <- cna_pos
    }
    
    
    var_genes <- data.frame(Genes,Positions)
    cna_gene<- data.frame(CNA_geneName,CNA_copyNumber,CNA_position)
    Ploidy <- rep(2,length(var_genes$Genes)) # initialize ploidy
    
    common <- intersect(var_genes$Genes, cna_gene$CNA_geneName) # find common Vars & CNAs
    position_var_gene <- var_genes[Genes==common,]
    pos_cna_gene <- cna_gene[CNA_geneName==common,]
    Genes_CN = data.frame(Genes,Ploidy)
    
    if (dim(position_var_gene[1]) > 0) {
        # if common variants exist
        Genes_CN = data.frame(Genes,Ploidy)
        Genes_CN[Genes_CN$Genes == common, 2]
        cn <- cna_gene[CNA_geneName == common, 2]
        Genes_CN[Genes_CN$Genes == common, 2] = as.numeric(as.character(cn))
    }
    
    
    
    # create temporary txt file with Gene AlleleFreq Depth and Ploidy for R_allfit.py
    # headers SNV	Allele_freq	Depth	Ploidy PathologicalPurity ComputationalPurity
    PathologicalPurity = rep(Pathological_Purity/100,length(Genes))
    ComputationalPurity = rep(Computational_Purity/100,length(Genes))
    temp_df <- data.frame(Genes_CN$Genes,100*Allele_Freqs,Depths,Genes_CN$Ploidy,PathologicalPurity,
                          ComputationalPurity) # allele freq in %
    colnames(temp_df) <- c("ID","Allele_freq","Depth","Ploidy","PathologicalPurity","ComputationalPurity")
    write.table(temp_df, file = "tabledata.txt", sep = "\t", dec = ".",
                row.names = FALSE, col.names = TRUE)
    
    pur <- py_run_file("All-FIT2.py")
    
    # access output file
    mystring <- read_file("tabledata_out.txt")
    
    allfit_pur <- as.numeric(gsub('^.*tabledata_out\t\\s*|\\s*\t.*$', '', mystring))
    allfit_pur_CI <- gsub('^.*tabledata_out\t*\t\\s*|\\s*".*$', '', mystring)
    allfit_pur_CI <- gsub('^.*\t\\s*|\\s*\n.*$', '', allfit_pur_CI)
    allfit_pur <- rep(100*allfit_pur,length(Genes))
    
    # Define Dataframe: 
    df_GeneralInfo <- data.frame(rep(ReportID,length(ShortVariants)),rep(SubmittedDiagnosis,length(ShortVariants)),
                                 rep(Gender, length(ShortVariants)),rep(DOB,length(ShortVariants)),
                                 rep(SpecSite,length(ShortVariants)),rep(CollDate,length(ShortVariants)),
                                 rep(ReceivedDate,length(ShortVariants)),rep(CountryOfOrigin,length(ShortVariants)),
                                 rep(Pathological_Purity,length(ShortVariants)),rep(Computational_Purity,length(ShortVariants)),
                                 allfit_pur,rep(TMB_Score,length(ShortVariants)),rep(MSS,length(ShortVariants)))
    colnames(df_GeneralInfo) <- c("Report_ID","Diagnosis","Gender","DOB","Specimen_Site","Collection_Date",
                                  "Recieved_Date","Country_Of_Origin","Pathological_Purity",
                                  "Computational_Purity","AllFIT_Purity","TMB_Score","MS_Status")
    
    df_ShortVariants <- data.frame(Genes_CN$Genes,Protein_Change,Allele_Freqs,Genes_CN$Ploidy,Depths,Positions,Strand,
                                   Mutations,Functional_Effect)
    colnames(df_ShortVariants) <- c("Gene","Protein_Change","Allele_Freq","Ploidy","Depth","Position",
                                    "Strand","Mutation","Functional_Effect")
    
    df <- cbind(df_GeneralInfo,df_ShortVariants)
    df <- df[order(Allele_Freqs),] 
    
    # extract model predictions
    pathModels <- read.delim("pathologicalPurity_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
    compModdels <- read.delim("computationalPurity_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
    allfitModels <- read.delim("AllFIT_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
    
    df_models <- data.frame(pathModels$Model,compModdels$Model,allfitModels$Model)
    colnames(df_models) <- c("Pathological Purity Model","Computational Purity Model","All-FIT Purity Model")
    
    df <- cbind(df,df_models)
    
    contents <- df
    
    return(contents)
}
## -------------------------------------------------------------------------------------------------------
# Read general file format  function
generalFileFormat <- function(inputFile){
    
    fileData <- read_excel(inputFile)
    # check input headers
    cn = colnames(fileData)
    required_headers = c("Sample_ID","Gene","VAF","Depth","Copy_Number")
    chk = is.element(required_headers,cn)
    if (isTRUE(all(chk))) {
        print("all headers are present")
    } else {
        stop("Required Fields Missing; Required Field headers: Sample_ID, Gene, VAF, Depth, Copy_Number")
    }
    
    # check if pathalogical and/or computational purities exist
    chk1 = is.element("Pathological_Purity", cn)
    chk2 = is.element("Computational_Purity", cn)
    
    if (isTRUE(all(chk1))) {
        print("pathological purity present")
        Pathological_Purity = fileData$Pathological_Purity
    } else {
        Pathological_Purity = rep(100*0.01,length(fileData$Gene))
        fileData = cbind(fileData,Pathological_Purity)
        Pathological_Purity = fileData$Pathological_Purity
    }
    
    if (isTRUE(all(chk2))) {
        print("computational purity present")
        Computational_Purity = fileData$Computational_Purity
    } else {
        Computational_Purity = rep(1,length(fileData$Gene))
        fileData = cbind(fileData,Computational_Purity)
        Computational_Purity = fileData$Computational_Purity
    }
    
    
    all_purities = list() 
    All_pathModels = data.frame()
    All_compModels = data.frame()
    All_allfitModels = data.frame()
    for (i in 1:length(unique(fileData$Sample_ID))) {
        tempID = unique(fileData$Sample_ID)[i]
        temp_df = filter(fileData, Sample_ID == tempID)
        
        # input to all-fit and lohgic
        temp_df <- data.frame(temp_df$Gene,temp_df$VAF,temp_df$Depth,temp_df$Copy_Number,temp_df$Pathological_Purity/100,temp_df$Computational_Purity) # allele freq in %
        colnames(temp_df) <- c("Gene","Allele_freq","Depth","Ploidy","PathologicalPurity","ComputationalPurity")
        write.table(temp_df, file = "tabledata.txt", sep = "\t" , dec = ".",
                    row.names = FALSE, col.names = TRUE)
        
        pur <- py_run_file("All-FIT2.py")
        
        # access output file
        mystring <- read_file("tabledata_out.txt")
        allfit_pur <- as.numeric(gsub('^.*tabledata_out\t\\s*|\\s*\t.*$', '', mystring))
        allfit_pur <- rep(100*allfit_pur,length(temp_df$Gene))
        
        
        #all_purities = rbind(all_purities,allfit_pur)
        all_purities <- append(all_purities, allfit_pur)
        
        # extract model predictions
        pathModels <- read.delim("pathologicalPurity_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
        compModdels <- read.delim("computationalPurity_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
        allfitModels <- read.delim("AllFIT_Models_tabledata_out.txt", header = TRUE, sep = "\t", dec = ".")
        
        # keep only gene name in SNV column
        pathModels$SNV <- sub("\\:.*","\\",pathModels$SNV)
        compModdels$SNV <- sub("\\:.*","\\",compModdels$SNV)
        allfitModels$SNV <- sub("\\:.*","\\",allfitModels$SNV)
        
        pathModels = pathModels[ order(match(pathModels$SNV, temp_df$Gene)), ]
        
        All_pathModels = rbind(All_pathModels,pathModels)
        All_compModels = rbind(All_compModels,compModdels)
        All_allfitModels = rbind(All_allfitModels,allfitModels)
    }
    all_purities = as.data.frame(unlist(all_purities))
    colnames(all_purities) <- "All_FIT_Purity"
    
    fileData = cbind(fileData,all_purities)
    # rearrange rows
    
    
    fileData_models <- data.frame(All_pathModels$Model,All_compModels$Model,All_allfitModels$Model)
    colnames(fileData_models) <- c("Pathological Purity Model","Computational Purity Model","All-FIT Purity Model")
    
    fileData = cbind(fileData,fileData_models)
    contents <- fileData
    
    return(contents)
}
## -------------------------------------------------------------------------------------------------------
foundationXML_VUS <- function(inputFile){
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    
    # ID
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    # VUSs
    VariantProperties <- gsub('^.*<VariantProperties>\\s*|\\s*</VariantProperties>.*$', '', xmlText)
    VariantProperties <- str_split(VariantProperties, " \n ")
    VariantProperties <- VariantProperties[[1]][] # all VUSs with attributes
    VUSs <- matrix(data=NA, nrow=length(VariantProperties), ncol=2)
    for (i in 1 : length(VariantProperties)) {
        line <- VariantProperties[i]
        gName <- gsub('^.*geneName=\"\\s*|\\s*\".*$', '', line)
        variantName <- gsub('^.*variantName=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(gName)>12) {
            gName = ""
            variantName = ""
        }
        
        VUSs[i,1] <- gName
        VUSs[i,2] <- variantName
    }
    
    # Define Dataframe: 
    df_VUS <- data.frame(VUSs)
    df_IDs <- data.frame(rep(ReportID,length(VUSs)))
    df_VUS <- cbind(df_IDs,df_VUS)
    colnames(df_VUS) <- c("Report_ID","Gene","Protein_Change")
    
    VUS <- df_VUS
    
    return(VUS)
}
## -------------------------------------------------------------------------------------------------------
foundationXML_CNA <- function(inputFile){
    
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    
    # ID
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    
    # Copy Number Alterations
    CopyNumberAlterations <- gsub('^.*<copy-number-alterations>\\s*|\\s*</copy-number-alterations>.*$', '', xmlText)
    CopyNumberAlterations <- str_split(CopyNumberAlterations, " </copy-number-alteration> ")
    CopyNumberAlterations<- CopyNumberAlterations[[1]][] # all CNAs with attributes
    
    CNA_geneName <- c()
    CNA_copyNumber <- c()
    CNA_position <- c()
    CNA_type <- c()
    CNA_status <- c()
    CNA_numberOfExons <- c()
    for (i in 1 : length(CopyNumberAlterations)) {
        line <- CopyNumberAlterations[i]
        cna_geneName <- gsub('^.*gene=\"\\s*|\\s*\".*$', '', line)
        
        cna_copynumber <- gsub('^.*copy-number=\"\\s*|\\s*".*$', '', line)
        
        cna_pos <- gsub('^.*position=\"\\s*|\\s*\".*$', '', line)
        
        cna_tp <- gsub('^.*type=\"\\s*|\\s*\".*$', '', line)
        
        cna_stat <- gsub('^.*status=\"\\s*|\\s*\".*$', '', line)
        
        cna_exons <- gsub('^.*number-of-exons=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(cna_copynumber)>5) {
            cna_geneName = ""
            cna_copynumber = ""
            cna_pos = ""
            cna_tp = ""
            cna_stat = ""
            cna_exons = ""
        }
        CNA_geneName[i] <- cna_geneName
        CNA_copyNumber[i] <- cna_copynumber
        CNA_position[i] <- cna_pos
        CNA_type[i] <- cna_tp
        CNA_status[i] <- cna_stat
        CNA_numberOfExons[i] <- cna_exons
        
    }
    df_IDs <- data.frame(rep(ReportID,length(CNA_geneName)))
    df_CNA <- data.frame(CNA_geneName,CNA_copyNumber,CNA_position,CNA_type,CNA_status,CNA_numberOfExons)
    df_CNA <- cbind(df_IDs,df_CNA)
    colnames(df_CNA) <- c("Report_ID","Gene","Copy_Number","Position","Type","Status","number_Exons")
    
    # Define Dataframe: 
    
    CN_Alterations <- df_CNA
    
    return(CN_Alterations)
}

## -------------------------------------------------------------------------------------------------------
foundationXML_Rearrangements <- function(inputFile){
    
    fileData <- xmlParse(inputFile)
    
    # Read the text from the file
    xmlText <- paste(readLines(inputFile), "\n", collapse="")
    
    # ID
    ReportID <- gsub('^.*<ReportId>\\s*|\\s*</ReportId>.*$', '', xmlText)
    
    # Rearrangements
    Rearrangements <- gsub('^.*<rearrangements>\\s*|\\s*</rearrangements>.*$', '', xmlText)
    Rearrangements <- str_split(Rearrangements, " </rearrangement> ")
    Rearrangements<- Rearrangements[[1]][] # all CNAs with attributes
    
    Rearrangement_targetedGene <- c()
    Rearrangement_description <- c()
    Rearrangement_inFrame <- c()
    Rearrangement_otherGene <- c()
    Rearrangement_Position1 <- c()
    Rearrangement_Position2 <- c()
    Rearrangement_status <- c()
    Rearrangement_type <- c()
    for (i in 1 : length(Rearrangements)) {
        line <- Rearrangements[i]
        ra_targetedGene <- gsub('^.*targeted-gene=\"\\s*|\\s*\".*$', '', line)
        
        ra_description <- gsub('^.*<rearrangement description=\"\\s*|\\s*\".*$', '', line)
        
        ra_inFrame <- gsub('^.*in-frame=\"\\s*|\\s*\".*$', '', line)
        
        ra_otherGene <- gsub('^.*other-gene=\"\\s*|\\s*\".*$', '', line)
        
        ra_pos1 <- gsub('^.*pos1=\"\\s*|\\s*\".*$', '', line)
        
        ra_pos2 <- gsub('^.*pos2=\"\\s*|\\s*\".*$', '', line)
        
        ra_status <- gsub('^.*status=\"\\s*|\\s*\".*$', '', line)
        
        ra_type <- gsub('^.*type=\"\\s*|\\s*\".*$', '', line)
        
        if (str_length(ra_targetedGene)>20) {
            ra_targetedGene = ""
            ra_description = ""
            ra_inFrame = ""
            ra_otherGene = ""
            ra_pos1 = ""
            ra_pos2 = ""
            ra_status = ""
            ra_type = ""
        }
        Rearrangement_targetedGene[i] <- ra_targetedGene
        Rearrangement_description[i] <- ra_description
        Rearrangement_inFrame[i] <- ra_inFrame
        Rearrangement_otherGene[i] <- ra_otherGene
        Rearrangement_Position1[i] <- ra_pos1
        Rearrangement_Position2[i] <- ra_pos2
        Rearrangement_status[i] <- ra_status
        Rearrangement_type[i] <- ra_type
        
    }
    df_IDs <- data.frame(rep(ReportID,length(Rearrangement_targetedGene)))
    df_RA <- data.frame(Rearrangement_targetedGene,Rearrangement_description,
                        Rearrangement_otherGene,Rearrangement_Position1,Rearrangement_Position2,
                        Rearrangement_inFrame,Rearrangement_status,Rearrangement_type)
    df_RA <- cbind(df_IDs,df_RA)
    colnames(df_RA) <- c("Report_ID","Tageted_Gene","Description","Other_Gene","Position1",
                         "Position2","In_Frame","Status","Type")
    
    # Define Dataframe: 
    
    Rearrangements <- df_RA
    
    
    return(Rearrangements)
}
############################################################################################################
# Define UI ----
ui <- fluidPage(theme = shinytheme("cosmo"),
                
                # App theme ----
                #shinythemes::themeSelector(),  # <--- Add this somewhere in the UI
                
                # Define the sidebar with one input
                sidebarPanel(width = 12,
                    selectInput("fileType", "File Type:", 
                                choices= c("General","Foundation_xml")),
                    hr(),
                    helpText(HTML("General: input file with minimum fields. <br/> Foundation xml: xml file specific to FoundationMedicine"))
                    ),
                textOutput("result"),
                
                # topbar layout
                topbar( fileInput("file1", "Choose File",
                                  multiple = FALSE)),
                # Horizontal line ----
                tags$hr(),
                
                mainPanel(
                    fluidRow(
                        column(12, tabsetPanel(
                            id = 'dataset',
                            # Output: Data file ----
                            tabPanel("Variants", DT::dataTableOutput("contents", height = "500px")),
                            tabPanel("VUS", DT::dataTableOutput("VUS", height = "500px")),
                            tabPanel("CN_Alterations", DT::dataTableOutput("CN_Alterations", height = "500px")),
                            tabPanel("Rearrangements", DT::dataTableOutput("Rearrangements", height = "500px")),
                            tabPanel("Plots",plotOutput("p1_histogramVAFs"))
                        ) # tabsetPanel end
                        ), # column end
                        #column(6, plotOutput("scatter1"))
                    ), # fluidRow end
                    
                    # horizontal line
                    tags$hr(),
                    
                    # plots for selcted row
                    fluidRow(
                        column(12,plotOutput("models_purities"))
                    )# fluid row end
                    
                    
                ) # mainPanel end
                
) # fluidPage end




# Define server logic to read selected file ----
server <- function(input, output,session) {
    
    inputData <- reactive({
        req(input$file1)
        inputFile <- input$file1$datapath
        
        if (req(input$fileType) == "General") {
            paste("Selected File type: ", input$fileType)
            # function related to general format
            contents <- generalFileFormat(inputFile)
        } else if (input$fileType == "Foundation_xml") {
            paste("Selected File type: Foundation_xml")
            
            # functions related to foundation xml
            contents <- foundationXML_general(inputFile)
        }
        
        
    })
    
    output$result <- renderText({
        paste("Selected File type: ", input$fileType)
    })
    
    
    output$contents <- DT::renderDataTable(
        datatable( data = inputData()
                   , extensions = 'Buttons',
                   selection = 'multiple'
                   , options = list( 
                       dom = "Blfrtip"
                       , buttons = 
                           list("copy", list(
                               extend = "collection"
                               , buttons = c("csv", "excel", "pdf")
                               , text = "Download"
                           ) ) # end of buttons customization
                   ) # end of options
        ) # end of datatables
    )
    
    
    
    # plot data 
    output$models_purities <- renderPlot({
        s = input$contents_rows_selected
        plotData <- data.frame(inputData()[s,c(14,16,17,18,9,10,11,23,24,25)])
        print(plotData)
        
        # models only
        p1<-ggplot(plotData, aes(x=Pathological_Purity, y=Allele_Freq*100)) +
            geom_point(aes(shape = Pathological.Purity.Model,  color = Pathological.Purity.Model)) + 
            geom_point(shape=22, fill="blue", color="blue", size=8,alpha = 1/10) +
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) +# germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + xlab("Pathological Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        p1 <- p1 + theme(legend.position = "right")
        
        p2<-ggplot(plotData, aes(x=Computational_Purity, y=Allele_Freq*100)) +
            geom_point(aes(shape = Computational.Purity.Model,  color = Computational.Purity.Model)) + 
            geom_point(shape=22, fill="blue", color="blue", size=8,alpha = 1/10) +
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) +# germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + xlab("Computational Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        p2 <- p2 + theme(legend.position = "right")
        
        p3<- ggplot(plotData, aes(x=AllFIT_Purity, y=Allele_Freq*100)) +
            geom_point(aes(shape = All.FIT.Purity.Model,  color = All.FIT.Purity.Model)) + 
            geom_point(shape=22, fill="blue", color="blue", size=8,alpha = 1/10) +
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) + # germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + xlab("All-FIT Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        p3 <- p3 + theme(legend.position = "right")
        
        p1|p2|p3
    })
    
    
    output$VUS <- DT::renderDataTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        inputFile <- input$file1$datapath
        VUS <- foundationXML_VUS(inputFile)
        
    },extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                   lengthMenu = list(c(10,25,50,-1),
                                     c(10,25,50,"All"))))
    
    output$CN_Alterations <- DT::renderDataTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        req(input$file1)
        inputFile <- input$file1$datapath
        CN_Alterations <- foundationXML_CNA(inputFile)
        
        
    },extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                   lengthMenu = list(c(10,25,50,-1),
                                     c(10,25,50,"All"))))
    
    output$Rearrangements <- DT::renderDataTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        req(input$file1)
        inputFile <- input$file1$datapath
        Rearrangements <- foundationXML_Rearrangements(inputFile)
        
        
    },extensions = 'Buttons',
    options = list(dom = 'Blfrtip',
                   buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                   lengthMenu = list(c(10,25,50,-1),
                                     c(10,25,50,"All"))))
    output$p1_histogramVAFs <- renderPlot({
        req(input$file1)
        inputFile <- input$file1$datapath
        contents <- foundationXML_general(inputFile)
        
        p1<-ggplot(contents, aes(x = Allele_Freq ))+
            geom_bar(aes(y = (..count..)/sum(..count..))) +
            scale_y_continuous(labels = percent) + theme_classic() +
            theme(legend.position="top") + xlim(c(0,1)) + ylim(c(0,1)) +
            xlab("Allele Frequencies") + ylab("Percentage (%)")
        
        p2 <- ggplot(contents,aes(x=Depth)) +
            geom_histogram(stat = "bin", binwidth = 50) + theme_classic() + 
            theme(legend.position = "top") + xlab("Depth") + ylab("Count")
        
        p3 <- ggplot(contents, aes(x=Pathological_Purity, y=Allele_Freq*100)) + 
            geom_point(aes(color = Gene)) +
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) +# germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + theme(legend.position = "top") + xlab("Pathological Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        
        p4 <- ggplot(contents, aes(x=Computational_Purity, y=Allele_Freq*100)) + 
            geom_point(aes(color = Gene)) +
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) +# germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + theme(legend.position = "top") + xlab("Computational Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        
        p5 <- ggplot(contents, aes(x=AllFIT_Purity, y=Allele_Freq*100)) + 
            geom_point(aes(color = Gene)) +
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) +# germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + theme(legend.position = "top") + xlab("All-FIT Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        
        (p1|p2)/(p3|p4|p5)
        
        
    })
    
    
    output$models <- renderPlot({
        s = input$x1_rows_selected
        # models only
        ggplot() + 
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color = '#00CCCC',
                       data = data.frame(x1=0,x2=100,y1=50,y2=50),curvature = 0) +# germline cn = 1
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#CC0033',
                       data = data.frame(x1=0,x2=100,y1=0,y2=100),curvature = 0.13) + # somatic LOH cnmut=1
            geom_abline(slope = 1, intercept = 0,color='#FF9900') + # somatic cn = 2
            geom_abline(slope = 0.5, intercept = 0,color='#66CCFF') +  # somatic cn = 1
            geom_abline(slope = 0.5, intercept = 50,color='#99FF99')+  # germline cn = 2
            geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2),color='#6699FF',
                       data = data.frame(x1=0,x2=100,y1=50,y2=100),curvature = 0.13) + # germline loh cn = 1
            theme_classic() + theme(legend.position = "top") + xlab("All-FIT Purity") + 
            ylab("Allele Frequency")  + xlim(c(0,100)) + ylim(c(0,100))
        
        
        
        
        
        
    })
}


# Run the app ----
shinyApp(ui = ui, server = server)


