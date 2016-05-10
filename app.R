# Emilia F Gan
# Student Number: 1138607
# Email: efgan@uw.edu
# Copyright E Gan 2016
# Submitted in partial fulfillment of the requirements
# of the degree Master of Science

# Principal Component Analysis (PCA) Analysis

# Opening a CSV File and reading into R
# Modeled after: http://shiny.rstudio.com/gallery/upload-file.html
# Downloading plots learned from: https://www.youtube.com/watch?v=LSnWGmVkB6A

library(shiny)
library(ggplot2) # for plotting color plots with legend
library(edgeR) # for DGE  
library(limma)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 9*1024^2)

#--------------------------------------------------- UI CODE --------------------------------------------------------------
ui <- fluidPage(
    titlePanel("PCA Analysis"),
    hr(),
    tabsetPanel(
        tabPanel("PCA Plot",
                 sidebarLayout(
                     
                     # These elements will display in the sidebar panel at the L of the page.
                     sidebarPanel(
                         # Allows user to select data file from local file system.
                         fileInput('file1', 'Choose file to upload',
                                   accept = c(
                                       'text/csv',
                                       'text/comma-separated-values',
                                       'text/tab-separated-values',
                                       'text/plain',
                                       '.csv',
                                       '.tsv'
                                   )
                         ),
                         tags$hr(),
                         
                         # Note: Plot will only display color-coded groups if the header "Type" is present
                         # when data file is formatted with samples displayed in rows.
                         radioButtons('sep', 'Separator:',
                                      c(Comma=',',
                                        Semicolon=';',
                                        Tab='\t'),
                                      ','),
                         tags$hr(),
                         
                         # Collects formatting data for proper analysis.
                         radioButtons('geneDir', 'Genes are listed in:',
                                      c(Rows = 'Row',
                                        Columns = 'Col'),
                                      'Row'),
                         
                         tags$hr(),
                         
                         
                         # Collects information on how to do PCA.
                         radioButtons('dir', 'Display PCA Results For:',
                                      c(Rows = 'R',
                                        Columns = 'C'),
                                      'C'),
                         
                         tags$hr(),
                         
                         # Determines whether counts file represents normalized or unnormalized data.
                         # Normalized data are not processed further prior to PCA.
                         # Unnormalized data are normalized by the app using the TMM method as recommended
                         # by Peter Linsley of BRI using function written by Kristen Dang.
                         radioButtons('norm', 'Apply TMM Normalization:',
                                      c(No = 'No_Norm',
                                        Yes = 'Norm'),
                                      'No_Norm'),
                         tags$hr(),
                         
                         # If no grouping data provided for samples, default color is BLACK.
                         # If grouping data are provided, but color plot output not desired,
                         # the color and color key can be suppressed by selecting 'Black' radion button.
                         radioButtons('col', 'Display points:',
                                      c(Black = 'BLK',
                                        Color = 'CLR'),
                                      'BLK'),
                         
                         # Option to have points on plot displayed with different symbols for each sample subtype.
                         conditionalPanel(
                             condition = "input.col == 'CLR'",
                             checkboxInput('Shapes', 'Shapes', FALSE)
                         ),
                         
                         tags$hr(),
                         
                         # Allows user to input a custom title for the plot.
                         # If left blank, no title will be displayed.
                         textInput('title', "Plot Title", ''),
                         
                         # Suppresses error messages about NULL while awaiting user selection.
                         # Credit to: Joe Cheng[R Studio] at
                         # https://groups.google.com/forum/#!topic/shiny-discuss/FyMGa2R_Mgs
                         tags$style(type="text/css",
                                    ".shiny-output-error { visibility: hidden; }",
                                    ".shiny-output-error:before { visibility: hidden; }"
                         ),
                         width = 3
                     ),
                     
                     # These elements will display in the main panel (R side of page).
                     mainPanel(
                         plotOutput('pcaPlot',click = "plot_click"),
                         conditionalPanel(
                             condition = "input.col == 'CLR'",
                             downloadButton('downloadColor', 'Download the color plot as a pdf', class = "buttCol")
                         ),
                         conditionalPanel(
                             condition = "input.col == 'BLK'",
                             downloadButton('downloadBlack', 'Download the black & white plot as a pdf', class = "buttCol")
                         ),
                         tags$head(tags$style(".buttCol{background-color:#faf0e6;} .buttCol{color: black;}")),
                         tags$hr(),
                         conditionalPanel(
                             condition = "input.col == 'CLR'",
                             textOutput('warning')
                         ),
                         verbatimTextOutput("info"),
                         tags$hr(),
                         tableOutput('table')
                     )  # end mainPanel
                 )  # end sidebarLayout
        ),  # end first tabPanel
        tabPanel("Help",
                 tags$h1("File Formatting Issues", style = "color:darkblue"),
                 tags$p(" It is anticipated that most issues that arise while trying to use this application will be related
                        to file formatting. A very specific file format is required. This format has several variations, 
                        depending on the type of data included in the file. The figures below illustrate accepatable file formats, 
                        as they would appear in a spreadsheet program (such as Excel) before being coverted into one of the accepted
                        file formats (.csv, .tsv, or .txt).", style="font-size:20px"),
                 tags$img(height = 414  , width = 671 , src = "Permitted.PNG"),
                 tags$p(" In the figure above, genes are arranged by column and samples by row. Transposing this arrangement,
                        so that genes are in rows, while samples are in columns is permitted and should not cause any issues 
                        in the functioning of the application.", style="font-size:20px")
                 )
                 )  # end tabsetPanel
                 )  # end fluidPage

#----------------------------------------------- SERVER CODE ----------------------------------------
server <- function(input,output){
    # Function to read data from a file, after checking file in not NULL.
    getData <- function(dataFile){
        inFile <- dataFile
        if (is.null(inFile)) {
            return(NULL)
        }
        else {
            read.csv(inFile$datapath, header = TRUE, sep = input$sep, quote = input$quote, row.names = 1, stringsAsFactors = FALSE)
        }
    }
    
    # Read in the data from the input file.
    myData <- reactive({getData(input$file1)})
    
    # Get data dimensions as 2 item list [#rows, rcols].
    dataDimensions <- reactive(dim(myData()))
    
    # Get information on whether sample data are in ROWS or COLUMNS.
    direction <- reactive({input$dir})
    
    #-----------------------------------------------------------------------------------------
    #           TASK ONE: FORMAT THE DATA
    #-----------------------------------------------------------------------------------------
    
    # Get the data into a "tidier" format with row names identified as such and
    # Sample/Gene identifiers separated from the numerical data table elements.
    
    # 1) Check that ROW TYPE DATA is present 
    # i.e. Last COLUMN contains TYPE data
    hasGroupDataInLastCol <- function(df){
        if(!is.null(df$Type)){
            return(TRUE)
        }
        return(FALSE)
    }
    
    # 2) Check that COLUMN TYPE DATA is present 
    # i.e. Last ROW contains TYPE data
    hasGroupDataInLastRow <- function(df){
        numRows = dataDimensions()[1]
        if(row.names(df)[numRows] == "Type"){
            return(TRUE)
        }
        return(FALSE)
    }
    
    # 3) Get GROUP TYPE information from last COLUMN
    # (check if present first using hasGroupDataInLastCol)
    # Need to specify to R that these should be considered a factor with levels.
    getGroupInfoLastCol <- function(df){
        dimensions <- dim(df)
        numRows <- dimensions[1]
        if(isTRUE(hasGroupDataInLastRow(df))) {
            return(as.factor(df$Type[1:numRows-1]))
        }
        else {
            return(as.factor(df$Type))
        }
    }
    
    # 4) Get GROUP TYPE information from last ROW
    # (check if present first using hasGroupDataInLastRow before calling this function)
    # Need to specify to R that these should be considered a factor with levels. 
    getGroupInfoLastRow <- function(df){  
        dimensions <- dim(df)
        numRows <- dimensions[1]
        numCols <- dimensions[2]
        
        if(isTRUE(hasGroupDataInLastCol(df))){
            lastRowInfo <- df[numRows,1:numCols-1]
        } else {
            lastRowInfo <- df[numRows,1:numCols]
        }
        
        lastRowInfo <- data.frame(t(lastRowInfo))
        
        return(t(lastRowInfo$Type)) 
    }
    
    # 5) Now, we're ready to do the actual reformatting:
    # Remove the first column from the data frame -- it should contain text (gene or sample names).
    # If it contains data, the input file is NOT PROPERLY FORMATTED!
    # At the end of this code block, only numerical values should be in the body of the data frame.
    # The data frame's column and row names should now be appropriately recognized as such by R.
    formattedData <- function(){
        dimensions <- dim(myData())
        numRows <- dimensions[1]
        numCols <- dimensions[2]
        
        formatted <- myData()
        # If there is GROUP data in the last column, it needs to be removed.
        if(isTRUE(hasGroupDataInLastCol(formatted)) == TRUE) {
            numCols = numCols - 1
        }
        # If there is group data in the last row, it needs to be removed.
        if(isTRUE(hasGroupDataInLastRow(formatted)) == TRUE){
            numRows = numRows - 1
        }
        # reformat
        formatted <- formatted[1:numRows, 1:numCols]
        # get the row names and store them as "r_names"
        r_names <- row.names(formatted)
        
        formatted <- lapply(formatted, as.numeric)
        formatted <- as.data.frame(formatted)
        
        # Add back the row name info as formal row names since they got lost 
        # in the two steps above.
        row.names(formatted) <- r_names
        
        return(formatted)
    }
    
    # Normalize the counts data using TMM, if user selects this option. 
    normalized <- function(){
        if(input$geneDir == 'Col') {
            counts_m <- t(formattedData())
        }
        else {
            counts_m <- formattedData()
        }
        dge <- DGEList(counts=counts_m)
        # can now use TMM to normalize the counts data
        # this function Calculates normalization factors to scale the raw library sizes
        # using the weighted trimmed mean of M-values (to the reference), 
        # where the weights are from the delta method on Binomial data
        data.dge = calcNormFactors(dge, method=c("TMM"), refColumn = NULL,logratioTrim = .3, sumTrim = 0.05, doWeighting=TRUE, Acutoff=-1e10, p=0.75)
        #### NEXT: Use this function from Kristen Dang ####
        # This function uses the library size from the datafile and the
        # normalization factors calculated to return a matrix of normalized counts
        # See DGEList documentation: default for lib.size = colSums(counts)
        # default norm.factors = rep(1, ncol(counts)) (i.e. a vector of '1's)
        return_TMM_matrix_noRound=function(inDGEObject){
            eff.libsize = inDGEObject$samples$lib.size * inDGEObject$samples$norm.factors
            tndata = 1e6*t(inDGEObject$counts)/eff.libsize
            return(t(tndata))
        }
        
        # Call the function, then convert matrix obtained to a data frame:
        data.tmm = return_TMM_matrix_noRound(data.dge) #### A function from Kristen Dang
        data.tmm = as.data.frame(data.tmm)
        
        if(input$geneDir == 'Col'){
            return(t(data.tmm))
        }
        
        return(data.tmm)
    }
    
    #-----------------------------------------------------------------------------------------
    #           TASK TWO: PERFORM PCA
    #-----------------------------------------------------------------------------------------
    
    # Perform PCA
    # If 'R' selected, remember to transpose the data.
    pcaTest <- function(){
        if(input$norm == 'Norm'){
            #forPCA <- as.data.frame(normalized())
            forPCA <- normalized()
        }
        else {
            forPCA <- formattedData()
        }
        
        if(input$dir == 'R'){
            pr.out <- prcomp(t(forPCA))
        }
        else {
            pr.out <- prcomp(forPCA)
        }
        return(pr.out)
    }
    
    # create a data frame with the sample names and PC1 and PC2 eigenvectors
    # Number of rows will vary depending on whether plotting samples or genes in PC1 v Pc2 plot. 
    dataframe <- reactive({
        if(input$geneDir == 'Row') {
            if(input$dir == 'R'){
                if(input$norm != 'Norm'){
                    Names <- row.names(formattedData())
                    
                }
                else {
                    Names <- row.names(normalized())
                }
            }
            else {
                if(input$norm != 'Norm'){
                    Names <- colnames(formattedData())
                }
                else {
                    Names <- colnames(normalized())
                }
            }
        }
        else {
            if(input$dir == 'R'){
                if(input$norm != 'Norm'){
                    Names <- row.names(formattedData())
                    
                }
                else {
                    Names <- row.names(normalized())
                }
            }
            else {
                if(input$norm != 'Norm'){
                    Names <- colnames(formattedData())
                }
                else {
                    Names <- colnames(normalized())
                }
            }
        }
        data.frame("PC_1" = pcaTest()$rotation[,1], "PC_2" = pcaTest()$rotation[,2], row.names = Names)
    })
    
    # Get the proportions of variance info for the axes labels
    proportions <- reactive({
        cumsum(pcaTest()$sdev^2 / sum(pcaTest()$sdev^2))
    })
    
    #-----------------------------------------------------------------------------------------
    #           TASK THREE: PLOT PCA RESULTS AND DOWNLOAD PLOT
    #-----------------------------------------------------------------------------------------
    
    # Get GROUP MEMBERSHIP information for COLOR-CODING the PC1 v PC2 Plot
    typeInfo <- function(){
        # If elements to show in plot are arranged in Rows and attributes in Columns:
        if(input$dir == 'R') {
            #   then Element Group data located in last column
            if(isTRUE(hasGroupDataInLastCol(myData()))){
                Groups = getGroupInfoLastCol(myData())
            }
            else{
                Groups = NULL
            }
        }
        # If elements to show in plot are arranged in Columns and attributes in Rows:
        else {
            #   then Element Group data in last row
            if(isTRUE(hasGroupDataInLastRow(myData()))){
                Groups = t(getGroupInfoLastRow(myData()))
            }
            else{
                Groups = NULL
            }
        }
        return(Groups)
    }
    
    # Display noted to user that color selections are only functional
    # when appropriate group type data is supplied.
    output$warning <- renderText({
        if(input$col == 'CLR' && is.null(typeInfo())){
            "COLOR options will be non-functional as no GROUPING DATA has been supplied."
        }
        else {""}
    }) 
    
    # Plot the 2 pricipal components. 
    plotInput = function(){
        axes <- proportions()
        xaxis <- axes[1]
        yaxis <- axes[2] - xaxis
        df = dataframe()
        
        if(is.null(typeInfo()) | input$col == "BLK") {
            print("TypeInfo seems to be null")
            plot(df$PC_1, df$PC_2, pch = 19, cex = 1.5,main = input$title, xlab = paste(
                "PC 1 (", round(xaxis,3)*100, "%)", sep=""), ylab = paste("PC 2 (", round(yaxis,3)*100, "%)", sep=""))
            
        } else {
            
            Key = typeInfo()
            
            if(input$Shapes == TRUE){
                ggplot(df, aes(x = df$PC_1, y = df$PC_2)) +
                    geom_point(mapping = aes(colour = Key, shape = Key), size = 5) + labs(x = paste("PC 1 (", round(xaxis,3)*100, "%)", sep=""), y = paste("PC 2 (", round(yaxis,3)*100, "%)", sep="")) +
                    labs(title = input$title) + theme_bw() 
            }
            else{
                ggplot(df, aes(x = df$PC_1, y = df$PC_2)) +
                    geom_point(mapping = aes(colour = Key), size = 5) + labs(x = paste("PC 1 (", round(xaxis,3)*100, "%)", sep=""), y = paste("PC 2 (", round(yaxis,3)*100, "%)", sep="")) +
                    labs(title = input$title) + theme_bw()  
            }
        }
    }
    
    output$pcaPlot <- renderPlot({
        plotInput()
    })
    
    # Display point data for points on plot close to location of mouse click 
    # If multiple points are located very close together, all will be listed.
    output$info <- renderPrint({
        nearPoints(dataframe(), input$plot_click, xvar = "PC_1" , yvar = "PC_2")
    })
    
    #--------------------------------------------------------------------------------
    # Download handler for ggplot2 plot.
    output$downloadColor <- downloadHandler(           
        filename = 'downloadedPlot.pdf',
        content = function(file) {
            device <- function(..., width, height) {
                grDevices::pdf(..., width = width, height = height)
            }
            ggsave(file, plot = plotInput(), device = device)
        }
    )
    
    # Download handler for basic R plot.
    output$downloadBlack <- downloadHandler( 
        filename = 'downloadedPlot.pdf',
        content = function(file) {
            pdf(file)
            plot(dataframe()$PC_1, dataframe()$PC_2, pch = 19, cex = 1.5,main = input$title, xlab = paste(
                "PC 1 (", round(proportions()[1],3)*100, "%)", sep=""), ylab = paste("PC 2 (", round((proportions()[2] - proportions()[1]),3)*100, "%)", sep=""))
            
            dev.off()
            
        }
    )
    
    
    # Display data table below plot.
    output$table <- renderTable(formattedData())
}

shinyApp(ui = ui, server = server)