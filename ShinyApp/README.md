# Instructions how to run ShinyApp

To run the shinyApp locally you should have R and the following packages installed: shiny, d3heatmap, DT, plotly, viridis, dplyr, rfigshare. Please download the files on your computer in a new directory and in R type: 

library(shiny)

runApp("path_to_directory_containing_RData_folder")

This app downloads data from figshare and therefore you should have a figshare to download the needed files (https://figshare.com/articles/chick_NC-GRN_RData/xxxxxxx and https://figshare.com/articles/chick_NC-GRN_images/xxxxxxx). You will be directed to authenticate the connection:

Use a local file ('.httr-oauth'), to cache OAuth access credentials between R sessions?

To which you should answer yes (1) and you will be directed to a figshare page to Log-in. Then you will have to allow access to ropensci and you will be redirected back to R. You should now have a working App. 

Alternatively, you can use RStudio to visualise the shinyApp: Please download the files onto your computer in a new directory and open the ui.R file with RStudio. Then click "Run App" which is in the top right corner of the ui.R script for the application to run. The same verification process will happen for figshare. 
