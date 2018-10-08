# Instructions how to run ShinyApp

To run the shinyApp locally you should have R and the following packages installed: shiny, d3heatmap, DT, ggplot2, viridis, dplyr, tidyr, igraph, GGally. Please download the files on your computer in a new directory and in R type: 

library(shiny)

runApp("path_to_directory_containing_RData_folder")

Alternatively, you can use RStudio to visualise the shinyApp: Please download the files onto your computer in a new directory and open the ui.R file with RStudio. Then click "Run App" which is in the top right corner of the ui.R script for the application to run.
