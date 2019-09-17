rm(list=ls())

# load packages ----------------------
library(d3heatmap)
library(DT)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(igraph)
library(GGally)

# Server ------------------------

shinyServer(function(input, output, session){
  
  # load necessary data ----------------------
  withProgress(message = 'Please wait',
               detail = 'Loading the files...', value = 0, {
                  load('RData/RNASEQ.RData')
                  incProgress(1/9)
                  load("RData/diffbind_data.RData")
                  incProgress(2/9)
                  load("RData/cor_genes.RData")
                  incProgress(3/9)
                  load("RData/wgcnaclusters.RData")
                  incProgress(4/9)
                  load("RData/wgcnagenes.RData")
                  incProgress(5/9)
                  load("RData/enhancer_data.RData")
                  incProgress(6/9)
                  load("RData/cluster_RN.RData")
                  incProgress(7/9)
                  load("RData/TFBS_motifs.RData")
                  incProgress(8/9)
                  load("RData/DiffBind_motifs.RData")
               }
  )
  
  
###  RNA-seq Analysis Part -------------------------------------------------------------
  updateSelectizeInput(session, 'Gene_Name', choices = unique(de_5to6$Geneid), selected="SOX10",server = TRUE)
  
  output$SelectedGene <- renderText({ 
    paste("Gene name: ", input$Gene_Name )
  })

  geneinfo <- reactive({
    if (input$Gene_Name=="") {
      return("") 
    } else {
      annotation <- de_5to6[which(de_5to6$Geneid==input$Gene_Name),1:3]
      return(annotation) 
    }
  })
  
  output$Ensembl_Geneid <- renderText({ 
    paste("Ensembl ID: ", paste(geneinfo()["EnsemblID"],collapse = ","))
  })
  
  output$Full_Gene_Name <- renderText({ 
      paste("Full name: ",  paste(geneinfo()["GeneName"],collapse = ","))
  })
  
  fpkms <- reactive({
    if (input$Gene_Name=="")
      return("") 
    fpkms_mygene <- filter(fpkm_genes,ensembl_ID==geneinfo()$EnsemblID)[,3:12]  %>% gather()
    fpkms_mygene$group <- ifelse(grepl("6ss", fpkms_mygene$key),"5-6ss",ifelse(grepl("10ss", fpkms_mygene$key),"8-10ss","Negative"))
    return(fpkms_mygene)
  })

  output$plot_fpkm <- renderPlot({
    if (input$Gene_Name=="")
      return("") 
    ggplot(fpkms()) + geom_boxplot(aes(factor(group), value), fill=c("green4","dodgerblue","grey")) + geom_point(aes(factor(group), value)) + coord_flip() +
        labs(y = "FPKM", x = "") + theme_classic()
      #plot_ly(fpkms(), x = ~value, y = ~group, color = ~group, type = "box", marker = list(color=c("green","blue","grey")))

  })

  # DE genes ----------------------------
  output$table_de <- renderTable({
    if(input$Gene_Name=="")
      return(NULL)
    # subset df
    my_5to6 <-  select(filter(de_5to6, de_5to6$Geneid == input$Gene_Name),c(EnsemblID,log2FoldChange,padj))
    my_8to10 <-  select(filter(de_8to10, de_8to10$Geneid == input$Gene_Name),c(EnsemblID,log2FoldChange,padj))
    my_5to6vs8to10 <- select(filter(de_5to6vs8to10, de_5to6vs8to10$Geneid == input$Gene_Name),c(EnsemblID,log2FoldChange,padj))
    all_together <- rbind(my_5to6,my_8to10,my_5to6vs8to10)
    all_together$Stage <- rep(c("Stage 5 to 6 vs neg", "Stage 8 to 10 vs neg","Stage 5 to 6 vs 8 to 10"), times=c(nrow(my_5to6), nrow(my_8to10),nrow(my_5to6vs8to10)))
    all_together <- all_together[,c(4,1:3)]
    all_together$Significant <- ifelse(all_together$padj<0.05,"Yes","No")
    all_together$padj <- format(all_together$padj, scientific=T)
    names(all_together)[4] <- "Padjusted"
    return(all_together)
  })
  ### WGCNA clusters ------------------------------------------------
  
  updateSelectizeInput(session, 'WGCNA_Choice', choices = c(paste0("Cluster",1:13),unique(wgcna_genes)), selected="Cluster1",server = TRUE)
  
  wgcnacluster <- reactive({
    if(input$WGCNA_Choice==""){
      return()
    } else {
      if(grepl("Cluster",input$WGCNA_Choice)) {
        mycluster <- gsub("Cluster",replacement = "", input$WGCNA_Choice)
        return(paste0("cl",mycluster))
      } else {
        genes_cl <- data.frame(Geneid=c(cl1$Geneid,cl2$Geneid,cl3$Geneid,cl4$Geneid,cl5$Geneid,cl6$Geneid,cl7$Geneid,cl8$Geneid,cl9$Geneid,cl10$Geneid,cl11$Geneid,cl12$Geneid, cl13$Geneid),
                               cluster=c(rep(c("cl1","cl2","cl3","cl4","cl5","cl6","cl7","cl8","cl9","cl10","cl11","cl12","cl13"), times=c(nrow(cl1),nrow(cl2),nrow(cl3),nrow(cl4),nrow(cl5),nrow(cl6),nrow(cl7),nrow(cl8),nrow(cl9),nrow(cl10),nrow(cl11),nrow(cl12), nrow(cl13)))))
        mycluster2 <- as.character(unique(genes_cl$cluster[genes_cl$Geneid==input$WGCNA_Choice]))
        return(mycluster2)
      }
    }
  })
  
  output$text_wgcna <- renderText({
    if(input$WGCNA_Choice==""){
      return("")
    } else {
      if (length(wgcnacluster())==1) {
        return(paste0("Showing the genes for cluster: Cluster", gsub("cl",replacement="", wgcnacluster()),"."))
      } else {
        return(paste0("This gene belongs to multiple clusters: ", paste("Cluster",gsub("cl",replacement="", wgcnacluster()),collapse=","), ". Showing only the results from the first of them. To view one of the others, choose it directly from the above options. "))
      }
    }
  })  
  
  output$wgcnaimage <- renderImage({
    if(input$WGCNA_Choice==""){
      # setting the ClusterI as a default so that no error comes up
      #wgcna_file <- normalizePath(file.path('./images/Cluster_I_WGCNA.png', sep=''))
      list(src = file.path('./images/Cluster_I_WGCNA.png', sep=''),contentType = 'image/png', width = 600,height = 600)
    } else {
      #wgcna_file <- normalizePath(file.path('./images', paste0('Cluster_', as.roman(gsub("cl",replacement="", wgcnacluster()[1])), '_WGCNA.png', sep='')))
      list(src = file.path('./images', paste0('Cluster_', as.roman(gsub("cl",replacement="", wgcnacluster()[1])), '_WGCNA.png', sep='')),contentType = 'image/png', width = 600,height =600)
    }
  }, deleteFile = FALSE)
  
  output$WGCNAtable <- DT::renderDataTable({
    if(input$WGCNA_Choice==""){
      return()
    } else {
      DT::datatable(get(wgcnacluster()[1]),escape=FALSE)
    }
  })      
  
### Regulatory elements ----------------------------
  updateSelectizeInput(session, 'Gene_Name_reg', choices = gene_names_2ndtab, selected="SOX10",server = TRUE)

# ATTN: Again here all the cluster 9 have been replaced by cluster 1 in principle but the variables remain the same to save me from changing/breaking everything  
  # Enhancers ----------------------------
  output$table_enh_clust <- renderTable({
    if(input$Gene_Name_reg=="")
      return(NULL)
    # subset enhancers clusters
    my_enh_cl3 <-  select(filter(enh_cl3, enh_cl3$gene_name == input$Gene_Name_reg),-gene_name)
    my_enh_cl4 <-  select(filter(enh_cl4, enh_cl4$gene_name == input$Gene_Name_reg),-gene_name)
    my_enh_cl9 <-  select(filter(enh_cl9, enh_cl9$gene_name == input$Gene_Name_reg),-gene_name)
    all_together_enh <- rbind(my_enh_cl3,my_enh_cl4,my_enh_cl9)
    all_together_enh$Cluster <- rep(c("Cluster 3", "Cluster 4", "Cluster 1"), times=c(nrow(my_enh_cl3), nrow(my_enh_cl4), nrow(my_enh_cl9)))
    return(all_together_enh)
  })

  # TFBS in enhancers ----------------------------

  mygene_df <- reactive({
    if (input$Gene_Name_reg=="")
      return("") 
    mygene1 <- filter(cluster3_motifs, gene_name==input$Gene_Name_reg) 
    #mygene2 <- filter(cluster4_motifs, gene_name==input$Gene_Name_reg) 
    mygene2 <- filter(cluster9_motifs, gene_name==input$Gene_Name_reg)
    if(all(c(nrow(mygene1), nrow(mygene2))!=0))  {
      mygene <- rbind(mygene1, mygene2)
    } else {
      myind <- which(c(nrow(mygene1), nrow(mygene2))!=0)
      if(length(myind)==1) mygene <- get(paste("mygene", myind[1], sep=""))
      if(length(myind)==0) return(NULL)
    }
    mygene$cluster <- rep(c("Cl3","Cl1"), times=c(nrow(mygene1), nrow(mygene2)) )
    row.names(mygene) <- paste(mygene$cluster, mygene$Chr, mygene$Start, mygene$End, sep="_")
    mygene <- select_if(mygene, not_all_na) %>% select(-c(gene_name, cluster)) 
    mygene <- mygene[,c(5:ncol(mygene))]
    if (nrow(mygene)<2) {
      mygene <- rbind(NA,mygene)
      row.names(mygene)[1] <- ""
    }
    return(mygene)
  })
  
  mygene_mat_df <- reactive({
    if (input$Gene_Name_reg=="")
      return("") 
    mygene_mat1 <- filter(motifs_3, gene_name==input$Gene_Name_reg)
    mygene_mat2 <- filter(motifs_9, gene_name==input$Gene_Name_reg) 
    if(all(c(nrow(mygene_mat1), nrow(mygene_mat2))!=0))  {
      mygene_mat <- rbind(mygene_mat1, mygene_mat2)
    } else {
      myind <- which(c(nrow(mygene_mat1), nrow(mygene_mat2))!=0)
      if(length(myind)==1) mygene_mat <- get(paste("mygene_mat", myind[1], sep=""))
      if(length(myind)==0) return(NULL)
    }
    mygene_mat$cluster <- rep(c("Cl3","Cl1"), times=c(nrow(mygene_mat1), nrow(mygene_mat2)) )
    row.names(mygene_mat) <- paste(mygene_mat$cluster, mygene_mat$Chr, mygene_mat$Start, mygene_mat$End, sep="_")
    mygene_mat <- select_if(mygene_mat, not_all_0) %>% select(-c(gene_name, cluster)) 
    mygene_mat <- mygene_mat[,c(5:ncol(mygene_mat))]
    if (nrow(mygene_mat)<2) {
      mygene_mat <- rbind(NA,mygene_mat)
      row.names(mygene_mat)[1] <- "" 
    }
    return(mygene_mat)
    })
    
  output$TFBS <- renderD3heatmap({
    if(any(input$Gene_Name_reg=="")){
      return()
    } else if (is.null(mygene_mat_df())) {
      return()
    } else {
      d3heatmap(mygene_mat_df(),dendrogram ="none",
            colors = "YlOrRd",
            show_grid = FALSE,
            anim_duration = 0,
            xaxis_font_size = "6px",
            yaxis_font_size = "8px",
            xaxis_height = 500,
#            yaxis_width = 150,
            cellnote=mygene_df())
    }
  })
  
  # ATAC-Seq - DiffBind----------------------------
  output$table_diffbind <- renderTable({
    if(input$Gene_Name_reg=="")
      return(NULL)
    # subset df
    my_diffbind_5to6 <-  select(filter(diffbind_5to6, diffbind_5to6$gene_name == input$Gene_Name_reg),-gene_name)
    my_diffbind_8to10 <-  select(filter(diffbind_8to10, diffbind_8to10$gene_name == input$Gene_Name_reg),-gene_name)
    all_together_diff <- rbind(my_diffbind_5to6,my_diffbind_8to10)
    all_together_diff$Stage <- rep(c("Stage 5 to 6 vs neg", "Stage 8 to 10 vs neg"), times=c(nrow(my_diffbind_5to6), nrow(my_diffbind_8to10)))
    all_together_diff$p.value <- format(all_together_diff$p.value, scientific=T)
    all_together_diff$FDR <- format(all_together_diff$FDR, scientific=T)
    return(all_together_diff)
  })

    # TFBS in DiffBind data ----------------------------
  
  mygene_diffbind_df <- reactive({
    if (input$Gene_Name_reg=="")
      return("") 
    mygene1_diff <- filter(diffbind_5to6_motifs, gene_name==input$Gene_Name_reg) 
    mygene2_diff <- filter(diffbind_8to10_motifs, gene_name==input$Gene_Name_reg)
    if(all(c(nrow(mygene1_diff), nrow(mygene2_diff))!=0))  {
      mygene_diff <- rbind(mygene1_diff, mygene2_diff)
    } else {
      myind <- which(c(nrow(mygene1_diff), nrow(mygene2_diff))!=0)
      if(length(myind)==1) mygene_diff <- get(paste("mygene", myind[1],"_diff", sep=""))
      if(length(myind)==0) return(NULL)
    }
    mygene_diff$group <- rep(c("Stage 5-6","Stage 8-10"), times=c(nrow(mygene1_diff), nrow(mygene2_diff)) )
    row.names(mygene_diff) <- paste(mygene_diff$group, mygene_diff$Chr, mygene_diff$Start, mygene_diff$End, sep="_")
    mygene_diff <- select_if(mygene_diff, not_all_na) %>% select(-c(gene_name, group)) 
    mygene_diff <- mygene_diff[,c(5:ncol(mygene_diff))]
    if (nrow(mygene_diff)<2) {
      mygene_diff <- rbind(NA,mygene_diff)
      row.names(mygene_diff)[1] <- ""
    }
    return(mygene_diff)
  })
  
  mygene_diff_mat_df <- reactive({
    if (input$Gene_Name_reg=="")
      return("") 
    mygene_diff_mat1 <- filter(TFBS_diffbind_5to6, gene_name==input$Gene_Name_reg)
    mygene_diff_mat2 <- filter(TFBS_diffbind_8to10, gene_name==input$Gene_Name_reg) 
    if(all(c(nrow(mygene_diff_mat1), nrow(mygene_diff_mat2))!=0))  {
      mygene_diff_mat <- rbind(mygene_diff_mat1, mygene_diff_mat2)
    } else {
      myind <- which(c(nrow(mygene_diff_mat1), nrow(mygene_diff_mat2))!=0)
      if(length(myind)==1) mygene_diff_mat <- get(paste("mygene_diff_mat", myind[1], sep=""))
      if(length(myind)==0) return(NULL)
    }
    mygene_diff_mat$group <- rep(c("Stage 5-6","Stage 8-10"), times=c(nrow(mygene_diff_mat1), nrow(mygene_diff_mat2)) )
    row.names(mygene_diff_mat) <- paste(mygene_diff_mat$group, mygene_diff_mat$Chr, mygene_diff_mat$Start, mygene_diff_mat$End, sep="_")
    mygene_diff_mat <- select_if(mygene_diff_mat, not_all_0) %>% select(-c(gene_name, group)) 
    mygene_diff_mat <- mygene_diff_mat[,c(5:ncol(mygene_diff_mat))]
    if (nrow(mygene_diff_mat)<2) {
      mygene_diff_mat <- rbind(NA,mygene_diff_mat)
      row.names(mygene_diff_mat)[1] <- "" 
    }
    return(mygene_diff_mat)
  })
  
  output$TFBS_diffbind <- renderD3heatmap({
    if(any(input$Gene_Name_reg=="")){
      return()
    } else if (is.null(mygene_diff_mat_df())) {
      return()
    } else {
      d3heatmap(mygene_diff_mat_df(),dendrogram ="none",
                colors = "YlOrRd",
                show_grid = FALSE,
                anim_duration = 0,
                xaxis_font_size = "6px",
                yaxis_font_size = "8px",
                xaxis_height = 500,
                #            yaxis_width = 150,
                cellnote=mygene_diffbind_df())
    }
  })

    # Regulatory networks ----------------------------
  
  updateSelectizeInput(session, "Gene_Name_RN", choices = unique(c(cluster3_RN$Target_Gene, cluster9_RN$Target_Gene)), selected="SOX10",server = TRUE)

  gene_RN_cl3 <- reactive({
    if (input$Gene_Name_RN=="")
      return("") 
    filtered_gene <- filter(cluster3_RN, Target_Gene==input$Gene_Name_RN) %>%
    group_by(Upstream_Input) %>% 
    summarize(n()) %>% 
    left_join(tf_info, by=c("Upstream_Input"="Gene"))
    return(filtered_gene)
  })
  
  gene_RN_cl9 <- reactive({
    if (input$Gene_Name_RN=="")
      return("") 
    filtered_gene <- filter(cluster9_RN, Target_Gene==input$Gene_Name_RN) %>%
      group_by(Upstream_Input) %>% 
      summarize(n()) %>% 
      left_join(tf_info, by=c("Upstream_Input"="Gene"))
    return(filtered_gene)
  })
  
  output$plot_rn_3 <- renderPlot({
    if (input$Gene_Name_RN=="")
      return("") 
    if (nrow(gene_RN_cl3())<1) {
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("This gene is not in Cluster 3."), 
           cex = 1.6, col = "black")
    } else{
      geneedges <- cbind(gene_RN_cl3()[,1], to= input$Gene_Name_RN, gene_RN_cl3()[,2], gene_RN_cl3()[,3])
      geneedges$to <- as.character(geneedges$to)
      names(geneedges)[c(1,3,4)] <- c("from","weight","type")
      if (any(geneedges$from==input$Gene_Name_RN)) geneedges <- geneedges[c(which(geneedges$from==input$Gene_Name_RN),seq(1, nrow(geneedges))[geneedges$from!=input$Gene_Name_RN]),]
      genenodes <- data.frame(id=unique(c(geneedges$from, geneedges$to)))
      
      net <- graph.data.frame(geneedges, genenodes, directed=T)
      E(net)$width <- geneedges$weight
      #E(net)$edge.lty <- ifelse(geneedges$type=="activator",1,2)
      l <- layout_with_gem(net) # or layout_with_dh?
      
      plot(net, edge.arrow.size=.4,  vertex.color="slategray1", vertex.size=2.8* max(nchar(geneedges$from)),  layout=l, 
           edge.color=c("khaki","gold","orange","red","darkred","sienna","tan4")[geneedges$weight],vertex.shape	="rectangle", cex.col="black",
           #edge.lty = rep(c(1,2), each=16))
           edge.lty = ifelse(geneedges$type=="activator",1,2))
      legend(x=-1.5, y=-1.1, paste("Occurences =",sort(unique(geneedges$weight))), lty=1, lwd=sort(unique(geneedges$weight)),
             col=unique(c("khaki","gold","orange","red","darkred","sienna","tan4")[sort(geneedges$weight)]), pt.cex=2, cex=.8, bty="n", ncol=1)
      legend(x=0.8, y=-1.1, c("activator","repressor"), lty=c(1,2), pt.cex=2, cex=.8, bty="n", ncol=1)
    }
  })
  
  output$plot_rn_9 <- renderPlot({
    if (input$Gene_Name_RN=="")
      return("") 
    if (nrow(gene_RN_cl9())<1) {
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, paste("This gene is not in 1."), 
           cex = 1.6, col = "black")
    } else{
    geneedges <- cbind(gene_RN_cl9()[,1], to= input$Gene_Name_RN, gene_RN_cl9()[,2], gene_RN_cl9()[,3])
    geneedges$to <- as.character(geneedges$to)
    names(geneedges)[c(1,3,4)] <- c("from","weight","type")
    if (any(geneedges$from==input$Gene_Name_RN)) geneedges <- geneedges[c(which(geneedges$from==input$Gene_Name_RN),seq(1, nrow(geneedges))[geneedges$from!=input$Gene_Name_RN]),]
    genenodes <- data.frame(id=unique(c(geneedges$from, geneedges$to)))
    
    net <- graph.data.frame(geneedges, genenodes, directed=T)
    E(net)$width <- geneedges$weight
    #E(net)$edge.lty <- ifelse(geneedges$type=="activator",1,2)
    l <- layout_with_gem(net) # or layout_with_dh?
    
    plot(net, edge.arrow.size=.4,  vertex.color="grey79", vertex.size=2.8* max(nchar(geneedges$from)),  layout=l, 
         edge.color=c("khaki","gold","orange","red","darkred","sienna","tan4")[geneedges$weight],vertex.shape	="rectangle", cex.col="black",
         #edge.lty = rep(c(1,2), each=16))
         edge.lty = ifelse(geneedges$type=="activator",1,2))
    legend(x=-1.5, y=-1.1, paste("Occurences =",sort(unique(geneedges$weight))), lty=1, lwd=sort(unique(geneedges$weight)),
           col=unique(c("khaki","gold","orange","red","darkred","sienna","tan4")[sort(geneedges$weight)]), pt.cex=2, cex=.8, bty="n", ncol=1)
    legend(x=0.8, y=-1.1, c("activator","repressor"), lty=c(1,2), pt.cex=2, cex=.8, bty="n", ncol=1)
    }
  })
  
    
### Coexpression correlations ------------------------------------------------

  updateSelectizeInput(session, 'Cor_Gene', choices = row.names(my_exprs), selected="SOX10",server = TRUE)

  cor_matrix <- reactive({
    if(any(input$Cor_Gene=="")){
      return()
    } else {
      cor(t(my_exprs),t(my_exprs[input$Cor_Gene,,drop=F]))
    }
  })

  correl_genes <- reactive({
    if (length(input$Cor_Gene)==1) {
      corgenes <- row.names(cor_matrix())[cor_matrix()>input$cor_thres]
    } else {
      corgenes <- input$Cor_Gene
    }
    return(corgenes)
  })

  output$text_cor <- renderText({
    if(any(input$Cor_Gene=="")){
      return("")
    } else {
      if (length(correl_genes())<2) {
        return("No correlated genes for these choices, please select another gene or a lower correlation threshold.")
      } else if (length(correl_genes())>500) {
        return("Too many correlated genes for these choices, please select another gene or a higher correlation threshold.")
      } else {
        return(NULL)
      }
    }
  })

  output$corgenelen <- reactive({
    length(input$Cor_Gene)==1
  })
  
  outputOptions(output, "corgenelen", suspendWhenHidden = FALSE) 
  
  output$heatmap_coexpr <- renderD3heatmap({
    if(any(input$Cor_Gene=="")){
      return()
    } else {
      # not working for a vector of one
      if (length(correl_genes())<2) {
        return(NULL)
      } else if (length(correl_genes())>500) {
        return(NULL)
      } else {
        myfont <- paste0(round(1/log10(length(correl_genes()))*10,0),"px")
        d3heatmap(my_exprs[correl_genes(),, drop=F],
                  dendrogram ="column",
                  colors = "YlOrRd",
                  show_grid = FALSE,
                  anim_duration = 0,
                  xaxis_font_size = "6px",
                  yaxis_font_size = myfont,
                  xaxis_height = 150,
                  yaxis_width = 150,
                  labCol = paste0("sample",1:ncol(my_exprs))
        )
      }
    }
  })

  
})
