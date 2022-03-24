#' @import knitr
#' @import kableExtra
#' @import RColorBrewer
#' @import ellipsis
#' @import tidyverse
#' @import tidyr
#' @import scales
#' @import ggpubr
#' @import stringr
#' @import lme4
#' @import rpart
#' @import data.table
#' @import vegan
#' @import factoextra
#' @import FactoMineR
#' @import dendextend
#' @import cluster
#' @import NbClust
#' @import fpc
#' @import ade4
#' @import BBmisc
#' @import clustree
#' @import plotly
#' @import networkD3
#' @import htmlwidgets
#' @import labdsv
#' @import ggfortify
#' @import corrplot
#' @import ggforce
#' @import concaveman
#' @import forcats
#' @import ggrepel
#' @import ggdendro
#' @import stats
#' @import rmarkdown
#' @import xfun
#' @import tinytex
#' @import withr
globalVariables(c("Dim.1" ,"Dim.2", "MED_stocks", "MeanDim1", "MeanDim2" ,"area", "as.dendrogram","catch_cluster",
                  "catch_cluster_stock", "catch_ship", "catch_ship_stock","clust_size", "cluster", "cluster_nr" ,"cmdscale",
                  "colorRampPalette", "cor", "cum_share",dist ,"hclust", "k", "landings" ,"landings_stock", "landkg" ,"major.area" ,"r",
                  "share_level","share_stock", "ship_ID", "size" ,"species", "stock", "total_landings", "weight_landed"))

# Section I - Functions #######
## 1) Transform base data ####

#' @title Transformation of basic catch data
#'
#' @description This function is used to transform stock-based fisheries catch data into a suitable format for the clustering procedure and cluster validation of the FleetSegmentation package.
#' It uses the share of the species on the total catch of each ship instead of the amount (or other seafood) of fish caught.
#' See example_catchdata for a example data frame format. Catch needs to be given in kilograms.
#' @param data The basic catch data frame.
#' @keywords transformation
#' @keywords data
#' @export catchdata_transformation
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
catchdata_transformation <- function(data){
  names(data)<-c("ship_ID","stock","landings")
  stocks_sum <- data %>%
    group_by(ship_ID, stock) %>%
    summarise(weight_landed = sum(landings)) %>%
    group_by(ship_ID)%>%
    mutate(share_stock= weight_landed/sum(weight_landed)) %>%
    ungroup()

  stocks_sum <- stocks_sum[,-ncol(stocks_sum)]
  stocks_wide <- spread(stocks_sum, key = "stock", value = weight_landed, fill = 0)
  stocks_wide <- as.data.frame(stocks_wide)
  rownames(stocks_wide) <- stocks_wide$ship_ID
  stocks_wide <- stocks_wide[,-c(1)]

  catchdata <- stocks_wide/rowSums(stocks_wide)
  if(ncol(catchdata) == 0){
    warning("Your transformed catchdata does not contain any columns. Are you sure, your original data included more than one stock and more than one vessel?")
  }
  return(catchdata)
}


## 2) Check number of clusters with table ####
#' @title Compute most important indices for best number of clusters in table format
#'
#' @description The fleet segmentation package uses the the average silhouettes, the Mantel test, the Davis-Bouldin index,
#' the SD-index and the Calinski-Harabasz index. This function gives the values of those indices for the given maximal number of clusters in a table, which can be printed in a basic or html format or stored as a data frame.
#' A modified (metric converted) Bray-Curtis distance matrix is computed from the input data, the clustering is performed as a hierarchical agglomerative clustering (HAC) using the average linkage link function.
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param max_clusternumber The maximum number of clusters to be expected. Defaults to 1 less than the number of ships in the catchdata-frame, up to a maximum of 15.
#' @param style The output style, defaults to `basic`, which prints a data frame in the console and can be stored. For a html-version, use `html`
#' @param distance The distance measure used. Defaults to modified (metric conversion) Bray-Curtis distance distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param method The link function used. Defaults to average linkage. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @keywords number of clusters
#' @keywords table
#' @export numberclust_table
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' numberclust_table(catchdata = catchdata,max_clusternumber = 15)
#' numberclust_table(catchdata = catchdata,max_clusternumber = 15,style = "html")
#' optclust_table <- numberclust_table(catchdata = catchdata,max_clusternumber = 15)
numberclust_table <- function(catchdata,max_clusternumber=ifelse(nrow(catchdata)<= 15,(nrow(catchdata)-1),15),style="basic", distance="jaccard", method="average") {
  # calculate distance matrix
  distmat <- vegdist(catchdata, method = distance)
  # perform clustering
  hc1_average <- hclust(distmat,method = method)
  # best number of clusters - silhouette
  silhouette <- fviz_nbclust(catchdata, hcut, method = "silhouette" ,k.max=max_clusternumber,linecolor = "#00AAAA") + theme_bw() + ggtitle("Average silhouettes")
  silhouette_df <- silhouette$data
  silhouette_df$clusters <- as.numeric(silhouette_df$clusters)
  silhouette_best <- as.numeric(which.max(silhouette_df$y))
  silhouette_index <- round(max(silhouette_df$y),3)
  # best number of clusters - mantel test
  #Mantel test
  # Function to compute a binary distance matrix from groups
  grpdist <- function(X)
  {
    gr <- as.data.frame(as.factor(X))
    distgr <- daisy(gr, "gower")
    distgr
  }


  kt <- data.frame(k=1:max_clusternumber, r=0)
  for (i in 2:max_clusternumber) {                                                 #for every row except first and last
    gr <- dendextend::cutree(hc1_average, i)                                    #cut the clustering into i number um clusters and name that gr
    distgr <- grpdist(gr)                                           #calculate the gower distance (pairwise dissimilaritys) of clustering with i number of clusters
    mt <- cor(distmat, distgr, method="pearson")          #Compare the result with the original distance matrix using pearsons correlation coefficient
    kt[i,2] <- mt                                                   #write the result into column2
  }
  mantel_best <- which.max(kt$r)
  mantel_index <- round(max(kt$r),3)
  # Davis-Bouldin index
  DB <- NbClust(data = catchdata, diss = distmat,distance = NULL, min.nc = 1,max.nc =
                  max_clusternumber, method = "average", index ="db")

  DB_Best <- DB$Best.nc[1]
  DB_Index <- round(DB$Best.nc[2],3)
  #SD index
  SD <- NbClust(data = catchdata,diss = distmat,distance = NULL, min.nc = 1,max.nc =
                  max_clusternumber, method = "average", index ="sdindex")

  SD_Best <- SD$Best.nc[1]
  SD_Index <- round(SD$Best.nc[2],3)
  #Pseudo F
  PsF <- NbClust(data = catchdata, diss = distmat,distance = NULL, min.nc = 1,max.nc =
                   max_clusternumber, method = "average", index ="ch")

  PsF_Best <- PsF$Best.nc[1]
  PsF_index <-round(PsF$Best.nc[2],3)

  #write table
  Indices <- c("Average silhouettes","Mantel test","Davis-Bouldin index","SD index","Calinski-Harabasz index")
  Opt_Number_clust <- c(silhouette_best,mantel_best,DB_Best,SD_Best,PsF_Best)
  Index_Values <- c(silhouette_index,mantel_index,DB_Index,SD_Index,PsF_index)

  table <- data.frame(Indices,Opt_Number_clust,Index_Values)
  names(table)<-c("Indices","Optimal number of clusters","Index value")
  if(style=="basic"){
    return(table)
  }
  if(style=="html"){
    (kbl(table,format = "html", align = "c",caption = "Optimal number of clusters",booktabs = T) %>%
       kable_styling(full_width = T,bootstrap_options = c("striped"),fixed_thead = T) %>%
       footnote(general = "The average silhouette score has proven to be most accurate in estimating the best number of clusters in fisheries catch data. The SD and Calinski-Harabasz index can take Inf values in case of few sample vessels.
                     If the scores of the indices are contradictory, visualise the indices via numberclust_plot() and trace the clustering with NumberClust_ClustTree()."))
  }
}

## 3) Check number of clusters with plot ####
#' @title Plot most important diagnostics and indices for best number of clusters
#'
#' @description This function uses a scree plot of the explained variation, the average silhouettes, the Mantel test and a dendrogram for estimating the best number of clusters for the clustering.
#' It creates a plotgrid with a scree plot (upper left), a plot of the average silhouette score including the maximum (upper right),
#' a plot of the Mantel test including the best result (lower left) and a dendrogram of the clustering (lower left). The dendrogram can either
#' be displayed without any further features or with the recommended cutting heights (setting dend_method to "range") or a cut at a defined linkage distance (setting dend_method to "cut").
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param max_clusternumber The maximum number of clusters to be expected. Defaults to 1 less than the number of ships in the catchdata-frame, up to a maximum of 15.
#' @param distance The distance measure used. Defaults to modified (metric-converted) Bray-Curtis distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param method The link function used. Defaults to average linkage. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param dend_method The style of the plotted dendrogram. "basic" returns a blank dendrogram, "range" a dendrogram with the recommended cutting heights depicted, "cut" enables the user to cut the dendrogram at a height of his choice. The resulting number of clusters will be shown in the plot.
#' @param dend_cut The height used to cut the dendrogram. Defaults to 0.75, which showed to be the appropriate cutting height for various fleet data sets.
#' @param range_min The lower border of the cutting range set for a dendrogram with the "range"- method
#' @param range_max The upper border of the cutting range set for a dendrogram with the "range"- method
#' @keywords number of clusters
#' @keywords plot
#' @keywords grid
#' @export numberclust_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' numberclust_plot(catchdata = catchdata,max_clusternumber = 15)
numberclust_plot <- function(catchdata,max_clusternumber=ifelse(nrow(catchdata)<= 15,(nrow(catchdata)-1),15), distance= "jaccard", method = "average",dend_method = "basic", dend_cut = 0.75, range_min = 0.4, range_max = 0.8) {
  # calculate distance matrix
  distmat <- vegdist(catchdata, method = distance)
  # perform clustering
  hc1_average <- hclust(distmat, method = method)
  # maximum number of clusters
  # elbow method plot
  elbow <- fviz_nbclust(catchdata, hcut, method = "wss", k.max = max_clusternumber,linecolor = "#00AAAA") + theme_bw() + ggtitle("Elbow method")
  # silhouette plot
  silhouette <- fviz_nbclust(catchdata, hcut, method = "silhouette" ,k.max=max_clusternumber,linecolor = "#00AAAA") + theme_bw() + ggtitle("Average silhouettes")
  #Mantel test
  # Function to compute a binary distance matrix from groups
  grpdist <- function(X)
  {
    gr <- as.data.frame(as.factor(X))
    distgr <- daisy(gr, "gower")
    distgr
  }


  kt <- data.frame(k=1:max_clusternumber, r=0)
  for (i in 2:max_clusternumber) {                                                 #for every row except first and last
    gr <- dendextend::cutree(hc1_average, i)                                    #cut the clustering into i number um clusters and name that gr
    distgr <- grpdist(gr)                                           #calculate the gower distance (pairwise dissimilaritys) of clustering with i number of clusters
    mt <- cor(distmat, distgr, method="pearson")          #Compare the result with the original distance matrix using pearsons correlation coefficient
    kt[i,2] <- mt                                                   #write the result into column2
  }
  k.best <- which.max(kt$r)
  mantel_plot <- ggplot(kt)+geom_point(aes(k,r),colour="#00AAAA")+
    geom_line(aes(k,r),colour="#00AAAA")+
    geom_vline(xintercept = k.best, linetype="dashed",colour="#00AAAA")+
    theme_bw()+
    labs(title = "Mantel test",y="Pearson's correlation",x="Number of clusters k")+
    scale_x_continuous(breaks = seq(1,max_clusternumber))

  #Dendrogramm#
  dend <- as.dendrogram(hc1_average)
  if(any(dend_cut > 1 | range_max > 1 | range_min > 1 | dend_cut < 0 | range_max < 0 | range_min < 0 )){
         stop("You have selected an invalid range limit or cutting height. Please make sure to select values between 0 and 1 for the arguments dend_cut, range_min, and range_max")
  }
  if(dend_method == "basic"){
    dend <- dend %>% set("branches_lwd", 0.7)
    ggdend <- as.ggdend(dend)
    dendro <- ggplot(ggdend,labels = F,theme = theme_minimal())+
      theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
      labs(y="Linkage distance")+
      scale_y_continuous(breaks = c(0,.25,.5,.75,1))
  }
  if(dend_method == "range"){
    dend <- dend %>% set("branches_lwd", 0.7)
    ggdend <- as.ggdend(dend)
    dend_clusters_min <- dendextend::cutree(dend,h = range_min)
    nr_cluster_min <- as.numeric(n_distinct(dend_clusters_min))
    dend_clusters_max <- dendextend::cutree(dend,h = range_max)
    nr_cluster_max <- as.numeric(n_distinct(dend_clusters_max))
    dend_data <- dendro_data(dend, type = "rectangle")
    dend_data_xpos_range <- max(dend_data$segments$x)*.7
    dend_data_ypos_range <- max(dend_data$segments$y)
    dend_label_range <- paste(nr_cluster_max,"to",nr_cluster_min,"clusters")
    dendro <- ggplot(ggdend,labels = F,theme = theme_minimal())+
      theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
      labs(y="Linkage distance")+
      scale_y_continuous(breaks = c(0,.25,.5,.75,1))+
      geom_hline(yintercept = (range_max - (range_max-range_min)/2), size=.8, color="#252525",linetype="dashed",alpha=.8)+
      geom_hline(yintercept = range_max, size=.8, color="#5f9ea0",linetype="dashed",alpha=.8)+
      geom_hline(yintercept = range_min, size=.8, color="#5f9ea0",linetype="dashed",alpha=.8)+
      geom_label(aes(dend_data_xpos_range,dend_data_ypos_range,label=dend_label_range),colour="black",fill="white",size=3,fontface="bold")
  }

  if(dend_method == "cut"){
    suppressWarnings(
      dend <- dend %>% set("branches_lwd", 0.7) %>%
        color_branches(h=as.numeric(dend_cut),col=c("#29505a","#03396c","#005b96","#6497b1","#b3cde0","#a3c1ad","#a0d6b4","#5f9ea0","#317873","#85d2b9","#5d9c94","#3b898b","#447a94","#3c6081")))

    ggdend <- as.ggdend(dend)
    dend_clusters <- dendextend::cutree(dend,h = as.numeric(dend_cut))
    nr_cluster <- as.numeric(n_distinct(dend_clusters))
    dend_data <- dendro_data(dend, type = "rectangle")
    dend_data_xpos <- max(dend_data$segments$x)*.7
    dend_data_ypos <- max(dend_data$segments$y)
    dend_label <- paste(nr_cluster,"clusters")
    dendro <-   suppressWarnings(
      ggplot(ggdend,labels = F,theme = theme_minimal())+
        theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
        labs(y="Linkage distance")+
        scale_y_continuous(breaks = c(0,.25,.5,.75,1))+
        geom_hline(yintercept = dend_cut, size=.8, color="#252525",linetype="dashed",alpha=.8)+
        geom_label(aes(dend_data_xpos,dend_data_ypos,label=dend_label),colour="black",fill="white",size=3,fontface="bold")
    )
  }


  # arrange plots
  ggarrange(elbow,silhouette,mantel_plot,dendro,nrow=2,ncol = 2, labels = c("A","B","C","D"), label.x = .9)
}

####  4) Plot clustering dendrogramm ####
#' @title Plot dendrogram to identify number of clusters
#'
#' @description This function plots a dendrogram of the clustering. The dendrogram can either be displayed without any further features
#' or with the recommended cutting heights (setting dend_method to "range") or a cut at a defined linkage distance (setting dend_method to "cut").
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param distance The distance measure used. Defaults to modified (metric-converted) Bray-Curtis distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param method The link function used. Defaults to average linkage. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param dend_method The style of the plotted dendrogram. "basic" returns a blank dendrogram, "range" a dendrogram with the recommended cutting heights depicted, "cut" enables the user to cut the dendrogram at a height of his choice. The resulting number of clusters will be shown in the plot.
#' @param dend_cut The height used to cut the dendrogram. Defaults to 0.75, which showed to be the appropriate cutting height for various fleet data sets.
#' @param range_min The lower border of the cutting range set for a dendrogram with the "range"- method
#' @param range_max The upper border of the cutting range set for a dendrogram with the "range"- method
#' @keywords number of clusters
#' @keywords plot
#' @keywords dendrogram
#' @export numberclust_dendrogram
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' numberclust_dendrogram(catchdata=catchdata) ### basic dendrogram
#' numberclust_dendrogram(catchdata=catchdata, dend_method= "range")
#' ### dendrogram with range of good cutting heights
#' numberclust_dendrogram(catchdata=catchdata, dend_method= "cut", dend_cut = 0.75)
#' ### dendrogram with cut at a linkage distance of 0.75 and the number of resulting clusters
numberclust_dendrogram <- function(catchdata, distance= "jaccard", method = "average",dend_method = "basic", dend_cut = 0.75, range_min = 0.4, range_max = 0.8){
  # calculate distance matrix
  distmat <- vegdist(catchdata, method = distance)
  # perform clustering
  hc1_average <- hclust(distmat, method = method)
  dend <- as.dendrogram(hc1_average)
  if(any(dend_cut > 1 | range_max > 1 | range_min > 1 | dend_cut < 0 | range_max < 0 | range_min < 0 )){
    stop("You have selected an invalid range limit or cutting height. Please make sure to select values between 0 and 1 for the arguments dend_cut, range_min, and range_max")
  }
  if(dend_method == "basic"){
    dend <- dend %>% set("branches_lwd", 0.7)
    ggdend <- as.ggdend(dend)
    dendro <- ggplot(ggdend,labels = F,theme = theme_minimal())+
      theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
      labs(y="Linkage distance")+
      scale_y_continuous(breaks = c(0,.25,.5,.75,1))
  }
  if(dend_method == "range"){
    dend <- dend %>% set("branches_lwd", 0.7)
    ggdend <- as.ggdend(dend)
    dend_clusters_min <- dendextend::cutree(dend,h = range_min)
    nr_cluster_min <- as.numeric(n_distinct(dend_clusters_min))
    dend_clusters_max <- dendextend::cutree(dend,h = range_max)
    nr_cluster_max <- as.numeric(n_distinct(dend_clusters_max))
    dend_data <- dendro_data(dend, type = "rectangle")
    dend_data_xpos_range <- max(dend_data$segments$x)*.7
    dend_data_ypos_range <- max(dend_data$segments$y)
    dend_label_range <- paste(nr_cluster_max,"to",nr_cluster_min,"clusters")
    dendro <- ggplot(ggdend,labels = F,theme = theme_minimal())+
      theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
      labs(y="Linkage distance")+
      scale_y_continuous(breaks = c(0,.25,.5,.75,1))+
      geom_hline(yintercept = (range_max - (range_max-range_min)/2), size=.8, color="#252525",linetype="dashed",alpha=.8)+
      geom_hline(yintercept = range_max, size=.8, color="#5f9ea0",linetype="dashed",alpha=.8)+
      geom_hline(yintercept = range_min, size=.8, color="#5f9ea0",linetype="dashed",alpha=.8)+
      geom_label(aes(dend_data_xpos_range,dend_data_ypos_range,label=dend_label_range),colour="black",fill="white",size=3,fontface="bold")
  }
  if(dend_method == "cut"){
    suppressWarnings(
      dend <- dend %>% set("branches_lwd", 0.7) %>%
        color_branches(h=as.numeric(dend_cut),col=c("#29505a","#03396c","#005b96","#6497b1","#b3cde0","#a3c1ad","#a0d6b4","#5f9ea0","#317873","#85d2b9","#5d9c94","#3b898b","#447a94","#3c6081")))

    ggdend <- as.ggdend(dend)
    dend_clusters <- dendextend::cutree(dend,h = as.numeric(dend_cut))
    nr_cluster <- as.numeric(n_distinct(dend_clusters))
    dend_data <- dendro_data(dend, type = "rectangle")
    dend_data_xpos <- max(dend_data$segments$x)*.7
    dend_data_ypos <- max(dend_data$segments$y)
    dend_label <- paste(nr_cluster,"clusters")
    dendro <- ggplot(ggdend,labels = F,theme = theme_minimal())+
      theme(axis.text.x = element_blank(),axis.title.x = element_blank())+
      labs(y="Linkage distance")+
      scale_y_continuous(breaks = c(0,.25,.5,.75,1))+
      geom_hline(yintercept = dend_cut, size=.8, color="#252525",linetype="dashed",alpha=.8)+
      geom_label(aes(dend_data_xpos,dend_data_ypos,label=dend_label),colour="black",fill="white",size=3,fontface="bold")
  }
  suppressWarnings(dendro)
}

####  5) Visualize clustering tree ####
#' @title Visualize clustering process with clustering tree
#'
#' @description This function creates a clustering tree with the clustree()-function from the eponymous package. It visualizes the clustering process by showing the splits
#' of the clusters in a tree plot. This a very useful method for identifying major segmentations of big groups in the data and ultimately deciding on how
#' many clusters to use.
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param max_clusternumber The maximum number of clusters to be expected. Defaults to 1 less than the number of ships in the catchdata-frame, up to a maximum of 15.
#' @param distance The distance measure used. Defaults to modified (metric-converted) Bray-Curtis distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param method The link function used. Defaults to average linkage. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @keywords number of clusters
#' @keywords clustree
#' @export numberclust_clustree
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' numberclust_clustree(catchdata = catchdata,max_clusternumber = 15)
numberclust_clustree <- function(catchdata,max_clusternumber=ifelse(nrow(catchdata)<= 15,(nrow(catchdata)-1),15), distance = "jaccard", method = "average") {
  # calculate distance matrix
  distmat <- vegdist(catchdata, method = distance)
  # perform clustering
  hc1_average <- hclust(distmat, method = method)
  tmp2 <- data.frame(dendextend::cutree(hc1_average, k=c(1:max_clusternumber)))
  colnames(tmp2) <- gsub("X","K",colnames(tmp2))
  tmp2$ship_ID <- rownames(tmp2)
  # Define the number of colors you want
  mycolors <- rev(colorRampPalette(brewer.pal(8, "YlGnBu"))(max_clusternumber))
  suppressWarnings(
    clustree(x = tmp2, prefix = "K", node_size_range = c(3,20), node_label = "size",node_label_nudge=-0.3,show_axis=T)+
      scale_colour_manual(values = mycolors)+theme(legend.position = "none")+
      scale_fill_manual(values = rep("lightblue",max_clusternumber))+
      labs(y="Number of clusters k")
  )
}

####  6) Perform clustering #####
#' @title Perform clustering
#'
#' @description This is the core function to perform the clustering of the catchdata. Use the number of suitable number of clusters estimated with numberclust_table() and numberclust_plot().
#' A modified (metric-converted) Bray-Curtis distance matrix is computed from the input data, the clustering is performed as a hierarchical agglomerative clustering (HAC) using the average linkage link function.
#' The function creates a new data frame, which can be printed or stored.
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param n_cluster The number of clusters to be generated. No default.
#' @param distance The distance measure used. Defaults to modified (metric-converted) Bray-Curtis distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @param method The link function used. Defaults to average linkage. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @keywords clustering
#' @export segmentation_clustering
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
segmentation_clustering <- function(catchdata,n_cluster, distance = "jaccard",method = "average") {
  # calculate distance matrix
  distmat <- vegdist(catchdata, method = distance)
  # perform clustering
  hc1_average <- hclust(distmat, method = method)
  # clustering
  catchdata$ship_ID <- rownames(catchdata)
  catchdata$cluster <- dendextend::cutree(hc1_average, k=n_cluster)
  catchdata <- catchdata[,c((ncol(catchdata)-1),ncol(catchdata))]
  catchdata$cluster <- sub("^", "cluster ", catchdata$cluster)
  catchdata$cluster <- as.factor(catchdata$cluster)
  rownames(catchdata) <- seq(1, length(catchdata$ship_ID))
  clustering <- data.frame(catchdata)
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels),ordered = T)
  return(clustering)
}

#### 7) Tabelize stock shares of clusters ####
#' @title Tabelize stock shares of clusters
#'
#' @description This function creates a table with the average shares of the stocks on the catch of each clusters vessel. The result is a data frame, which can be printed or stored as an object.
#' Alternatively, an html-table can be created by setting style = "html".
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param style The style, in which the result is printed. Defaults to "basic". "html" produces an html-table.
#' @keywords clustering
#' @keywords table
#' @keywords stockshares
#' @export clustering_stockshares_table
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' clustering_stockshares_table(data = stockdata,clustering = clustering)
#' clustering_stockshares_table(data = stockdata,clustering = clustering,style = "html")
#' segments_stockshares <- clustering_stockshares_table(data = stockdata,clustering = clustering)
clustering_stockshares_table <- function(data,clustering, style="basic"){
  names(data) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  #Calculate spec shares
  data <- data %>%
    group_by(ship_ID) %>%
    left_join(clustering,by="ship_ID") %>%
    group_by(cluster,stock) %>%
    mutate(catch_cluster_stock =sum(landings)) %>%
    group_by(cluster) %>%
    mutate(catch_cluster =sum(landings)) %>%
    group_by(cluster,stock) %>%
    summarise(share_stock = catch_cluster_stock/catch_cluster*100) %>%
    unique()%>%
    arrange(cluster,desc(share_stock)) %>%
    ungroup()

  if(style=="basic"){
    return(data)
  }
  if(style=="html"){
    data$cluster <- as.character(data$cluster)
    data <- dplyr::filter(data, share_stock >= 5)
    data$share_stock <- round(data$share_stock, digits = 2)
    print(kbl(data[,2:3],format = "html", align = "c",caption = "Catch composition of clusters",
              col.names=c("ICES stock", "Share of stock on total catch [%]")) %>%
            pack_rows(index = table(fct_inorder(data$cluster))) %>%
            kable_styling(full_width = T,bootstrap_options = c("striped"),fixed_thead = T) %>%
            footnote(general = "Only catch shares larger than 5% are depicted."))

  }
}

##### 8) Plot stock shares of clusters ####
#' @title Plot stock shares of clusters
#'
#' @description This is function creates an overview barplot of the average shares of stocks on the catch of each clusters vessels.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param min_share The minimum average percentage share a stock has to have to be labelled in the plot. Defaults to 5\%.
#' @param label_wrap Indicates the number of characters per line before a line break in the stock labels. Defaults to 6.
#' @param display_cluster_size Indicates, whether the number of vessels in each cluster should be displayed in the plot. Defaults to FALSE.
#' @keywords clustering
#' @keywords plot
#' @keywords stockshares
#' @export clustering_stockshares_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' clustering_stockshares_plot(data = stockdata,clustering = clustering)
#' clustering_stockshares_plot(data = stockdata,clustering = clustering,
#' min_share=10,label_wrap=10, display_cluster_size=TRUE)
clustering_stockshares_plot <- function(data,clustering, min_share=5,label_wrap=6, display_cluster_size=F){
  names(data) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  #Calculate spec shares
  data <- data %>%
    group_by(ship_ID) %>%
    left_join(clustering,by="ship_ID") %>%
    group_by(cluster,stock) %>%
    mutate(catch_cluster_stock =sum(landings)) %>%
    group_by(cluster) %>%
    mutate(catch_cluster =sum(landings)) %>%
    group_by(cluster,stock) %>%
    summarise(share_stock = catch_cluster_stock/catch_cluster) %>%
    unique()%>%
    ungroup() %>%
    arrange(cluster,desc(share_stock))

  data$label <- data$stock
  data$label <- as.character(data$label)
  data$label[data$share_stock < (min_share/100)] <- " "
  data$label <- as.factor(data$label)

  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  data$cluster <- factor(data$cluster, levels = rev(clusterlevels))
  data <- data %>%
    arrange(share_stock,cluster)
  data$share_level <- "minimal"
  data$share_level[data$share_stock > .1] <- "low"
  data$share_level[data$share_stock > .25] <- "medium"
  data$share_level[data$share_stock > .5] <- "high"
  data$share_level[data$share_stock > .75] <- "very high"
  data$share_level <- factor(data$share_level, levels = c("very high","high","medium","low","minimal"))
  data <- data %>%
    group_by(cluster) %>%
    mutate(cum_share=cumsum(share_stock)) %>%
    ungroup()%>%
    arrange(cluster,share_level,share_stock)
  data$label <- gsub("_"," ",data$label)

  n_vessel <- clustering %>% group_by(cluster) %>% summarise(n_vessel = n_distinct(ship_ID))
  n_vessel$cluster <- factor(n_vessel$cluster, levels = rev(clusterlevels))

  if(display_cluster_size==F){
    stock_plot <- ggplot()+
      geom_col(data=data, aes(cluster, share_stock,fill=share_level),position = "fill", colour="black")+
      scale_fill_manual(values = rev(c("#e7f3f6","#c4d9df","#62919c","#29505a","#04172d")),labels=c("very high (> 75%)","high (50-75%)","medium (25-50%)","low (10-25%)","minimal (< 10%)"),guide = guide_legend(reverse=TRUE))+
      geom_label(data=data[data$share_stock > (min_share/100),], aes(cluster, (cum_share-0.5*share_stock),label=str_wrap(label,width = label_wrap)),size=4, color="black",fill="white")+
      theme_bw() +
      theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold",size=12),
            legend.position = "bottom",legend.title = element_blank(),
            plot.margin=unit(c(5.5,30,5.5,5.5),"points"), legend.text = element_text(size = 10),axis.text.y = element_text(size=8))+
      scale_y_continuous(labels=percent, expand = c(0, 0),breaks = c(.25,.50,.75,1),limits = c(0,1))+
      scale_x_discrete(labels=rev(seq(1,clust_number,1)))+
      coord_flip()
  }
  if(display_cluster_size==T){
    stock_plot <- ggplot()+
      geom_col(data=data, aes(cluster, share_stock,fill=share_level),position = "fill", colour="black")+
      scale_fill_manual(values = rev(c("#e7f3f6","#c4d9df","#62919c","#29505a","#04172d")),labels=c("very high (> 75%)","high (50-75%)","medium (25-50%)","low (10-25%)","minimal (< 10%)"),guide = guide_legend(reverse=TRUE))+
      geom_label(data=data[data$share_stock > (min_share/100),], aes(cluster, (cum_share-0.5*share_stock),label=str_wrap(label,width = label_wrap)),size=4, color="black",fill="white")+
      geom_text(data=n_vessel, aes(cluster,1.1,label=n_vessel),size=4,show.legend = F)+
      theme_minimal() +
      theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold",size=12),
            legend.position = "bottom",legend.title = element_blank(),
            plot.margin=unit(c(5.5,30,5.5,5.5),"points"), legend.text = element_text(size = 10),
            axis.text.y = element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_y_continuous(labels=percent, expand = c(0, 0),breaks = c(.25,.50,.75,1),limits = c(0,1.2))+
      scale_x_discrete(labels=rev(seq(1,clust_number,1)))+
      coord_flip()
    }
  stock_plot
}


##### 9) Plot assemblage shares of clusters ####
#' @title Plot assemblage shares of clusters
#'
#' @description This is function creates an overview barplot of the average shares of stocks on the catch of each clusters vessels.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param min_share The minimum average percentage share a stock has to have to be labelled in the plot. Defaults to 5\%.
#' @param label_wrap Indicates the number of characters per line before a line break in the stock labels. Defaults to 6.
#' @param display_cluster_size Indicates, whether the number of vessels in each cluster should be displayed in the plot. Defaults to FALSE.
#' @keywords clustering
#' @keywords plot
#' @keywords assemblage
#' @export clustering_assemblageshares_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' clustering_assemblageshares_plot(data = stockdata,clustering = clustering)
#' clustering_assemblageshares_plot(data = stockdata,clustering = clustering,
#' min_share=10,label_wrap=10, display_cluster_size=TRUE)
clustering_assemblageshares_plot <- function(data,clustering, min_share=5,label_wrap=6, display_cluster_size=F){
  names(data) <- c("ship_ID", "stock", "landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels, levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))

  assemblage_red <- assemblage %>%
    dplyr::select(species_code,target_assemblage_code,target_assemblage)
suppressMessages(
  data <- data %>%
    mutate(species_code = toupper(sub("\\..*", "", stock))) %>%
    mutate(species_code = toupper(sub("\\-.*", "", species_code))) %>%
    left_join(assemblage_red) %>%
    mutate(target_assemblage_code = ifelse(species_code == "NEP", "CRU",target_assemblage_code),
           target_assemblage = ifelse(species_code == "NEP", "Crustaceans",target_assemblage),
           target_assemblage_code = replace_na(target_assemblage_code,"UNK"),
           target_assemblage = replace_na(target_assemblage, "unknown")) %>%
    left_join(clustering) %>%
    group_by(cluster,target_assemblage_code) %>%
    summarise(catch_cluster = sum(landings)) %>%
    group_by(cluster) %>%
    mutate(share_assemblage = catch_cluster/sum(catch_cluster)) %>%
    unique() %>% ungroup() %>% arrange(cluster, desc(share_assemblage))
  )

  data$label <- data$target_assemblage_code
  data$label <- as.character(data$label)
  data$label[data$share_assemblage < (min_share/100)] <- " "
  data$label <- as.factor(data$label)
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels, levels)
  }
  data$cluster <- factor(data$cluster, levels = rev(clusterlevels))
  data <- data %>% arrange(share_assemblage, cluster)

  data$share_level <- "minimal"
  data$share_level[data$share_assemblage > 0.1] <- "low"
  data$share_level[data$share_assemblage > 0.25] <- "medium"
  data$share_level[data$share_assemblage > 0.5] <- "high"
  data$share_level[data$share_assemblage > 0.75] <- "very high"
  data$share_level <- factor(data$share_level, levels = c("very high",
                                                          "high", "medium", "low", "minimal"))
  data <- data %>% group_by(cluster) %>% mutate(cum_share = cumsum(share_assemblage)) %>%
    ungroup() %>% arrange(cluster, share_level, share_assemblage)
  data$label <- gsub("_", " ", data$label)
  n_vessel <- clustering %>% group_by(cluster) %>% summarise(n_vessel = n_distinct(ship_ID))
  n_vessel$cluster <- factor(n_vessel$cluster, levels = rev(clusterlevels))

  n_vessel <- clustering %>% group_by(cluster) %>% summarise(n_vessel = n_distinct(ship_ID))
  n_vessel$cluster <- factor(n_vessel$cluster, levels = rev(clusterlevels))

  if(display_cluster_size==F){
    assemblage_plot <- ggplot() + geom_col(data = data, aes(cluster,
                                                       share_assemblage, fill = share_level), position = "fill",
                                      colour = "black") + scale_fill_manual(values = rev(c("#ef745c", "#c15955",
                                                                                           "#923e4d", "#632345", "#34073d")), labels = c("very high (> 75%)",
                                                                                                                                         "high (50-75%)", "medium (25-50%)", "low (10-25%)",
                                                                                                                                         "minimal (< 10%)"), guide = guide_legend(reverse = TRUE)) +
      geom_label(data = data[data$share_assemblage > (min_share/100),], aes(cluster, (cum_share - 0.5 * share_assemblage),label = str_wrap(label, width = label_wrap)),
      size = 4, color = "black", fill = "white") +
      theme_bw() + theme(axis.title.x = element_blank(),
                         axis.title.y = element_text(face = "bold", size = 12),
                         legend.position = "bottom", legend.title = element_blank(),
                         plot.margin = unit(c(5.5, 30, 5.5, 5.5), "points"),
                         legend.text = element_text(size = 10), axis.text.y = element_text(size = 8)) +
      scale_y_continuous(labels = percent, expand = c(0,
                                                      0), breaks = c(0.25, 0.5, 0.75, 1), limits = c(0,
                                                                                                     1)) +
      scale_x_discrete(labels = rev(seq(1, clust_number,1))) +
      coord_flip()
  }
  if(display_cluster_size==T){
    assemblage_plot <- ggplot() + geom_col(data = data, aes(cluster,
                                                            share_assemblage, fill = share_level), position = "fill",
                                           colour = "black") +
      scale_fill_manual(values = rev(c("#ef745c", "#c15955","#923e4d", "#632345", "#34073d")), labels = c("very high (> 75%)","high (50-75%)", "medium (25-50%)", "low (10-25%)",                                                                                                                              "minimal (< 10%)"), guide = guide_legend(reverse = TRUE)) +
      geom_label(data = data[data$share_assemblage > (min_share/100),], aes(cluster, (cum_share - 0.5 * share_assemblage),
             label = str_wrap(label, width = label_wrap)),
      size = 4, color = "black", fill = "white") +
      geom_text(data=n_vessel, aes(cluster,1.1,label=n_vessel),size=4,show.legend = F)+
      theme_minimal() +
      theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold",size=12),
            legend.position = "bottom",legend.title = element_blank(),
            plot.margin=unit(c(5.5,30,5.5,5.5),"points"), legend.text = element_text(size = 10),
            axis.text.y = element_text(size=8),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      scale_y_continuous(labels=percent, expand = c(0, 0),breaks = c(.25,.50,.75,1),limits = c(0,1.2))+
      scale_x_discrete(labels=rev(seq(1,clust_number,1)))+
      coord_flip()

    }
  assemblage_plot

}



##### 10) Plot stock shares of single cluster ####
#' @title Plot stock shares of clusters
#'
#' @description This is function creates an overview barplot of the average shares of stocks on the catch of all the vessels in one cluster.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param min_share The minimum average percentage share a stock has to have to be labelled in the plot. Defaults to 5\%.
#' @param cluster.number The number of the cluster of which vessel stock shares should be displayed.
#' @param label_wrap Indicates the number of characters per line before a line break in the stock labels. Defaults to 15.
#' @keywords clustering
#' @keywords plot
#' @keywords stockshares
#' @export single_cluster_stockshares
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' single_cluster_stockshares(data = stockdata,clustering = clustering, cluster.number=1)
single_cluster_stockshares <- function(data,clustering, min_share=5,cluster.number,label_wrap=15){
  names(data) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  #Calculate spec shares
  data <- data %>%
    group_by(ship_ID) %>%
    left_join(clustering,by="ship_ID") %>%
    dplyr::filter(cluster == levels(cluster)[cluster.number]) %>%
    group_by(ship_ID,stock) %>%
    mutate(catch_ship_stock =sum(landings)) %>%
    group_by(ship_ID) %>%
    mutate(catch_ship =sum(landings)) %>%
    group_by(ship_ID,stock) %>%
    summarise(share_stock = catch_ship_stock/catch_ship) %>%
    unique()%>%
    ungroup() %>%
    arrange(ship_ID,desc(share_stock))

  data$label <- data$stock
  data$label <- as.character(data$label)
  data$label[data$share_stock < (min_share/100)] <- " "
  data$label <- as.factor(data$label)

  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))

  data <- data %>%
    arrange(share_stock,ship_ID)
  data$share_level <- "minimal"
  data$share_level[data$share_stock > .1] <- "low"
  data$share_level[data$share_stock > .25] <- "medium"
  data$share_level[data$share_stock > .5] <- "high"
  data$share_level[data$share_stock > .75] <- "very high"
  data$share_level <- factor(data$share_level, levels = c("very high","high","medium","low","minimal"))
  data <- data %>%
    group_by(ship_ID) %>%
    mutate(cum_share=cumsum(share_stock)) %>%
    ungroup()%>%
    arrange(ship_ID,share_level,share_stock)
  data$label <- gsub("_"," ",data$label)

  stock_plot <- ggplot()+
    geom_col(data=data, aes(ship_ID, share_stock,fill=share_level),position = "fill", colour="black")+
    scale_fill_manual(values = rev(c("#e7f3f6","#c4d9df","#62919c","#29505a","#04172d")),labels=c("very high (> 75%)","high (50-75%)","medium (25-50%)","low (10-25%)","minimal (< 10%)"),guide = guide_legend(reverse=TRUE))+
    geom_label(data=data[data$share_stock > (min_share/100),], aes(ship_ID, (cum_share-0.5*share_stock),label=str_wrap(label,width = label_wrap)),size=4, color="black",fill="white")+
    theme_bw() +
    theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold",size=12),
          legend.position = "bottom",legend.title = element_blank(),
          plot.margin=unit(c(5.5,30,5.5,5.5),"points"), legend.text = element_text(size = 10),axis.text.y = element_text(size=8))+
    scale_y_continuous(labels=percent, expand = c(0, 0),breaks = c(.25,.50,.75,1),limits = c(0,1))+
    coord_flip()
  stock_plot
}


##### 11) Plot number of ships in clusters ####
#' @title Plot number of ships in clusters
#'
#' @description This is function creates an overview plot of the number of ships in each cluster.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @keywords clustering
#' @keywords plot
#' @keywords size
#' @export cluster_size_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' cluster_size_plot(clustering = clustering)
cluster_size_plot <- function(clustering){
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  #plot size of clusters
  clust_size_red <- clustering %>%
    group_by(cluster)%>%
    summarise(size = n_distinct(ship_ID))
  clust_size_plot_ymax <- max(clust_size_red$size)+max(clust_size_red$size)*0.15
  clust_number <- n_distinct(clust_size_red$cluster)
  clust_size_red$cluster <- factor(clust_size_red$cluster, levels = c(clusterlevels))

  cluster_size_plot <- ggplot(clust_size_red, aes(cluster,size))+
    geom_col(colour="black",fill="#04172d",alpha=0.8)+
    geom_text(aes(label=size), vjust=-0.5, size=4, fontface="bold")+
    labs(y="Number of ships in cluster",x="cluster")+
    scale_x_discrete(labels=seq(1,clust_number,1))+
    theme_bw()+
    theme(axis.title.y =  element_text(face="bold",size=12),axis.title.x = element_blank())+
    scale_y_continuous(limits = c(0,clust_size_plot_ymax))
  cluster_size_plot
}

##### 12) Plot length of ships in clusters ####
#' @title Plot length of ships in clusters
#'
#' @description This is function creates an overview mixed dot- and boxplot of the length of the ships in each cluster. The length can be given in cm or m, the function will auto-transform from cm to m.
#' A boxplot will only be drawn for clusters containing more than 5 ships.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param shiplength A data frame containing the length of the ships clustered.
#' @keywords clustering
#' @keywords plot
#' @keywords shiplength
#' @export shiplength_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' shiplength_plot(clustering = clustering,shiplength = example_lengthdata)
shiplength_plot <- function(clustering,shiplength){
  names(shiplength)<-c("ship_ID","loa")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  cluster_length <- left_join(clustering,shiplength,by="ship_ID")
  cluster_length$unit <- ifelse(mean(cluster_length$loa >= 500),"cm","m")
  cluster_length$length <- ifelse(cluster_length$unit=="cm",cluster_length$loa/100,cluster_length$loa)

  cluster_length <- cluster_length %>%
    group_by(cluster)%>%
    mutate(clust_size= n_distinct(ship_ID))
  cluster_length$cluster <- factor(cluster_length$cluster, levels = c(clusterlevels))
  cluster_length <- cluster_length %>% arrange(cluster)

  clust_length_plot <- ggplot(data = cluster_length, aes(cluster,length))+
    geom_point(alpha=0)+
    geom_boxplot(data= dplyr::filter(cluster_length, clust_size > 5), aes(cluster, length),colour="black",fill="#29505a",alpha=.8)+
    geom_point(data= dplyr::filter(cluster_length, clust_size <= 5), size=5,aes(cluster, length),colour="black",fill="#29505a",alpha=.8,shape=21)+
    labs(y="length [m]",x="cluster")+
    theme_bw()+
    theme(axis.title =  element_text(face="bold",size=12))+
    scale_x_discrete(labels = c(1:clust_number))+
    scale_y_continuous(limits = c(0,(max(cluster_length$length)*1.2)))
  clust_length_plot
}


#### 13) Plot catch of single ships in clusters ####
#' @title Plot of catch of single ships in clusters
#'
#' @description This is function creates an overview mixed box- and dotplot of the catch of single ships in each cluster. A boxplot will only be drawn for clusters containing more than 5 ships.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @keywords clustering
#' @keywords plot
#' @keywords catch
#' @export singleship_catch_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' singleship_catch_plot(data = stockdata,clustering = clustering)
singleship_catch_plot <- function(data, clustering){
  names(data) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  singleships_tons <- data %>%
    left_join(clustering,by="ship_ID") %>%
    group_by(ship_ID,cluster) %>%
    summarise(total_landings = sum(landings)) %>%
    group_by(cluster)%>%
    mutate(clust_size=n_distinct(ship_ID))%>%
    ungroup()

  clust_number <- n_distinct(singleships_tons$cluster)
  singleships_tons$cluster <- factor(singleships_tons$cluster, levels = c(clusterlevels))
  singleships_tons <- singleships_tons %>% arrange(cluster)

  singleships_tons_plot <- ggplot(singleships_tons, aes(cluster, total_landings/1000))+
    geom_point(alpha=0)+
    geom_boxplot(data= dplyr::filter(singleships_tons,clust_size > 5), aes(cluster, total_landings/1000),fill="#62919c",alpha=.8)+
    geom_point(data= dplyr::filter(singleships_tons,clust_size <= 5), aes(cluster, total_landings/1000),shape=21, size=3,fill="#62919c",alpha=.8)+
    labs(y="Annual catch / ship [t]",x="cluster")+
    scale_x_discrete(labels = c(1:clust_number))+
    theme_bw()+
    theme(axis.title =  element_text(face="bold",size=12))+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE), limits = c(0,(max(singleships_tons$total_landings)/1000*1.2)))

  singleships_tons_plot
}

#### 14) Catch of all ships in cluster ####
#' @title Plot of total catch of ships in clusters
#'
#' @description This is function creates an overview barplot of the total catch of all ships in each cluster.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @keywords clustering
#' @keywords plot
#' @keywords catch
#' @export clustercatch_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' clustercatch_plot(data = stockdata,clustering = clustering)
clustercatch_plot <- function(data, clustering){
  names(data) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  cluster_tons <- data %>%
    left_join(clustering,by="ship_ID") %>%
    group_by(cluster) %>%
    summarise(total_landings = sum(landings)) %>%
    ungroup()

  clust_number <- n_distinct(cluster_tons$cluster)
  cluster_tons$cluster <- factor(cluster_tons$cluster, levels = c(clusterlevels))
  cluster_tons <- cluster_tons %>% arrange(cluster)

  cluster_tons_plot <- ggplot(cluster_tons, aes(cluster, total_landings/1000))+
    geom_col(fill="#62919c",alpha=.8,colour="black")+
    labs(y="Annual catch / cluster [t]",x="cluster")+
    scale_x_discrete(labels = c(1:clust_number))+
    theme_bw()+
    theme(axis.title =  element_text(face="bold",size=12))+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))
  cluster_tons_plot
}

##### 15) Grid of shiplength and catch plots ####
#' @title Plotgrid of number of ships, ship length, catch of single ships and total catch of ships in clusters
#'
#' @description This is function creates an overview plot grid of
#' 1) A barplot of the number of ships in each cluster
#' 2) A mixed dot- and boxplot of the length of ships in the clusters
#' 3) A mixed dot- and boxplot of the catch of single ships in the clusters.
#' 4) A barplot of the total catch of all ships in each cluster.
#' Boxplots will only be drawn for clusters containing more than 5 ships.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param shiplength A data frame containing the length of the ships clustered.
#' @keywords clustering
#' @keywords plot
#' @keywords grid
#' @export clustering_plotgrid
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' clustering_plotgrid(data = stockdata,clustering = clustering,shiplength =example_lengthdata)
clustering_plotgrid <- function(data,clustering,shiplength){
  dataframe <- data
  colnames(shiplength)<-c("ship_ID","loa")
  colnames(dataframe) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = (clusterlevels))
  #plot size of clusters
  clust_size_red <- clustering %>%
    group_by(cluster)%>%
    summarise(size = n_distinct(ship_ID))
  clust_size_plot_ymax <- max(clust_size_red$size)+max(clust_size_red$size)*0.15
  clust_number <- n_distinct(clust_size_red$cluster)
  clust_size_red$cluster <-factor(clust_size_red$cluster, levels = (clusterlevels))
  cluster_size_plot <- ggplot(clust_size_red, aes(cluster,size))+
    geom_col(colour="black",fill="#04172d",alpha=0.8)+
    geom_text(aes(label=size), vjust=-0.5, size=4, fontface="bold")+
    labs(y="Number of ships in cluster",x="cluster", title = " ")+
    scale_x_discrete(labels=seq(1,clust_number,1))+
    theme_bw()+
    theme(axis.title.y =  element_text(face="bold",size=10),axis.title.x = element_blank())+
    scale_y_continuous(limits = c(0,clust_size_plot_ymax))
  #plot length of ships in clusters
  cluster_length <- left_join(clustering,shiplength, by="ship_ID") %>%
    group_by(cluster) %>%
    mutate(clust_size = n_distinct(ship_ID)) %>%
    ungroup()
  cluster_length$unit <- ifelse(mean(cluster_length$loa >= 500),"cm","m")
  cluster_length$length <- ifelse(cluster_length$unit=="cm",cluster_length$loa/100,cluster_length$loa)
  clust_number <- n_distinct(cluster_length$cluster)
  cluster_length$cluster <- factor(cluster_length$cluster, levels = c(clusterlevels),ordered = T)
  cluster_length <- cluster_length %>% arrange(cluster)

  cluster_length_plot <- ggplot(data = cluster_length, aes(cluster,length))+
    geom_point(alpha=0)+
    geom_boxplot(data= dplyr::filter(cluster_length, clust_size > 5), aes(cluster, length),colour="black",fill="#29505a",alpha=.8)+
    geom_point(data= dplyr::filter(cluster_length, clust_size <= 5), size=3,aes(cluster, length),colour="black",fill="#29505a",alpha=.8,shape=21)+
    labs(y="length [m]",x="cluster", title = " ")+
    theme_bw()+
    theme(axis.title.y =  element_text(face="bold",size=10),axis.title.x = element_blank())+
    scale_x_discrete(labels = c(1:clust_number))+
    scale_y_continuous(limits = c(0,(max(cluster_length$length)*1.2)))
  #plot catch of single ships in clusters
  singleships_tons <- dataframe %>%
    left_join(clustering,by="ship_ID") %>%
    group_by(ship_ID,cluster) %>%
    summarise(total_landings = sum(landings)) %>%
    group_by(cluster)%>%
    mutate(clust_size=n_distinct(ship_ID))%>%
    ungroup()

  clust_number <- n_distinct(singleships_tons$cluster)
  singleships_tons$cluster <- factor(singleships_tons$cluster, levels = c(clusterlevels),ordered = T)
  singleships_tons <- singleships_tons %>% arrange(cluster)

  singleships_tons_plot <- ggplot(singleships_tons, aes(cluster, total_landings/1000))+
    geom_point(alpha=0)+
    geom_boxplot(data= dplyr::filter(singleships_tons,clust_size > 5), aes(cluster, total_landings/1000),fill="#62919c",alpha=.8)+
    geom_point(data= dplyr::filter(singleships_tons,clust_size <= 5), aes(cluster, total_landings/1000),shape=21, size=3,fill="#62919c",alpha=.8)+
    labs(y="Annual catch / ship [t]",x="cluster", title = " ")+
    scale_x_discrete(labels = c(1:clust_number))+
    theme_bw()+
    theme(axis.title =  element_text(face="bold",size=10))+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE), limits = c(0,(max(singleships_tons$total_landings)/1000*1.2)))
  #Plot catch of cluster
  cluster_tons <- dataframe %>%
    left_join(clustering,by="ship_ID") %>%
    group_by(cluster) %>%
    summarise(total_landings = sum(landings)) %>%
    ungroup()

  clust_number <- n_distinct(cluster_tons$cluster)
  cluster_tons$cluster <- factor(cluster_tons$cluster, levels = c(clusterlevels),ordered = T)
  cluster_tons <- cluster_tons %>% arrange(cluster)
  cluster_tons_plot <- ggplot(cluster_tons, aes(cluster, total_landings/1000))+
    geom_col(fill="#62919c",alpha=.8,colour="black")+
    labs(y="Annual catch / cluster [t]",x="cluster", title = " ")+
    scale_x_discrete(labels = c(1:clust_number))+
    theme_bw()+
    theme(axis.title =  element_text(face="bold",size=10))+
    scale_y_continuous(labels=function(x) format(x, big.mark = ",", decimal.mark = ".", scientific = FALSE))
  cluster_tons_plot
  # grid
  ggarrange(cluster_size_plot,cluster_length_plot,singleships_tons_plot,cluster_tons_plot,nrow = 2,ncol = 2,align = "hv", labels = c("A","B","C","D"), label.x = .9)
}


#### 16) Tabelize HHI of catch ####
#' @title Table of HHI of  overall catch of clusters.
#'
#' @description This function creates an overview table of the Herfindahl-Hirschmann Index (HHI) of the catch clusters. The HHI is equivalent to the Simpson Index,
#' it is a measure of degree of concentration/consolidation. In case of fisheries catch data, it indicates, whether a fishery is of rather targeted or mixed nature.
#' The HHI can take values between 0 and 1. An index value below 0.25 indicates a highly mixed fishery, between 0.25 and 0.5 a rather mixed fishery, between 0.5 and 0.75 a rather targeted fishery and above 0.75 a highly targeted fishery.
#' The result of the function is a data frame, which can be printed or stored as an object.
#' Alternatively, an html-table can be created by setting style = "html".
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param style The style, in which the result is printed. Defaults to "basic". "html" produces an html-table.
#' @keywords clustering
#' @keywords plot
#' @keywords HHI
#' @export HHI_table
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' HHI_table(data = stockdata,clustering = clustering)
HHI_table <-  function(data, clustering,style="basic"){
  dataframe <- data
  names(dataframe) <- c("ship_ID","stock","landings")
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))

  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  HHI_stocks <- dataframe %>%
    group_by(ship_ID) %>%
    mutate(share_stock = landings/sum(landings))%>%
    summarise(HHI = sum(share_stock^2)) %>%
    left_join(clustering, by="ship_ID")%>%
    group_by(cluster)%>%
    summarise(across(HHI,list(median = ~median(.x, na.rm = TRUE), min = ~min(.x, na.rm=T),max = ~max(.x, na.rm=T)))) %>%
    ungroup()

  HHI_clusters <- left_join(dataframe, clustering, by="ship_ID") %>%
    group_by(cluster, stock) %>%
    summarise(landings_stock = sum(landings)) %>%
    group_by(cluster) %>%
    mutate(share_stock = landings_stock/sum(landings_stock))%>%
    summarise(HHI_clust = sum(share_stock^2)) %>%
    ungroup()

  HHI <- left_join(HHI_stocks,HHI_clusters, by="cluster")
  HHI$cluster <- factor(HHI$cluster, levels = (clusterlevels))
  HHI[,2:5] <- round(HHI[,2:5], digits = 2)
  names(HHI) <- c("cluster","HHI (median)","HHI (min)","HHI (max)", "HHI (cluster)")

  if(style=="basic"){
    return(HHI)
  }
  if(style=="html"){
    print(kbl(HHI[,2:5],format = "html", align = "c",caption = "Catch HHI of ships in clusters") %>%
            pack_rows(index = table(fct_inorder(HHI$cluster))) %>%
            kable_styling(full_width = T,bootstrap_options = c("striped","hover","responsive"),fixed_thead = T))
  }

}
#### 17) Plot HHI of catch ####
#' @title Mixed box- and dotplot of Herfindahl-Hirschmann Index (HHI) of catch of single ships and overall catch of clusters.
#'
#' @description This function creates a mixed box- and dotplot of the HHI of the catch of 1) single ships and 2) clusters. The HHI is equivalent to the Simpson Index,
#' it is a measure of degree of concentration/consolidation. In case of fisheries catch data, it indicates, whether a fishery is of rather targeted or mixed nature.
#' The HHI can take values between 0 and 1. An index value below 0.25 indicates a highly mixed fishery, between 0.25 and 0.5 a rather mixed fishery, between 0.5 and 0.75 a rather targeted fishery and above 0.75 a highly targeted fishery.
#' @param data The original, untransformed data that was used for the clustering.
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @keywords clustering
#' @keywords plot
#' @keywords HHI
#' @export HHI_plot
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' HHI_plot(data = stockdata,clustering = clustering)
HHI_plot <- function(data, clustering){
  dataframe <- data
  names(dataframe) <- c("ship_ID","stock","landings")
  HHI_stocks <- dataframe %>%
    group_by(ship_ID) %>%
    mutate(share_stock = landings/sum(landings))%>%
    summarise(HHI_stocks = sum(share_stock^2)) %>%
    left_join(clustering, by="ship_ID")

  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))

  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  HHI_clusters <- left_join(dataframe, clustering, by="ship_ID") %>%
    group_by(cluster, stock) %>%
    summarise(landings_stock = sum(landings)) %>%
    group_by(cluster) %>%
    mutate(share_stock = landings_stock/sum(landings_stock))%>%
    summarise(HHI_clusters = sum(share_stock^2)) %>%
    ungroup()

  HHI_clusters$cluster <- factor(HHI_clusters$cluster, levels = rev(clusterlevels))

  HHI_stocks <- HHI_stocks %>% group_by(cluster) %>% mutate(clust_size = n_distinct(ship_ID))

  HHI_stocks$cluster <- factor(HHI_stocks$cluster, levels = rev(clusterlevels))
  HHI_clusters$cluster <- factor(HHI_clusters$cluster, levels = rev(clusterlevels))

  HHI_plot <- ggplot()+
    geom_point(data=HHI_clusters, aes(cluster, HHI_clusters),alpha=0)+
    geom_point(data=dplyr::filter(HHI_stocks,clust_size <= 5), aes(cluster, HHI_stocks),fill="#29505a",size=3,shape=21, colour="black",alpha=.8)+
    geom_boxplot(data=dplyr::filter(HHI_stocks,clust_size > 5), aes(cluster, HHI_stocks),fill="#29505a", colour="black",alpha=.6)+
    geom_point(data=HHI_clusters, aes(cluster, HHI_clusters),fill="#c4d9df",size=5,shape=21, colour="black",alpha=.9)+
    theme_bw() +
    theme(axis.title.y = element_blank(),axis.text.y = element_text(face="bold",size=12),
          legend.position = "right",legend.title = element_blank())+
    labs(y="HHI")+
    scale_y_continuous(breaks = c(0,.25,.50,.75,1),limits = c(0,1),labels = c(0,.25,.50,.75,1))+
    coord_flip()
  HHI_plot
}

#### 18) MDS of clustering ####
#' @title Multi-dimensional scaling (MDS) of the clustering
#'
#' @description This is function creates an MDS of the clustering result. The MDS can be either 2-dimensional (which is the default setting), or 3-dimensional.
#' 3-dimensional MDS are harder to interpret, but due to the nature of compositional catch data, 2-dimensional MDS often have a poor goodness of fit (GoF) and have to be treated with caution.
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param dim The dimensions of the MDS. Use `2` for a 2-dimensional, classic MDS and `3` for a 3-dimensional MDS.
#' @param GoF Display goodness of fit in the MDS plot. Defaults to TRUE
#' @param distance The distance measure used. Defaults to modified (metric conversion) Bray-Curtis distance distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @keywords clustering
#' @keywords MDS
#' @export clustering_MDS
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' clustering_MDS(catchdata = catchdata,clustering = clustering, GoF=TRUE)
#' clustering_MDS(catchdata = catchdata,clustering = clustering,dim = 3)
clustering_MDS <- function(catchdata,clustering, dim=2,GoF=T, distance="jaccard"){
  clust_number <- as.numeric(n_distinct(as.character(clustering$cluster)))
  clusterlevels <- c()
  for (x in 1:clust_number) {
    levels <- as.vector(c(paste("cluster", x)))
    clusterlevels <- c(clusterlevels,levels)
  }
  clustering$cluster <- factor(clustering$cluster, levels = clusterlevels)
  catchdata_clustering <- catchdata
  catchdata_clustering$ship_ID <- rownames(catchdata_clustering)
  catchdata_clustering <- left_join(catchdata_clustering,clustering,by="ship_ID")
  suppressWarnings(
    mds <- catchdata %>%
      vegdist(method = distance) %>%
      cmdscale() %>%
      as_tibble() %>%
      mutate(cluster = catchdata_clustering$cluster)%>%
      mutate(ship_ID = catchdata_clustering$ship_ID)
  )
  mdspalette <- c("#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590","#277da1",
                  "#360568","#5b2a86","#7785ac","#9ac6c5","#a5e6ba","#dbebc0","#779cab","#627c85","#35524a","#0a2342",
                  "#e54b4b","#f46197","#efc3f5","#e28413","#ba9593","#cf0e12","#757761","#51bbfe","#94ecbe","#003459",
                  "#eefc57","#e2a0ff","#32de8a","#e87ea1","#53f4ff")

  colnames(mds) <- c("Dim.1", "Dim.2","cluster","ship_ID")

  # Plot MDS
  mds_midpoints_cluster <- mds %>%
    group_by(cluster)%>%
    mutate(MeanDim1 = mean(Dim.1)) %>%
    mutate(MeanDim2 = mean(Dim.2))
  mds_midpoints_cluster <- unique(mds_midpoints_cluster[,c("MeanDim1","MeanDim2","cluster")])
  mds_midpoints_cluster$cluster_nr <- as.numeric(gsub("cluster ", "", mds_midpoints_cluster$cluster))

  # get goodness of fit
  fit <- suppressWarnings(cmdscale(vegdist(catchdata,distance="jaccard"),T, k=2)) # k is the number of dim
  GOF <- fit$GOF[1]
  GOF_x_pos <- max(mds$Dim.1*.7)
  GOF_y_pos <- max(mds$Dim.2*.9)
  GOF_label <- paste("GoF =",round(GOF,digits = 2))

  if(dim==2){
    TwoD_MDS <- ggplot(mds) +
      geom_point(aes(Dim.1,Dim.2,colour=cluster),size=2)+
      geom_mark_hull(aes(Dim.1,Dim.2,fill=cluster, colour=cluster),concavity = 5,expand=0,radius=0,alpha=0.15)+
      geom_label_repel(data=mds_midpoints_cluster,aes(MeanDim1, MeanDim2, label=cluster_nr, fill=cluster),colour="white",show.legend = F,fontface="bold")+
      scale_fill_manual(values = mdspalette)+
      scale_colour_manual(values = mdspalette)+
      labs(x="Dimension 1", y="Dimension 2")+
      theme_bw()+
      theme(legend.position = "none")
    if(GoF==T){
      return(TwoD_MDS + geom_label(data=tibble(),aes(GOF_x_pos,GOF_y_pos,label=GOF_label),colour="black",size=4,fontface="bold", alpha=.5))
    }
    else{
      return(TwoD_MDS)
    }
  }

  if(dim==3){
    suppressWarnings(
      mds_3d <- catchdata %>%
        vegdist(method = distance) %>%
        cmdscale(k = 3) %>%
        data.frame() %>%
        mutate(cluster = catchdata_clustering$cluster)%>%
        mutate(eunr = catchdata_clustering$ship_ID)
    )
    suppressWarnings(colnames(mds_3d) <- c("Dim.1", "Dim.2","Dim.3","cluster","ship_ID"))
    options(warn = -1)
    return(suppressWarnings(plot_ly(data = mds_3d,x=~Dim.1, y=~Dim.2, z=~Dim.3, type="scatter3d",mode="markers",color = ~cluster,colors = mdspalette)))
  }
}

#### 19) MDS of cluster assemblages ####
#' @title Multi-dimensional scaling (MDS) of the assemblage-based cluster catches
#'
#' @description This is function creates an MDS of the assemblage-based catches of the clustering result. The MDS can be either 2-dimensional (which is the default setting), or 3-dimensional.
#' 3-dimensional MDS are harder to interpret, but due to the nature of compositional catch data, 2-dimensional MDS often have a poor goodness of fit (GoF) and have to be treated with caution.
#' @param data The original, untransformed data that was used for the clustering.
#' @param catchdata The transformed catchdata created with catchdata_transformation()
#' @param clustering The result of the clustering procedure, stored as a data frame.
#' @param dim The dimensions of the MDS. Use `2` for a 2-dimensional, classic MDS and `3` for a 3-dimensional MDS.
#' @param GoF Display goodness of fit in the MDS plot. Defaults to TRUE
#' @param distance The distance measure used. Defaults to modified (metric conversion) Bray-Curtis distance distance. CAUTION! The clustering approach for the fleet segmentation is designed to work with modified (metric-converted) Bray-Curtis distance and the average linkage method! Changing either of them is not advised!
#' @keywords clustering
#' @keywords MDS
#' @export cluster_assemblages_MDS
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' catchdata <- catchdata_transformation(data = stockdata)
#' clustering <- segmentation_clustering(catchdata = catchdata,n_cluster = 6)
#' cluster_assemblages_MDS(data =data, catchdata = catchdata,clustering = clustering, GoF=TRUE)
cluster_assemblages_MDS <- function(data,catchdata,clustering, interactive=F,GoF=T, distance="jaccard"){

  assemblage <- data %>%
    mutate(species_code = toupper(sub("\\..*", "", stock))) %>%
    mutate(species_code = toupper(sub("\\-.*", "", species_code))) %>%
    left_join(assemblage_red) %>%
    mutate(target_assemblage_code = ifelse(species_code == "NEP", "CRU",target_assemblage_code),
           target_assemblage = ifelse(species_code == "NEP", "Crustaceans",target_assemblage),
           target_assemblage_code = replace_na(target_assemblage_code,"UNK"),
           target_assemblage = replace_na(target_assemblage, "unknown")) %>%
    left_join(clustering) %>%
    group_by(cluster,target_assemblage_code) %>%
    summarise(catch_cluster = sum(landings)) %>%
    group_by(cluster) %>%
    mutate(share_assemblage = catch_cluster/sum(catch_cluster)) %>%
    ungroup()

  assemblage_matrix <- assemblage %>%
    dplyr::select(cluster,target_assemblage_code,share_assemblage) %>%
    pivot_wider(names_from = target_assemblage_code,values_from = share_assemblage,values_fill = 0)%>%
    ungroup()

  assemblage_table <- as.data.frame(assemblage_matrix)
  assemblage_table <- assemblage_table[,-1]
  rownames(assemblage_table) <- assemblage_matrix$cluster

  mds.assemblage <- cmdscale(vegdist(x = assemblage_table,method = distance)) %>% as.data.frame()
  mds.assemblage$names <- rownames(mds.assemblage)

  # get goodness of fit
  fit <- suppressWarnings(cmdscale(vegdist(catchdata,method=distance),T, k=2)) # k is the number of dim
  GOF <- fit$GOF[1]
  GOF_x_pos <- max(mds.assemblage[1])*.7
  GOF_y_pos <- max(mds.assemblage[2])*.9
  GOF_label <- paste("GoF =",round(GOF,digits = 2))

  ### mds
  assemblage_mds <- ggplot(mds.assemblage, aes(V1, V2, label=names)) +
    geom_point(colour="black",fill="#B1D0E8", size=4,shape=22) +
    geom_text(colour="blue", check_overlap = TRUE, size=2.5,
              hjust = "center", vjust = "bottom", nudge_x = 0, nudge_y = 0.025) +
    theme_minimal()+
    theme(axis.title = element_blank())

  if(GoF==T){
    return(assemblage_mds + geom_label(data=tibble(),aes(GOF_x_pos,GOF_y_pos,label=GOF_label),colour="black",size=4,fontface="bold", alpha=.5))
  }
  else{
    return(assemblage_mds)
  }
}


#### 20) Assigning ICES-Stocks ####
#' @title Function to assign ICES stocks to data.
#'
#' @description This function assigns ICES stocks to catch data based on the caught species and the FAO fishing area where it was caught.
#' Due to limits in available data, especially on a fine spatial resolution, not all ICES-stocks can be taken into account. Therefore, stocks of Nephrops norvegicus
#' and Ammodytes spp. are based on EU definitions, not on ICES definitions. Additionally, the function automatically creates stocks based on species name and FAO area for species
#' which are caught in larger quantities (> 10000 kg) but are not part of defined ICES-stock. This feature can be shut off if the user decides not to apply the procedure and
#' use only defined ICES-stocks.
#' @param data A basic catchdata frame consisting of 4 columns:
#' 1) A unique ID for each ship.
#' 2) The 3-letter species code for the caught species as specified by ICES.
#' 3) The full area code of the division where the fish was caught. The area can be an FAO area, e.g. "27.4.b" for Central North Sea or "27.5.b.1" for Faroe Plateau.
#' If the fishing area is one of the Mediterranean Management area, please indicate "GSA", e.g. "GSA 17".
#' 4) The catch weight of the respective species in the respective area, in kg.
#' @param reduce Indicates, whether or not the resulting data frame is reduced to Ship ID, ICES stock and catch weight. Defaults to TRUE.
#' Turned to FALSE, the resulting data frame will resemble the original data with the ICES-Stock column added. This format is used for control and correction purposes,
#' but not for further calculations.
#' @param auto.generate Indicates whether or stocks should be automatically generated for
#' a) species, which are not assessed in ICES-Stocks or
#' b) assessed by ICES, but caught out of stock-managed areas.
#' The automatically generated stocks comprise the species name and the FAO area.
#' The relevant quantity is defined by the argument threshold.auto.generate
#' @param threshold.auto.generate Threshold of automatic generation of ICES stocks. Only relevant if auto.generate = T. Defaults to 100.
#' @param min.share The minimal share a stock has to have on at least one vessels catch to be included in the stock dataframe. Defaults to 0, i.e. every stock is retained by default.
#' @keywords Stocks
#' @export assign_stocks
#' @examples
#' data <- example_catchdata
#' stockdata <- assign_stocks(data=data)
#' stockdata <- assign_stocks(data=data, reduce=FALSE, auto.generate=FALSE)
assign_stocks <- function(data, reduce=T, auto.generate=T,threshold.auto.generate=100,min.share=0){
  dataframe <- data
  # Code for assigning stocks
  names(dataframe) <- c("ship_ID","species","area","landkg")
  dataframe$stock <- "Bycatch/Unknown"
  dataframe$area <- gsub(pattern = ",",replacement = ".",x = dataframe$area)
  dataframe$area <- gsub(pattern = ";",replacement = ".",x = dataframe$area)

  ### The following is Bernies Code for the ICES-Stocks #
  # A) Ostsee ----

  dataframe$stock[dataframe$species=="HER" & dataframe$area=="27.3.a.21" | dataframe$species=="HER" & dataframe$area=="27.3.c.22" | dataframe$species=="HER" & dataframe$area=="27.3.d.24"
                        |  dataframe$species=="HER" & dataframe$area=="27.3.b.23" |  dataframe$species=="HER" & dataframe$area=="27.3.b" ]<-"her.27.20-24"
  # Skagerrak ist hier ausgelassen, da dies meisten Hering in 27.3.a.20 dem Nordseebestand dataframeugerechnet werden!!!
  dataframe$stock[dataframe$species=="HER" & dataframe$area=="27.3.d" |dataframe$species=="HER" & dataframe$area=="27.3.d.25" | dataframe$species=="HER" & dataframe$area=="27.3.d.26" | dataframe$species=="HER" & dataframe$area=="27.3.d.27"
                       |  dataframe$species=="HER" & dataframe$area=="27.3.d.28.2" |  dataframe$species=="HER" & dataframe$area=="27.3.d.29" |  dataframe$species=="HER" & dataframe$area=="27.3.d.32"]<-"her.27.25-2932"
  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.3.d.25" | dataframe$species=="COD" & dataframe$area=="27.3.d.26" | dataframe$species=="COD" & dataframe$area=="27.3.d.27"
                       |  dataframe$species=="COD" & dataframe$area=="27.3.d.28.2" |  dataframe$species=="COD" & dataframe$area=="27.3.d.29" |  dataframe$species=="COD" & dataframe$area=="27.3.d.30"
                       |  dataframe$species=="COD" & dataframe$area=="27.3.d.31" |  dataframe$species=="COD" & dataframe$area=="27.3.d.32"]<-"cod.27.24-32"
  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.3.c.22" | dataframe$species=="COD" & dataframe$area=="27.3.b.23" | dataframe$species=="COD" & dataframe$area=="27.3.d.24"| dataframe$species=="COD" & dataframe$area=="27.3.d"]<-"cod.27.22-24"
  dataframe$stock[dataframe$species=="SPR" & dataframe$area=="27.3.c.22" | dataframe$species=="SPR" & dataframe$area=="27.3.b.23" | dataframe$species=="SPR" & dataframe$area=="27.3.d.24" | dataframe$species=="SPR" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="SPR" & dataframe$area=="27.3.d.26" | dataframe$species=="SPR" & dataframe$area=="27.3.d.27" | dataframe$species=="SPR" & dataframe$area=="27.3.d.28.2" | dataframe$species=="SPR" & dataframe$area=="27.3.d.29"
                       | dataframe$species=="SPR" & dataframe$area=="27.3.d.30" | dataframe$species=="SPR" & dataframe$area=="27.3.d.31" | dataframe$species=="SPR" & dataframe$area=="27.3.d.32"| dataframe$species=="SPR" & dataframe$area=="27.3.d"]<-"spr.27.22-32"
  dataframe$stock[dataframe$species=="PLE" & dataframe$area=="27.3.a.21" | dataframe$species=="PLE" & dataframe$area=="27.3.c.22" | dataframe$species=="PLE" & dataframe$area=="27.3.b.23"]<-"ple.27.21-23"
  dataframe$stock[dataframe$species=="SOL" & dataframe$area=="27.3.a.n" | dataframe$species=="SOL" & dataframe$area=="27.3.a" | dataframe$species=="SOL" & dataframe$area=="27.3.a.21" | dataframe$species=="SOL" & dataframe$area=="27.3.c.22"
                       | dataframe$species=="SOL" & dataframe$area=="27.3.b.23"| dataframe$species=="SOL" & dataframe$area=="27.3.d.24" | dataframe$species=="SOL" & dataframe$area=="27.3.d.26" | dataframe$species=="SOL" & dataframe$area=="27.3.d.25" | dataframe$species=="SOL" & dataframe$area=="27.3.a.20"]<-"sol.27.20-24"

  dataframe$stock[dataframe$species=="DAB" & dataframe$area=="27.3.c.22" | dataframe$species=="DAB" & dataframe$area=="27.3.b.23" | dataframe$species=="DAB" & dataframe$area=="27.3.d.24" | dataframe$species=="DAB" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="DAB" & dataframe$area=="27.3.d.26" | dataframe$species=="DAB" & dataframe$area=="27.3.d.27" | dataframe$species=="DAB" & dataframe$area=="27.3.d.28.2" | dataframe$species=="DAB" & dataframe$area=="27.3.d.29"
                       | dataframe$species=="DAB" & dataframe$area=="27.3.d.30" | dataframe$species=="DAB" & dataframe$area=="27.3.d.31" | dataframe$species=="DAB" & dataframe$area=="27.3.d.32" | dataframe$species=="DAB" & dataframe$area=="27.3.d"]<-"dab.27.22-32"
  dataframe$stock[dataframe$species=="FLE" & dataframe$area=="27.3.b.23" | dataframe$species=="FLE" & dataframe$area=="27.3.c.22"]<-"fle.27.2223"
  dataframe$stock[dataframe$species=="FLE" & dataframe$area=="27.3.d.24" | dataframe$species=="FLE" & dataframe$area=="27.3.d.25"]<-"fle.27.2425"
  dataframe$stock[dataframe$species=="TUR" & dataframe$area=="27.3.c.22" | dataframe$species=="TUR" & dataframe$area=="27.3.b.23" | dataframe$species=="TUR" & dataframe$area=="27.3.d.24" | dataframe$species=="TUR" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="TUR" & dataframe$area=="27.3.d.26" | dataframe$species=="TUR" & dataframe$area=="27.3.d.27" | dataframe$species=="TUR" & dataframe$area=="27.3.d.28.2" | dataframe$species=="TUR" & dataframe$area=="27.3.d.29"
                       | dataframe$species=="TUR" & dataframe$area=="27.3.d.30" | dataframe$species=="TUR" & dataframe$area=="27.3.d.31" | dataframe$species=="TUR" & dataframe$area=="27.3.d.32"| dataframe$species=="TUR" & dataframe$area=="27.3.d"]<-"tur.27.22-32"
  dataframe$stock[dataframe$species=="PLE" & dataframe$area=="27.3.d.24" | dataframe$species=="PLE" & dataframe$area=="27.3.d.25" | dataframe$species=="PLE" & dataframe$area=="27.3.d.26"
                       | dataframe$species=="PLE" & dataframe$area=="27.3.d.27" | dataframe$species=="PLE" & dataframe$area=="27.3.d.28.2" | dataframe$species=="PLE" & dataframe$area=="27.3.d.29"
                       | dataframe$species=="PLE" & dataframe$area=="27.3.d.30" | dataframe$species=="PLE" & dataframe$area=="27.3.d" | dataframe$species=="PLE" & dataframe$area=="27.3.d.31" | dataframe$species=="PLE" & dataframe$area=="27.3.d.32"]<-"ple.27.24-32"
  dataframe$stock[dataframe$species=="BLL" & dataframe$area=="27.3.c.22" | dataframe$species=="BLL" & dataframe$area=="27.3.b.23" | dataframe$species=="BLL" & dataframe$area=="27.3.d.24" | dataframe$species=="BLL" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="BLL" & dataframe$area=="27.3.d.26" | dataframe$species=="BLL" & dataframe$area=="27.3.d.27" | dataframe$species=="BLL" & dataframe$area=="27.3.d.28.2" | dataframe$species=="BLL" & dataframe$area=="27.3.d.29"
                       | dataframe$species=="BLL" & dataframe$area=="27.3.d.30" | dataframe$species=="BLL" & dataframe$area=="27.3.d.31" | dataframe$species=="BLL" & dataframe$area=="27.3.d.32"]<-"bll.27.22-32"
  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.3.a.21" | dataframe$species=="COD" & dataframe$area=="27.3.a"]<-"cod.27.21"
  #dataframe$stock[dataframe$species=="SPR" & dataframe$area=="27.3.a.n" | dataframe$species=="SPR" & dataframe$area=="27.3.a.21"]<-"spr.27.3a"
  dataframe$stock[dataframe$species=="WHG" & dataframe$area=="27.3.a.n" | dataframe$species=="WHG" & dataframe$area=="27.3.a" | dataframe$species=="WHG" & dataframe$area=="27.3.a.21" | dataframe$species=="WHG" & dataframe$area=="27.3.a.20"]<-"whg.27.3a"
  dataframe$stock[dataframe$species=="RNG" & dataframe$area=="27.3.a.n" | dataframe$species=="RNG" & dataframe$area=="27.3.a.21" | dataframe$species=="RNG" & dataframe$area=="27.3.a.20"]<-"rng.27.3a"

  # B) Nordsee ----

  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.4.a" | dataframe$species=="COD" & dataframe$area=="27.4.b" | dataframe$species=="COD" & dataframe$area=="27.4.c"
                       | dataframe$species=="COD" & dataframe$area=="27.7.d" | dataframe$species=="COD" & dataframe$area=="27.3.a.n" | dataframe$species=="COD" & dataframe$area=="27.3.a.20"]<-"cod.27.47d20"
  dataframe$stock[dataframe$species=="HAD" & dataframe$area=="27.4.a" | dataframe$species=="HAD" & dataframe$area=="27.4.b" | dataframe$species=="HAD" & dataframe$area=="27.4.c"
                       | dataframe$species=="HAD" & dataframe$area=="27.6.a" | dataframe$species=="HAD" & dataframe$area=="27.3.a.n" | dataframe$species=="HAD" & dataframe$area=="27.3.a" | dataframe$species=="HAD" & dataframe$area=="27.3.c.22" |dataframe$species=="HAD" & dataframe$area=="27.3.a.20"]<-"had.27.46a20"
  dataframe$stock[dataframe$species=="HER" & dataframe$area=="27.4.a" | dataframe$species=="HER" & dataframe$area=="27.4.b" | dataframe$species=="HER" & dataframe$area=="27.4.c"
                       | dataframe$species=="HER" & dataframe$area=="27.7.d" | dataframe$species=="HER" & dataframe$area=="27.3.a.n" | dataframe$species=="HER" & dataframe$area=="27.3.a.20"]<-"her.27.3a47d"
  ## OBACHT! Eigentlich m?ssten auch die Hering aus
  # 27.3.a.21 in diese Kategorie bdataframew. einige aus diesem Bestand in den "her-3a22"-Bestand, da es im Skagerrak und Kattegat starke Vermischungen gibt. Norbert hat jeweils den
  # "Verteilungsschl?ssel", wie die F?nge welchem Bestand dataframeugeordnet werden. F?r den Flottenbericht belasse ich es aber bei dieser Unsch?rfe!
  dataframe$stock[dataframe$species=="JAX" & dataframe$area=="27.4.b" | dataframe$species=="JAX" & dataframe$area=="27.4.c" | dataframe$species=="JAX" & dataframe$area=="27.3.a.n"
                       | dataframe$species=="JAX" & dataframe$area=="27.3.a.21" | dataframe$species=="JAX" & dataframe$area=="27.7.d" | dataframe$species=="JAX" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="JAX" & dataframe$area=="27.4.b" | dataframe$species=="JAX" & dataframe$area=="27.4.c" | dataframe$species=="JAX" & dataframe$area=="27.3.a.n"
                       | dataframe$species=="JAX" & dataframe$area=="27.3.a.21" | dataframe$species=="JAX" & dataframe$area=="27.7.d"| dataframe$species=="JAX" & dataframe$area=="27.3.c.22"
                       |dataframe$species=="HOM" & dataframe$area=="27.4.b" | dataframe$species=="HOM" & dataframe$area=="27.4.c" | dataframe$species=="HOM" & dataframe$area=="27.3.a.n"
                       | dataframe$species=="HOM" & dataframe$area=="27.3.a.21" | dataframe$species=="HOM" & dataframe$area=="27.7.d" | dataframe$species=="HOM" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="HOM" & dataframe$area=="27.4.b" | dataframe$species=="HOM" & dataframe$area=="27.4.c" | dataframe$species=="HOM" & dataframe$area=="27.3.a.n"
                       | dataframe$species=="HOM" & dataframe$area=="27.3.a.21" | dataframe$species=="HOM" & dataframe$area=="27.7.d" | dataframe$species=="HOM" & dataframe$area=="27.3.c.22" ]<-"hom.27.3a4bc7d"
  dataframe$stock[dataframe$species=="PLE" & dataframe$area=="27.4.a" | dataframe$species=="PLE" & dataframe$area=="27.4.b" | dataframe$species=="PLE" & dataframe$area=="27.4.c"
                       | dataframe$species=="PLE" & dataframe$area=="27.3.a.n" | dataframe$species=="PLE" & dataframe$area=="27.3.a.20" | dataframe$species=="PLE" & dataframe$area=="27.3.a"]<-"ple.27.420"
  dataframe$stock[dataframe$species=="PLE" & dataframe$area=="27.7.d"]<-"ple.27.7d"
  dataframe$stock[dataframe$species=="DAB" & dataframe$area=="27.4.a" | dataframe$species=="DAB" & dataframe$area=="27.4.b" | dataframe$species=="DAB" & dataframe$area=="27.4.c"
                       | dataframe$species=="DAB" & dataframe$area=="27.3.a.n" | dataframe$species=="DAB" & dataframe$area=="27.3.a.21" | dataframe$species=="DAB" & dataframe$area=="27.3.a.20" | dataframe$species=="DAB" & dataframe$area=="27.3.a"]<-"dab.27.3a4"
  dataframe$stock[dataframe$species=="BLL" & dataframe$area=="27.4.a" | dataframe$species=="BLL" & dataframe$area=="27.4.b" | dataframe$species=="BLL" & dataframe$area=="27.4.c"
                       | dataframe$species=="BLL" & dataframe$area=="27.3.a.n" | dataframe$species=="BLL" & dataframe$area=="27.3.a" | dataframe$species=="BLL" & dataframe$area=="27.3.a.21" | dataframe$species=="BLL" & dataframe$area=="27.7.d"
                       | dataframe$species=="BLL" & dataframe$area=="7E" | dataframe$species=="BLL" & dataframe$area=="27.3.a.20"]<-"bll.27.3a47de"
  dataframe$stock[dataframe$species=="ANF" & dataframe$area=="27.4.a" |dataframe$species=="MON" & dataframe$area=="27.4.a" | dataframe$species=="ANF" & dataframe$area=="27.4.b" | dataframe$species=="ANF" & dataframe$area=="27.4.c"
                       | dataframe$species=="ANF" & dataframe$area=="27.6.a" | dataframe$species=="ANF" & dataframe$area=="6B"| dataframe$species=="ANF" & dataframe$area=="27.6.b"| dataframe$species=="ANF" & dataframe$area=="27.6.b.2"| dataframe$species=="ANF" & dataframe$area=="27.3.a.n"
                       | dataframe$species=="ANF" & dataframe$area=="27.3.a" | dataframe$species=="ANF" & dataframe$area=="27.3.a.21" | dataframe$species=="ANF" & dataframe$area=="27.3.a.20" | dataframe$species=="MON" & dataframe$area=="27.4.b" ]<-"anf.27.3a46"
  dataframe$stock[dataframe$species=="BSS" & dataframe$area=="27.4.b" | dataframe$species=="BSS" & dataframe$area=="27.4.c" | dataframe$species=="BSS" & dataframe$area=="7A"
                       | dataframe$species=="BSS" & dataframe$area=="27.7.d" | dataframe$species=="BSS" & dataframe$area=="7E" | dataframe$species=="BSS" & dataframe$area=="7F"
                       | dataframe$species=="BSS" & dataframe$area=="7G" | dataframe$species=="BSS" & dataframe$area=="7H" | dataframe$species=="BSS" & dataframe$area=="27.3.c.22"]<-"bss.27.4bc7ad-h"
  dataframe$stock[dataframe$species=="FLE" & dataframe$area=="27.4.a" | dataframe$species=="FLE" & dataframe$area=="27.4.b" | dataframe$species=="FLE" & dataframe$area=="27.4.c"
                       | dataframe$species=="FLE" & dataframe$area=="27.3.a.n" | dataframe$species=="FLE" & dataframe$area=="27.3.a.21" | dataframe$species=="FLE" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="FLE" & dataframe$area=="27.3.a" | dataframe$species=="FLE" & dataframe$area=="27.3.d.26" | dataframe$species=="FLE" & dataframe$area=="27.3.d"]<-"fle.27.3a4"
  dataframe$stock[dataframe$species=="MEG" & dataframe$area=="27.4.a" | dataframe$species=="MEG" & dataframe$area=="6A" | dataframe$species=="MEG" & dataframe$area=="27.2.a.2"| dataframe$species=="MEG" & dataframe$area=="27.2.a"
                       | dataframe$species=="MEG" & dataframe$area=="27.4.b" | dataframe$species=="MEG" & dataframe$area=="27.3.a"| dataframe$species=="MEG" & dataframe$area=="27.3.a.n"]<-"meg.27.4a6a"
  dataframe$stock[dataframe$species=="NOP" & dataframe$area=="27.4.a" | dataframe$species=="NOP" & dataframe$area=="27.4.b" | dataframe$species=="NOP" & dataframe$area=="27.4.c"
                       | dataframe$species=="NOP" & dataframe$area=="27.3.a.n" | dataframe$species=="NOP" & dataframe$area=="27.3.a.21" | dataframe$species=="NOP" & dataframe$area=="27.3.a.20" ]<-"nop.27.34"
  dataframe$stock[dataframe$species=="POK" & dataframe$area=="27.4.a" | dataframe$species=="POK" & dataframe$area=="27.4.b" | dataframe$species=="POK" & dataframe$area=="27.4.c"
                       | dataframe$species=="POK" & dataframe$area=="27.6.b" |dataframe$species=="POK" & dataframe$area=="27.6.a" | dataframe$species=="POK" & dataframe$area=="6B" | dataframe$species=="POK" & dataframe$area=="27.3.a.n"
                       | dataframe$species=="POK" & dataframe$area=="27.3.a.21" | dataframe$species=="POK" & dataframe$area=="27.3.a.20" | dataframe$species=="POK" & dataframe$area=="27.3.a"
                       | dataframe$species=="POK" & dataframe$area=="27.3.d.24" | dataframe$species=="POK" & dataframe$area=="27.3.c.22"]<-"pok.27.3a46"
  dataframe$stock[dataframe$species=="SOL" & dataframe$area=="27.4.a" | dataframe$species=="SOL" & dataframe$area=="27.4.b" | dataframe$species=="SOL" & dataframe$area=="27.4.c"]<-"sol.27.4"
  dataframe$stock[dataframe$species=="SPR" & dataframe$area=="27.4.a" | dataframe$species=="SPR" & dataframe$area=="27.4.b" | dataframe$species=="SPR" & dataframe$area=="27.4.c"
                       | dataframe$species=="SPR" & dataframe$area=="27.3.a.21" | dataframe$species=="SPR" & dataframe$area=="27.3.a.20" | dataframe$species=="SPR" & dataframe$area=="27.3.a.n" | dataframe$species=="SPR" & dataframe$area=="27.3.a"]<-"spr.27.3a4"
  dataframe$stock[dataframe$species=="WHG" & dataframe$area=="27.4.a" | dataframe$species=="WHG" & dataframe$area=="27.4.b" | dataframe$species=="WHG" & dataframe$area=="27.4.c"
                       | dataframe$species=="WHG" & dataframe$area=="27.7.d"]<-"whg.27.47d"
  dataframe$stock[dataframe$species=="POL" & dataframe$area=="27.4.a" | dataframe$species=="POL" & dataframe$area=="27.4.b" | dataframe$species=="POL" & dataframe$area=="27.4.c"
                       | dataframe$species=="POL" & dataframe$area=="27.3.a.n" | dataframe$species=="POL" & dataframe$area=="27.3.a.21" | dataframe$species=="POL" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="POL" & dataframe$area=="27.3.a"| dataframe$species=="POL" & dataframe$area=="27.3.c.22"| dataframe$species=="POL" & dataframe$area=="27.1.b"
                       | dataframe$species=="POL" & dataframe$area=="27.2.a.2"| dataframe$species=="POL" & dataframe$area=="27.2.a"]<-"pol.27.3a4"
  dataframe$stock[dataframe$species=="TUR" & dataframe$area=="27.4.a" | dataframe$species=="TUR" & dataframe$area=="27.4.b" | dataframe$species=="TUR" & dataframe$area=="27.4.c"]<-"tur.27.4"
  dataframe$stock[dataframe$species=="TUR" & dataframe$area=="27.3.a.21" | dataframe$species=="TUR" & dataframe$area=="27.3.a" |dataframe$species=="TUR" & dataframe$area=="27.3.a.20" | dataframe$species=="TUR" & dataframe$area=="27.3.a.n"]<-"tur.27.3a"

  dataframe$stock[dataframe$species=="LEM" & dataframe$area=="27.4.a" | dataframe$species=="LEM" & dataframe$area=="27.4.b" | dataframe$species=="LEM" & dataframe$area=="27.4.c"
                       | dataframe$species=="LEM" & dataframe$area=="27.7.d" | dataframe$species=="LEM" & dataframe$area=="27.3.a.n" | dataframe$species=="LEM" & dataframe$area=="27.3.a.21"
                       | dataframe$species=="LEM" & dataframe$area=="27.3.a"| dataframe$species=="LEM" & dataframe$area=="27.3.a.20" | dataframe$species=="LEM" & dataframe$area=="27.3.a" | dataframe$species=="LEM" & dataframe$area=="27.2.a.2"| dataframe$species=="LEM" & dataframe$area=="27.3.c.22"]<-"lem.27.3a47d"

  dataframe$stock[dataframe$species=="WIT" & dataframe$area=="27.4.a" | dataframe$species=="WIT" & dataframe$area=="27.4.b" | dataframe$species=="WIT" & dataframe$area=="27.4.c"
                       | dataframe$species=="WIT" & dataframe$area=="27.7.d" | dataframe$species=="WIT" & dataframe$area=="27.3.a.n" | dataframe$species=="WIT" & dataframe$area=="27.3.a.21"
                       | dataframe$species=="WIT" & dataframe$area=="27.3.a.20" | dataframe$species=="WIT" & dataframe$area=="27.3.a" | dataframe$species=="WIT" & dataframe$area=="27.3.c.22"
                       | dataframe$species == "WIT" & dataframe$area=="27.3.d.24" ]<-"wit.27.3a47d"
  dataframe$stock[dataframe$species=="GUG" & dataframe$area=="27.4.a" | dataframe$species=="GUG" & dataframe$area=="27.4.b" | dataframe$species=="GUG" & dataframe$area=="27.4.c"
                       | dataframe$species=="GUG" & dataframe$area=="27.7.d" | dataframe$species=="GUG" & dataframe$area=="27.3.a.n" | dataframe$species=="GUG" & dataframe$area=="27.3.a.21"
                       | dataframe$species=="GUG" & dataframe$area=="27.3.a.20"]<-"gug.27.3a47d"

  # C) Groenland und Nordostarktis ----
  GREandNEARC <- c("27.1.","27.1.a" ,"27.1","27.2.a.1","27.2.b.1","27.2.a.2","27.2.b.2","27.2.b","27.1.b" ,"27.2.a")
  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.1." |dataframe$species=="COD" & dataframe$area=="27.1.a" | dataframe$species=="COD" & dataframe$area=="27.1" | dataframe$species=="COD" & dataframe$area=="27.1.b" | dataframe$species=="COD" & dataframe$area=="27.2.a"
                       | dataframe$species=="COD" & dataframe$area=="27.2.a.1" | dataframe$species=="COD" & dataframe$area=="27.2.a.2" | dataframe$species=="COD" & dataframe$area=="27.2.b"
                       | dataframe$species=="COD" & dataframe$area=="27.2.b.1" | dataframe$species=="COD" & dataframe$area=="27.2.b.2"]<-"cod.27.1-2"
  dataframe$stock[dataframe$species=="HAD" & dataframe$area=="27.1." | dataframe$species=="HAD" & dataframe$area=="27.1.a" | dataframe$species=="HAD" & dataframe$area=="27.1.b" | dataframe$species=="HAD" & dataframe$area=="27.2.a"
                       | dataframe$species=="HAD" & dataframe$area=="27.2.a.1" | dataframe$species=="HAD" & dataframe$area=="27.2.a.2" | dataframe$species=="HAD" & dataframe$area=="27.2.b"
                       | dataframe$species=="HAD" & dataframe$area=="27.2.b.1" | dataframe$species=="HAD" & dataframe$area=="27.2.b.2"]<-"had.27.1-2"
  dataframe$stock[dataframe$species=="POK" & dataframe$area=="27.1.a" | dataframe$species=="POK" & dataframe$area=="27.1.b" | dataframe$species=="POK" & dataframe$area=="27.2.a"
                       | dataframe$species=="POK" & dataframe$area=="27.2.a.1" | dataframe$species=="POK" & dataframe$area=="27.2.a.2" | dataframe$species=="POK" & dataframe$area=="27.2.b"
                       | dataframe$species=="POK" & dataframe$area=="27.2.b.1" | dataframe$species=="POK" & dataframe$area=="27.2.b.2"]<-"pok.27.1-2"
  dataframe$stock[dataframe$species=="GHL" & dataframe$area=="27.1." |dataframe$species=="GHL" & dataframe$area=="27.1.a" | dataframe$species=="GHL" & dataframe$area=="27.1.b" | dataframe$species=="GHL" & dataframe$area=="27.2.a"
                       | dataframe$species=="GHL" & dataframe$area=="27.2.a.1" | dataframe$species=="GHL" & dataframe$area=="27.2.a.2" | dataframe$species=="GHL" & dataframe$area=="27.2.b"
                       | dataframe$species=="GHL" & dataframe$area=="27.2.b.1" | dataframe$species=="GHL" & dataframe$area=="27.2.b.2"]<-"ghl.27.1-2"
  dataframe$stock[dataframe$species=="LIN" & dataframe$area=="27.1.a" | dataframe$species=="LIN" & dataframe$area=="27.1.b" | dataframe$species=="LIN" & dataframe$area=="27.2.a"
                       | dataframe$species=="LIN" & dataframe$area=="27.2.a.1" | dataframe$species=="LIN" & dataframe$area=="27.2.a.2" | dataframe$species=="LIN" & dataframe$area=="27.2.b"
                       | dataframe$species=="LIN" & dataframe$area=="27.2.b.1" | dataframe$species=="LIN" & dataframe$area=="27.2.b.2"]<-"lin.27.1-2"

  dataframe$stock[dataframe$species=="REG" & dataframe$area=="27.1.a" | dataframe$species=="REG" & dataframe$area=="27.1.b" | dataframe$species=="REG" & dataframe$area=="27.2.a"
                       | dataframe$species=="REG" & dataframe$area=="27.2.a.1" | dataframe$species=="REG" & dataframe$area=="27.2.a.2" | dataframe$species=="REG" & dataframe$area=="27.2.b"
                       | dataframe$species=="REG" & dataframe$area=="27.2.b.1" | dataframe$species=="REG" & dataframe$area=="27.2.b.2"]<-"reg.27.1-2"
  dataframe$stock[dataframe$species=="USK" & dataframe$area=="27.1.a" |dataframe$species=="USK" & dataframe$area=="27.1.c" | dataframe$species=="USK" & dataframe$area=="27.1.b" | dataframe$species=="USK" & dataframe$area=="27.2.a"
                       | dataframe$species=="USK" & dataframe$area=="27.2.a.1" | dataframe$species=="USK" & dataframe$area=="27.2.a.2" | dataframe$species=="USK" & dataframe$area=="27.2.b"
                       | dataframe$species=="USK" & dataframe$area=="27.2.b.1" | dataframe$species=="USK" & dataframe$area=="27.2.b.2"]<-"usk.27.1-2"

  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.14.a" | dataframe$species=="COD" & dataframe$area=="27.14.b"
                       | dataframe$species=="COD" & dataframe$area=="27.14.b.1"| dataframe$species=="COD" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="COD" & dataframe$area=="21.1.F"]<-"cod.2127.1f14"
  dataframe$stock[dataframe$species=="USK" & dataframe$area=="27.14.a" | dataframe$species=="USK" & dataframe$area=="27.14.b"
                       | dataframe$species=="USK" & dataframe$area=="27.14.b.1"| dataframe$species=="USK" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="USK" & dataframe$area=="27.5.a" | dataframe$species=="USK" & dataframe$area=="27.5.a.1"
                       | dataframe$species=="USK" & dataframe$area=="27.5.a.2"]<-"usk.27.5a14"
  dataframe$stock[dataframe$species=="GHL" & dataframe$area=="27.14.a" | dataframe$species=="GHL" & dataframe$area=="27.14.b" | dataframe$species=="GHL" & dataframe$area=="27.14.b.1"| dataframe$species=="GHL" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="GHL" & dataframe$area=="27.5.a" | dataframe$species=="GHL" & dataframe$area=="27.5.a.1" | dataframe$species=="GHL" & dataframe$area=="27.5.a.2"
                       | dataframe$species=="GHL" & dataframe$area=="27.5.b" | dataframe$species=="GHL" & dataframe$area=="27.5.b.1" | dataframe$species=="GHL" & dataframe$area=="27.5.b.1.a"
                       | dataframe$species=="GHL" & dataframe$area=="27.5.b.1.b" | dataframe$species=="GHL" & dataframe$area=="27.5.b.2" | dataframe$species=="GHL" & dataframe$area=="27.6.a"
                       | dataframe$species=="GHL" & dataframe$area=="27.6.b" | dataframe$species=="GHL" & dataframe$area=="27.6.b"| dataframe$species=="GHL" & dataframe$area=="27.6.b.1" | dataframe$species=="GHL" & dataframe$area=="27.6.b.2"
                       | dataframe$species=="GHL" & dataframe$area=="27.12.a" | dataframe$species=="GHL" & dataframe$area=="27.12.a.1" | dataframe$species=="GHL" & dataframe$area=="27.12.a.2"
                       | dataframe$species=="GHL" & dataframe$area=="27.12.a.3" | dataframe$species=="GHL" & dataframe$area=="27.12.a.4" | dataframe$species=="GHL" & dataframe$area=="27.12.b"
                       | dataframe$species=="GHL" & dataframe$area=="27.12.c"] <- "ghl.27.561214"
  dataframe$stock[dataframe$species=="GHL" & dataframe$area=="21.1.F" | dataframe$species=="GHL" & dataframe$area=="21.1.E" | dataframe$species=="GHL" & dataframe$area=="21.1.D"
                       | dataframe$species=="GHL" & dataframe$area=="21.1.C" | dataframe$species=="GHL" & dataframe$area=="21.1.B" | dataframe$species=="GHL" & dataframe$area=="21.1.A"]<-"ghl-NAFO"
  dataframe$stock[dataframe$species=="BLI" & dataframe$area=="27.14.a" | dataframe$species=="BLI" & dataframe$area=="27.14.b" | dataframe$species=="BLI" & dataframe$area=="27.14.b.1"| dataframe$species=="BLI" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="BLI" & dataframe$area=="27.5.a" | dataframe$species=="BLI" & dataframe$area=="27.5.a.1" | dataframe$species=="BLI" & dataframe$area=="27.5.a.2" ]<-"bli.27.5a14"
  dataframe$stock[dataframe$species=="ARU" & dataframe$area=="27.14.a" | dataframe$species=="ARU" & dataframe$area=="27.14.b" | dataframe$species=="ARU" & dataframe$area=="27.14.b.1"| dataframe$species=="ARU" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="ARU" & dataframe$area=="27.5.a" | dataframe$species=="ARU" & dataframe$area=="27.5.a.1" | dataframe$species=="ARU" & dataframe$area=="27.5.a.2" ]<-"aru.27.5a14" # Argentina silus (Greater silver smelt)
  dataframe$stock[dataframe$species=="REG" & dataframe$area=="27.14.a" | dataframe$species=="REG" & dataframe$area=="27.14.b" | dataframe$species=="REG" & dataframe$area=="27.14.b.1"| dataframe$species=="REG" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="REG" & dataframe$area=="27.5.a" | dataframe$species=="REG" & dataframe$area=="27.5.a.1" | dataframe$species=="REG" & dataframe$area=="27.5.a.2"
                       | dataframe$species=="REG" & dataframe$area=="27.5.b" | dataframe$species=="REG" & dataframe$area=="27.5.b.1" | dataframe$species=="REG" & dataframe$area=="27.5.b.1.a"
                       | dataframe$species=="REG" & dataframe$area=="27.5.b.1.b" | dataframe$species=="REG" & dataframe$area=="27.5.b.2" | dataframe$species=="REG" & dataframe$area=="27.6.a"
                       | dataframe$species=="REG" & dataframe$area=="27.6.b" | dataframe$species=="REG" & dataframe$area=="27.6.b" | dataframe$species=="REG" & dataframe$area=="27.6.b.1" | dataframe$species=="REG" & dataframe$area=="27.6.b.2"
                       | dataframe$species=="REG" & dataframe$area=="27.12.a" | dataframe$species=="REG" & dataframe$area=="27.12.a.1" | dataframe$species=="REG" & dataframe$area=="27.12.a.2"
                       | dataframe$species=="REG" & dataframe$area=="27.12.a.3" | dataframe$species=="REG" & dataframe$area=="27.12.a.4" | dataframe$species=="REG" & dataframe$area=="27.12.b"
                       | dataframe$species=="REG" & dataframe$area=="27.12.c"] <- "reg.27.561214"

  dataframe$stock[dataframe$species=="REB" & dataframe$area=="27.1.a" | dataframe$species=="REB" & dataframe$area=="27.1.b" | dataframe$species=="REB" & dataframe$area=="27.2.a"
                       | dataframe$species=="REB" & dataframe$area=="27.2.a.1" | dataframe$species=="REB" & dataframe$area=="27.2.a.2" | dataframe$species=="REB" & dataframe$area=="27.2.b"
                       | dataframe$species=="REB" & dataframe$area=="27.2.b.1" | dataframe$species=="REB" & dataframe$area=="27.2.b.2"]<-"reb.27.1-2"

  ## Hier definier ich noch einmal die Fänge in der Banane als pelagische Fänge (falls dies bei DTS statt TM auftauchen sollte!)
  dataframe$stock[dataframe$species=="REB" & dataframe$area=="27.2.a.1" | dataframe$species=="REB" & dataframe$area=="27.2.b.1"]<-"reb.27.1-2 pel"



  ## Falls Fänge von RED im Pelagial sind:
  dataframe$stock[dataframe$species=="RED" & dataframe$area=="27.2.a.1" | dataframe$species=="RED" & dataframe$area=="27.2.b.1"]<-"reb.27.1-2 pel"

  ##Fänge in der Irmingersee (pelagisch)
  dataframe$stock[dataframe$species=="RED" & dataframe$area=="27.14.b.1" | dataframe$species=="REB" & dataframe$area=="27.14.b.1"]<-"reb.2127.dp"

  ## Demersaler Grönlandbestand
  dataframe$stock[dataframe$species=="REB" & dataframe$area=="27.14.b" | dataframe$species=="REB" & dataframe$area=="27.14.b.2"]<-"reb.27.14b dem"

  # D) Weit verbreitete Bestaende + weitere ----

  dataframe$stock[dataframe$species=="HER" & dataframe$area=="27.1.a" | dataframe$species=="HER" & dataframe$area=="27.1.b" | dataframe$species=="HER" & dataframe$area=="27.2.a"
                       | dataframe$species=="HER" & dataframe$area=="27.2.a.1" | dataframe$species=="HER" & dataframe$area=="27.2.a.2" | dataframe$species=="HER" & dataframe$area=="27.2.b"
                       | dataframe$species=="HER" & dataframe$area=="27.2.b.1" | dataframe$species=="HER" & dataframe$area=="27.2.b.2" | dataframe$species=="HER" & dataframe$area=="27.14.a"
                       | dataframe$species=="HER" & dataframe$area=="27.5.a" | dataframe$species=="HER" & dataframe$area=="27.5.a.1" | dataframe$species=="HER" & dataframe$area=="27.5.a.2"
                       | dataframe$species=="HER" & dataframe$area=="27.5.b" | dataframe$species=="HER" & dataframe$area=="27.5.b.1" | dataframe$species=="HER" & dataframe$area=="27.5.b.1.a"
                       | dataframe$species=="HER" & dataframe$area=="27.5.b.1.b" | dataframe$species=="HER" & dataframe$area=="27.5.b.2" | dataframe$species=="HER" & dataframe$area=="27.4.a"]<-"her.27.1-24a514a"

  dataframe$stock[dataframe$species=="MAC" & dataframe$area=="27.1.a" | dataframe$species=="MAC" & dataframe$area=="27.1.b" | dataframe$species=="MAC" & dataframe$area=="27.2.a"
                       | dataframe$species=="MAC" & dataframe$area=="27.2.a.1" | dataframe$species=="MAC" & dataframe$area=="27.2.a.2" | dataframe$species=="MAC" & dataframe$area=="27.2.b"
                       | dataframe$species=="MAC" & dataframe$area=="27.2.b.1" | dataframe$species=="MAC" & dataframe$area=="27.2.b.2" | dataframe$species=="MAC" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="MAC" & dataframe$area=="27.3.a.n" | dataframe$species=="MAC" & dataframe$area=="27.3.a.21" | dataframe$species=="MAC" & dataframe$area=="27.3.c.22"
                       | dataframe$species=="MAC" & dataframe$area=="27.3.b.23" | dataframe$species=="MAC" & dataframe$area=="27.3.d.24"| dataframe$species=="MAC" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="MAC" & dataframe$area=="27.3.d.26"| dataframe$species=="MAC" & dataframe$area=="27.3.d.27"| dataframe$species=="MAC" & dataframe$area=="27.3.d.28.2"
                       | dataframe$species=="MAC" & dataframe$area=="27.3.d.29"| dataframe$species=="MAC" & dataframe$area=="27.3.d.30"| dataframe$species=="MAC" & dataframe$area=="27.3.d.31"
                       | dataframe$species=="MAC" & dataframe$area=="27.3.d.32"| dataframe$species=="MAC" & dataframe$area=="27.4.a"| dataframe$species=="MAC" & dataframe$area=="27.4.b"
                       | dataframe$species=="MAC" & dataframe$area=="27.4.c" | dataframe$species=="MAC" & dataframe$area=="27.5.a" | dataframe$species=="MAC" & dataframe$area=="27.5.a.1"
                       | dataframe$species=="MAC" & dataframe$area=="27.5.a.2" | dataframe$species=="MAC" & dataframe$area=="27.5.b" | dataframe$species=="MAC" & dataframe$area=="27.5.b.1"
                       | dataframe$species=="MAC" & dataframe$area=="27.5.b.1.a" | dataframe$species=="MAC" & dataframe$area=="27.5.b.1.b" | dataframe$species=="MAC" & dataframe$area=="27.5.b.2"
                       | dataframe$species=="MAC" & dataframe$area=="27.6.a" | dataframe$species=="MAC" & dataframe$area=="27.6.b"| dataframe$species=="MAC" & dataframe$area=="27.6.b.1" | dataframe$species=="MAC" & dataframe$area=="27.6.b.2"
                       | dataframe$species=="MAC" & dataframe$area=="27.6.a" | dataframe$species=="MAC" & dataframe$area=="27.7.a" | dataframe$species=="MAC" & dataframe$area=="27.7.b"
                       | dataframe$species=="MAC" & dataframe$area=="27.7.c" | dataframe$species=="MAC" & dataframe$area=="27.7.c.1" | dataframe$species=="MAC" & dataframe$area=="27.7.c.2" | dataframe$species=="MAC" & dataframe$area=="27.7.d"
                       | dataframe$species=="MAC" & dataframe$area=="27.7.e" | dataframe$species=="MAC" & dataframe$area=="27.7.f" | dataframe$species=="MAC" & dataframe$area=="27.7.g"
                       | dataframe$species=="MAC" & dataframe$area=="27.7.h" | dataframe$species=="MAC" & dataframe$area=="27.7.j.1" | dataframe$species=="MAC" & dataframe$area=="27.7.j.2"
                       | dataframe$species=="MAC" & dataframe$area=="27.7.k" | dataframe$species=="MAC" & dataframe$area=="27.7.k.1" | dataframe$species=="MAC" & dataframe$area=="27.7.k.2"  | dataframe$species=="MAC" & dataframe$area=="27.8.a"
                       | dataframe$species=="MAC" & dataframe$area=="27.8.b" | dataframe$species=="MAC" & dataframe$area=="27.8.c" | dataframe$species=="MAC" & dataframe$area=="27.8.d.1"
                       | dataframe$species=="MAC" & dataframe$area=="27.8.d.2" | dataframe$species=="MAC" & dataframe$area=="27.8.e"| dataframe$species=="MAC" & dataframe$area=="27.8.e.1" | dataframe$species=="MAC" & dataframe$area=="27.8.e.2"
                       | dataframe$species=="MAC" & dataframe$area=="27.9.a" | dataframe$species=="MAC" & dataframe$area=="27.14.a" | dataframe$species=="MAC" & dataframe$area=="27.14.b"
                       | dataframe$species=="MAC" & dataframe$area=="27.14.b.1"| dataframe$species=="MAC" & dataframe$area=="27.14.b.2" | dataframe$species=="MAC" & dataframe$area=="27.3.a"
                       | dataframe$species=="MAC" & dataframe$area=="27.8.d"]<-"mac.27.nea"

  dataframe$stock[dataframe$species=="JAX" & dataframe$area=="27.8.a" | dataframe$species=="JAX" & dataframe$area=="27.8.b" | dataframe$species=="JAX" & dataframe$area=="27.8.c"
                       | dataframe$species=="JAX" & dataframe$area=="27.8.d" | dataframe$species=="JAX" & dataframe$area=="27.8.d.1" | dataframe$species=="JAX" & dataframe$area=="27.8.d.2" | dataframe$species=="JAX" & dataframe$area=="27.8.e.1"
                       | dataframe$species=="JAX" & dataframe$area=="27.8.e.2" | dataframe$species=="JAX" & dataframe$area=="27.2.a" | dataframe$species=="JAX" & dataframe$area=="27.2.a.1"
                       | dataframe$species=="JAX" & dataframe$area=="27.2.a.2" | dataframe$species=="JAX" & dataframe$area=="27.4.a" | dataframe$species=="JAX" & dataframe$area=="27.5.b"
                       | dataframe$species=="JAX" & dataframe$area=="27.5.b.1" | dataframe$species=="JAX" & dataframe$area=="27.5.b.1.a" | dataframe$species=="JAX" & dataframe$area=="27.5.b.1.b"
                       | dataframe$species=="JAX" & dataframe$area=="27.5.b.2" | dataframe$species=="JAX" & dataframe$area=="27.6.a" | dataframe$species=="JAX" & dataframe$area=="27.7.a"
                       | dataframe$species=="JAX" & dataframe$area=="27.7.b" | dataframe$species=="JAX" & dataframe$area=="27.7.c" | dataframe$species=="JAX" & dataframe$area=="27.7.c.1" | dataframe$species=="JAX" & dataframe$area=="27.7.c.2"
                       | dataframe$species=="JAX" & dataframe$area=="27.7.e" | dataframe$species=="JAX" & dataframe$area=="27.7.f" | dataframe$species=="JAX" & dataframe$area=="27.7.g"
                       | dataframe$species=="JAX" & dataframe$area=="27.7.h" | dataframe$species=="JAX" & dataframe$area=="27.7.j.1" | dataframe$species=="JAX" & dataframe$area=="27.7.j.2"
                       | dataframe$species=="JAX" & dataframe$area=="27.7.k" |dataframe$species=="JAX" & dataframe$area=="27.7.k.1" | dataframe$species=="JAX" & dataframe$area=="27.7.k.2"
                       | dataframe$species=="HOM" & dataframe$area=="27.8.a" | dataframe$species=="HOM" & dataframe$area=="27.8.b" | dataframe$species=="HOM" & dataframe$area=="27.8.c"
                       | dataframe$species=="HOM" & dataframe$area=="27.8.d.1" | dataframe$species=="HOM" & dataframe$area=="27.8.d.2" | dataframe$species=="HOM" & dataframe$area=="27.8.e.1"
                       | dataframe$species=="HOM" & dataframe$area=="27.8.e.2" | dataframe$species=="HOM" & dataframe$area=="27.2.a" | dataframe$species=="HOM" & dataframe$area=="27.2.a.1"
                       | dataframe$species=="HOM" & dataframe$area=="27.2.a.2" | dataframe$species=="HOM" & dataframe$area=="27.4.a" | dataframe$species=="HOM" & dataframe$area=="27.5.b"
                       | dataframe$species=="HOM" & dataframe$area=="27.5.b.1" | dataframe$species=="HOM" & dataframe$area=="27.5.b.1.a" | dataframe$species=="HOM" & dataframe$area=="27.5.b.1.b"
                       | dataframe$species=="HOM" & dataframe$area=="27.5.b.2" | dataframe$species=="HOM" & dataframe$area=="27.6.a" | dataframe$species=="HOM" & dataframe$area=="27.7.a"
                       | dataframe$species=="HOM" & dataframe$area=="27.7.b" | dataframe$species=="HOM" & dataframe$area=="27.7.c.1" | dataframe$species=="HOM" & dataframe$area=="27.7.c.2"
                       | dataframe$species=="HOM" & dataframe$area=="27.7.e" | dataframe$species=="HOM" & dataframe$area=="27.7.f" | dataframe$species=="HOM" & dataframe$area=="27.7.g"
                       | dataframe$species=="HOM" & dataframe$area=="27.7.h" | dataframe$species=="HOM" & dataframe$area=="27.7.j.1" | dataframe$species=="HOM" & dataframe$area=="27.7.j.2"
                       | dataframe$species=="HOM" & dataframe$area=="27.7.k.1" | dataframe$species=="HOM" & dataframe$area=="27.7.k.2" ]<-"hom.27.2a4a5b6a7a-ce-k8"

  dataframe$stock[dataframe$species=="WHB" & dataframe$area=="27.1.a" | dataframe$species=="WHB" & dataframe$area=="27.1.b" | dataframe$species=="WHB" & dataframe$area=="27.2.a"
                       | dataframe$species=="WHB" & dataframe$area=="27.2.a.1" | dataframe$species=="WHB" & dataframe$area=="27.2.a.2" | dataframe$species=="WHB" & dataframe$area=="27.2.b"
                       | dataframe$species=="WHB" & dataframe$area=="27.2.b.1" | dataframe$species=="WHB" & dataframe$area=="27.2.b.2" | dataframe$species=="WHB" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="WHB" & dataframe$area=="27.3.a.n" | dataframe$species=="WHB" & dataframe$area=="27.3.a" | dataframe$species=="WHB" & dataframe$area=="27.3.a.21" | dataframe$species=="WHB" & dataframe$area=="27.3.c.22"
                       | dataframe$species=="WHB" & dataframe$area=="27.3.b.23" | dataframe$species=="WHB" & dataframe$area=="27.3.d.24"| dataframe$species=="WHB" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="WHB" & dataframe$area=="27.3.d.26"| dataframe$species=="WHB" & dataframe$area=="27.3.d.27"| dataframe$species=="WHB" & dataframe$area=="27.3.d.28.2"
                       | dataframe$species=="WHB" & dataframe$area=="27.3.d.29"| dataframe$species=="WHB" & dataframe$area=="27.3.d.30"| dataframe$species=="WHB" & dataframe$area=="27.3.d.31"
                       | dataframe$species=="WHB" & dataframe$area=="27.3.d.32"| dataframe$species=="WHB" & dataframe$area=="27.4.a"| dataframe$species=="WHB" & dataframe$area=="27.4.b"
                       | dataframe$species=="WHB" & dataframe$area=="27.4.c" | dataframe$species=="WHB" & dataframe$area=="27.5.a" | dataframe$species=="WHB" & dataframe$area=="27.5.a.1"
                       | dataframe$species=="WHB" & dataframe$area=="27.5.a.2" | dataframe$species=="WHB" & dataframe$area=="27.5.b" | dataframe$species=="WHB" & dataframe$area=="27.5.b.1"
                       | dataframe$species=="WHB" & dataframe$area=="27.5.b.1.a" | dataframe$species=="WHB" & dataframe$area=="27.5.b.1.b" | dataframe$species=="WHB" & dataframe$area=="27.5.b.2"
                       | dataframe$species=="WHB" & dataframe$area=="27.6.b.1" | dataframe$species=="WHB" & dataframe$area=="27.6.b"| dataframe$species=="WHB" & dataframe$area=="27.6.b.2"
                       | dataframe$species=="WHB" & dataframe$area=="27.6.a" | dataframe$species=="WHB" & dataframe$area=="27.7.a" | dataframe$species=="WHB" & dataframe$area=="27.7.b"
                       | dataframe$species=="WHB" & dataframe$area=="27.7.c" | dataframe$species=="WHB" & dataframe$area=="27.7.c.1" | dataframe$species=="WHB" & dataframe$area=="27.7.c.2" | dataframe$species=="WHB" & dataframe$area=="27.7.d"
                       | dataframe$species=="WHB" & dataframe$area=="27.7.e" | dataframe$species=="WHB" & dataframe$area=="27.7.f" | dataframe$species=="WHB" & dataframe$area=="27.7.g"
                       | dataframe$species=="WHB" & dataframe$area=="27.7.h" | dataframe$species=="WHB" & dataframe$area=="27.7.j.1" | dataframe$species=="WHB" & dataframe$area=="27.7.j.2"
                       | dataframe$species=="WHB" & dataframe$area=="27.7.k" | dataframe$species=="WHB" & dataframe$area=="27.7.k.1" | dataframe$species=="WHB" & dataframe$area=="27.7.k.2"  | dataframe$species=="WHB" & dataframe$area=="27.8.a"
                       | dataframe$species=="WHB" & dataframe$area=="27.8.b" | dataframe$species=="WHB" & dataframe$area=="27.8.c" | dataframe$species=="WHB" & dataframe$area=="27.8.d.1"
                       | dataframe$species=="WHB" & dataframe$area=="27.8.d.2" | dataframe$species=="WHB" & dataframe$area=="27.8.e.1" | dataframe$species=="WHB" & dataframe$area=="27.8.e.2"
                       | dataframe$species=="WHB" & dataframe$area=="27.9.a"  | dataframe$species=="WHB" & dataframe$area=="27.9.b.1" | dataframe$species=="WHB" & dataframe$area=="27.9.b.2"
                       | dataframe$species=="WHB" & dataframe$area=="27.12.a" | dataframe$species=="WHB" & dataframe$area=="27.12.a.1"  | dataframe$species=="WHB" & dataframe$area=="27.12.a.2" | dataframe$species=="WHB" & dataframe$area=="27.12.a.3"
                       | dataframe$species=="WHB" & dataframe$area=="27.12.a.4"  | dataframe$species=="WHB" & dataframe$area=="27.12.b" | dataframe$species=="WHB" & dataframe$area=="27.12.c"
                       | dataframe$species=="WHB" & dataframe$area=="27.14.a" | dataframe$species=="WHB" & dataframe$area=="27.14.b" | dataframe$species=="WHB" & dataframe$area=="27.14.b.1"
                       | dataframe$species=="WHB" & dataframe$area=="27.14.b.2" ]<-"whb.27.1-91214"

  dataframe$stock[dataframe$species=="HKE" & dataframe$area=="27.3.a.n" | dataframe$species=="HKE" & dataframe$area=="27.3.a.21" | dataframe$species=="HKE" & dataframe$area=="27.4.a"
                       | dataframe$species=="HKE" & dataframe$area=="27.4.b" | dataframe$species=="HKE" & dataframe$area=="27.4.c" | dataframe$species=="HKE" & dataframe$area=="27.6.a"
                       | dataframe$species=="HKE" & dataframe$area=="27.6.b" | dataframe$species=="HKE" & dataframe$area=="27.6.b.1" | dataframe$species=="HKE" & dataframe$area=="27.6.b.2" | dataframe$species=="HKE" & dataframe$area=="27.7.a"
                       | dataframe$species=="HKE" & dataframe$area=="27.7.b" | dataframe$species=="HKE" & dataframe$area=="27.7.c" | dataframe$species=="HKE" & dataframe$area=="27.7.c.1" | dataframe$species=="HKE" & dataframe$area=="27.7.c.2"
                       | dataframe$species=="HKE" & dataframe$area=="27.7.d" | dataframe$species=="HKE" & dataframe$area=="27.7.e" | dataframe$species=="HKE" & dataframe$area=="27.7.f"
                       | dataframe$species=="HKE" & dataframe$area=="27.7.g" | dataframe$species=="HKE" & dataframe$area=="27.7.h" | dataframe$species=="HKE" & dataframe$area=="27.7.j.1"
                       | dataframe$species=="HKE" & dataframe$area=="27.7.j.2" | dataframe$species=="HKE" & dataframe$area=="27.7.k.1" | dataframe$species=="HKE" & dataframe$area=="27.7.k.2"
                       | dataframe$species=="HKE" & dataframe$area=="27.8.a" | dataframe$species=="HKE" & dataframe$area=="27.8.b" | dataframe$species=="HKE" & dataframe$area=="27.8.d.1"
                       | dataframe$species=="HKE" & dataframe$area=="27.8.d.2" | dataframe$species=="HKE" & dataframe$area=="27.3.a.20" | dataframe$species=="HKE" & dataframe$area=="27.3.a"
                       | dataframe$species=="HKE" & dataframe$area=="27.2.a.2" | dataframe$species=="HKE" & dataframe$area=="27.3.d.24" | dataframe$species=="HKE" & dataframe$area=="27.3.c.22" ]<-"hke.27.3a46-8abd"

  dataframe$stock[dataframe$species=="HER" & dataframe$area=="27.7.a" | dataframe$species=="HER" & dataframe$area=="27.7.g" | dataframe$species=="HER" & dataframe$area=="27.7.h"
                       | dataframe$species=="HER" & dataframe$area=="27.7.j.1" | dataframe$species=="HER" & dataframe$area=="27.7.j.2" | dataframe$species=="HER" & dataframe$area=="27.7.k.1"
                       | dataframe$species=="HER" & dataframe$area=="27.7.k.2"] <-"her.27.irls"

  dataframe$stock[dataframe$species=="HER" & dataframe$area=="27.6.a" | dataframe$species=="HER" & dataframe$area=="27.7.b" | dataframe$species=="HER" & dataframe$area=="27.7.c.1"
                       | dataframe$species=="HER" & dataframe$area=="27.7.c.2"]<-"her.27.6a7bc"

  dataframe$stock[dataframe$species=="ARU" & dataframe$area=="27.3.a.n" | dataframe$species=="ARU" & dataframe$area=="27.3.a.21" | dataframe$species=="ARU" & dataframe$area=="27.4.a"
                       | dataframe$species=="ARU" & dataframe$area=="27.4.b" | dataframe$species=="ARU" & dataframe$area=="27.4.c" |  dataframe$species=="ARU" & dataframe$area=="27.1.a"
                       | dataframe$species=="ARU" & dataframe$area=="27.1.b" | dataframe$species=="ARU" & dataframe$area=="27.2.a" | dataframe$species=="ARU" & dataframe$area=="27.2.a.1"
                       | dataframe$species=="ARU" & dataframe$area=="27.2.a.2" | dataframe$species=="ARU" & dataframe$area=="27.2.b" | dataframe$species=="ARU" & dataframe$area=="27.2.b.1"
                       | dataframe$species=="ARU" & dataframe$area=="27.2.b.2" | dataframe$species=="ARU" & dataframe$area=="27.3.a.20"]<-"aru.27.123a4" # Argentina silus (Greater silver smelt)

  dataframe$stock[dataframe$species=="ARU" & dataframe$area=="27.5.b" | dataframe$species=="ARU" & dataframe$area=="27.5.b.1" | dataframe$species=="ARU" & dataframe$area=="27.5.b.1.a"
                       | dataframe$species=="ARU" & dataframe$area=="27.5.b.1.b" | dataframe$species=="ARU" & dataframe$area=="27.5.b.2"
                       | dataframe$species=="ARU" & dataframe$area=="27.6.a"]<-"aru.27.5b6a" # Argentina silus (Greater silver smelt)

  dataframe$stock[dataframe$species=="PIL" & dataframe$area=="27.7.a" | dataframe$species=="PIL" & dataframe$area=="27.7.b"
                       | dataframe$species=="PIL" & dataframe$area=="27.7.c.1" | dataframe$species=="PIL" & dataframe$area=="27.7.c.2" | dataframe$species=="PIL" & dataframe$area=="27.7.d"
                       | dataframe$species=="PIL" & dataframe$area=="27.7.e" | dataframe$species=="PIL" & dataframe$area=="27.7.f" | dataframe$species=="PIL" & dataframe$area=="27.7.g"
                       | dataframe$species=="PIL" & dataframe$area=="27.7.h" | dataframe$species=="PIL" & dataframe$area=="27.7.j.1" | dataframe$species=="PIL" & dataframe$area=="27.7.j.2"
                       | dataframe$species=="PIL" & dataframe$area=="27.7.k.1" | dataframe$species=="PIL" & dataframe$area=="27.7.k.2"]<-"pil.27.7"

  dataframe$stock[dataframe$species=="PIL" & dataframe$area=="27.8.d.1" | dataframe$species=="PIL" & dataframe$area=="27.8.d.2" |
                         dataframe$species=="PIL" & dataframe$area=="27.8.a" | dataframe$species=="PIL" & dataframe$area=="27.8.b"]<-"pil.27.8abd"

  dataframe$stock[dataframe$species=="PIL" & dataframe$area=="27.8.c" | dataframe$species=="PIL" & dataframe$area=="27.9.a" ]<-"pil.27.8c9a"

  dataframe$stock[dataframe$species=="BLI" & dataframe$area=="27.3.a.n" | dataframe$species=="BLI" & dataframe$area=="27.3.a.21" | dataframe$species=="BLI" & dataframe$area=="27.4.a"
                       | dataframe$species=="BLI" & dataframe$area=="27.1.a" | dataframe$species=="BLI" & dataframe$area=="27.1.b" | dataframe$species=="BLI" & dataframe$area=="27.2.a"
                       | dataframe$species=="BLI" & dataframe$area=="27.2.a.1" | dataframe$species=="BLI" & dataframe$area=="27.2.a.2" | dataframe$species=="BLI" & dataframe$area=="27.2.b"
                       | dataframe$species=="BLI" & dataframe$area=="27.2.b.1" | dataframe$species=="BLI" & dataframe$area=="27.2.b.2" | dataframe$species=="BLI" & dataframe$area=="27.8.a"
                       | dataframe$species=="BLI" & dataframe$area=="27.8.b" | dataframe$species=="BLI" & dataframe$area=="27.8.c" | dataframe$species=="BLI" & dataframe$area=="27.8.d.1"
                       | dataframe$species=="BLI" & dataframe$area=="27.8.d.2" | dataframe$species=="BLI" & dataframe$area=="27.8.e.1" | dataframe$species=="BLI" & dataframe$area=="27.8.e.2"
                       | dataframe$species=="BLI" & dataframe$area=="27.9.a"  | dataframe$species=="BLI" & dataframe$area=="27.9.b.1" | dataframe$species=="BLI" & dataframe$area=="27.9.b.2"
                       | dataframe$species=="BLI" & dataframe$area=="27.12.a.1"  | dataframe$species=="BLI" & dataframe$area=="27.12.a.2" | dataframe$species=="BLI" & dataframe$area=="27.12.a.3"
                       | dataframe$species=="BLI" & dataframe$area=="27.12.a.4"  | dataframe$species=="BLI" & dataframe$area=="27.12.b" | dataframe$species=="BLI" & dataframe$area=="27.12.c"
                       | dataframe$species=="BLI" & dataframe$area=="27.3.a.20" | dataframe$species=="BLI" & dataframe$area=="27.4.b" | dataframe$species=="BLI" & dataframe$area=="21.1.c" ]<-"bli.27.nea"

  dataframe$stock[dataframe$species=="BLI" & dataframe$area %in%
                         unique(grep("27.7|27.5.b|27.6",dataframe$area,value = T))] <- "bli.27.5b67"

  dataframe$stock[ dataframe$species=="LIN" & dataframe$area=="27.6.b.1" | dataframe$species=="LIN" & dataframe$area=="27.6.b" | dataframe$species=="LIN" & dataframe$area=="27.6.b.2" | dataframe$species=="LIN" & dataframe$area=="27.3.a.20"
                        | dataframe$species=="LIN" & dataframe$area=="27.6.a" | dataframe$species=="LIN" & dataframe$area=="27.7.a" | dataframe$species=="LIN" & dataframe$area=="27.7.b"
                        | dataframe$species=="LIN" & dataframe$area=="27.7.c" | dataframe$species=="LIN" & dataframe$area=="27.7.c.1" | dataframe$species=="LIN" & dataframe$area=="27.7.c.2" | dataframe$species=="LIN" & dataframe$area=="27.7.d"
                        | dataframe$species=="LIN" & dataframe$area=="27.7.e" | dataframe$species=="LIN" & dataframe$area=="27.7.f" | dataframe$species=="LIN" & dataframe$area=="27.7.g"
                        | dataframe$species=="LIN" & dataframe$area=="27.7.h" | dataframe$species=="LIN" & dataframe$area=="27.7.j.1" | dataframe$species=="LIN" & dataframe$area=="27.7.j.2"
                        | dataframe$species=="LIN" & dataframe$area=="27.7.k" |dataframe$species=="LIN" & dataframe$area=="27.7.k.1" | dataframe$species=="LIN" & dataframe$area=="27.7.k.2"  | dataframe$species=="LIN" & dataframe$area=="27.8.a"
                        | dataframe$species=="LIN" & dataframe$area=="27.8.b" | dataframe$species=="LIN" & dataframe$area=="27.8.c" | dataframe$species=="LIN" & dataframe$area=="27.8.d.1"
                        | dataframe$species=="LIN" & dataframe$area=="27.8.d.2" | dataframe$species=="LIN" & dataframe$area=="27.8.e.1" | dataframe$species=="LIN" & dataframe$area=="27.8.e.2"
                        | dataframe$species=="LIN" & dataframe$area=="27.9.a"  | dataframe$species=="LIN" & dataframe$area=="27.9.b.1" | dataframe$species=="LIN" & dataframe$area=="27.9.b.2"
                        | dataframe$species=="LIN" & dataframe$area=="27.3.a.n" | dataframe$species=="LIN" & dataframe$area=="27.3.a.21" | dataframe$species=="LIN" & dataframe$area=="27.4.a"
                        | dataframe$species=="LIN" & dataframe$area=="27.12.a.1"  | dataframe$species=="LIN" & dataframe$area=="27.12.a.2" | dataframe$species=="LIN" & dataframe$area=="27.12.a.3"
                        | dataframe$species=="LIN" & dataframe$area=="27.12.a.4"  | dataframe$species=="LIN" & dataframe$area=="27.12.b" | dataframe$species=="LIN" & dataframe$area=="27.12.c"
                        | dataframe$species=="LIN" & dataframe$area=="27.14.a" | dataframe$species=="LIN" & dataframe$area=="27.14.b" | dataframe$species=="LIN" & dataframe$area=="27.14.b.1"
                        | dataframe$species=="LIN" & dataframe$area=="27.14.b.2" | dataframe$species=="LIN" & dataframe$area=="27.3.a"| dataframe$species=="LIN" & dataframe$area=="27.3.c.22"
                        | dataframe$species=="LIN" & dataframe$area=="27.4.b" | dataframe$species=="LIN" & dataframe$area=="27.4.c"]<-"lin.27.3a4a6-91214"

  dataframe$stock[dataframe$species=="USK" & dataframe$area=="27.3.a.n" | dataframe$species=="USK" & dataframe$area=="27.3.a" | dataframe$species=="USK" & dataframe$area=="27.3.a.21" | dataframe$species=="USK" & dataframe$area=="27.4.a"
                       | dataframe$species=="USK" & dataframe$area=="27.4.b" | dataframe$species=="USK" & dataframe$area=="27.4.c" | dataframe$species=="USK" & dataframe$area=="27.6.a"
                       | dataframe$species=="USK" & dataframe$area=="27.5.b" | dataframe$species=="USK" & dataframe$area=="27.5.b.1" | dataframe$species=="USK" & dataframe$area=="27.5.b.1.a"
                       | dataframe$species=="USK" & dataframe$area=="27.5.b.1.b" | dataframe$species=="USK" & dataframe$area=="27.5.b.2" | dataframe$species=="USK" & dataframe$area=="27.7.a"
                       | dataframe$species=="USK" & dataframe$area=="27.7.b" | dataframe$species=="USK" & dataframe$area=="27.7.c.1" | dataframe$species=="USK" & dataframe$area=="27.7.c.2"
                       | dataframe$species=="USK" & dataframe$area=="27.7.d" | dataframe$species=="USK" & dataframe$area=="27.7.e" | dataframe$species=="USK" & dataframe$area=="27.7.f"
                       | dataframe$species=="USK" & dataframe$area=="27.7.g" | dataframe$species=="USK" & dataframe$area=="27.7.h" | dataframe$species=="USK" & dataframe$area=="27.7.j.1"
                       | dataframe$species=="USK" & dataframe$area=="27.7.j.2" | dataframe$species=="USK" & dataframe$area=="27.7.k.1" | dataframe$species=="USK" & dataframe$area=="27.7.k.2"
                       | dataframe$species=="USK" & dataframe$area=="27.8.a" | dataframe$species=="USK" & dataframe$area=="27.8.b" | dataframe$species=="USK" & dataframe$area=="27.8.c"
                       | dataframe$species=="USK" & dataframe$area=="27.8.d.1" | dataframe$species=="USK" & dataframe$area=="27.8.d.2" | dataframe$species=="USK" & dataframe$area=="27.8.e.1"
                       | dataframe$species=="USK" & dataframe$area=="27.8.e.2" | dataframe$species=="USK" & dataframe$area=="27.9.a"  | dataframe$species=="USK" & dataframe$area=="27.9.b.1"
                       | dataframe$species=="USK" & dataframe$area=="27.9.b.2"  | dataframe$species=="USK" & dataframe$area=="27.12.b"
                       | dataframe$species=="USK" & dataframe$area=="27.3.a.20" ]<-"usk.27.3a45b6a7-912b"

  dataframe$stock[dataframe$species=="RNG" & dataframe$area=="27.1." |dataframe$species=="RNG" & dataframe$area=="27.1.a" | dataframe$species=="RNG" & dataframe$area=="27.1.b" | dataframe$species=="RNG" & dataframe$area=="27.2.a"
                       | dataframe$species=="RNG" & dataframe$area=="27.2.a.1" | dataframe$species=="RNG" & dataframe$area=="27.2.a.2" | dataframe$species=="RNG" & dataframe$area=="27.2.b"
                       | dataframe$species=="RNG" & dataframe$area=="27.2.b.1" | dataframe$species=="RNG" & dataframe$area=="27.2.b.2" | dataframe$species=="RNG" & dataframe$area=="27.4.a"
                       | dataframe$species=="RNG" & dataframe$area=="27.4.b" | dataframe$species=="RNG" & dataframe$area=="27.4.c" | dataframe$species=="RNG" & dataframe$area=="27.5.a"
                       | dataframe$species=="RNG" & dataframe$area=="27.8.a" | dataframe$species=="RNG" & dataframe$area=="27.8.b" | dataframe$species=="RNG" & dataframe$area=="27.8.c"
                       | dataframe$species=="RNG" & dataframe$area=="27.8.d.1" | dataframe$species=="RNG" & dataframe$area=="27.8.d.2" | dataframe$species=="RNG" & dataframe$area=="27.8.e.1"
                       | dataframe$species=="RNG" & dataframe$area=="27.8.e.2" | dataframe$species=="RNG" & dataframe$area=="27.9.a"  | dataframe$species=="RNG" & dataframe$area=="27.9.b.1"
                       | dataframe$species=="RNG" & dataframe$area=="27.9.b.2" | dataframe$species=="RNG" & dataframe$area=="27.14.a" | dataframe$species=="RNG" & dataframe$area=="27.14.b.2"
                       | dataframe$species=="RNG" & dataframe$area=="27.5.a.2" | dataframe$species=="RNG" & dataframe$area=="27.14.b.1"| dataframe$species=="RNG" & dataframe$area=="27.14.b"]<-"rng.27.1245a8914ab"

  dataframe$stock[dataframe$species=="CAP" & dataframe$area=="27.2.a" | dataframe$species=="CAP" & dataframe$area=="27.2.a.1" | dataframe$species=="CAP" & dataframe$area=="27.2.a.2"
                       | dataframe$species=="CAP" & dataframe$area=="27.14.a" | dataframe$species=="CAP" & dataframe$area=="27.14.b" | dataframe$species=="CAP" & dataframe$area=="27.14.b.1"
                       | dataframe$species=="CAP" & dataframe$area=="27.14.b.2" | dataframe$species=="CAP" & dataframe$area=="27.5.a" | dataframe$species=="CAP" & dataframe$area=="27.5.a.1"
                       | dataframe$species=="CAP" & dataframe$area=="27.5.a.2" | dataframe$species=="CAP" & dataframe$area=="27.5.b" | dataframe$species=="CAP" & dataframe$area=="27.5.b.1"
                       | dataframe$species=="CAP" & dataframe$area=="27.5.b.1.a" | dataframe$species=="CAP" & dataframe$area=="27.5.b.1.b" | dataframe$species=="CAP" & dataframe$area=="27.5.b.2"
                       ]<-"cap.27.2a514"

  dataframe$stock[dataframe$species=="GUR" & dataframe$area=="27.3.a.n" | dataframe$species=="GUR" & dataframe$area=="27.3.a.21" | dataframe$species=="GUR" & dataframe$area=="27.3.c.22"
                       | dataframe$species=="GUR" & dataframe$area=="27.3.b.23" | dataframe$species=="GUR" & dataframe$area=="27.3.d.24"| dataframe$species=="GUR" & dataframe$area=="27.3.d.25"
                       | dataframe$species=="GUR" & dataframe$area=="27.3.d.26"| dataframe$species=="GUR" & dataframe$area=="27.3.d.27"| dataframe$species=="GUR" & dataframe$area=="27.3.d.28.2"
                       | dataframe$species=="GUR" & dataframe$area=="27.3.d.29"| dataframe$species=="GUR" & dataframe$area=="27.3.d.30"| dataframe$species=="GUR" & dataframe$area=="27.3.d.31"
                       | dataframe$species=="GUR" & dataframe$area=="27.3.d.32"| dataframe$species=="GUR" & dataframe$area=="27.4.a"| dataframe$species=="GUR" & dataframe$area=="27.4.b"
                       | dataframe$species=="GUR" & dataframe$area=="27.4.c" | dataframe$species=="GUR" & dataframe$area=="27.5.a" | dataframe$species=="GUR" & dataframe$area=="27.5.a.1"
                       | dataframe$species=="GUR" & dataframe$area=="27.5.a.2" | dataframe$species=="GUR" & dataframe$area=="27.5.b" | dataframe$species=="GUR" & dataframe$area=="27.5.b.1"
                       | dataframe$species=="GUR" & dataframe$area=="27.5.b.1.a" | dataframe$species=="GUR" & dataframe$area=="27.5.b.1.b" | dataframe$species=="GUR" & dataframe$area=="27.5.b.2"
                       | dataframe$species=="GUR" & dataframe$area=="27.6.b"| dataframe$species=="GUR" & dataframe$area=="27.6.b.1" | dataframe$species=="GUR" & dataframe$area=="27.6.b.2" | dataframe$species=="GUR" & dataframe$area=="27.3.a.20"
                       | dataframe$species=="GUR" & dataframe$area=="27.6.a" | dataframe$species=="GUR" & dataframe$area=="27.7.a" | dataframe$species=="GUR" & dataframe$area=="27.7.b"
                       | dataframe$species=="GUR" & dataframe$area=="27.7.c.1" | dataframe$species=="GUR" & dataframe$area=="27.7.c.2" | dataframe$species=="GUR" & dataframe$area=="27.7.d"
                       | dataframe$species=="GUR" & dataframe$area=="27.7.e" | dataframe$species=="GUR" & dataframe$area=="27.7.f" | dataframe$species=="GUR" & dataframe$area=="27.7.g"
                       | dataframe$species=="GUR" & dataframe$area=="27.7.h" | dataframe$species=="GUR" & dataframe$area=="27.7.j.1" | dataframe$species=="GUR" & dataframe$area=="27.7.j.2"
                       | dataframe$species=="GUR" & dataframe$area=="27.7.k.1" | dataframe$species=="GUR" & dataframe$area=="27.7.k.2"  | dataframe$species=="GUR" & dataframe$area=="27.8.a"
                       | dataframe$species=="GUR" & dataframe$area=="27.8.b" | dataframe$species=="GUR" & dataframe$area=="27.8.c" | dataframe$species=="GUR" & dataframe$area=="27.8.d.1"
                       | dataframe$species=="GUR" & dataframe$area=="27.8.d.2" | dataframe$species=="GUR" & dataframe$area=="27.8.e.1" | dataframe$species=="GUR" & dataframe$area=="27.8.e.2"
                       ]<-"gur.27.3-8"


  ### Rectangles are unknown, therefore for NEP and SAN, the EU TAC Boundaries are used

  # nep.3a
  dataframe$stock[dataframe$species=="NEP" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.3.a")))]<-"nep.3a"

  # nep.EU2a4 and # nep.NOR4
  dataframe$stock[dataframe$species=="NEP" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.2.a|27.4")))]<-"nep.2a4"

  # nep.5bc6
  dataframe$stock[dataframe$species=="NEP" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.5.b|27.5.c|27.6")))]<-"nep.5bc6"

  # nep.7
  dataframe$stock[dataframe$species=="NEP" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.7")))]<-"nep.7"
  # nep.8abde
  dataframe$stock[dataframe$species == "NEP" & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[c(1:2,4:5)],collapse = "|"),dataframe$area,value = T))]<-"nep.8abde"
  # nep.8c
  dataframe$stock[dataframe$species=="NEP" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.8.c")))]<-"nep.8c"

  ## SAN in 2a, 3a and 4
  dataframe$stock[dataframe$species=="SAN" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.2.a|27.3.a|27.4")))]<-"san.27.2a3a4"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # san.27.6a
  # Sandeel (Ammodytes spp.) in Division 6.a (West of Scotland)
  dataframe$stock[dataframe$species %in% c("SAN") & dataframe$area %in%
                         unique(grep(paste0("27.6.a"),dataframe$area,value = T))]<-"san.27.6a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # san.sa.6
  # Sandeel (Ammodytes spp.) in subdivisions 20-22, Sandeel Area 6 (Kattegat)
  dataframe$stock[dataframe$species %in% c("SAN") & dataframe$area %in%
                         unique(grep(paste0("27.3.a.",20:22,collapse = "|"),dataframe$area,value = T))]<-"san.sa.6"


  # E) Elasmobranchier ----
  # First, all unknown rays (they may be overwritten later)
  ray_specs <- c("RJC","JAI","RJM","RJN","RJR","RJH","BHY","JAR", "RJI", "RJF", "JRS", "RJE", "RJU", "RJY", "JAY" ,"JAM", "JRG", "RJK", "BMU"
                 ,"JAZ" ,"BZB", "RJL", "BZM", "DPV", "RJS","RJT", "RJO","RRW", "BAM" ,"BEA", "BYF", "JAF", "JRB" ,"SRR" ,"RJA","JRM", "RJB", "JFV"
                 ,"RJG", "JAD" ,"JAL", "JAT", "RJD" ,"JAK" ,"JRD" ,"RFS" ,"JRC", "JRX", "JAV", "JRR","RHJ", "BYD" ,"BVS" ,"BVU")

  # raj.27.1012
  # Rays and skates (Rajidae) (mainly thornback ray (Raja clavata)) in subareas 10 and 12 (Azores grounds and north of Azores)
  dataframe$stock[dataframe$species %in% c(ray_specs) & dataframe$area %in%
                         unique(grep("27.10|27.12",dataframe$area,value = T))] <- "raj.27.1012"

  # raj.27.89a
  # Rays and skates (Rajidae) in Subarea 8 and Division 9.a (Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c(ray_specs) & dataframe$area %in%
                         unique(grep("27.8|27.9.a",dataframe$area,value = T))] <- "raj.27.89a"

  # "raj.27.3a47d"###
  # Rays and skates (Rajidae) in Subarea 4 and in divisions 3.a and 7.d (North Sea, Skagerrak, Kattegat, and eastern English Channel)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species %in% c(ray_specs) & dataframe$area %in%
                         unique(grep("27.4|27.3.a|27.7.d",dataframe$area,value = T))]<-"raj.27.3a47d"

  # "raj.27.67a-ce-h"
  # Rays and skates (Rajidae) in Subarea 6 and divisions 7.a-c and 7.e-h (Rockall and West of Scotland, southern Celtic Seas, western English Channel)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species %in% c(ray_specs) & dataframe$area %in%
                         unique(grep("27.6",dataframe$area,value = T)) |
                         dataframe$species %in% c(ray_specs) & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[c(1:3,5:8)],collapse = "|"),dataframe$area,value = T))]<-"raj.27.67a-ce-h"


  dataframe$stock[dataframe$species=="DGS"] <-"dgs.27.nea" # Dornhai

  dataframe$stock[dataframe$species=="RJH" & dataframe$area=="27.4.a" | dataframe$species=="RJH" & dataframe$area=="27.6.b.1" | dataframe$species=="RJH" & dataframe$area=="27.6.b"
                       | dataframe$species=="RJH" & dataframe$area=="27.6.b.2" | dataframe$species=="RJH" & dataframe$area=="27.6.a"] <- "rjh.27.4a6" ##Blonde ray

  dataframe$stock[dataframe$species=="RJH" & dataframe$area=="27.4.c" | dataframe$species=="RJH" & dataframe$area=="27.4.b" | dataframe$species=="RJH" & dataframe$area=="27.7.d" ] <- "rjh.27.4c7d" ##Blonde ray

  dataframe$stock[dataframe$species=="RJN" & dataframe$area=="27.3.a.n" | dataframe$species=="RJN" & dataframe$area=="27.3.a.21" | dataframe$species=="RJN" & dataframe$area=="27.4.a"
                       | dataframe$species=="RJN" & dataframe$area=="27.4.b" | dataframe$species=="RJN" & dataframe$area=="27.4.c" | dataframe$species=="RJN" & dataframe$area=="27.3.a.20" ] <- "rjn.27.3a4" ## Kuckucksrochen

  dataframe$stock[dataframe$species=="SYC" & dataframe$area=="27.3.a.n" | dataframe$species=="SYC" & dataframe$area=="27.3.a.21" | dataframe$species=="SYC" & dataframe$area=="27.4.a"
                       | dataframe$species=="SYC" & dataframe$area=="27.4.b" | dataframe$species=="SYC" & dataframe$area=="27.4.c" | dataframe$species=="SYC" & dataframe$area=="27.7.d"
                       | dataframe$species=="SYC" & dataframe$area=="27.3.a.20"] <- "syc.27.3a47d" ## Kleingefleckter Katdataframeenhai

  dataframe$stock[dataframe$species=="RJC" & dataframe$area=="27.3.a.n" | dataframe$species=="RJC" & dataframe$area=="27.3.a.21" | dataframe$species=="RJC" & dataframe$area=="27.4.a"
                       | dataframe$species=="RJC" & dataframe$area=="27.4.b" | dataframe$species=="RJC" & dataframe$area=="27.4.c" | dataframe$species=="RJC" & dataframe$area=="27.7.d"
                       | dataframe$species=="RJC" & dataframe$area=="27.3.a.20"] <- "rjc.27.3a47d" ## Nagelrochen



  #First - ANF #
  # in 27.2.a.2, 27.6.b.2, 27.7.c.2, 27.7.k.2, 27.1.b, 27.3.a
  #27.3.a an 27.6.b belong to stock anf.27.3a46, and go there
  # This is an official ICES-stock #
  #  ANF in 7 and 8.a-b and 8.d is an own stock
  dataframe$stock[dataframe$species=="ANF" & dataframe$area %in%
                         unique(grep("27.7|27.8.a|27.8.b|27.8.d",dataframe$area,value = T))]<-"anf.27.78abd"


  # The rest is "northern" anglerfish
  # This is an official ICES-stock #
  dataframe$stock[dataframe$species=="ANF" & dataframe$area=="27.2.a.2" |dataframe$species=="ANF" & dataframe$area=="27.2.a.1" |
                         dataframe$species=="ANF" & dataframe$area=="27.2.a" |dataframe$species=="ANF" & dataframe$area=="27.1.b" |
                         dataframe$species=="MON" & dataframe$area=="27.2.a.2" | dataframe$species=="MON" & dataframe$area=="27.2.a.1"
                       | dataframe$species=="MON" & dataframe$area=="27.2.a"]<-"anf.27.1-2"


  # in 27.7 it belongs to bli.27.5b67
  # This is an official ICES-stock #
  dataframe$stock[dataframe$species=="ANF" & dataframe$area %in%
                         unique(grep("27.6|27.5.b|27.7",dataframe$area,value = T))]<-"bli.27.5b67"


  # Next - COD
  # in 27.3.a and 27.5.a.2
  # First is added to stock cod.27.21, second to
  # Second is an own stock (Iceland grounds)
  dataframe$stock[dataframe$species=="COD" & dataframe$area=="27.5.a.2" ]<-"cod.27.5a"


  # Next is CYO ###
  # Portuguese dogfish
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="CYO"]<-"cyo.27.nea"

  # Next - DAB
  # In 27.3.a
  # This belongs to the stock dab.27.3a4


  # Next - FLE
  # In 27.3.a
  # This belongs to the stock fle.27.3a4


  # in 27.3.d.28.2 it belongs to bwq.27.2628
  dataframe$stock[dataframe$species=="FLE" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.3.d.28")))]<-"bwq.27.2628"


  # Next is GUQ ###
  # This is leafscale gulper shark
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="GUQ"]<-"guq.27.nea"

  # Next - HAD
  # In 27.3.a and 27.3.c.22
  # Treat as part orf North sea/ skagerrag stock had.27.46a20
  # In 27.14.b.2
  # Treat as part of Iceland/Greenland Stock had.27.5.a
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="HAD" & dataframe$area=="27.5.a"]<-"had.27.5a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="HAD" & dataframe$area=="27.5.b"]<-"had.27.5b"



  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="HAD" & dataframe$area %in%
                         unique(grep(paste0("27.7",letters[2:11],collapse = "|"),dataframe$area,value = T))]<-"had.27.7b-k"


  # Next - HKE
  # In 27.3.a, 27.2.a.2 and 27.3.d.24
  # These are all part of hke.27.3a46-8abd
  dataframe$stock[dataframe$species=="HKE" & dataframe$area %in%
                         unique(grep("27.5|27.2|27.3|27.7",dataframe$area,value = T))] <- "hke.27.3a46-8abd"



  # Next - LIN
  # In 27.3.a, 27.3.c.22, 27.4.b und 27.4.c
  # This is ling - They all belong to lin.27.3a4a6-91214

  # LIN in Iceland are separate stocks
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="LIN" & dataframe$area %in%
                         c("27.5.a")  ]<-"lin.27.5a"
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="LIN" & dataframe$area %in%
                         c("27.5.b")  ]<-"lin.27.5b"

  # Next - NOP
  # This is Norway Pout
  # 27.6.a West of Scotland, Celtic Sea
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="NOP" & dataframe$area=="27.6.a"]<-"nop.27.6a"

  # Next - ORY
  # This is orange roughy
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="ORY"]<-"ory.27.nea"


  # Next - POK
  # In 27.3.a, 27.3.d.24, 27.3.c.22, 27.14.b.2
  # The first three a part of pok.27.3a46
  # The last is most likely part of the Iceland stock pok.27.5a
  dataframe$stock[dataframe$species=="POK" & dataframe$area=="27.14.b.2" ]<-"pok.27.5a
  "
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="POK" & dataframe$area=="27.5.b" ]<-"pok.27.5b"

  # Next - POL
  # IN 27.3.a, 27.2.a.2, 27.3.c.22, 27.1.b
  # They are all assumed to belong to the same stock pol.27.3a4

  # in 27.6 and 27.7, they are a separate stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="POL" & dataframe$area %in%
                         unique(grep("27.6|27.7",dataframe$area,value = T)) ]<-"pol.27.67"


  # Next - POR
  # This is porbeagle (Lamna nasus)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="POR"]<-"por.27.nea"

  # Next - RJH
  # This is blonde ray and belongs to elasmobranchs
  # In 27.4.b
  # There is a stock for 27.4.a and for 27.4.c.
  # They are all caught pretty south and added to south-western stock rjh.27.4c7d

  # in 27.7 - this is celtic seas stock
  # Roundnose grenadier (Coryphaenoides rupestris) in subareas 6-7 and divisions 5.b and 12.b (Celtic Seas and the English Channel, Faroes grounds, and western Hatton Bank)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="RNG" & dataframe$area %in%
                         unique(grep("27.6|27.7|27.5.b|27.12.b",dataframe$area,value = T))]<-"rng.27.5b6712b"


  # Next - SOL
  # In 27.3.a and 27.3.d.25
  # Part of sol.27.20-24

  # Next - SPR
  # In 27.3.a
  # Part of spr.27.3a4

  # Next - USK
  # In 27.3.a
  # Part of usk.27.3a45b6a7-912b

  # In 27.1
  dataframe$stock[dataframe$species=="USK" & dataframe$area %in%
                         GREandNEARC]<-"usk.27.1-2"


  #
  # Next - WHB
  # In 27.3.a
  # Part of whb.27.1-91214nea

  # Next - WHG
  # In 27.3.d.24, 27.3.c.22, 27.3.a, 27.2.a.2, 27.3.d.25
  # 27.3.a part of whg.27.3a

  # West of scotland is an official stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="WHG" & dataframe$area=="27.6.a" ]<-"whg.27.6a"

  # rockall is an official stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="WHG" & dataframe$area=="27.6.b" ]<-"whg.27.6b"

  # southern celtic sea and western channel is an official stock
  # Whiting (Merlangius merlangus) in divisions 7.b-c and 7.e-k (southern Celtic Seas and eastern English Channel)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="WHG" & dataframe$area %in%
                         unique(grep(paste0("27.7",letters[c(2:3,5:11)],collapse = "|"),dataframe$area,value = T))]<-"whg.27.7b-ce-k"


  # Next - ANE
  # in 27.7.e
  # They most likely belong to the ices stock ane.27.8
  dataframe$stock[dataframe$species=="ANE" & dataframe$area %in%
                    unique(grep(paste0("27.8"),dataframe$area,value = T))]<-"ane.27.8"

  # Next - ARG, ARY and ARU
  # These species are very hard to tell apart and will be treated together

  dataframe$stock[dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.3.a.n" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.3.a.21" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.4.a"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.4.b" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.4.c" |  dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.1.a"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.1.b" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.2.a" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.2.a.1"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.2.a.2" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.2.b" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.2.b.1"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.2.b.2" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.3.a.20"]<-"aru.27.123a4" # Argentina silus (Greater silver smelt)

  dataframe$stock[dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.5.b" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.5.b.1" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.5.b.1.a"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.5.b.1.b" | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.5.b.2"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.6.a"]<-"aru.27.5b6a" # Argentina silus (Greater silver smelt)

  dataframe$stock[dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.14.b.2"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.14.b.1"
                       | dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area=="27.14.b"] <- "aru.27.5a14"

  dataframe$stock[dataframe$species%in% c("ARU", "ARG", "ARY") & dataframe$area %in%
                         unique(grep("27.6.b|27.7|27.8|27.9|27.10|27.12",dataframe$area,value = T))] <- "aru.27.6b7-1012"

  # Next BOR
  #In 27.6.a and 27.7.c.2
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="BOR" & dataframe$area %in%
                         unique(grep("27.6|27.7|27.8",dataframe$area,value = T))] <- "boc.27.6-8"

  # Next COD
  # This is cod
  # lots of stocks
  dataframe$stock[dataframe$species == "COD" & dataframe$area =="27.5.a"]<-"cod.27.5a"

  dataframe$stock[dataframe$species == "COD" & dataframe$area %in% c("27.5.b","27.5.b.1.a","27.5.b.1.b")] <- "cod.27.5b1"

  dataframe$stock[dataframe$species == "COD" & dataframe$area %in% c("27.5.a","27.5.a.1","27.5.a.2")] <- "cod.27.5a"

  dataframe$stock[dataframe$species == "COD" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.6.a")))] <- "cod.27.6a"

  dataframe$stock[dataframe$species == "COD" & dataframe$area %in%
                         unique(str_subset(dataframe$area,c("27.6.b")))] <- "cod.27.6b"


  dataframe$stock[dataframe$species == "COD" & dataframe$area %in%
                         unique(grep("21.1.a|21.1.b|21.1.c|21.1.d|21.1.e",dataframe$area,value = T))] <- "cod.21.1a-e"

  dataframe$stock[dataframe$species == "COD" & dataframe$area %in%
                         unique(grep("21.1.f|27.14",dataframe$area,value = T))] <- "cod.2127.1f14"

  # Next is ELE
  #This is european eal#
  # ICES treats it as an official stock ele.2737.nea
  # So I will do the same
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="ELE"  %in%
                         unique(grep("27.",dataframe$area,value = T))]<-"ele.2737.nea"

  # Next is MUL and MUX
  # These are mullets
  # 27.3.d.24 27.4.b    27.4.c    27.3.c.22
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species %in% c("MUL", "MUX") & dataframe$area %in%
                         unique(grep("27.3|27.4|27.7.d",dataframe$area,value = T))]<-"mur.27.3a47d"

  # Saithe on iceland grounds is an own stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "POK" & dataframe$area=="27.5.a"]<-"pok.27.5a"

  # Saithe in the northeast is an own stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "POK" & dataframe$area %in%
                         GREandNEARC]<-"pok.27.1-2"


  # Next are RED
  # This is unkwon redfish #
  # in 27.4.a   27.2.a.2 27.4.b
  # They most likely belong to reg.27.561214 (caught in mixed demersal fishery)
  dataframe$stock[dataframe$species=="RED"]<-"reg.27.561214"

  # Next is RHG
  # This is roundhead grenadier
  # Official Ices Stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="RHG" & dataframe$area %in%
                         unique(grep("27.",dataframe$area,value = T))]<-"rhg.27.nea"


  # Next is RJF
  # This is shagreen ray
  # Stock in areas 6-7
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="RJF" & dataframe$area %in%
                         unique(grep("27.7|27.6",dataframe$area,value = T))]<-"rjf.27.67"




  # Next is RJR
  # This is starry ray (caught in northern waters 27.2.b.2 and 27.1.b)
  # There is an official ices stock rjr.27.23a4
  # All catches will be considered as this stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="RJR" & dataframe$area %in%
                         unique(grep("27.3.a|27.2|27.4",dataframe$area,value = T))]<-"rjr.27.23a4"
  # Next is SAL
  # This is salmon
  # 27.3.d.24 27.3.c.22 27.4.b
  # In baltic, this is an official stock (sal.27.22-31)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "SAL" & dataframe$area %in%
                         unique(grep(paste0("27.3.c.",22:31,collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species == "SAL" & dataframe$area %in%
                         unique(grep(paste0("27.3.d.",22:31,collapse = "|"),dataframe$area,value = T))]<-"sal.27.22-31"

  # In North Sea, it belongs to the atlantic stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "SAL" & dataframe$area %in%
                         unique(grep(paste0("27.",4:14,collapse = "|"),dataframe$area,value = T))]<-"sal.neac.all"



  # Next is SBR
  # This is blackspot (=red) seabream
  # Official stock in 27.6-8
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "SBR" & dataframe$area %in%
                         unique(grep("27.7|27.6|27.8",dataframe$area,value = T))] <- "sbr.27.6-8"

  dataframe$stock[dataframe$species == "SBR" & dataframe$area %in%
                         unique(str_subset(dataframe$area,"27.10")) ] <- "sbr.27.10"

  dataframe$stock[dataframe$species == "SBR" & dataframe$area %in%
                         unique(str_subset(dataframe$area,"27.9")) ] <- "sbr.27.9"



  # Next is SCK
  # This is kitefin shark
  # Official stock in 27.nea
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "SCK" & dataframe$area %in%
                         unique(str_subset(dataframe$area,"27")) ] <- "sck.27.nea"


  # Next is SPR
  # Sprat
  # English Channel (27.7.d und 27.7.e)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species=="SPR" & dataframe$area=="27.7.d" |
                         dataframe$species=="SPR" & dataframe$area=="27.7.e"]<-"spr.27.7de"


  # Next is TRO and TRS
  # This is trout
  # They are valuable bycatch in the baltic and have their own stock trs.27.22-32
  #27.3.d.24 27.3.c.22 27.4.b
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species %in% c("TRS","TRO") & dataframe$area %in%
                         unique(grep(paste0("27.3.c.",22:32,collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("TRS","TRO") & dataframe$area %in%
                         unique(grep(paste0("27.3.d.",22:32,collapse = "|"),dataframe$area,value = T))]<-"trs.27.22-32"


  # Next is WHG
  # In Celtic Sea (27.7b-ce-k)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "WHG" & dataframe$area=="27.7.e"]<-"whg.27.7b-ce-k"

  # Next is WHG
  # In Celtic Sea (27.7b-ce-k)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "WHG" & dataframe$area=="27.7.e"]<-"whg.27.7b-ce-k"


  ## Now, at the remaining ICES-Stocks
  # agn.27.nea
  # Angel shark (Squatina squatina) in subareas 1-10, 12 and 14 (the Northeast Atlantic and adjacent waters)
  ## THIS IS AN OFFICIAL ICES STOCK ###
  dataframe$stock[dataframe$species == "ANG" & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"agn.27.nea"

  # alf.27.nea
  ## THIS IS AN OFFICIAL ICES STOCK ###
  # Alfonsinos (Beryx spp.) in subareas 1-10, 12 and 14 (the Northeast Atlantic and adjacent waters)
  dataframe$stock[dataframe$species %in% c("ALF","BRX","BXD","BYS") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"alf.27.nea"

  # ane.27.9a
  ## THIS IS AN OFFICIAL ICES STOCK ###
  # Anchovy (Engraulis encrasicolus) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("ANE") & dataframe$area %in%
                         unique(str_subset(dataframe$area,"27.9.a"))]<-"ane.27.9a"

  # ank.27.78abd
  ## THIS IS AN OFFICIAL ICES STOCK ###
  # Black-bellied anglerfish (Lophius budegassa) in Subarea 7 and divisions 8.a-b and 8.d (Celtic Seas, Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("ANK") & dataframe$area %in%
                         unique(grep(paste0("27.7|27.8.a|8.b|8.c|8.d"),dataframe$area,value = T))]<-"ank.27.78abd"


  # ank.27.8c9a
  ## THIS IS AN OFFICIAL ICES STOCK ###
  # Black-bellied anglerfish (Lophius budegassa) in divisions 8.c and 9.a (Cantabrian Sea, Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("ANK") & dataframe$area %in%
                         unique(grep(paste0("27.8.c|9.a"),dataframe$area,value = T))]<-"ank.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # bsf.27.nea
  #Black scabbardfish (Aphanopus carbo) in subareas 1, 2, 4-8, 10, and 14, and divisions 3.a, 9.a, and 12.b (Northeast Atlantic and Arctic Ocean)
  dataframe$stock[dataframe$species %in% c("BSF") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1,2,4:8,10,14),collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("BSF") & dataframe$area %in%
                         unique(grep(paste0("27.3.a|9.a|27.12.b"),dataframe$area,value = T)) ]<-"bsf.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # bsk.27.nea
  # Basking shark (Cetorhinus maximus) in Subareas 1-10, 12 and 14 (Northeast Atlantic and adjacent waters)
  dataframe$stock[dataframe$species %in% c("BSK") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"bsk.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #bss.27.6a7bj
  # Seabass (Dicentrarchus labrax) in divisions 6.a, 7.b, and 7.j (West of Scotland,  West of Ireland, eastern part of southwest of Ireland)
  dataframe$stock[dataframe$species %in% c("BSS") & dataframe$area %in%
                         unique(grep(paste0("6.a|7.b|7.j"),dataframe$area,value = T))]<-"bss.27.6a7bj"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #bss.27.8ab
  # Seabass (Dicentrarchus labrax) in divisions 8.a-b (northern and central Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("BSS") & dataframe$area %in%
                         unique(grep(paste0("8.a|8.b"),dataframe$area,value = T))]<-"bss.27.8ab"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  bss.27.8c9a
  # Seabass (Dicentrarchus labrax) in divisions 8.c and 9.a (southern Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("BSS") & dataframe$area %in%
                         unique(grep(paste0("8.c|9.a"),dataframe$area,value = T))]<-"bss.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # bwp.27.2729-32
  # Baltic flounder (Platichthys solemdali) in subdivisions 27 and 29â€“32 (northern central and northern Baltic Sea)
  dataframe$stock[dataframe$species %in% c("BWP","FLE") & dataframe$area %in%
                         unique(grep(paste0(".d.",c(27,29:32),collapse = "|"),dataframe$area,value = T))]<-"bss.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # bwq.27.2425
  # Flounder (Platichthys spp) in subdivisions 24 and 25 (west of Bornholm and southwestern central Baltic)
  dataframe$stock[dataframe$species %in% c("BWP","FLE") & dataframe$area %in%
                         unique(grep(paste0(".d.",c(24:25),collapse = "|"),dataframe$area,value = T))]<-"bwq.27.2425"

  # cap.27.1-2
  # Capelin (Mallotus villosus) in subareas 1 and 2 (Northeast Arctic), excluding Division 2.a west of 5Â°W (Barents Sea capelin)
  dataframe$stock[dataframe$species %in% c("CAP") & dataframe$area %in%
                         GREandNEARC]<-"cap.27.1-2"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.21.1
  # Cod (Gadus morhua) in NAFO Subarea 1, inshore (West Greenland cod)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         unique(grep(paste0("21.1"),dataframe$area,value = T))]<-"cod.21.1"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.27.1-2
  # coast Cod (Gadus morhua) in subareas 1 and 2 (Norwegian coastal waters cod)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         GREandNEARC]<-"cod.27.1-2"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.27.21 Cod
  # (Gadus morhua) in Subdivision 21 (Kattegat)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         ("27.3.a.21")]<-"cod.27.21"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.27.5b2
  # Cod (Gadus morhua) in Subdivision 5.b.2 (Faroe Bank)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         ("27.5.b.2")]<-"cod.27.5b2"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.27.6b
  # Cod (Gadus morhua) in Division 6.b (Rockall)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         unique(grep(paste0("27.6.b"),dataframe$area,value = T))]<-"cod.27.6b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.27.7a
  # Cod (Gadus morhua) in Division 7.a (Irish Sea)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         unique(grep(paste0("27.7.a"),dataframe$area,value = T))]<-"cod.27.7a"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # cod.27.7e-k
  # Cod (Gadus morhua) in divisions 7.e-k (eastern English Channel and southern Celtic Seas)
  dataframe$stock[dataframe$species %in% c("COD") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[5:11],collapse = "|"),dataframe$area,value = T))]<-"cod.27.7e-k"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # gag.27.nea
  # Tope (Galeorhinus galeus) in subareas 1-10, 12 and 14 (the Northeast Atlantic and adjacent waters)
  dataframe$stock[dataframe$species == "GAG" & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"gag.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # gfb.27.nea
  # gfb.27.nea Greater forkbeard (Phycis blennoides) in subareas 1-10, 12 and 14 (the Northeast Atlantic and adjacent waters)
  dataframe$stock[dataframe$species %in% c("GFB","FOR","FOX") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"gfb.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # gur.27.3-8
  # Red gurnard (Chelidonichthys cuculus) in subareas 3-8 (Northeast Atlantic)
  dataframe$stock[dataframe$species %in% c("GUR") & dataframe$area %in%
                         unique(grep(paste0("27.",c(3:8),collapse = "|"),dataframe$area,value = T))]<-"gur.27.3"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # had.27.6b
  # Haddock (Melanogrammus aeglefinus) in Division 6.b (Rockall)
  dataframe$stock[dataframe$species %in% c("HAD") & dataframe$area %in%
                         unique(grep(paste0("27.6.b"),dataframe$area,value = T))]<-"had.27.6b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # had.27.7a
  #  Haddock (Melanogrammus aeglefinus) in Division 7.a (Irish Sea)
  dataframe$stock[dataframe$species %in% c("HAD") & dataframe$area %in%
                         unique(grep(paste0("27.7.a"),dataframe$area,value = T))]<-"had.27.7a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # her.27.28
  # Herring (Clupea harengus) in Subdivision 28.1 (Gulf of Riga)
  dataframe$stock[dataframe$species %in% c("HER") & dataframe$area %in%
                         unique(grep(paste0("27.3.d.28"),dataframe$area,value = T))]<-"her.27.28"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # her.27.3031
  # Herring (Clupea harengus) in Subdivisions 30 and 31 (Gulf of Bothnia)
  dataframe$stock[dataframe$species %in% c("HER") & dataframe$area %in%
                         unique(grep(paste0("27.3.d.30|27.3.d.31"),dataframe$area,value = T))]<-"her.27.3031"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # her.27.5a
  # Herring (Clupea harengus) in Division 5.a, summer-spawning herring (Iceland grounds)
  dataframe$stock[dataframe$species %in% c("HER") & dataframe$area %in%
                         unique(grep(paste0("27.5.a"),dataframe$area,value = T))]<-"her.27.5a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # her.27.nirs
  #Herring (Clupea harengus) in Division 7.a North of 52Â°30â€™N (Irish Sea)
  dataframe$stock[dataframe$species %in% c("HER") & dataframe$area %in%
                         unique(grep(paste0("27.7.a"),dataframe$area,value = T))]<-"her.27.nirs"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # her.27.nirs
  #Herring (Clupea harengus) in Division 7.a North of 52Â°30â€™N (Irish Sea)
  dataframe$stock[dataframe$species %in% c("HER") & dataframe$area %in%
                         unique(grep(paste0("27.7.a"),dataframe$area,value = T))]<-"her.27.nirs"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # hke.27.8c9a
  # Hake (Merluccius merluccius) in divisions 8.c and 9.a, Southern stock (Cantabrian Sea and  Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("HKE") & dataframe$area %in%
                         unique(grep(paste0("27.8.c|27.9.a"),dataframe$area,value = T))]<-"hke.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # hom.27.9a
  # Horse mackerel (Trachurus trachurus) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("HOM","JAX") & dataframe$area %in%
                         unique(grep(paste0("27.9.a"),dataframe$area,value = T))]<-"hom.27.9a"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # jaa.27.10a2
  # Blue jack mackerel (Trachurus picturatus) in Subdivision 10.a.2 (Azores grounds)
  dataframe$stock[dataframe$species %in% c("JAA") & dataframe$area %in%
                         unique(grep(paste0("10.a.2"),dataframe$area,value = T))]<-"jaa.27.10a2"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ldb.27.7b-k8abd
  # Four-spot megrim (Lepidorhombus boscii) in divisions 7.b-k, 8.a-b, and 8.d (west and southwest of Ireland, Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("LDB") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[2:11],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("LDB") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[c(1:2,4)],collapse = "|"),dataframe$area,value = T))]<-"ldb.27.7b-k8abd"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ldb.27.8c9a
  # Four-spot megrim (Lepidorhombus boscii) in divisions 8.c and 9.a (southern Bay of Biscay and Atlantic Iberian waters East)
  dataframe$stock[dataframe$species %in% c("LDB") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[3],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("LDB") & dataframe$area %in%
                         unique(grep(paste0("27.9.",letters[c(1)],collapse = "|"),dataframe$area,value = T))]<-"ldb.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # lez.27.4a6a
  # Megrim (Lepidorhombus spp.) in divisions 4.a and 6.a (northern North Sea, West of Scotland)
  dataframe$stock[dataframe$species %in% c("LEZ","MEG","LDB") & dataframe$area %in%
                         unique(grep(paste0("27.4.",letters[1],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("LEZ","MEG") & dataframe$area %in%
                         unique(grep(paste0("27.6.",letters[c(1)],collapse = "|"),dataframe$area,value = T))]<-"lez.27.4a6a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # lez.27.6b
  # Megrim (Lepidorhombus spp.) in Division 6.b (Rockall)
  dataframe$stock[dataframe$species %in% c("LEZ","MEG","LDB") & dataframe$area %in%
                         unique(grep(paste0("27.6.",letters[c(2)],collapse = "|"),dataframe$area,value = T))]<-"lez.27.6b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # meg.27.7b-k8abd
  # Megrim (Lepidorhombus whiffiagonis) in divisions 7.b-k, 8.a-b, and 8.d (west and southwest of Ireland, Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("LEZ","MEG") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[2:11],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("LEZ","MEG") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[c(1:2,4)],collapse = "|"),dataframe$area,value = T))]<-"meg.27.7b-k8abd"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # meg.27.8c9a
  # Megrim (Lepidorhombus whiffiagonis) in divisions 8.c and 9.a (Cantabrian Sea and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("LEZ","MEG") & dataframe$area %in%
                         unique(grep(paste0("27.9.",letters[1],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("LEZ","MEG") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[3],collapse = "|"),dataframe$area,value = T))]<-"meg.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # mon.27.78abd
  # White anglerfish (Lophius piscatorius) in Subarea 7 and divisions 8.a-b and 8.d (Celtic Seas, Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("MON","ANF") & dataframe$area %in%
                         unique(grep(paste0("27.7.",collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("MON","ANF") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[c(1:2,4)],collapse = "|"),dataframe$area,value = T))]<-"mon.27.78abd"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # mon.27.8c9a
  # White anglerfish (Lophius piscatorius) in divisions 8.c and 9.a (Cantabrian Sea and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("MON","ANF") & dataframe$area %in%
                         unique(grep(paste0("27.9.",letters[c(1)],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("MON","ANF") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[c(3)],collapse = "|"),dataframe$area,value = T))]<-"mon.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # mur.27.67a-ce-k89a
  # Striped red mullet (Mullus surmuletus) in subareas 6 and 8, and divisions 7.a-c, 7.e-k, and 9.a (North Sea, Bay of Biscay, southern Celtic Seas, and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("MUR") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[c(1:3,5:11)],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("MUR") & dataframe$area %in%
                         unique(grep("27.6|27.8|27.9.a",dataframe$area,value = T))]<-"mur.27.67a-ce-k89a"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # nop.27.3a4
  # Norway pout (Trisopterus esmarkii) in Subarea 4 and Division 3.a (North Sea, Skagerrak and Kattegat)
  dataframe$stock[dataframe$species %in% c("NOP") & dataframe$area %in%
                         unique(grep("27.4|27.3.a",dataframe$area,value = T))]<-"nop.27.3a4"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # pil.27.8c9a
  # Sardine (Sardina pilchardus) in divisions 8.c and 9.a (Cantabrian Sea and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("PIL") & dataframe$area %in%
                         unique(grep("27.8.c|27.9.a",dataframe$area,value = T))]<-"pil.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.7a
  # Plaice (Pleuronectes platessa) in Division 7.a (Irish Sea)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep("27.7.a",dataframe$area,value = T))]<-"ple.27.7a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.7bc
  # Plaice (Pleuronectes platessa) in divisions 7.b-c (West of Ireland)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep("27.7.b|27.7.c",dataframe$area,value = T))]<-"ple.27.7bc"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.7d
  # Plaice (Pleuronectes platessa) in Division 7.d (eastern English Channel)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep("27.7.d",dataframe$area,value = T))]<-"ple.27.7d"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.7e
  # Plaice (Pleuronectes platessa) in Division 7.e (western English Channel)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep("27.7.e",dataframe$area,value = T))]<-"ple.27.7e"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.7fg
  # Plaice (Pleuronectes platessa) in divisions 7.f and 7.g (Bristol Channel, Celtic Sea)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep("27.7.f|27.7.g",dataframe$area,value = T))]<-"ple.27.7fg"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.7fg
  # ple.27.7h-k Plaice (Pleuronectes platessa) in divisions 7.h-k (Celtic Sea South, southwest of Ireland)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[c(8:11)],collapse = "|"),dataframe$area,value = T))]<-"ple.27.7fg"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # ple.27.89a
  # Plaice (Pleuronectes platessa) in Subarea 8 and Division 9.a (Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("PLE") & dataframe$area %in%
                         unique(grep("27.8|27.9.a",dataframe$area,value = T))]<-"ple.27.89a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # pol.27.89a
  # Pollack (Pollachius pollachius) in Subarea 8 and Division 9.a (Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("POL") & dataframe$area %in%
                         unique(grep("27.8|27.9.a",dataframe$area,value = T))]<-"pol.27.89a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # pra.27.1-2
  # Northern shrimp (Pandalus borealis) in subareas 1 and 2 (Northeast Arctic)
  dataframe$stock[dataframe$species %in% c("PRA") & dataframe$area %in%
                         GREandNEARC]<-"pra.27.1-2"

  ## WARNING: PRA in 4a cant be distinguished, so it will be treated as one stock
  ## THIS IS AN OFFICIAL ICES STOCK ###
  # pra.27.3a4a
  # Northern shrimp (Pandalus borealis) in divisions 3.a and 4.a East and West
  dataframe$stock[dataframe$species %in% c("PRA") & dataframe$area %in%
                         unique(grep("27.3.a|27.4.a",dataframe$area,value = T))]<-"pra.27.3a4a"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # reb.27.14b
  # Beaked redfish  (Sebastes mentella) in Division 14.b, demersal (Southeast Greenland)
  dataframe$stock[dataframe$species %in% c("RED","REB") & dataframe$area %in%
                         unique(grep("27.14.b.2",dataframe$area,value = T))] <- "reb.27.14b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # reb.27.5a14
  # Beaked redfish  (Sebastes mentella) in Subarea 14 and Division 5.a, Icelandic slope stock (East of Greenland, Iceland grounds)
  dataframe$stock[dataframe$species %in% c("RED","REB") & dataframe$area %in%
                         unique(grep("27.14.b.1|27.14.a|27.5.a",dataframe$area,value = T)) |
                         dataframe$species %in% c("RED","REB") & dataframe$area %in%
                         c("27.14.b")] <- "reb.27.5a14"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rja.27.nea
  # White skate (Rostroraja alba) in subareas 1-10, 12 and 14 (the Northeast Atlantic and adjacent waters)
  dataframe$stock[dataframe$species %in% c("RJA") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"rja.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjb.27.3a4
  # Common skate complex (Blue skate (Dipturus batis) and flapper skate (Dipturus intermedius) in Subarea 4 and Division 3.a (North Sea, Skagerrak and Kattegat)
  dataframe$stock[dataframe$species %in% c("RJB") & dataframe$area %in%
                         unique(grep(paste0("27.4|27.3.a",collapse = "|"),dataframe$area,value = T))]<-"rja.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjb.27.67a-ce-k
  # Common skate complex (Blue skate (Dipturus batis) and flapper skate (Dipturus intermedius) in Subarea 6 and divisions 7.aâ€“c and 7.eâ€“k (Celtic Seas and western English Channel)
  dataframe$stock[dataframe$species %in% c("RJB") & dataframe$area %in%
                         unique(grep(paste0("27.6"),dataframe$area,value = T))|
                         dataframe$species %in% c("RJB") & dataframe$area %in%
                         unique(grep(paste0("27.7",letters[c(1:3,5:11)],collapse = "|"),dataframe$area,value = T))]<-"rjb.27.67a-ce-k"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjb.27.89a
  # Common skate complex (Blue skate (Dipturus batis) and flapper skate (Dipturus intermedius) in Subarea 8 and Division 9.a (Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("RJB") & dataframe$area %in%
                         unique(grep(paste0("27.8|27.9.a"),dataframe$area,value = T))]<-"rjb.27.89a"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  #rjc.27.6
  # Thornback ray (Raja clavata) in Subarea 6 (West of Scotland)
  dataframe$stock[dataframe$species %in% c("RJC") & dataframe$area %in%
                         unique(grep(paste0("27.6"),dataframe$area,value = T))]<-"rjc.27.6"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #rjc.27.7afg
  # Thornback ray (Raja clavata) in Subarea 6 (West of Scotland)
  dataframe$stock[dataframe$species %in% c("RJC") & dataframe$area %in%
                         unique(grep(paste0("27.7",letters[c(1,6:7)],collapse = "|"),dataframe$area,value = T))]<-"rjc.27.7afg"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjc.27.7e
  # Thornback ray (Raja clavata) in Division 7.e (western English Channel)
  dataframe$stock[dataframe$species %in% c("RJC") & dataframe$area %in%
                         unique(grep(paste0("27.7.e"),dataframe$area,value = T))]<-"rjc.27.7e"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjc.27.8
  # rjc.27.8 Thornback ray (Raja clavata) in Subarea 8 (Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("RJC") & dataframe$area %in%
                         unique(grep(paste0("27.8"),dataframe$area,value = T))]<-"rjc.27.8"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjc.27.9a
  # Thornback ray (Raja clavata) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("RJC") & dataframe$area %in%
                         unique(grep(paste0("27.9.a"),dataframe$area,value = T))]<-"rjc.27.9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rje.27.7de #
  # Small-eyed ray (Raja microocellata) in divisions 7.d and 7.e (English Channel)
  dataframe$stock[dataframe$species %in% c("RJE") & dataframe$area %in%
                         unique(grep(paste0("27.7.d|27.7.e"),dataframe$area,value = T))]<-"rje.27.7de"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rje.27.7fg
  # Small-eyed ray (Raja microocellata) in divisions 7.f and 7.g (Bristol Channel, Celtic Sea North)
  dataframe$stock[dataframe$species %in% c("RJE") & dataframe$area %in%
                         unique(grep(paste0("27.7.f|27.7.g"),dataframe$area,value = T))]<-"rje.27.7fg"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjh.27.4a6
  # Blonde ray (Raja brachyura) in Subarea 6 and Division 4.a (North Sea and West of Scotland)
  dataframe$stock[dataframe$species %in% c("RJH") & dataframe$area %in%
                         unique(grep(paste0("27.6|27.4.a"),dataframe$area,value = T))]<-"rjh.27.4a6"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjh.27.7afg
  # Blonde ray (Raja brachyura) in divisions 7.a and 7.f-g (Irish Sea, Bristol Channel, Celtic Sea North)
  dataframe$stock[dataframe$species %in% c("RJH") & dataframe$area %in%
                         unique(grep(paste0("27.7.a|27.7.f|27.7.g"),dataframe$area,value = T))]<-"rjh.27.7afg"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjh.27.7e
  # Blonde ray (Raja brachyura) in Division 7.e (western English Channel)
  dataframe$stock[dataframe$species %in% c("RJH") & dataframe$area %in%
                         unique(grep(paste0("27.7.e"),dataframe$area,value = T))]<-"rjh.27.7e"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  rjh.27.9a
  # Blonde ray (Raja brachyura) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("RJH") & dataframe$area %in%
                         unique(grep(paste0("27.9.a"),dataframe$area,value = T))]<-"rjh.27.9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rji.27.67
  # Sandy ray (Leucoraja circularis) in subareas 6-7 (West of Scotland, southern Celtic Seas, English Channel)
  dataframe$stock[dataframe$species %in% c("RJI") & dataframe$area %in%
                         unique(grep(paste0("27.6|27.7"),dataframe$area,value = T))]<-"rji.27.67"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjm.27.3a47d
  # Spotted ray (Raja montagui) in Subarea 4 and Divisions 3.a and 7.d (North Sea, Skagerrak, Kattegat, and eastern English Channel)
  dataframe$stock[dataframe$species %in% c("RJM") & dataframe$area %in%
                         unique(grep(paste0("27.4|27.3.a|27.7.d"),dataframe$area,value = T))]<-"rjm.27.3a47d"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjm.27.67bj
  # Spotted ray (Raja montagui) in Subarea 6 and divisions 7.b and 7.j (West of Scotland, west and southwest of Ireland)
  dataframe$stock[dataframe$species %in% c("RJM") & dataframe$area %in%
                         unique(grep(paste0("27.7.b|27.7.j."),dataframe$area,value = T))]<-"rjm.27.67bj"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  rjm.27.7ae-h
  # Spotted ray (Raja montagui) in divisions 7.a and 7.e-h (southern Celtic Seas and western English Channel)
  dataframe$stock[dataframe$species %in% c("RJM") & dataframe$area %in%
                         unique(grep(paste0("27.7",letters[c(1,5:8)],collapse = "|"),dataframe$area,value = T))]<-"rjm.27.7ae-h"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #rjm.27.8
  # Spotted ray (Raja montagui) in  Subarea 8 (Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("RJM") & dataframe$area %in%
                         unique(grep(paste0("27.8"),dataframe$area,value = T))]<-"rjm.27.8"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjm.27.9a
  # Spotted ray (Raja montagui) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("RJM") & dataframe$area %in%
                         unique(grep(paste0("27.9.a"),dataframe$area,value = T))]<-"rjm.27.9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjn.27.678abd
  # Cuckoo ray (Leucoraja naevus) in subareas 6-7 and divisions 8.a-b and 8.d (West of Scotland, southern Celtic Seas, and western English Channel, Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("RJN") & dataframe$area %in%
                         unique(grep(paste0("27.6|27.7|27.8.a|27.8.b|27.8.d"),dataframe$area,value = T))]<-"rjn.27.678abd"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  rjn.27.8c
  # Cuckoo ray (Leucoraja naevus) in Division 8.c (Cantabrian Sea)
  dataframe$stock[dataframe$species %in% c("RJN") & dataframe$area %in%
                         unique(grep(paste0("27.8.c"),dataframe$area,value = T))]<-"rjn.27.8c"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rjn.27.9a
  # Cuckoo ray (Leucoraja naevus) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("RJN") & dataframe$area %in%
                         unique(grep(paste0("27.9.a"),dataframe$area,value = T))]<-"rjn.27.9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rju.27.7bj
  # Undulate ray (Raja undulata) in divisions 7.b and 7.j (west and southwest of Ireland)
  dataframe$stock[dataframe$species %in% c("RJU") & dataframe$area %in%
                         unique(grep(paste0("27.7.b|27.7.j"),dataframe$area,value = T))]<-"rju.27.7bj"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  rju.27.7de
  # Undulate ray (Raja undulata) in divisions 7.d and 7.e (English Channel)
  dataframe$stock[dataframe$species %in% c("RJU") & dataframe$area %in%
                         unique(grep(paste0("27.7.d|27.7.e"),dataframe$area,value = T))]<-"rju.27.7de"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rju.27.8ab
  # Undulate ray (Raja undulata) in divisions 8.a-b (northern and central Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("RJU") & dataframe$area %in%
                         unique(grep(paste0("27.8.a|27.8.b"),dataframe$area,value = T))]<-"rju.27.8ab"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rju.27.8c
  # Undulate ray (Raja undulata) in Division 8.c (Cantabrian Sea)
  dataframe$stock[dataframe$species %in% c("RJU") & dataframe$area %in%
                         unique(grep(paste0("27.8.c"),dataframe$area,value = T))]<-"rju.27.8c"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rju.27.9a
  # Undulate ray (Raja undulata) in Division 9.a (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("RJU") & dataframe$area %in%
                         unique(grep(paste0("27.9.a"),dataframe$area,value = T))]<-"rju.27.9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # rng.27.5a10b12ac14b
  # Roundnose grenadier (Coryphaenoides rupestris) in Divisions 10.b and 12.c, and Subdivisions 12.a.1, 14.b.1, and 5.a.1 (Oceanic Northeast Atlantic and northern Reykjanes Ridge)
  dataframe$stock[dataframe$species %in% c("RNG") & dataframe$area %in%
                         unique(grep(paste0("27.10.b|27.12.c|27.12.a.1|27.14.b.1|27.5.a.1"),dataframe$area,value = T))]<-"rng.27.5a10b12ac14b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sal.27.32
  # Salmon (Salmo salar) in Subdivision 32 (Gulf of Finland)
  dataframe$stock[dataframe$species %in% c("SAL") & dataframe$area %in%
                         unique(grep(paste0("27.3.d.32"),dataframe$area,value = T))]<-"sal.27.32"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sal.wgc.all
  # Salmon (Salmo salar) in Subarea 14 and NAFO division 1 (east and west of Greenland)
  dataframe$stock[dataframe$species %in% c("SAL") & dataframe$area %in%
                         unique(grep(paste0("27.14|21.1"),dataframe$area,value = T))]<-"sal.wgc.all"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sbr.27.10
  # Blackspot seabream (Pagellus bogaraveo) in Subarea 10 (Azores grounds)
  dataframe$stock[dataframe$species %in% c("SBR") & dataframe$area %in%
                         unique(grep(paste0("27.10."),dataframe$area,value = T))]<-"sbr.27.10"


  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sbr.27.9
  # Blackspot seabream (Pagellus bogaraveo) in Subarea 9 (Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("SBR") & dataframe$area %in%
                         unique(grep(paste0("27.9."),dataframe$area,value = T))]<-"sbr.27.9"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sdv.27.nea
  # Smooth-hound (Mustelus spp.) in subareas 1-10, 12 and 14 (the Northeast Atlantic and adjacent waters)
  dataframe$stock[dataframe$species %in% c("SDV") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:10,12,14),collapse = "|"),dataframe$area,value = T))]<-"sdv.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  seh.27.1
  # Harp seals (Pagophilus groenlandicus) in Subarea 1 (Barents and White sea stock)
  dataframe$stock[dataframe$species %in% c("SEH") & dataframe$area %in%
                         c("27.1.","27.1.a","27.1" ,"27.1.b")]<-"seh.27.1"

  # This cant be distinguished by area code alone. All SEH caught in 27.1 will be considered part of the Barents and White Sea stock

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # seh.27.125a14
  # Harp seals (Pagophilus groenlandicus) in subareas 1, 2 and 14 and Division 5.a (Greenland Sea stock)
  dataframe$stock[dataframe$species %in% c("SEH") & dataframe$area %in%
                         c(  "27.2.a.1", "27.2.b.1", "27.2.a.2","27.2.b.2","27.2.b","27.2.a" )|
                         dataframe$species %in% c("SEH") & dataframe$area %in%
                         unique(grep(paste0("27.14|27.5.a"),dataframe$area,value = T))]<-"seh.27.125a14"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sez.27.2514
  # Hooded seals (Cystophora cristata) in subareas 2, 5 and 14 (Greenland Sea stock)
  dataframe$stock[dataframe$species %in% c("SEZ") & dataframe$area %in%
                         unique(grep(paste0("27.14|27.5|27.2"),dataframe$area,value = T))]<-"sez.27.2514"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sho.27.67
  # Black-mouth dogfish (Galeus melastomus) in subareas 6 and 7 (West of Scotland, southern Celtic Seas, and English Channel)
  dataframe$stock[dataframe$species %in% c("SHO") & dataframe$area %in%
                         unique(grep(paste0("27.6|27.7"),dataframe$area,value = T))]<-"sho.27.67"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sho.27.89a
  # Black-mouth dogfish (Galeus melastomus) in Subarea 8 and Division 9.a (Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("SHO") & dataframe$area %in%
                         unique(grep(paste0("27.8|27.9.a"),dataframe$area,value = T))]<-"sho.27.89a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.7a
  # Sole (Solea solea) in Division 7.a (Irish Sea)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.7.a"),dataframe$area,value = T))]<-"sol.27.7a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.7bc
  # Sole (Solea solea) in divisions 7.b and 7.c (West of Ireland)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.7.b|27.7.c"),dataframe$area,value = T))]<-"sol.27.7bc"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.7d
  # Sole (Solea solea) in Division 7.d (eastern English Channel)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.7.d"),dataframe$area,value = T))]<-"sol.27.7d"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.7e
  # Sole (Solea solea) in Division 7.e (western English Channel)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.7.e"),dataframe$area,value = T))]<-"sol.27.7e"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.7fg
  # Sole (Solea solea) in divisions 7.f and 7.g (Bristol Channel, Celtic Sea)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.7.f|27.7.g"),dataframe$area,value = T))]<-"sol.27.7fg"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.7h-k
  # Sole (Solea solea) in Divisions 7.h-k (Celtic Sea South, southwest of Ireland)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[8:11],collapse = "|"),dataframe$area,value = T))]<-"sol.27.7h-k"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.8ab
  # Sole (Solea solea) in divisions 8.a-b (northern and central Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.8.a|27.8.b"),dataframe$area,value = T))]<-"sol.27.8ab"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # sol.27.8c9a
  # Sole (Solea solea) in divisions 8.c and 9.a (Cantabrian Sea and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("SOL") & dataframe$area %in%
                         unique(grep(paste0("27.8.c|27.9.a"),dataframe$area,value = T))]<-"sol.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  spr.27.67a-cf-k
  # Sprat (Sprattus sprattus) in Subarea 6 and Divisions 7.a-c and 7.f-k (West of Scotland, southern Celtic Seas)
  dataframe$stock[dataframe$species %in% c("SPR") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[c(1:3,6:11)],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("SPR") & dataframe$area %in%
                         unique(grep(paste0("27.6"),dataframe$area,value = T))]<-"spr.27.67a-cf-k"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # syc.27.3a47d
  # Lesser spotted dogfish (Scyliorhinus canicula) in Subarea 4 and divisions 3.a and 7.d (North Sea, Skagerrak and Kattegat, eastern English Channel)
  dataframe$stock[dataframe$species %in% c("SYC") & dataframe$area %in%
                         unique(grep(paste0("27.4|27.3.a|27.7.d"),dataframe$area,value = T))]<-"syc.27.3a47d"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # syc.27.67a-ce-j
  # Lesser spotted dogfish (Scyliorhinus canicula) in Subarea 6 and divisions 7.a-c and 7.e-j  (West of Scotland, Irish Sea, southern Celtic Seas)
  dataframe$stock[dataframe$species %in% c("SYC") & dataframe$area %in%
                         unique(grep(paste0("27.7.",letters[c(1:3,6:11)],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("SPR") & dataframe$area %in%
                         unique(grep(paste0("27.6"),dataframe$area,value = T))]<-"spr.27.67a-cf-k"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # syc.27.8abd
  # Lesser spotted dogfish (Scyliorhinus canicula) in divisions 8.a-b and 8.d (Bay of Biscay)
  dataframe$stock[dataframe$species %in% c("SYC") & dataframe$area %in%
                         unique(grep(paste0("27.8.",letters[c(1:2,4)],collapse = "|"),dataframe$area,value = T))]<-"syc.27.8abd"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # syc.27.8c9a
  # Lesser spotted dogfish (Scyliorhinus canicula) in divisions 8.c and 9.a (Cantabrian Sea and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("SYC") & dataframe$area %in%
                         unique(grep(paste0("27.8.c|27.9.a"),dataframe$area,value = T))]<-"syc.27.8c9a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # syt.27.67
  # Greater-spotted dogfish (Scyliorhinus stellaris) in subareas 6 and 7 (West of Scotland, southern Celtic Sea, and the English Channel)
  dataframe$stock[dataframe$species %in% c("SYT") & dataframe$area %in%
                         unique(grep(paste0("27.6|27.7"),dataframe$area,value = T))]<-"syt.27.67"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # thr.27.nea
  # Thresher sharks (Alopias spp.) in Subareas 10, 12, Divisions 7.c-k, 8.d-e, and Subdivisions 5.b.1, 9.b.1, 14.b.1 (Northeast Atlantic)
  dataframe$stock[dataframe$species %in% c("THR") & dataframe$area %in%
                         unique(grep(paste0("27.10|27.12|27.5.b.1|27.9.b.1|27.14.b.1"),dataframe$area,value = T)) |
                         dataframe$species %in% c("THR") & dataframe$area %in%
                         unique(grep(paste0("27.7", letters[3:11],collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("THR") & dataframe$area %in%
                         unique(grep(paste0("27.8", letters[4:5],collapse = "|"),dataframe$area,value = T))]<-"thr.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # tsu.27.nea
  # Roughsnout grenadier (Trachyrincus scabrus) in subareas 1-2, 4-8, 10, 12, 14 and Division 3a (Northeast Atlantic and Arctic Ocean)
  dataframe$stock[dataframe$species %in% c("TSU") & dataframe$area %in%
                         unique(grep(paste0("27.",c(1:2,4:8,10,12,14),collapse = "|"),dataframe$area,value = T)) |
                         dataframe$species %in% c("TSU") & dataframe$area %in%
                         unique(grep(paste0("27.3.a"),dataframe$area,value = T))] <- "tsu.27.nea"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # usk.27.12ac
  # Tusk (Brosme brosme) in Subarea 12, excluding Division 12.b (southern Mid-Atlantic Ridge)
  dataframe$stock[dataframe$species %in% c("USK") & dataframe$area %in%
                         unique(grep(paste0("27.12"),dataframe$area,value = T)) &
                         dataframe$area != "27.12.b"]<-"usk.27.12ac"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # usk.27.6b
  # Tusk (Brosme brosme) in Division 6.b (Rockall)
  dataframe$stock[dataframe$species %in% c("USK") & dataframe$area %in%
                         unique(grep(paste0("27.6.b"),dataframe$area,value = T))]<-"usk.27.6b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # whg.27.6b
  # Whiting (Merlangius merlangus) in Division 6.b (Rockall)
  dataframe$stock[dataframe$species %in% c("WHG") & dataframe$area %in%
                         unique(grep(paste0("27.6.b"),dataframe$area,value = T))]<-"whg.27.6b"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  # whg.27.7a
  # Whiting (Merlangius merlangus) in Division 7.a (Irish Sea)
  dataframe$stock[dataframe$species %in% c("WHG") & dataframe$area %in%
                         unique(grep(paste0("27.7.a"),dataframe$area,value = T))]<-"whg.27.7a"

  ## THIS IS AN OFFICIAL ICES STOCK ###
  #  whg.27.89a
  # Whiting (Merlangius merlangus) in Subarea 8 and Division 9.a (Bay of Biscay and Atlantic Iberian waters)
  dataframe$stock[dataframe$species %in% c("WHG") & dataframe$area %in%
                         unique(grep(paste0("27.9.a|27.8"),dataframe$area,value = T))]<-"whg.27.89a"

  ### TUNA & LPF stocks ####
  ## First is BFT - Bluefin Tuna
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  #  BFT-E
  #  Bluefin Tuna in Mediterranean and Eastern Atlantic
  dataframe$stock[dataframe$species %in% c("BFT") & dataframe$area %in%
                    unique(grep(paste0("27.|37.|34.|47.|GSA"),dataframe$area,value = T))]<-"BFT-E"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  #  BFT-W
  #  Bluefin Tuna in Western Atlantic
  dataframe$stock[dataframe$species %in% c("BFT") & dataframe$area %in%
                    unique(grep(paste0("21.|31.|41."),dataframe$area,value = T))]<-"BFT-W"


  ## Next is ALB - Albacore ##
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## ALB-M
  ## Medditeranean Albacore
  dataframe$stock[dataframe$species %in% c("ALB") & dataframe$area %in%
                    unique(grep(paste0("GSA|37."),dataframe$area,value = T))]<-"ALB-M"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## ALB-N
  ## Northern Albacore
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Albacore stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("ALB") & dataframe$area %in%
                    unique(grep(paste0("27.|31.|21.|34.1|34.2|34.3.1|34.3.2|34.4.2"),dataframe$area,value = T))]<-"ALB-N"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## ALB-S
  ## Southern Albacore
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Albacore stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("ALB") & dataframe$area %in%
                    unique(grep(paste0("34.4.1|34.3.3|34.3.4|34.3.5|34.3.6|41.|47.|48.6"),dataframe$area,value = T))]<-"ALB-S"

  ### Next is BET - Bigeye Tuna
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## BET-A
  ## Atlantic Bigeye Tuna
  dataframe$stock[dataframe$species %in% c("BET") & dataframe$area %in%
                    unique(grep(paste0("27.|21.|31.|34.|37.|41.|47.|48."),dataframe$area,value = T))]<-"BET-A"

  ### Next is YFT - Yellowfin Tuna
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## YFT-A
  ## Atlantic Yellowfin Tuna
  dataframe$stock[dataframe$species %in% c("YFT") & dataframe$area %in%
                    unique(grep(paste0("27.|21.|31.|34.|37.|41.|47.|48."),dataframe$area,value = T))]<-"YFT-A"

  ## Next is SKJ - Skipjack Tuna
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SKJ-E
  ## Eastern Skipjack Tuna
  dataframe$stock[dataframe$species %in% c("SKJ") & dataframe$area %in%
                    unique(grep(paste0("27.|37.|34.|47."),dataframe$area,value = T))]<-"SKJ-E"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SKJ-W
  ## Western Skipjack Tuna
  dataframe$stock[dataframe$species %in% c("SKJ") & dataframe$area %in%
                    unique(grep(paste0("21.|31.|41."),dataframe$area,value = T))]<-"SKJ-W"

  ## Next is SWO - This is swordfish
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SWO-M
  ## Medditeranean Swordfish
  dataframe$stock[dataframe$species %in% c("SWO") & dataframe$area %in%
                    unique(grep(paste0("GSA|37."),dataframe$area,value = T))]<-"SWO-M"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SWO-N
  ## Northern Swordfish
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the swordfish stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("SWO") & dataframe$area %in%
                    unique(grep(paste0("27.|31.|21.|34.1|34.2|34.3.1|34.3.2|34.4.2"),dataframe$area,value = T))]<-"SWO-N"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SWO-S
  ## Southern Swordfish
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Swordfish stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("WO") & dataframe$area %in%
                    unique(grep(paste0("34.4.1|34.3.3|34.3.4|34.3.5|34.3.6|41.|47.|48.6"),dataframe$area,value = T))]<-"SWO-S"


  ## Next is BUM - This is Blue Marlin
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## BUM-A
  ## Atlantic Blue Marlin
  dataframe$stock[dataframe$species %in% c("BUM") & dataframe$area %in%
                    unique(grep(paste0("27.|21.|31.|34.|37.|41.|47.|48."),dataframe$area,value = T))]<-"BUM-A"

  ## Next is WHM - This is White Marlin
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## WHM-A
  ## Atlantic White Marlin
  dataframe$stock[dataframe$species %in% c("WHM") & dataframe$area %in%
                    unique(grep(paste0("27.|21.|31.|34.|37.|41.|47.|48."),dataframe$area,value = T))]<-"WHM-A"

  ## Next is SAI - This is Sailfish
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SAI-E
  ## Eastern Sailfish
  dataframe$stock[dataframe$species %in% c("SAI") & dataframe$area %in%
                    unique(grep(paste0("27.|37.|34.|47."),dataframe$area,value = T))]<-"SAI-E"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SAI-W
  ## Western Sailfish
  dataframe$stock[dataframe$species %in% c("SAI") & dataframe$area %in%
                    unique(grep(paste0("21.|31.|41."),dataframe$area,value = T))]<-"SAI-W"


  ## Next is SPF - This is Longbill Spearfish
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SPF-E
  ## Eastern Longbill Spearfish
  dataframe$stock[dataframe$species %in% c("SPF") & dataframe$area %in%
                    unique(grep(paste0("27.|37.|34.|47."),dataframe$area,value = T))]<-"SPF-E"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SPF-W
  ## Western Longbill Spearfish
  dataframe$stock[dataframe$species %in% c("SPF") & dataframe$area %in%
                    unique(grep(paste0("21.|31.|41."),dataframe$area,value = T))]<-"SPF-W"

  ## Next is BSH - Blue shark ##
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## BSH-N
  ## Northern Blue shark
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Blue shark stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("BSH") & dataframe$area %in%
                    unique(grep(paste0("27.|31.|21.|34.1|34.2|34.3.1|34.3.2|34.4.2|GSA|37."),dataframe$area,value = T))]<-"BSH-N"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## BSH-S
  ## Southern Blue shark
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Blue shark stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("BSH") & dataframe$area %in%
                    unique(grep(paste0("34.4.1|34.3.3|34.3.4|34.3.5|34.3.6|41.|47.|48.6"),dataframe$area,value = T))]<-"BSH-S"

  ## Next is SMA - Shortfin Mako ##
  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SMA-N
  ## Northern Shortfin Mako
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Blue shark stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("SMA") & dataframe$area %in%
                    unique(grep(paste0("27.|31.|21.|34.1|34.2|34.3.1|34.3.2|34.4.2|GSA|37."),dataframe$area,value = T))]<-"SMA-N"

  ## THIS IS AN OFFICIAL ICCAT STOCK ###
  ## SMA-S
  ## Southern Shortfin Mako
  ## CAUTION - AN ASSUMPTION IS MADE HERE - THE FAO-AREAS DONT FULLY MATCH THE ICES AREAS
  ## The North-South line of the Shortfin Mako stock is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe$stock[dataframe$species %in% c("SMA") & dataframe$area %in%
                    unique(grep(paste0("34.4.1|34.3.3|34.3.4|34.3.5|34.3.6|41.|47.|48.6"),dataframe$area,value = T))]<-"SMA-S"

  ### Join MED stocks ####
  MED_stocks_wdf <- MED_stocks
  names(MED_stocks_wdf) <- c(c("species","area","stock"))
  dataframe_nostock <- dataframe[,-which(names(dataframe) %in% c("stock"))]

  dataframe_med <- left_join(dataframe_nostock,MED_stocks_wdf)
  dataframe_med <- dataframe_med[!is.na(dataframe_med$stock),]

  if(nrow(dataframe_med)>0){
    dataframe_bind <- anti_join(dataframe,dataframe_med, by=c("ship_ID","species","area"))
    dataframe_join <- rbind(dataframe_bind,dataframe_med)
  } else{
    dataframe_join <- dataframe
  }

  ## Include the ICCAT tuna and shark bycatch species
  ICCAT_Bycatch <- c("ASM","BAU","BBM","BEP","BIL","BIP","BKJ","BLF","BLM","BLT","BON",
                     "BOP","BRS","BUK","CER","CHY","COM","DBM","DOT","FRI","FRZ","GUT",
                     "KAK","KAW","KGM","KGX","KOS","LEB","LOT","LTA","MAC","MAS","MAW",
                     "MLS","MOS","MSP","NPH","PAP","QUM","RSP","SBF","SFA","SHM","SIE",
                     "SLT","SSM","SSP","STS","TUN","TUS","TUX","WAH","AGN","ALV","API",
                     "ASK","BRO","BSK","BTH","CCA","CCB","CCE","CCG","CCL","CCN","CCO",
                     "CCP","CCR","CCS","CCT","CCV","CFB","CPL","CPU","CTI","CVX","CYO",
                     "CYY","DCA","DGH","DGS","DGX","DGZ","DNA","DOP","DUS","ETP","ETR",
                     "ETX","FAL","GAG","GAU","GNC","GSK","GUP","GUQ","HEI","HXC","ISB",
                     "LMA","LMP","LOO","MAK","MAN","MPO","MRB","MSK","MTR","MYL","NGB",
                     "NTC","OCS","ODH","OXN","OXY","PLS","POR","PSK","PTM","QUB","QUC",
                     "QUL","RDA","RDC","RFL","RHA","RHN","RHR","RHT","RHZ","RJF","RMB",
                     "RMH","RMJ","RMM","RMO","RMT","RSK","SBL","SCK","SCL","SDP","SDS",
                     "SDV","SHB","SHL","SHO","SHX","SKH","SMD","SOR","SPJ","SPK","SPL",
                     "SPN","SPY","SPZ","SSQ","STT","SUA","SUT","SYC","SYR","SYT","SYX",
                     "THR","TIG","TRK","TTO","WHS","MAE","JDY","RBY","RGL","RGI","RMN",
                     "MYM","MYO","RPP","RPM","MRM","RTB")

  ## Mediterranean Bycatch stocks
  dataframe_join$stock[dataframe_join$stock == "Bycatch/Unknown" &
                         dataframe_join$species %in% ICCAT_Bycatch &
                         dataframe_join$area %in%
                                 unique(grep(paste0("37.|GAS"),dataframe_join$area,value = T))] <-
    paste0(dataframe_join$species[dataframe_join$stock == "Bycatch/Unknown" &
                                    dataframe_join$species %in% ICCAT_Bycatch &
                                    dataframe_join$area %in%
                                            unique(grep(paste0("37.|GAS"),dataframe_join$area,value = T))],"-MD")


  ## Atlantic Northwest Bycatch stocks
  dataframe_join$stock[dataframe_join$stock == "Bycatch/Unknown" &
                         dataframe_join$species %in% ICCAT_Bycatch &
                         dataframe_join$area %in%
                         unique(grep(paste0("31.|21."),dataframe_join$area,value = T))] <-
    paste0(dataframe_join$species[dataframe_join$stock == "Bycatch/Unknown" &
                                    dataframe_join$species %in% ICCAT_Bycatch &
                                    dataframe_join$area %in%
                                    unique(grep(paste0("21.|31."),dataframe_join$area,value = T))],"-NW")

  ## Atlantic Southwest Bycatch stocks
  ## Caution, an assumption is made - Areas 48.2 and 48.4 are only partly within the ICCAT SW Area, but the majority of those areas is. Therefore, they are considered a part of ICCAT-SW
  dataframe_join$stock[dataframe_join$stock == "Bycatch/Unknown" &
                         dataframe_join$species %in% ICCAT_Bycatch &
                         dataframe_join$area %in%
                         unique(grep(paste0("41.|48.2|48.3|48.4"),dataframe_join$area,value = T))] <-
    paste0(dataframe_join$species[dataframe_join$stock == "Bycatch/Unknown" &
                                    dataframe_join$species %in% ICCAT_Bycatch &
                                    dataframe_join$area %in%
                                    unique(grep(paste0("41.|48.2|48.3|48.4"),dataframe_join$area,value = T))],"-SW")

  ## Atlantic Northeast Bycatch stocks
  ## The North-South line of ICCAT is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe_join$stock[dataframe_join$stock == "Bycatch/Unknown" &
                         dataframe_join$species %in% ICCAT_Bycatch &
                         dataframe_join$area %in%
                         unique(grep(paste0("27.|31.|21.|34.1|34.2|34.3.1|34.3.2|34.4.2"),dataframe_join$area,value = T))] <-
    paste0(dataframe_join$species[dataframe_join$stock == "Bycatch/Unknown" &
                                    dataframe_join$species %in% ICCAT_Bycatch &
                                    dataframe_join$area %in%
                                    unique(grep(paste0("27.|31.|21.|34.1|34.2|34.3.1|34.3.2|34.4.2"),dataframe_join$area,value = T))],"-NE")

  ## Atlantic Southeast Bycatch stocks
  ## The North-South line of ICCAT is at 5°N - FAO areas will be assigned according to where the majority of the area is
  dataframe_join$stock[dataframe_join$stock == "Bycatch/Unknown" &
                         dataframe_join$species %in% ICCAT_Bycatch &
                         dataframe_join$area %in%
                         unique(grep(paste0("34.4.1|34.3.3|34.3.4|34.3.5|34.3.6|41.|47.|48.6"),dataframe_join$area,value = T))] <-
    paste0(dataframe_join$species[dataframe_join$stock == "Bycatch/Unknown" &
                                    dataframe_join$species %in% ICCAT_Bycatch &
                                    dataframe_join$area %in%
                                    unique(grep(paste0("34.4.1|34.3.3|34.3.4|34.3.5|34.3.6|41.|47.|48.6"),dataframe_join$area,value = T))],"-SE")


  if(reduce==T & auto.generate==F){
    dataframe_red <- dataframe_join %>% group_by(ship_ID, stock) %>% summarise(landings = sum(landkg)) %>% ungroup()
  }
  if(reduce==F & auto.generate==F){
    dataframe_red <- dataframe_join %>% group_by(ship_ID,species,area, stock) %>% summarise(landings = sum(landkg)) %>% ungroup()
  }
  if(reduce==T & auto.generate==T){
    dataframe_red_I <- dataframe_join
    dataframe_red_I$major.area <- NA
    dataframe_red_I$major.area[dataframe_red_I$area %in% unique(grep(paste0("GSA"),dataframe_red_I$area,value = T))] <-strtrim(dataframe_red_I$area[dataframe_red_I$area %in% unique(grep(paste0("GSA"),dataframe_red_I$area,value = T))], 6)
    dataframe_red_I$major.area <- ifelse(grepl("GSA", dataframe_red_I$major.area), gsub(" ", ".", dataframe_red_I$major.area), dataframe_red_I$major.area)
    dataframe_red_I$major.area <- ifelse(grepl("^18|^21|^61|^71|^81|^88|^57|^58|^51|^47|^48|^41|^34|^31|^21|^87|^77|^67", dataframe_red_I$area),
                                         str_sub(dataframe_red_I$area[is.na(dataframe_red_I$major.area)], 1,2),
                                         dataframe_red_I$major.area)
    dataframe_red_I$major.area[is.na(dataframe_red_I$major.area)] <-strtrim(dataframe_red_I$area[is.na(dataframe_red_I$major.area)], 4)

    dataframe_red_I <- dataframe_red_I %>%
      group_by(ship_ID,species,area,major.area, stock) %>%
      summarise(landings = sum(landkg)) %>%
      group_by(species,major.area, stock) %>%
      mutate(total_landings = sum(landings)) %>%
      ungroup()
    unknowns <- dplyr::filter(dataframe_red_I, stock == "Bycatch/Unknown" & total_landings >= threshold.auto.generate)
    unknowns$species_lower <- tolower(unknowns$species)
    unknowns$area_red <- NA
    unknowns$area_red[unknowns$area %in% unique(grep(paste0("GSA"),unknowns$area,value = T))] <-strtrim(unknowns$area[unknowns$area %in% unique(grep(paste0("GSA"),unknowns$area,value = T))], 6)
    unknowns$area_red <- ifelse(grepl("GSA", unknowns$area_red), gsub(" ", ".", unknowns$area_red), unknowns$area_red)
    unknowns$area_red <- ifelse(grepl("^18|^21|^61|^71|^81|^88|^57|^58|^51|^47|^48|^41|^34|^31|^21|^87|^77|^67", unknowns$area),
                                  str_sub(unknowns$area[is.na(unknowns$area_red)], 1,2),
                                  unknowns$area_red)
    unknowns$area_red[is.na(unknowns$area_red)] <-strtrim(unknowns$area[is.na(unknowns$area_red)], 4)
    unknowns_mod <- unite(unknowns, stock, c("species_lower","area_red"),sep = ".",remove = F)

    dataframe_red_II <- anti_join(dataframe_red_I,unknowns_mod,by=c("species","area"))

    unknowns_mod_red <- unknowns_mod %>% group_by(ship_ID,stock) %>% summarise(landings=sum(landings)) %>% ungroup()
    dataframe_red_III <- dataframe_red_II %>% group_by(ship_ID,stock) %>% summarise(landings=sum(landings)) %>% ungroup()

    dataframe_red <- rbind(unknowns_mod_red,dataframe_red_III)
  }
  if(reduce==F & auto.generate==T){
    dataframe_red_I <- dataframe_join
    dataframe_red_I$major.area <- NA
    dataframe_red_I$major.area[dataframe_red_I$area %in% unique(grep(paste0("GSA"),dataframe_red_I$area,value = T))] <-strtrim(dataframe_red_I$area[dataframe_red_I$area %in% unique(grep(paste0("GSA"),dataframe_red_I$area,value = T))], 6)
    dataframe_red_I$major.area <- ifelse(grepl("GSA", dataframe_red_I$major.area), gsub(" ", ".", dataframe_red_I$major.area), dataframe_red_I$major.area)
    dataframe_red_I$major.area <- ifelse(grepl("^18|^21|^61|^71|^81|^88|^57|^58|^51|^47|^48|^41|^34|^31|^21|^87|^77|^67", dataframe_red_I$area),
                                         str_sub(dataframe_red_I$area[is.na(dataframe_red_I$major.area)], 1,2),
                                         dataframe_red_I$major.area)
    dataframe_red_I$major.area[is.na(dataframe_red_I$major.area)] <-strtrim(dataframe_red_I$area[is.na(dataframe_red_I$major.area)], 4)

    dataframe_red_I <- dataframe_red_I %>%
      group_by(ship_ID,species,area,major.area, stock) %>%
      summarise(landings = sum(landkg)) %>%
      group_by(species,major.area, stock) %>%
      mutate(total_landings = sum(landings)) %>%
      ungroup()
    unknowns <- dplyr::filter(dataframe_red_I, stock == "Bycatch/Unknown" & total_landings >= threshold.auto.generate)
    unknowns$species_lower <- tolower(unknowns$species)
    unknowns$area_red <- NA
    unknowns$area_red[unknowns$area %in% unique(grep(paste0("GSA"),unknowns$area,value = T))] <-strtrim(unknowns$area[unknowns$area %in% unique(grep(paste0("GSA"),unknowns$area,value = T))], 6)
    unknowns$area_red <- ifelse(grepl("GSA", unknowns$area_red), gsub(" ", ".", unknowns$area_red), unknowns$area_red)
    unknowns$area_red <- ifelse(grepl("^18|^21|^61|^71|^81|^88|^57|^58|^51|^47|^48|^41|^34|^31|^21|^87|^77|^67", unknowns$area),
                                         str_sub(unknowns$area[is.na(unknowns$area_red)], 1,2),
                                  unknowns$area_red)
    unknowns$area_red[is.na(unknowns$area_red)] <-strtrim(unknowns$area[is.na(unknowns$area_red)], 4)
    unknowns_mod <- unite(unknowns, stock, c("species_lower","area_red"),sep = ".",remove = F)

    dataframe_red_II <- anti_join(dataframe_red_I,unknowns_mod,by=c("species","area"))

    unknowns_mod_red <- unknowns_mod %>% group_by(ship_ID,species,area, stock) %>% summarise(landings=sum(landings)) %>% ungroup()
    dataframe_red_III <- dataframe_red_II %>% group_by(ship_ID,species,area, stock) %>% summarise(landings=sum(landings)) %>% ungroup()

    dataframe_red <- rbind(unknowns_mod_red,dataframe_red_III)
  }

  dataframe_red_minimals <- dataframe_red %>%
    group_by(ship_ID,stock) %>%
    summarise(landings=sum(landings)) %>%
    group_by(ship_ID) %>%
    mutate(share_stock=landings/sum(landings)) %>%
    dplyr::filter(share_stock >= (min.share/100))

  dataframe_red <- dplyr::filter(dataframe_red, stock %in% dataframe_red_minimals$stock)


  suppressWarnings(return(dataframe_red))
}



#### 21) Data preparation function ####
#' @title Function to automatically separate a fleet data frame by gear class and assign ICES and ICCAT stocks.
#'
#' @description This function separates given fleet data by gear class and creates separate, conveniently named and ready-to-use data frames.
#' It also assigns ICES stocks to catch data based on the caught species and the FAO fishing area where it was caught.
#' Due to limits in available data, especially on a fine spatial resolution, not all ICES-stocks can be taken into account. Therefore, stocks of Nephrops norvegicus
#' and Ammodytes spp. are based on EU definitions, not on ICES definitions. Additionally, the function automatically creates stocks based on species name and FAO area for species
#' which are caught in larger quantities (> 100 kg) but are not part of defined ICES-stock. This feature can be shut off if the user decides not to apply the procedure and
#' use only defined ICES-stocks.
#' @param fleetdata The underlying fleet data frame. With the following arguments, the columns are specified
#' @param vessel_ID The unquoted name of the column containing the vessel identifier
#' @param shiplength The unquoted name of the column containing the length of the vessel in m
#' @param gear The unquoted name of the column containing the code for the gear class
#' @param species The unquoted name of the column containing the three-digit code for the species caught
#' @param area The unquoted name of the column containing the code for the area in ICES/FAO coding (e.g., 27.4.a or GSA 16)
#' @param catch The unquoted name of the column containing the catch in kg
#' @param reduce Indicates, whether or not the resulting data frame is reduced to Ship ID, ICES stock and catch weight. Defaults to TRUE.
#' Turned to FALSE, the resulting data frame will resemble the original data with the ICES-Stock column added. This format is used for control and correction purposes,
#' but not for further calculations.
#' @param auto.generate Indicates whether or stocks should be automatically generated for
#' a) species, which are not assessed in ICES-Stocks or
#' b) assessed by ICES, but caught out of stock-managed areas.
#' The automatically generated stocks comprise the species name and the FAO area.
#' The relevant quantity is defined by the argument threshold.auto.generate
#' @param threshold.auto.generate Threshold of automatic generation of ICES stocks. Only relevant if auto.generate = T. Defaults to 100.
#' @param min.share The minimal share a stock has to have on at least one vessels catch to be included in the stock dataframe. Defaults to 0, i.e. every stock is retained by default.
#' @keywords data preparation
#' @export segmentation_datapreparation
#' @examples
#' segmentation_datapreparation(data = fleet_exmpl,
#'                               vessel_ID = ship_ID,
#'                               shiplength = loa,
#'                               gear = main_gear,
#'                               species = spec_code,
#'                               area = fao_area,
#'                               catch = landkg)

segmentation_datapreparation <- function(fleetdata,vessel_ID,shiplength, gear,species,area,catch,
                                         reduce=T,auto.generate=T,threshold.auto.generate=100, min.share=0){
  data <- fleetdata
  shiplength <<- data %>%
    dplyr::select({{vessel_ID}},{{shiplength}})  %>%
    rename(vessel_ID = {{vessel_ID}}, shiplength = {{shiplength}}) %>%
    unique()

  data$gear <- as.factor(data$gear)

  suppressMessages(
    data_red <- data %>%
      group_by({{vessel_ID}},{{gear}},{{species}},{{area}}) %>%
      summarise(landings = sum({{catch}})) %>%
      ungroup())
  suppressMessages(
    for (i in seq_along(levels(data$gear)))   {
      data_split <- data_red %>%
        group_split({{gear}},.keep = F)
      temp <- data_split[[i]]
      tempII <- assign_stocks(temp, reduce = reduce, auto.generate = auto.generate,
                              threshold.auto.generate = threshold.auto.generate, min.share = min.share)
      assign(levels(data$gear)[i],tempII,.GlobalEnv)
      remove(temp,tempII)
    }
  )
  cat(str_to_title(numbers_to_words(n_distinct(data$gear))),"gear class dataframes and a shiplength dataframe have been created.")

}

