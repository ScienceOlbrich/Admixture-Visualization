library(gridExtra)
# TAKE CARE OF REQUIRED LIBRARIES -----------------------------------------
 list.of.packages <- c("grid", "gridExtra", "tidyverse","tibble","magrittr")
 new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
 if(length(new.packages)) install.packages(new.packages)

# ROH FUNCTIONS -----------------------------------------------------------
 
#' Prepare Plink-ROH-output for plotting
#' @param: ROH.list, Plink output - data.table::fread(input="EGY_ROH.hom", header=TRUE)
#' @param: pop.info, Sample to Group mapping with two columns
#' @return: classified ROHs long format - colnames: "SAMPLE", "class", "MB", "POPULATION"
 prepROHPlot <- function(ROH.list, pop.info, MB.threshold=0.0) {
   require(tibble)
   require(dplyr)
   ## 1. add 'MB' column to the list
   ROH.list$MB <- ROH.list[[grep("KB",colnames(ROH.list), ignore.case=TRUE)]] / 1000
   ## drop everything that is smaller than 'MB.threshold' - DEFAULT: 0.0
   #ROH.list <- ROH.list[ROH.list$MB >= MB.threshold, ]
   ROH.list$MB[ROH.list$MB < MB.threshold] <- 0
   ## 2. add 'class' column to the list
   ROH.list$class <- factor(sapply(ROH.list$MB, function(value) ifelse(value <= 0.155, "short", ifelse(value <= 1.606, "medium", "long"))), levels=c("short","medium","long","all"))
   ## 3. equalize colnames in lists
   names(ROH.list)[1:2] <- c("SAMPLE","POPULATION")
   ## get sum of ROHs per individual
   df <- dplyr::left_join(aggregate(MB ~ SAMPLE + class, ROH.list, sum), pop.info, by="SAMPLE")
   ## get sum of all classes
   tmp <- dplyr::left_join(data.frame("SAMPLE"=unique(df$SAMPLE), "class"="all", stringsAsFactors=FALSE), aggregate(MB ~ SAMPLE + POPULATION, df, sum), by="SAMPLE")
   ## join rest of the inforamtion into it
   tmp <- dplyr::left_join(tmp, pop.info, by=c("SAMPLE","POPULATION"))
   
   tmp <- tmp[,names(df)]
   ## join frame for final output
   df <- rbind.data.frame(df, tmp)
   
   return(df)
 }
 
#' Plot classified ROHs for populations - prints two *.pdf for classiificatiion and frequency
#' @param: plot.df, output of prepROHPlot()
#' @param: ptitle, part of output-filename
#' @param: box.plot, TRUE to construct box-plots instead of violin-plots
#' @param: p.value, TRUE to calculate p.values
#' @return: list containing plots 'classification' and 'frequency'
plotROH <- function(plot.df, ptitle="untitled", box.plot=FALSE, p.value=FALSE) {
  library(ggpubr)
  o <- function(x) {
  subset(x, x == max(x) | x == min(x))
}

f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
  # define things to compare - we'll go with groups aka 'IID' for now
  comp.select <- "POPULATION"
  comp.text <- "Populations"
  c.vector <- RColorBrewer::brewer.pal(length(unique(plot.df[[eval(comp.select)]])), "Dark2")
  my_comparisons <- combn(unique(plot.df[[eval(comp.select)]]), m=2, function(x) c(x), simplify = FALSE)
  if(box.plot) {
    p_ROH_all <- ggplot(plot.df, aes(x = plot.df[,eval(comp.select)], y = MB)) + 
      geom_boxplot(notch = FALSE) +
      stat_summary(fun.data=f, geom="boxplot") + 
  	  stat_summary(fun = o, geom="point") +
  	  stat_boxplot(geom='errorbar',coef=1.5) + # coef: whiskers at 1.5 IQR
      #scale_y_continuous(trans='log2') + 
      geom_jitter(aes(color=plot.df[,eval(comp.select)]), alpha=0.25, position=position_jitter(0.2)) +
      scale_color_manual(values = c.vector) +
      facet_grid(rows=vars(class), scales="fixed", space="fixed") +
      # Remove x axis title
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      theme(axis.title.x = element_blank()) + 
      theme(legend.position="bottom") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      labs(colour = "Population", y = "Total length of ROHs (Mb)") + scale_y_continuous(trans = "log2")
    #facet_wrap(~ class, ncol=2, nrow=2) +
  } else {
    ## 1. facet_grid plot all the groups
    p_ROH_all <- ggpubr::ggviolin(plot.df, x = eval(comp.select), y = "MB", trim = T,add.params = list(fill = "white"),
                                  fill = eval(comp.select), add = c("boxplot"), palette = c.vector, alpha = 0.65,
                                  ylab = "Total length of ROHs (Mb)", xlab = comp.text,
                                  title = "All ROH comparison") +
      #scale_y_continuous(trans='log2') +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1)) +
      theme(axis.title.x = element_blank()) + 
      theme(legend.position="bottom") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      labs(fill = comp.text) + 
      facet_grid(rows=vars(class), scales="fixed", space="fixed") 
    #facet_wrap(~ class, ncol=2, nrow=2) +
  }
  
  ## 2. add statistical comparison
  if(p.value) {
    p_ROH_all <- p_ROH_all + ggpubr::stat_compare_means(comparisons = my_comparisons) + # Add pairwise comparisons p-value
      ggpubr::stat_compare_means() 
  }
  
  ### plot frequency
  pops <- unique(plot.df[,eval(comp.select)])
  popfreq.list <- list()
  for( pop.idx in 1:length(pops) ) {
    ## 1. select population to plot
    tmp <- plot.df[which(plot.df[,eval(comp.select)] %in% pops[pop.idx] & !(plot.df[,"class"] %in% "all")),]
    ## 2. get number of samples in population
    sample.scale <- length(unique(tmp$SAMPLE))
    
    popfreq.list[[eval(pops[pop.idx])]] <- ggplot(tmp, aes(MB, color=eval(comp.select))) +
      geom_histogram(aes(y=..count../ sample.scale), alpha=0.75, bins=100) +
      scale_color_manual(values = c.vector[pop.idx]) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=0, hjust=1)) +
      theme(axis.title.x = element_blank()) +
      theme(axis.title.y = element_blank()) + 
      theme(legend.position = "none") + 
      labs(title=pops[pop.idx])
    
    
  }
  
  n <- length(popfreq.list)
  nCol <- floor(sqrt(n))
  pdf(file=paste("ROH_frequency_",ptitle,".pdf",sep=""), width=6 * nCol, height=length(unique(plot.df$class))*4 / nCol, onefile=FALSE)
  popfreq.plot <- grid.arrange(grobs=popfreq.list, ncol=nCol, left="ROH frequency scaled by number of samples")
  dev.off()
  
  pdf(file=paste("ROH_classification_",ptitle,".pdf",sep=""), width=6, height=length(unique(plot.df$class))*6, onefile=FALSE)
  print(p_ROH_all)
  dev.off()
  
  ## prepare return
  p.list <- list()
  p.list[["classification"]] <- p_ROH_all
  p.list[["frequency"]] <- popfreq.plot
  return(p.list)
}


# ADMIXTURE FUNCTIONS -----------------------------------------------------

#' Select most important clusters of a Population in ADMIXTURE output
#' @param: plot_obj, Cluster List from ADMIXTURE
#' @param: aggCol, column to stratify by
#' @param: keeper, name of population to keep and evaluate by
#' @param: cutoff, significance cutoff
#' @param: select by max (TRUE) or by value I=[0,1]
#' @return: Cluster List sorted by importance of clusters containing only selected populations
selectPopulation <- function(plot_obj, aggCol="POPULATION_DATASET", keeper="Egyptian_EGYPTWGS", cutoff=0.05, max=TRUE) {
  ### DEBUG ##
  # plot_obj <- admix_plot
  ### DEBUG ###
  
  ## select by population cluster maximum
  #eval_population <- plot_obj[2614:dim(plot_obj)[1],]
  eval_population <- plot_obj[which(plot_obj[,aggCol] == eval(keeper)),]
  clusters = grep("Cluster", names(eval_population), ignore.case = TRUE)
  cluster2keep <- names(colMeans(eval_population[,clusters])[colMeans(eval_population[,clusters]) > cutoff])
  ### le DEBUG
  print(paste("With a cutoff of: ",cutoff," the following clusters remain: ",
              paste(cluster2keep, collapse=", "), sep=""))
  
  #idx <- c(cluster2keep, names(plot_obj)[which(!names(plot_obj) %in% names(plot_obj[,..clusters]))])
  #plot_obj <- plot_obj[,..idx]
  clusters = grep("Cluster", names(plot_obj), ignore.case = TRUE)
  tmp <- aggregate(list(plot_obj[,clusters]),  by=list(unlist(plot_obj[,aggCol])), mean)
  
  ##tmp <- aggregate(list(plot_obj[,..cluster2keep]),  by=list(unlist(plot_obj[,..aggCol])), mean)
  ## me being pedantic
  names(tmp)[1] <- eval(aggCol)
  ## drop to keep
  tmp <- tmp[which(!tmp[,aggCol] %in% keeper),]
  
  if(max == TRUE) {
    ## le DEBUG
    print("Selecting for populations with max value per cluster.")
    ## select max. in each remaining column
    #clusters <- clusters + 1
    #toKeep <- sapply(clusters, function(idx) which.max( tmp[,idx]))
    toKeep <- sapply(cluster2keep, function(idx) which.max( tmp[,idx]))
    
    ## extract names and put it together again
    toKeep <- c(tmp[toKeep, aggCol], keeper) 
  } else if(class(max) == "numeric") {
    ## le DEBUG
    print(paste("Selecting for populations with mean >= ",max," in either cluster.", sep=""))
    #clusters <- clusters + 1
    #toKeep <- unique(unlist(sapply(clusters, function(idx) which(tmp[,idx] >= max))))
    toKeep <- unique(unlist(sapply(cluster2keep, function(idx) which(tmp[,idx] >= max))))
    
    toKeep <- c(tmp[toKeep, aggCol], keeper) 
  }
  group.select <- sapply(1:dim(plot_obj)[1],function(x) ifelse(plot_obj[x,aggCol] %in% toKeep, TRUE,FALSE))
  plot_obj <- plot_obj[group.select,]
  
  return(plot_obj)
}

#' construct various types of admixture barplots
#' @param: admix.obj, Cluster List from ADMIXTURE
#' @param: aggCol, column to stratify by
#' @param: p.title, output filename and plot-title
#' @param: grid, TRUE use facet_grid and stratify by 'aggCol' else use facet_wrap with 'aggCol'
#' @param: plot, TRUE plot to file, else just return plot object
#' @param: colorScheme, TRUE use colors provided by cols variable: 
#' you probably want to adjust this section to  your liking
#' @param: cols, provide two colors for gradient coloring
#' @return: Cluster List sorted by importance of clusters containing only selected populations
admixtureBarPlot <- function(admix.obj, aggCol="POPULATION", p.title="noTitle", grid=TRUE, plot=TRUE, colorScheme=TRUE, cols=c("green","blue")) {
  require(ggplot2)
  ### DEBUG ###
  # admix.obj <- admix_plot
  ### DEBUG ###
  
  getCols = c(names(admix.obj)[grep("Cluster", names(admix.obj), ignore.case = TRUE)],"SAMPLE",eval(aggCol))
  admix.obj <- admix.obj[,getCols]
  
  admix.obj.long <- reshape2::melt(admix.obj, id.vars=c("SAMPLE",eval(aggCol)))
  admix.obj.long$value <- as.numeric(admix.obj.long$value)
  ## Adjust facet labels
  admix.obj.long[[eval(aggCol)]] <- factor(admix.obj.long[[eval(aggCol)]])
  ## resort cluster levels
  admix.obj.long$variable <- factor(admix.obj.long$variable, 
                                    levels = c("Cluster_4","Cluster_22","Cluster_7","Cluster_2","Cluster_6",
                                               "Cluster_21","Cluster_17","Cluster_9","Cluster_15","Cluster_13","Cluster_20",
                                               "Cluster_3","Cluster_1","Cluster_18","Cluster_5","Cluster_14","Cluster_24",
                                               "Cluster_10","Cluster_11","Cluster_16","Cluster_19","Cluster_8","Cluster_23",
                                               "Cluster_12"))
  
  if(colorScheme) {
    pal = grDevices::colorRampPalette(cols)
    cols = pal(length(unique(admix.obj.long$variable)))
    cols[1] <- "#FF7F0E"
    cols[2] <- "#2CA02C"
    cols[3] <- "#1F77B4"
    cols[4] <- "#F8DC0B"
    #cols[5] <- "#930000"
  } else {
    ## do nothing here
    cols[1] <- "#FF7F0E"
    cols[2] <- "#2CA02C"
    cols[3] <- "#1F77B4"
    cols[4] <- "#F8DC0B"
    #cols[5] <- "#930000"
  }
  # set #cols for facet wrap
  n <- length(unique(admix.obj.long[[eval(aggCol)]]))
  nCol <- floor(sqrt(n))
  
  if(grid == TRUE) {     ## plot in grid mode
    # Plot admixture barplot 
    admix.bar = ggplot(data=admix.obj.long, aes(x=SAMPLE, y=value, fill=variable))+
      geom_bar(stat = "identity")+
      scale_y_continuous(expand = c(0,0))+
      facet_grid(as.formula(paste("~", aggCol)), scales="free", space="free") +
      geom_bar(stat = "identity", width = 0.9) +
      scale_fill_manual(values = cols)+
      ylab("Admixture proportion") +
      # xlab("Individual")+
      theme(axis.text.x = element_blank(), # = element_text(angle = 90, hjust=1),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            strip.text = element_text(colour="black", size=12),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            legend.position = "none", #"top",
            legend.title = element_blank(),
            legend.text = element_text(size = 12))
  } else { ## plot in facet wrap mode
    # Plot admixture barplot 
    admix.bar = ggplot(data=admix.obj.long, aes(x=SAMPLE, y=value, fill=variable))+
      geom_bar(stat = "identity")+
      scale_y_continuous(expand = c(0,0) )+ #,limits=c(0,1)) +
      geom_bar(stat = "identity", width = 0.9) +
      facet_wrap(as.formula(paste("~", aggCol)), scales = "free_x", ncol = nCol) +
      #facet_wrap(~POPULATION_DATASET, scales = "free", ncol = nCol) +
      scale_fill_manual(values = cols)+
      ylab("Admixture proportion")+
      # xlab("Individual")+
      theme(axis.text.x = element_blank(), # = element_text(angle = 90, hjust=1),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            strip.text = element_text(colour="black", size=12),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            legend.position = "none", #"top",
            legend.title = element_blank(),
            legend.text = element_text(size = 12))
  }
  
  # if(plot == TRUE & grid == FALSE) {
  #   pdf(file=paste("Admixture",p.title,"barplot.pdf",sep="_"), width=12 * nCol, height=12 * (n/nCol), onefile=FALSE)
  #   print(admix.bar)
  #   dev.off()
  # } else if(plot == TRUE & grid == TRUE) {
  #   pdf(file=paste("Admixture",p.title,"barplot_grid.pdf",sep="_"), width=12 * n, height=12, onefile=FALSE)
  #   print(admix.bar)
  #   dev.off()
  # }
  
  if(plot == TRUE & grid == FALSE) {
    pdf(file=paste("Admixture",p.title,"barplot.pdf",sep="_"), width=70.7, height=100, onefile=FALSE)
    print(admix.bar)
    dev.off()
  } else if(plot == TRUE & grid == TRUE) {
    pdf(file=paste("Admixture",p.title,"barplot_grid.pdf",sep="_"), width=70.7, height=100, onefile=FALSE)
    print(admix.bar)
    dev.off()
  }
  
  return(admix.bar)
}

#' construct admixture pie-charts for all populations and save to list for further usage
#' @param: admix.obj, Cluster List from ADMIXTURE to plot pies for
#' @param: colorScheme, TRUE use colors provided by cols variable: 
#' you probably want to adjust this section to  your liking
#' @param: cols, provide two colors for gradient coloring
#' @param: radius, how big are the pies supposed to be, may be changed later but is convenient here
#' @return: list of pie-plot-objects
admixturePieHelper <- function(admix.obj, colorScheme=TRUE, cols=c("green","blue"), radius=3) {
  ### DEBUG ###
  # admix.obj <- admix_plot_important
  ### DEBUG ###
  # Define a function to plot pie charts using ggplot for each site
  pie_charts = function(admix_df, site, cols){
    # admix_df = dataframe in long format of admixture proportions per site 
    # site = string 
    # cols = vector of colours of length(clusters)
    ggplot(data = subset(admix_df, Group.1 == site),
           aes(x = "", y = value, fill = variable))+
      #geom_bar(width = 1, stat = "identity", colour = "#4f4e4d", show.legend = FALSE) +
      geom_bar(width = 2, stat = "identity", colour = "#919191", size = 0, show.legend = FALSE) +
      #geom_bar(width = 1, stat = "identity", show.legend = FALSE) +
      
      coord_polar(theta = "y", start=0) +
      scale_fill_manual(values = cols) +
      theme(legend.position="none") +
      theme_void()
  }
  
  ## aggregate cluster info
  clusters = grep("Cluster", names(admix.obj), ignore.case = TRUE)
  avg_admix = aggregate(admix.obj[, clusters], list(admix.obj$POPULATION_DATASET), mean)
  
  # Order alphabetically by site
  avg_admix = avg_admix[order(as.character(avg_admix$Group.1)), ]
  # Convert dataframe from wide to long format
  avg_admix = reshape2::melt(avg_admix, id.vars = "Group.1")
  avg_admix$variable <- factor(avg_admix$variable, 
                               levels = c("Cluster_4","Cluster_22","Cluster_7","Cluster_2","Cluster_6",
                                          "Cluster_21","Cluster_17","Cluster_9","Cluster_15","Cluster_13","Cluster_20",
                                          "Cluster_3","Cluster_1","Cluster_18","Cluster_5","Cluster_14","Cluster_24",
                                          "Cluster_10","Cluster_11","Cluster_16","Cluster_19","Cluster_8","Cluster_23",
                                          "Cluster_12"))
  avg_admix <- avg_admix[order(avg_admix$variable),]
  
  if(colorScheme) {
    pal = colorRampPalette(cols)
    cols = pal(length(unique(avg_admix$variable)))
  } else {
    ## do nothing here
  }
  # Apply function to all sites using for loop
  subsites <- unique(as.character(avg_admix$Group.1))
  pies = list()
  for (i in subsites){
    pies[[i]] = pie_charts(admix_df = avg_admix, site = i, cols = cols) 
  }
  
  r.frame <- list()
  r.frame[["pies"]] <- pies
  r.frame[["avg_admix"]] <- avg_admix
  r.frame[["radius"]] <- radius
  
  return(r.frame)
  
}

#' Given a Set of Points and Box sizes, find locations
#' Written by @zachp, updated by @slowkow - modified by Michael Olbrich
#' is actually not used ATM
findboxes <- function(
  df, xcol, ycol,
  boxsize,
  point_padding_x, point_padding_y,
  xlim, ylim,
  force = 1e-7, maxiter = 20000
) {
  
  ### DEBUG ###
  # df <- coords
  # xcol <- "Lon"
  # ycol <- "Lat"
  # boxsize <- "radius"
  # point_padding_x <- 0
  # point_padding_y <- 0
  # xlim = c(-7, 86)
  # ylim = c(-30, 62)
  ### DEBUG ###
  
  df <- as.data.frame(df)
  # x and y posiitons as a dataframe
  posdf <- df[c(xcol, ycol, boxsize)]
  
  #returnd a df where columns are points
  boxdf <- apply(posdf,1,function(row) { xval <- row[xcol]
  yval <- row[ycol]
  return(c(xval, 
           yval, 
           xval + (row[boxsize] * 2), 
           yval + (row[boxsize] * 2)))})                                       
  
  # columns are x1,y1,x2,y2
  posdf <- df[c(xcol, ycol)]
  boxmatrix <- as.matrix(t(boxdf))
  
  moved <- ggrepel:::repel_boxes(
    data_points = as.matrix(posdf),
    point_padding_x = point_padding_x,
    point_padding_y = point_padding_y,
    boxes = boxmatrix,
    xlim = xlim,
    ylim = ylim,
    hjust = 0.1,
    vjust = 0.1,
    force = force,
    maxiter = maxiter
  )
  
  posdf <- df[c("POPULATION_DATASET","POPULATION", boxsize, xcol, ycol)]
  
  
  finaldf <- cbind(posdf, moved)
  names(finaldf) <- c("POPULATION_DATASET","POPULATION","radius", "x1", "y1", "x2", "y2")
  return(finaldf)
}

#' put pies on the map
#' either use findboxes() or supply with coordinates in variable df1
#' @param: df1, containing two sets of coordinates as x1,y1 and x2,y2 
#' denoting coordinates to point to and the actual plotting coordinates (x2,y2)
mapPies <- function(locations, pies, force=1e-06, maxiter=2000, df1=NULL) {
  ### DEBUG ###
  # locations <- pop_loc
  # pies <- list(pies_unimportant,pies_remaining,pies_top5,pies_egypt)
  # df1 <- set.df1
  ### DEBUG ###
  
  if(is.null(names(pies))) {
    tmp.p <- NULL
    tmp.avg <- NULL
    tmp.rad <- NULL
    coords <- NULL
    subsites <- NULL
    for(idx in 1:length(pies)) {
      tmp.p <- c(tmp.p, pies[[idx]][["pies"]])
      tmp.avg <- rbind.data.frame(tmp.avg, pies[[idx]][["avg_admix"]])
      tmp.rad <- c(tmp.rad, pies[[idx]][["radius"]])
      
      tmp.subsites <- unique(as.character(pies[[idx]][["avg_admix"]]$Group.1))
      tmp.coords = locations[locations$POPULATION_DATASET %in% tmp.subsites, ]
      tmp.coords = tmp.coords[order(tmp.coords$POPULATION_DATASET), ] 
      tmp.coords$radius <- pies[[idx]][["radius"]]
      
      coords <- rbind.data.frame(coords, tmp.coords)
      subsites <- c(subsites, tmp.subsites)
    }
    
    pies[["pies"]] <- tmp.p
    pies[["avg_admix"]] <- tmp.avg
    
  } else {
    subsites <- unique(as.character(pies[["avg_admix"]]$Group.1))
    coords = locations[locations$POPULATION_DATASET %in% subsites, ]
    coords = coords[order(coords$POPULATION_DATASET), ] 
    coords$radius <- pies[["radius"]]
  }
  
  # df1 <- findboxes(coords,xcol = "Lon", ycol = "Lat", boxsize = "radius",point_padding_x = 0,
  #                   point_padding_y = 0,force = force, maxiter = maxiter,xlim = c(-15, 85),ylim = c(-10, 60))
  
  print(df1)
  coord.list = list()
  for (i in subsites){
    #coord.list[[i]] = c(subset(coords, POPULATION_DATASET == i)$Lon, subset(coords, POPULATION_DATASET == i)$Lat)
    coord.list[[i]] = c(subset(df1, POPULATION_DATASET == i)$x2, subset(df1, POPULATION_DATASET == i)$y2, subset(df1, POPULATION_DATASET == i)$radius)
  }
  #coord.list$radius <- coord.list$radius / 2
  # Convert ggplot pie charts to annotation_custom layers
  pies.ac = list()
  for (i in subsites){
    pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies$pies[[i]]),
                                     xmin = coord.list[[i]][[1]] - coord.list[[i]][[3]],
                                     xmax = coord.list[[i]][[1]] + coord.list[[i]][[3]],
                                     ymin = coord.list[[i]][[2]] - coord.list[[i]][[3]],
                                     ymax = coord.list[[i]][[2]] + coord.list[[i]][[3]])
  }
  
  ret <- list()
  ret[["pies"]] <- pies.ac
  ret[["df1"]] <- df1
  return(ret)
}

#' get the 'rnaturalearth' map of the planet with given resolution and map-limits
#' @param: m.size, one of "small, medium, large"
#' @param: limits, 4-element vector min-longitude, max-longitude, min-latitude, max-latitude
#' @param: grid, TRUE return map with fine grid for debugging and adjustments, else no grid
#' @return: map grob for further use
getWorldMap <- function(m.size="medium", limits=FALSE, grid=FALSE) {
  world <- ne_countries(scale = eval(m.size), returnclass = "sf")
  ## landmarks and countries
  world_points<- st_centroid(world)
  world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))
  if(!limits) {
    if(!grid) {
      basemap <- ggplot(data = world) + 
        geom_sf(fill= "antiquewhite") + 
        geom_text(data= world_points, aes(x=X, y=Y, label=name), color = gray(.5), fontface = "italic", size=4, check_overlap = TRUE) + 
        annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 2) + 
        annotation_scale(location = "bl", width_hint = 0.5) + 
        annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
        #coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97), expand = FALSE) + 
        xlab("Longitude") + ylab("Latitude") + 
        ggtitle("Map of Earth-Population-Admixture") + 
        theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue"))
    } else {
      basemap <- ggplot(data = world) + 
        geom_sf(fill= "antiquewhite") + 
        geom_text(data= world_points, aes(x=X, y=Y, label=name), color = gray(.5), fontface = "italic", size=4, check_overlap = TRUE) + 
        annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 2) + 
        annotation_scale(location = "bl", width_hint = 0.5) + 
        annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
        #coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97), expand = FALSE) + 
        xlab("Longitude") + ylab("Latitude") + 
        ggtitle("Map of Earth-Population-Admixture") + 
        theme(panel.grid.major = element_line(color = gray(0),linetype = "solid", size = 0.5),
              panel.grid.minor = element_line(color = gray(0),linetype = "solid", size = 0.5),
              panel.ontop = TRUE,
              panel.background = element_rect(fill = NA))      
    }
  } else {
    if(!grid) {
      basemap <- ggplot(data = world) + 
        geom_sf(fill= "antiquewhite") + 
        geom_text(data= world_points, aes(x=X, y=Y, label=name), color = gray(.5), fontface = "italic", size=4, check_overlap = TRUE) + 
        annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 2) + 
        annotation_scale(location = "bl", width_hint = 0.5) + 
        annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
        coord_sf(xlim = c(limits[1], limits[2]), ylim = c(limits[3], limits[4]), expand = FALSE) + 
        xlab("Longitude") + ylab("Latitude") + 
        ggtitle("Map of Earth-Population-Admixture") + 
        theme(panel.grid.major = element_line(color = gray(.5),linetype = "dashed", size = 0.5),
              panel.grid.minor = element_line(color = gray(.5),linetype = "dashed", size = 0.5),
              panel.background = element_rect(fill = "aliceblue"))
    } else {
      basemap <- ggplot(data = world) + 
        geom_sf(fill= "antiquewhite") + 
        geom_text(data= world_points, aes(x=X, y=Y, label=name), color = gray(.5), fontface = "italic", size=4, check_overlap = TRUE) + 
        annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", fontface = "italic", color = "grey22", size = 2) + 
        annotation_scale(location = "bl", width_hint = 0.5) + 
        annotation_north_arrow(location = "bl", which_north = "true", pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"), style = north_arrow_fancy_orienteering) + 
        coord_sf(xlim = c(limits[1], limits[2]), ylim = c(limits[3], limits[4]), expand = FALSE) + 
        xlab("Longitude") + ylab("Latitude") + 
        ggtitle("Map of Earth-Population-Admixture") + 
        ## Code to produce grid-lines for adjustments
        theme(panel.grid.major = element_line(color = gray(0),linetype = "solid", size = 0.5),
              panel.grid.minor = element_line(color = gray(0),linetype = "solid", size = 0.5),
              panel.ontop = TRUE,
              panel.background = element_rect(fill = NA))
    }
    
  }
  return(basemap)
}


# PCA FUNCTIONS -----------------------------------------------------------

#' Plot 2D-PCA combinations for axes 1 to 4 from prcomp() output
#' @param: pca.evec, eigenvector output of prcomp()
#' @param: pca.eval, eigenvalue output of prcomp()
#' @param: pop.info, metadata to plot and stratify by
#' @param: comp.select, column to stratify by
#' @param: p.title, plot/file title
#' @return: p.list, list with plot-objects named by axis combination, e.g. p12 for axes 1 and 2
plotPopPCA <- function(pca.evec, pca.eval, pop.info, comp.select="POPULATION", p.title="noTitle") {
  require(dplyr); require(ggplot2); require(scatterplot3d)
  
  ### DEBUG ###
  # pca.evec <- pca.vec
  # pca.eval <- pca.val
  # pop.info <- pop_info
  # comp.select="CONTINENT"
  # p.title="EGYPTGSA_BERGSTROEM_1000G"
  ### DEBUG ###
  
  ## set some values for plotting
  p.size <- 3
  p.alpha <- 0.35
  
  names(pca.evec) <- c("SAMPLE",paste("PC",1:(dim(pca.evec)[2]-2),sep=""),"pheno")
  ## calculate percentage of variation in PCs
  pca.eval <- pca.eval[1:(dim(pca.evec)[2]-2)]
  percentVar <- data.frame("axis.name"=paste("PC",1:(dim(pca.evec)[2]-2),sep=""),
                           "axis.percent"=NA, "axis.title"=NA)
  percentVar$axis.percent <- round(100*pca.eval/sum(pca.eval),2)
  percentVar$axis.title <- sapply(1:dim(percentVar)[1], function(idx) paste(percentVar$axis.name[idx]," (",percentVar$axis.percent[idx],"%)",sep=""))
  
  ## join with metadata
  plot.df <- dplyr::left_join(pca.evec, pop.info, by="SAMPLE") 
  ## order by region
  plot.df$CONTINENT[which(plot.df$CONTINENT %in% "Sub Saharan Africa")] <- "Sub-Saharan Africa"
  plot.df <- plot.df[order(plot.df$CONTINENT, decreasing=FALSE),]
  plot.df[,eval(comp.select)] <- factor(plot.df[,eval(comp.select)], levels=unique(plot.df[,eval(comp.select)]))
  
  additional <- plot.df[which(plot.df$POPULATION %in% "Egyptian"),]
  #plot.df <- plot.df[which(!(plot.df$POPULATION %in% "Egyptian")),]
  
  plot.df$alpha <- p.alpha
  plot.df$alpha[which(plot.df$POPULATION %in% "Egyptian")] <- 1
  plot.df$shape <- "normal"
  plot.df$shape[which(plot.df$POPULATION %in% "Egyptian")] <- "special"
  plot.df$shape2 <- 16
  plot.df$shape2[which(plot.df$POPULATION %in% "Egyptian")] <- 1
  
  c.names <- c("America","Asia","Egypt","Europe","Middle East","North Africa","Oceania","South Asia","Sub-Saharan Africa")
  cols <- c("#FFBFD4","#82CBFF","#000000","#2CA02C","#FF7F0E","#1F77B4","#800080","#930000","#F8DC0B")
  names(cols) <- c.names
  #p.colors <- plot.df[!duplicated(plot.df[,c('SUBREGION')]),c("SUBREGION","CONTINENT")]
  p.colors <- plot.df[!duplicated(plot.df[,c('CONTINENT')]),c("SUBREGION","CONTINENT")]
  
  p.colors$color <- sapply(p.colors$CONTINENT, function(idx) cols[idx])
  
  nRow <- 2 #round(sqrt(length(unique(plot.df[,eval(comp.select)]))))
  
  
  
  #### 2D-Projection
  ### PCA12
  p12 <- ggplot(plot.df,aes(x=PC1,y=PC2,color=plot.df[,eval(comp.select)])) + 
    geom_point(alpha=plot.df$alpha, size=p.size, shape=plot.df$shape2) +
    #scale_shape_manual(guide = 'none', values = c("normal" = 16, "special" = 1)) + 
    scale_color_manual(values = p.colors$color) +
    geom_point(data=additional, alpha=p.alpha, color="#000000",size=p.size, shape=1, show.legend = FALSE) + 
    xlab(percentVar$axis.title[1]) + ylab(percentVar$axis.title[2]) +
    ggtitle("WGS") + 
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    #facet_wrap(~ AFRICA_REGION, ncol=3, nrow=3) + 
    #theme(axis.title.x = element_blank()) + 
    theme(legend.position="bottom") +
    theme(axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank()) +
    labs(colour = comp.select) + 
    guides(color=guide_legend(nrow=nRow, byrow=TRUE)) +
    theme(legend.position="bottom", legend.title=element_blank())
  
  ### PCA13
  p13 <- ggplot(plot.df,aes(x=PC1,y=PC3,color=plot.df[,eval(comp.select)],)) + 
    geom_point(alpha=p.alpha, size=p.size) +
    xlab(percentVar$axis.title[1]) + ylab(percentVar$axis.title[3]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    theme(legend.position="bottom") +
    labs(colour = comp.select) + 
    theme(legend.position="bottom") + guides(color=guide_legend(nrow=nRow, byrow=TRUE))
  
  ### PCA14
  p14 <- ggplot(plot.df,aes(x=PC1,y=PC4,color=plot.df[,eval(comp.select)],)) + 
    geom_point(alpha=p.alpha, size=p.size) +
    xlab(percentVar$axis.title[1]) + ylab(percentVar$axis.title[4]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    theme(legend.position="bottom") +
    labs(colour = comp.select) + 
    theme(legend.position="bottom") + guides(color=guide_legend(nrow=nRow, byrow=TRUE))
  
  ### PCA23
  p23 <- ggplot(plot.df,aes(x=PC2,y=PC3,color=plot.df[,eval(comp.select)],)) + 
    geom_point(alpha=p.alpha, size=p.size) +
    xlab(percentVar$axis.title[2]) + ylab(percentVar$axis.title[3]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    theme(legend.position="bottom") +
    labs(colour = comp.select) + 
    theme(legend.position="bottom") + guides(color=guide_legend(nrow=nRow, byrow=TRUE))
  
  ### PCA24
  p24 <- ggplot(plot.df,aes(x=PC2,y=PC4,color=plot.df[,eval(comp.select)],)) + 
    geom_point(alpha=p.alpha, size=p.size) +
    xlab(percentVar$axis.title[2]) + ylab(percentVar$axis.title[4]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    theme(legend.position="bottom") +
    labs(colour = comp.select) + 
    theme(legend.position="bottom") + guides(color=guide_legend(nrow=nRow, byrow=TRUE))
  
  ### PCA34
  p34 <- ggplot(plot.df,aes(x=PC3,y=PC4,color=plot.df[,eval(comp.select)],)) + 
    geom_point(alpha=plot.df$alpha, size=p.size) +
    scale_shape_manual(guide = 'none', values = c("normal" = 16, "special" = 1)) + 
    scale_color_manual(values = p.colors$color) +
    geom_point(data=additional, alpha=p.alpha, color="#000000",size=p.size, shape=1, show.legend = FALSE) + 
    xlab(percentVar$axis.title[3]) + ylab(percentVar$axis.title[4]) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=0, hjust=1)) +
    #facet_wrap(~ AFRICA_REGION, ncol=3, nrow=3) + 
    #theme(axis.title.x = element_blank()) + 
    theme(legend.position="bottom") +
    theme(axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank()) +
    labs(colour = comp.select) + 
    guides(color=guide_legend(nrow=nRow, byrow=TRUE)) +
    theme(legend.position="bottom", legend.title=element_blank())
  # 
  
  
  pdf(file=paste("PCA_by_",comp.select,"_",p.title,".pdf",sep=""), width=12, height=14, onefile=TRUE)
  print(p12)
  print(p13)
  print(p14)
  print(p23)
  print(p24)
  print(p34)
  dev.off()
  
  p.list <- list()
  p.list[["pca12"]] <- p12
  p.list[["pca13"]] <- p13
  p.list[["pca14"]] <- p14
  p.list[["pca23"]] <- p23
  p.list[["pca24"]] <- p24
  p.list[["pca34"]] <- p34
  return(p.list)
}

#' Plot interactive 3D-PCA for axes 1 to 3 from prcomp() output
#' @param: pca.evec, eigenvector output of prcomp()
#' @param: pca.eval, eigenvalue output of prcomp()
#' @param: pop.info, metadata to plot and stratify by
#' @param: comp.select, column to stratify by
#' @param: p.title, plot/file title
#' @return: interactive plot object
plotPopPCA3D <- function(pca.evec, pca.eval, pop.info, comp.select="POPULATION", p.title="noTitle") {
  require(dplyr); require(ggplot2); require(plotly)
  
  ### DEBUG ###
  # pca.evec <- pca.vec
  # pca.eval <- pca.val
  # pop.info <- pop_info
  # comp.select="SUBREGION"
  ### DEBUG ###
  
  ## set some values for plotting
  p.size <- 3
  p.alpha <- 0.35
  
  names(pca.evec) <- c("SAMPLE",paste("PC",1:(dim(pca.evec)[2]-2),sep=""),"pheno")
  ## calculate percentage of variation in PCs
  pca.eval <- pca.eval[1:(dim(pca.evec)[2]-2)]
  percentVar <- data.frame("axis.name"=paste("PC",1:(dim(pca.evec)[2]-2),sep=""),
                           "axis.percent"=NA, "axis.title"=NA)
  percentVar$axis.percent <- round(100*pca.eval/sum(pca.eval),2)
  percentVar$axis.title <- sapply(1:dim(percentVar)[1], function(idx) paste(percentVar$axis.name[idx]," (",percentVar$axis.percent[idx],"%)",sep=""))
  
  ## join with metadata
  plot.df <- dplyr::left_join(pca.evec, pop.info, by="SAMPLE") 
  ## order by region
  plot.df <- plot.df[order(plot.df$CONTINENT, decreasing=FALSE),]
  plot.df[,eval(comp.select)] <- factor(plot.df[,eval(comp.select)], levels=unique(plot.df[,eval(comp.select)]))
  
  nRow <- round(sqrt(length(unique(plot.df[,eval(comp.select)]))))
  # c.vector <- rainbow(length(unique(plot.df[[eval(comp.select)]])))
  # names(c.vector) <- unique(plot.df[[eval(comp.select)]])
  
  c.names <- c("America","Asia","Egypt","Europe","Middle East","North Africa","Oceania","South Asia","Sub Saharan Africa")
  cols <- c("#FFBFD4","#82CBFF","#919191","#2CA02C","#FF7F0E","#1F77B4","#800080","#930000","#F8DC0B")
  names(cols) <- c.names
  
  test <- plot.df[!duplicated(plot.df[,c('SUBREGION')]),c("SUBREGION","CONTINENT")]
  test$p.color <- sapply(test$CONTINENT, function(idx) cols[idx])
  
  #test <- plot.df[!duplicated(plot.df[,c('POPULATION')]),c("POPULATION","CONTINENT")]
  cols.dark <- colorspace::darken(cols, -0.6)
  names(cols.dark) <- c.names
  scales::show_col(cols.dark)
  for(idx in unique(test$CONTINENT)) {
    pal = colorRampPalette(c(cols[idx],cols.dark[idx]))
    tmp.len <- length(which(test$CONTINENT %in% idx))
    tmp.col <- pal(tmp.len + 1)
    test$p.color[which(test$CONTINENT %in% idx)] <- tmp.col[1:tmp.len]
  }
  
  test$p.color[which(test$CONTINENT %in% "Egypt")] <- "black"  #"#919191"
  
  comp.select <- "SUBREGION"
  scales::show_col(test$p.color)
  
  ## PLOTTING from here
  fig <- plot_ly(plot.df, x = ~PC1, y = ~PC2, z = ~PC3, type="scatter3d", mode="markers", 
                 color = ~plot.df[[eval(comp.select)]], colors = test$p.color,
                 symbols = "p", size = 3, alpha=0.95)
  
  
  fig <- fig %>% layout(title = paste(""),#paste(p.title),
                        scene = list(xaxis = list(title = percentVar$axis.title[1],
                                                  gridcolor = 'rgb(255, 255, 255)',
                                                  #range = c(min(plot.df$PC1), max(plot.df$PC1)),
                                                  #type = 'log',
                                                  zerolinewidth = 1,
                                                  ticklen = 5,
                                                  gridwidth = 2),
                                     yaxis = list(title = percentVar$axis.title[2],
                                                  gridcolor = 'rgb(255, 255, 255)',
                                                  #range = c(min(plot.df$PC1), max(plot.df$PC1)),
                                                  #type='log',
                                                  zerolinewidth = 1,
                                                  ticklen = 5,
                                                  gridwith = 2),
                                     zaxis = list(title = percentVar$axis.title[3],
                                                  gridcolor = 'rgb(255, 255, 255)',
                                                  #type = 'log',
                                                  zerolinewidth = 1,
                                                  ticklen = 5,
                                                  gridwith = 2)),
                        paper_bgcolor = 'rgb(243, 243, 243)',
                        plot_bgcolor = 'rgb(243, 243, 243)')
  
  return(fig)
  
}

#' Plot fixed 3D-PCA for axes 1 to 3 from prcomp() output with given rotation and view-angle
#' @param: pca.evec, eigenvector output of prcomp()
#' @param: pca.eval, eigenvalue output of prcomp()
#' @param: pop.info, metadata to plot and stratify by
#' @param: comp.select, column to stratify by
#' @param: p.title, plot/file title
#' @param: phi, view-angle
#' @param: theta, rotation around central axis
#' @return: saves plot to file
plotPopPCA3DFixed <- function(pca.evec, pca.eval, pop.info, comp.select="POPULATION", p.title="noTitle", phi=0, theta=0) {
  require(dplyr); require(ggplot2); require(plot3D)
  library("grid")
  library("ggplotify")
  ### DEBUG ###
  # pca.evec <- pca.vec
  # pca.eval <- pca.val
  # pop.info <- pop_info
  # comp.select <- "SUBREGION"
  # p.title <- "DEVELOPMENT"
  # phi <- 0
  # theta <- 0
  ### DEBUG ###
  
  ## set some values for plotting
  pnt.size <- 1.5
  pnt.alpha <- 0.5
  pnt.type <- 16
  
  ## calculate percentage of variation in PCs and prep axis titles
  names(pca.evec) <- c("SAMPLE",paste("PC",1:(dim(pca.evec)[2]-2),sep=""),"pheno")
  pca.eval <- pca.eval[1:(dim(pca.evec)[2]-2)]
  percentVar <- data.frame("axis.name"=paste("PC",1:(dim(pca.evec)[2]-2),sep=""),
                           "axis.percent"=NA, "axis.title"=NA)
  percentVar$axis.percent <- round(100*pca.eval/sum(pca.eval),2)
  percentVar$axis.title <- sapply(1:dim(percentVar)[1], function(idx) paste(percentVar$axis.name[idx]," (",percentVar$axis.percent[idx],"%)",sep=""))
  
  ## join with metadata and order by region
  plot.df <- dplyr::left_join(pca.evec, pop.info, by="SAMPLE") 
  plot.df$CONTINENT[which(plot.df$CONTINENT %in% "Sub Saharan Africa")] <- "Sub-Saharan\nAfrica"
  plot.df <- plot.df[order(plot.df$CONTINENT, decreasing=FALSE),]
  plot.df[,eval(comp.select)] <- factor(plot.df[,eval(comp.select)], levels=unique(plot.df[,eval(comp.select)]))
  ## point type and sizes
  plot.df$pch <- pnt.type
  plot.df$pch[which(plot.df$POPULATION %in% "Egyptian")] <- 16
  plot.df$cex <- pnt.size
  plot.df$cex[which(plot.df$POPULATION %in% "Egyptian")] <- 3
  
  c.names <- c("America","Asia","Egypt","Europe","Middle East","North Africa","Oceania","South Asia","Sub-Saharan\nAfrica")
  cols <- c("#FFBFD4","#82CBFF","#919191","#2CA02C","#FF7F0E","#1F77B4","#800080","#930000","#F8DC0B")
  names(cols) <- c.names
  p.colors <- plot.df[!duplicated(plot.df[,c('SUBREGION')]),c("SUBREGION","CONTINENT")]
  p.colors$color <- sapply(p.colors$CONTINENT, function(idx) cols[idx])
  
  
  #p.colors <- plot.df[!duplicated(plot.df[,c('POPULATION')]),c("POPULATION","CONTINENT")]
  # cols.dark <- colorspace::darken(cols, -0.9)
  # names(cols.dark) <- c.names
  # scales::show_col(cols.dark)
  # for(idx in unique(p.colors$CONTINENT)) {
  #   pal = colorRampPalette(c(cols[idx],cols.dark[idx]))
  #   tmp.len <- length(which(p.colors$CONTINENT %in% idx))
  #   tmp.col <- pal(tmp.len + 1)
  #   p.colors$p.color[which(p.colors$CONTINENT %in% idx)] <- tmp.col[1:tmp.len]
  # }
  
  ## PLOTTING from here
  pdf(file=paste("PCA",p.title,"fixed.pdf",sep="_"), width=6 , height=6, onefile=FALSE)
  par(mar=c(2.1, 2.1, 2.1, 2.1), xpd=TRUE)
  pca.plot <- with(plot.df, scatter3D(x = PC1, y = PC2, z = PC3, image=TRUE, legend=FALSE,
                                      ## Colors
                                      colvar = as.integer(SUBREGION), col = p.colors$color, alpha=0.5, #bg="#919191",
                                      ## Points
                                      pch = plot.df$pch, cex = plot.df$cex, 
                                      ## Grid
                                      bty = "u",  
                                      col.panel =gray(0.98), expand =1, 
                                      col.grid = "darkgray",
                                      #ticktype = "detailed", 
                                      ## View Angle
                                      phi = phi, theta = theta,
                                      ## Labels
                                      xlab = percentVar$axis.title[1], ylab = percentVar$axis.title[2], 
                                      zlab = percentVar$axis.title[3],  #main = "Admixture PCA",
                                      ## Legend
                                      colkey = FALSE
                                      # colkey = list(at = c(0.5:20), side = 1,
                                      #               addlines = TRUE, length = 1, width = 0.5,
                                      #               labels = p.colors$SUBREGION, las=2)
                                      
                                      # colkey = list(length = 0.2, width = 0.4, shift = 0.15,
                                      #               cex.axis = 0.8, cex.clab = 0.85), lighting = TRUE, lphi = 90,
                                      # clab = c("height","m")
  )) 
  
  p.colors$CONTINENT[which(p.colors$CONTINENT %in% "Sub Saharan Africa")] <- "Sub Saharan\nAfrica"
  
  legend(-0.25, 0.37, legend=unique(p.colors$CONTINENT), bg="white", ncol=2, pch=16, pt.cex=1.5,
         col=unique(p.colors$color), cex=0.8)
  
  dev.off() 
  
}


