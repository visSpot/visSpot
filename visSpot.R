library('SPOT')
library('smoof')
library('ggplot2')
library('scales')
library('patchwork')
library('grid') #Required for multiplot-function in 3rd plot function  
library("GGally")
library('reshape2')
library('xtable')
library('tikzDevice')

visInit <- function(resAgg, x, fun, lower, upper, control, opt){
  
  if(length(resAgg) == 0){index <- 1:2} 
  else{index <- 1:(length(resAgg)+1)}

  lbl.lst <- list()
  
  name <- 'x'
  lbl.lst[[name]] <- deparse(substitute(x))
  name <- 'fun'
  lbl.lst[[name]] <- deparse(substitute(fun))
  name <- 'lower' 
  lbl.lst[[name]] <- deparse(substitute(lower))  
  name <- 'upper'
  lbl.lst[[name]] <- deparse(substitute(upper)) 
  name <- 'control'
  lbl.lst[[name]] <- deparse(substitute(control))
  name <- 'opt'
  lbl.lst[[name]] <- deparse(substitute(opt))
  
  name <- 'overview'
  resAgg[[name]][[tail(index, n=1) - 1]] <- lbl.lst
  
  if(missing(opt) && length(x) > 1){
    assign(
      paste0('res', 
             toString(tail(index, n=1)-1), 
             sep=''
             ),
      
      spot(x, fun, lower, upper, control)
    )
    opt = ''
    
  } else if (length(x) > 1){
    assign(
      paste0('res', 
             toString(tail(index, n=1)-1), 
             sep=''
             ),
      
      spot(x, fun, lower, upper, control, opt)
    )
  } else if (length(x) < 2){
    print("Dimension must greater 1")
  }
  
  name <-  paste0('res', 
                  toString(tail(index, n=1)-1), 
                  sep=''
                  )
  
  resAgg[[name]] <- get(paste0('res', 
                               toString(tail(index, n=1)-1), 
                               sep=''
                               )
                        )  
  
  assign('resAgg', 
         resAgg, 
         envir = .GlobalEnv
         )
}
visInputTable <- function(resAgg){
  
  spot.overview <- data.frame(matrix(unlist(resAgg[[1]]), 
                                     nrow=length(resAgg[[1]]), 
                                     byrow=T))
  
  colnames(spot.overview) <- c('x', 'fun', 'lower', 'upper', 'control', 'opt')
  assign('spot.overview', spot.overview)
  
}
visProgEval <- function(resAgg) {
  
  #https://stackoverflow.com/questions/23901907/
  #create-a-log-sequence-across-multiple-orders-of-magnitude
  
  x_main <- c(-1:-1 %o% 10^(7:0))
  y_main <- 0
  z_main <- c(1:1 %o% 10^(0:7))
  
  xyz_main <- c(x_main, y_main, z_main)
  
  x_minj <- c(-4:-1 %o% 10^(7:0))/4
  y_minj <- 0
  z_minj <- c(1:4 %o% 10^(0:7))/4
  
  xyz_minj <- c(x_minj, y_minj, z_minj)
  
  
  index <- 1:length(resAgg) 
  
  for(i in head(index, -1)){
    assign(paste0('res', 
                  toString(i)
                  ), 
           
           resAgg[[paste0('res', 
                          toString(i), 
                          sep='')]])
  }
  
  for(i in head(index, -1)){
    
    assign("n",
           1: nrow(get(paste0("res", 
                              toString(i), 
                              sep=''))$y)
           )
    
    assign("y", 
           get(paste0("res", 
                      toString(i), 
                      sep=''))$y)
    
    assign("inst", 
           rep(toString(i), 
               times = nrow(get(paste0("res",
                                       toString(i), 
                                       sep=''))$y)
               )
           )
    
    assign(paste0('df',
                  toString(i), 
                  sep=''), 
           data.frame(n,
                      y, 
                      inst)
           )
    
    if (i == 1){
      assign('df', 
             data.frame(n,
                        y, 
                        inst)
             )
      
      df <- df1
      
    } else {
      df <- rbind(df, 
                  get(paste0('df', 
                             toString(i), 
                             sep='')
                      )
                  )
    }
  }
  
  asinh_trans <- trans_new(name = 'asinh', 
                           transform = function(x) asinh(x), 
                           inverse = function(x) sinh(x)
                           )
  
  g <- ggplot()
  
  for(i in head(index, -1)) {
    g <- g + geom_line(data = get(paste0('df',
                                         toString(i), 
                                         sep='')
                                  ),
                       aes(n, 
                           y, 
                           color = inst
                           ), 
                       alpha=0.75, 
                       size = 0.25)
  }
  
  plotNy <- g + 
    scale_y_continuous(trans = asinh_trans, 
                       breaks= xyz_main, 
                       minor_breaks = xyz_minj,
                       labels = function(x) formatC(x, 
                                                    format = "e", 
                                                    digits = 2)
    )+ 
    theme_light() +
    theme(axis.text =element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust = 1), 
          text = element_text(size=6),
          legend.text  = element_text(size = 3),
          legend.key.size = unit(0.25, "lines"))+
    guides(shape = guide_legend(override.aes = list(size = 0.75)
                                )
           )+
    guides(color = guide_legend(override.aes = list(size = 0.75)
                                )
           )
  
  
  plotSBP <- ggplot(df, 
                    aes(x = inst, 
                        y = y, 
                        color = inst)
  )+
    
    geom_boxplot(size=0.25,
                 outlier.shape = 0.75,
                 outlier.color = "black",
                 outlier.size  = 0.25,
                 alpha = 0.01)+
    
    geom_jitter(alpha = 0.5, 
                width=0.2,
                size=0.25) + 
    
    scale_y_continuous(trans = asinh_trans, 
                       breaks= xyz_main, 
                       minor_breaks = xyz_minj,
                       labels = function(x) formatC(x, 
                                                    format = "e", 
                                                    digits = 2
                                                    )
    )+
    
    theme_light() + 
    
    labs(x = 'Spot run instance', fill = "Spot run instance")+
    
    theme(axis.text =element_text(size=8),
          axis.text.x = element_text(angle = 90, hjust = 1), 
          text = element_text(size=6),
          legend.text  = element_text(size = 3),
          legend.key.size = unit(0.25, "lines")
          )+
    guides(shape = guide_legend(override.aes = list(size = 0.75)))+
    guides(color = guide_legend(override.aes = list(size = 0.75)))
  
  assign('visProgEvalout',
         (plotNy + plotSBP)+ 
           plot_annotation(tag_levels = 'A'))
}
visPath2D <- function(resAgg) {
  
  index <- 1:length(resAgg)
  p <- list()
  
  for(i in head(index, -1)){
    
    df <- data.frame(x = resAgg[[paste0('res', 
                                        toString(i),
                                        sep='')
                                 ]][[3]][,1],
                     y = resAgg[[paste0('res', 
                                        toString(i),
                                        sep='')
                                 ]][[3]][,2]
                     )
    
    p[[i]] <- ggplot(df, aes(x, y)
                     )+
      
      scale_x_continuous(
        breaks = pretty_breaks(),
        expand = c(0.1, 0.1),
        labels = function(x) format(x, 
                                    scientific = TRUE
                                    )
      ) +
      
      scale_y_continuous(
        breaks= pretty_breaks(),
        expand = c(0.1, 0.1),
        labels = function(x) format(x, 
                                    scientific = TRUE
                                    )
      ) + 
      
      geom_segment(
        aes(xend=c(tail(x, n=-1), NA), 
            yend=c(tail(y, n=-1), NA)),
        arrow=arrow(length=unit(0.05,"inches"),
                    type = "closed",
                    angle=7.5
                    ),
        size=0.2,
        alpha=0.5
      ) +
      
      geom_point(
        data = df[which.min(resAgg[[paste0('res', 
                                           toString(i),
                                           sep='')
                                    ]][[4]][,1]),1:2], 
        
                    shape=21, color = "darkred", 
                    fill ="red", 
                    stroke = 0.5, 
                    size = 1, 
                    alpha=0.5
      ) +
      
      theme_light() +
      
      ggtitle(paste0('Spot run instance: ', 
                     toString(i), 
                     sept='')
              )+
      
      theme(plot.title = element_text(size = 9),
            axis.text.x = element_text(angle = 90, hjust = 1),
      ) + 
      labs(x = toString(paste0('x',
                               toString(1), 
                               sep='')
                        ), 
           y = toString(paste0('x',
                               toString(2), 
                               sep='')
                        )
           )
  }
  
  for(i in head(index, -1)){
    for(j in 1:(ncol(resAgg[[paste0('res', 
                                    toString(i),
                                    sep='')
                             ]][[3]])-1)){
      
      df <- data.frame(x = resAgg[[paste0('res',
                                          toString(i),
                                          sep='')
                                   ]][[3]][,j], 
                       y = resAgg[[paste0('res',
                                          toString(i),
                                          sep='')
                                   ]][[3]][,(j+1)]
                       )
      
      p[[i]][[j]] <- ggplot(df, 
                            aes(x, y)
                            )+
        
        scale_x_continuous(
          breaks= pretty_breaks(),
          expand = c(0.1, 0.1),
          labels = function(x) format(x, 
                                      scientific = TRUE)
        ) +
        
        scale_y_continuous(
          breaks= pretty_breaks(),
          expand =c(0.1, 0.1),
          labels = function(x) format(x, 
                                      scientific = TRUE)
        ) + 
        
        geom_segment(
          aes(xend=c(tail(x, n=-1), NA), 
              yend=c(tail(y, n=-1), NA)
              ),
          arrow=arrow(length=unit(0.05,"inches"),
                      type = "closed",
                      angle=7.5),
          size=0.2,
          alpha=0.5
        ) +
        
        geom_point(
          data = df[which.min(resAgg[[paste0('res', 
                                             toString(i),
                                             sep='')]][[4]][,1]),1:2],
          shape=21, color = "darkred", 
          fill ="red", 
          stroke = 0.5, 
          size = 1, 
          alpha=0.5
        ) +
        
        theme_light() +
        ggtitle(paste0('Spot run instance: ', toString(i), sept='')
                )+
        theme(plot.title = element_text(size = 9),
              axis.text.x = element_text(angle = 90, hjust = 1),
              text = element_text(size=6))+
        
        labs( x = toString(paste0('x',
                                  toString(j), 
                                  sep='')
                           ), 
              y = toString(paste0('x',
                                  toString(j+1), 
                                  sep='')
                           )
              )
      
    }
  }
  plot.res <- list()
  
  for(i in head(index, -1)){
    for(j in 1:(ncol(resAgg[[paste0('res',
                                    toString(i),
                                    sep='')]][[3]])-1)){
      
      if(i == 1 && j == 1) {plot.res <- p[[i]][[j]]}
      else {plot.res <- plot.res +  p[[i]][[j]]}
    }
  }
  
  plot.res <- plot.res+  
                plot_layout(ncol=2) +
                plot_annotation(tag_levels = 'A')
  
  assign('visPath2Dout', plot.res)                          
}
visDensMat <- function(resAgg){
  
  my_fn2 <- function(data, mapping, method="loess", ...){
    
    p <- ggplot(data, mapping = mapping) +  
      
      geom_point(
        shape=21,
        colour="black", 
        fill ="grey",
        stroke = 0.2,
        size = 0.5, 
        alpha = 0.75
      ) + 
      
      geom_smooth(
        method=loess, 
        size = 0.4, 
        alpha = 0.2,
        formula = y ~ x,
        colour="black", 
        fill ="grey",
        se = TRUE,
        ...) +  
      
      scale_x_continuous(
        breaks= pretty_breaks(),
        labels = function(x) format(x, scientific = TRUE)
      ) +
      
      scale_y_continuous(
        breaks= pretty_breaks(),
        labels = function(x) format(x, scientific = TRUE)
      ) + 
      
      theme_light()  
    
    p
    
  }
  my_dens_hist <- function(data, mapping, ...){
    
    p <- ggplot(data, mapping = mapping) +  
      
      geom_histogram(
        aes(y = ..density..), 
        bins = 8, 
        color="black", 
        fill="white"
      ) +
      
      scale_x_continuous(
        breaks= pretty_breaks(),
        labels = function(x) format(x, scientific = TRUE)
      ) +
      
      scale_y_continuous(
        breaks= pretty_breaks(),
        labels = function(x) format(x, scientific = TRUE)
      ) +
      
      scale_color_grey() + 
      
      theme_classic() 
    
    p
  }
  
  index <- 1:length(resAgg)
  pc <- list()
  
  for(i in head(index, -1)){
    df <- data.frame( x = resAgg[[paste0('res', 
                                         toString(i),
                                         sep='')
                                  ]][[3]],
            
                      y = resAgg[[paste0('res', 
                                         toString(i),
                                         sep='')
                                  ]][[4]])
    
    pc[[i]] <- ggpairs(
      df,
      lower = list(continuous= my_fn2),
      upper = list(continuous = "blank"),
      diag = list(continuous = my_dens_hist),
      axisLabels = c("show"),
      title = paste0('Spot run instance: ', toString(i))
    ) +
      
      theme(
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 5),
        plot.title = element_text(size = 9),
        legend.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_rect(fill = "white")
      ) 
  }
  
# multiplot from R cook-book 
# (Source: 
#http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
 
   multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)
                           ),
                       ncol = cols, 
                       nrow = ceiling(numPlots/cols)
                       )
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                 ncol(layout
                                                      )
                                                 )
                            )
                   )
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, 
                                        arr.ind = TRUE)
                                  )
        
        print(plots[[i]], 
              vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  plotout <- multiplot(plotlist = pc,
                       cols = 1)
  
  assign('visDensMatout', plotout) 
} 
visMeanSd2D  <- function(resAgg){
  
  index <- 1:length(resAgg)
  
  p <- list()
  pSd <- list()
  pMean <- list()
  
  for(i in head(index, -1)){
    
    fit <- list()
    predictionList <- list()
    
    fit <- buildKriging(resAgg[[paste0('res', 
                                       toString(i))
                                ]][['x']][,c(1,2)], 
                        resAgg[[paste0('res', 
                                       toString(i))
                                ]][['y']])
    
    fit$target <- c("y","s","ei")
    predictionList <- predict(fit, resAgg[[paste0('res', 
                                                  toString(i))
                                           ]][['x']][,c(1,2)]) 
    
    dfSd <- data.frame( x = resAgg[[paste0('res', 
                                           toString(i),
                                           sep='')
                                    ]][[3]][,1], 
                        
                        y = resAgg[[paste0('res', 
                                           toString(i),
                                           sep='')
                                    ]][[3]][,2], 
                        
                        z = predictionList[['s']])
    
    dfMean <- data.frame( x = resAgg[[paste0('res', 
                                             toString(i),
                                             sep='')
                                      ]][[3]][,1], 
                          
                          y = resAgg[[paste0('res', 
                                             toString(i),
                                             sep='')
                                      ]][[3]][,2], 
                          
                          z = predictionList[['y']])
    
    pSd[[i]] <- ggplot(dfSd, 
                       aes(x, y, z)
                       )+ 
      stat_summary_2d(aes(z = z)
                      )+
      scale_fill_gradientn(
        colours = terrain.colors(10), 
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2)
      )+  
      scale_x_continuous(
        breaks= pretty_breaks(), 
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2)
      )+
      scale_y_continuous(
        breaks= pretty_breaks(), 
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2)
      )+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90, 
                                       hjust = 1), 
            text = element_text(size=6),
            plot.title = element_text(size = 7),
            legend.title = element_text(size = 3), 
            legend.text  = element_text(size = 3),
            legend.key.size = unit(0.25, 
                                   "lines")
            )+
      guides(shape = guide_legend(override.aes = list(size = 0.75)
                                  )
             )+
      guides(color = guide_legend(override.aes = list(size = 0.75)
                                  )
             )+
      ggtitle(paste0('Spot run instance: ', toString(i), ' , sd')
              )+ 
      xlab(paste0('x', toString(1), sep='')
           )+ 
      ylab(paste0('x', toString(2), sep='')
           )
    
    pMean[[i]] <- ggplot(dfMean, 
                         aes(x, y, z)
                         )+ 
      stat_summary_2d(aes(z = asinh(z)))+
      scale_fill_gradientn(
        colours = terrain.colors(10), 
        labels = function(x) formatC(sinh(x), 
                                     format = "e", 
                                     digits = 2)
      )+  
      scale_x_continuous(
        breaks= pretty_breaks(), 
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2)
      )+
      scale_y_continuous(
        breaks= pretty_breaks(), 
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2)
      )+
      theme_light()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            text = element_text(size=6),
            plot.title = element_text(size = 7),
            legend.title = element_text(size = 3), 
            legend.text  = element_text(size = 3),
            legend.key.size = unit(0.25, "lines")
            )+
      guides(shape = guide_legend(override.aes = list(size = 0.75)
                                  )
             )+
      guides(color = guide_legend(override.aes = list(size = 0.75)
                                  )
             )+
      ggtitle(paste0('Spot run instance: ', toString(i), ' , mean')
              ) + 
      xlab(paste0('x', 
                  toString(1), 
                  sep='')
           ) + 
      ylab(paste0('x', 
                  toString(2), 
                  sep='')
           )
    
  }
  
  for(i in head(index, -1)){
    for(j in 1:(ncol(resAgg[[paste0('res', 
                                    toString(i),
                                    sep='')]][[3]])-1)){
      
      fit <- list()
      predictionList <- list()
      
      fit <- buildKriging(resAgg[[paste0('res', 
                                         toString(i))
                                  ]][['x']][,c(j,(j+1))], 
                          
                          resAgg[[paste0('res', 
                                         toString(i))
                                  ]][['y']])
      
      fit$target <- c("y","s","ei")
      predictionList <- predict(fit, resAgg[[paste0('res', 
                                                    toString(i))
                                             ]][['x']][,c(j,(j+1))]) 
      
      dfSd <- data.frame( x = resAgg[[paste0('res', 
                                             toString(i),
                                             sep='')
                                      ]][[3]][,j], 
                          y = resAgg[[paste0('res', 
                                             toString(i),
                                             sep='')
                                      ]][[3]][,(j+1)], 
                          z = predictionList[['s']])
      
      dfMean <- data.frame( x = resAgg[[paste0('res', 
                                               toString(i),
                                               sep='')
                                        ]][[3]][,j], 
                            y = resAgg[[paste0('res', 
                                               toString(i),
                                               sep='')
                                        ]][[3]][,(j+1)], 
                            z = predictionList[['y']])
      
      pSd[[i]][[j]] <- ggplot(dfSd, aes(x, y, z))+ 
        stat_summary_2d(aes(z = z))+
        scale_fill_gradientn(
          colours = terrain.colors(10), 
          labels = function(x) formatC(x, 
                                       format = "e", 
                                       digits = 2)
        )+  
        scale_x_continuous(
          breaks= pretty_breaks(), 
          labels = function(x) formatC(x, 
                                       format = "e", 
                                       digits = 2)
        )+
        scale_y_continuous(
          breaks= pretty_breaks(), 
          labels = function(x) formatC(x, 
                                       format = "e", 
                                       digits = 2)
        )+
        theme_light()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1), 
              text = element_text(size=6),
              plot.title = element_text(size = 7),
              legend.title = element_text(size = 3), 
              legend.text  = element_text(size = 3),
              legend.key.size = unit(0.25, "lines")
              )+
        guides(shape = guide_legend(override.aes = list(size = 0.75)
                                    )
               )+
        guides(color = guide_legend(override.aes = list(size = 0.75)
                                    )
               )+
        ggtitle(paste0('Spot run instance: ', 
                       toString(i), 
                       ' , sd')) + 
        xlab(paste0('x', 
                    toString(j), 
                    sep='')
             ) + 
        ylab(paste0('x', 
                    toString(j+1), 
                    sep='')
             )
      
      pMean[[i]][[j]] <- ggplot(dfMean, aes(x, y, z))+ 
        stat_summary_2d(aes(z = asinh(z)))+
        scale_fill_gradientn(
          colours = terrain.colors(10), 
          labels = function(x) formatC(sinh(x), 
                                       format = "e", 
                                       digits = 2)
        )+  
        scale_x_continuous(
          breaks= pretty_breaks(), 
          labels = function(x) formatC(x, 
                                       format = "e", 
                                       digits = 2)
        )+
        scale_y_continuous(
          breaks= pretty_breaks(), 
          labels = function(x) formatC(x, 
                                       format = "e", 
                                       digits = 2)
        )+
        theme_light()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1), 
              text = element_text(size=6),
              plot.title = element_text(size = 7),
              legend.title = element_text(size = 3), 
              legend.text  = element_text(size = 3),
              legend.key.size = unit(0.25, "lines"))+
        guides(shape = guide_legend(override.aes = list(size = 0.75)))+
        guides(color = guide_legend(override.aes = list(size = 0.75)))+
        ggtitle(paste0('Spot run instance: ', toString(i), ' , mean')) + 
        xlab(paste0('x', toString(j), sep='')) + 
        ylab(paste0('x', toString(j+1), sep=''))
    }
  }
  
  visMeanSdMap.res <- list()
  
  for(i in head(index, -1)){
    for(j in 1:(ncol(resAgg[[paste0('res', 
                                    toString(i),
                                    sep='')
                             ]][[3]])-1)){
      
      if(i == 1 && j == 1) {
        c <- 1
        p[[c]] <- (pMean[[i]][[j]] + pSd[[i]][[j]])}
      else {
        c <- c + 1
        p[[c]] <- (pMean[[i]][[j]] + pSd[[i]][[j]])}
    }
  }
  
  for(i in 1:length(p)){
    if(i == 1){visMeanSd2D.res <- p[[i]]}
    else{visMeanSd2D.res <- visMeanSd2D.res / p[[i]]}
    
  }
  
  if(length(p) == 1){visMeanSd2D.res <- visMeanSd2D.res+ 
                                        plot_layout(ncol=2)+ 
                                        plot_annotation(tag_levels = 'A')}
  
  else {visMeanSd2D.res <- visMeanSd2D.res+ 
                              plot_layout(ncol=1)+ 
                              plot_annotation(tag_levels = 'A')}
                
  visMeanSd2D.res
  
  assign('visMeanSd2Dout', visMeanSd2D.res)
}
visContour2D <- function(resAgg){
  
  createPlotModelList  <- function(resAgg){
    
    Model.list <- list()
    
    createPredictedMap <- function(resAgg, i, j){
      i <- i  
      j <- j
      
      f <- evaluateModel(buildKriging(resAgg[[paste0('res', 
                                                     toString(i))
                                              ]][['x']][,c(1,2)], 
                                      
                                      resAgg[[paste0('res', 
                                                     toString(i))
                                              ]][['y']], 
                                      control = list(algTheta = optimLHD)))
      
      x <- seq(min(resAgg[[paste0('res', 
                                  toString(i))
                           ]][['x']][,j]), 
               max(resAgg[[paste0('res', 
                                  toString(i))
                           ]][['x']][,j]), length = 25) 
      
      y <- seq(min(resAgg[[paste0('res', 
                                  toString(i))
                           ]][['x']][,j+1]), 
               max(resAgg[[paste0('res', 
                                  toString(i))
                           ]][['x']][,j+1]), length = 25)
      
      fn <- function(a,b){
        f(cbind(a,b))
      }	
      names(x) <- x                         
      names(y) <- y                        
      
      z.model <- outer(x, y, fn)
      
      z.melt<-melt(data = z.model, 
                   id.vars = "ID", 
                   measure.vars = c("x", "y", "z")) 
      names(z.melt) <- c("x", "y", "z")
      
      z.melt
    }
    
    index <- 1:length(resAgg)
    c <- 1
    for(i in head(index, -1)){
      for(j in 1:(ncol(resAgg[[paste0('res', 
                                      toString(i),
                                      sep='')
                               ]][[3]])-1)){
        
        Model.list[[c]] <- createPredictedMap(resAgg, i = i, j=j)
        
        c <- c+1
      }
    }
    
    Model.list
    
  }
  
  createPlotModelVar <- function(resAgg){
    
    dfVar.list <- list()
    
    createVarCatch <- function(resAgg, i, j){
      
      i <- i 
      j <- j
      
      df <- data.frame(x = resAgg[[paste0('res', 
                                          toString(i),
                                          sep='')
                                   ]][[3]][,j], 
                       
                       y = resAgg[[paste0('res', 
                                          toString(i),
                                          sep='')
                                   ]][[3]][,(j+1)])
      df
    }
    
    index <- 1:length(resAgg)
    c <- 1
    for(i in head(index, -1)){
      for(j in 1:(ncol(resAgg[[paste0('res', 
                                      toString(i),
                                      sep='')
                               ]][[3]])-1)){
        
        dfVar.list[[c]] <- createVarCatch(resAgg, i = i, j=j)
        
        c <- c+1
      }
    }
    
    dfVar.list
    
  }
  
  createPlotModelLabel <- function(resAgg){
    
    PlotModelLabel.list <- list()
    
    createLabelCatch <- function(resAgg, i, j){
      
      i <- i 
      j <- j
      
      Label.list <- list()
      Label.list <-c(paste0('Spot run instance: ', 
                            toString(i), 
                            sep=''), 
                     
                     paste0('x', 
                            toString(j), 
                            sep=''), 
                     
                     paste0('x', 
                            toString(j+1), 
                            sep='')
                     )
      
      Label.list
    }
    
    index <- 1:length(resAgg)
    c <- 1
    for(i in head(index, -1)){
      for(j in 1:(ncol(resAgg[[paste0('res', 
                                      toString(i),
                                      sep='')
                               ]][[3]])-1)){
        
        PlotModelLabel.list[[c]] <- createLabelCatch(resAgg, 
                                                     i = i, 
                                                     j=j
                                                     )
        c <- c+1
      }
    }
    
    PlotModelLabel.list
    
  }
  
  Model.list <- list()
  dfVar.list <- list()
  PlotModelLabel.list <- list()
  
  #Model.list <- createPlotModelList(resAgg)
  Model.list <- createPlotModelList(resAgg)
  
  dfVar.list <- createPlotModelVar(resAgg)
  PlotModelLabel.list <- createPlotModelLabel(resAgg)
  
  
  index <- 1:length(resAgg)
  
  c <- 1
  
  for(i in head(index, -1)){
    for(j in 1:(ncol(resAgg[[paste0('res', 
                                    toString(i),
                                    sep='')
                             ]][[3]])-1))
    {c <- c+1}
  }
  
  p <- list()
  
  for(i in 1:(c-1)){
    
    p[[i]] <- ggplot(Model.list[[i]], 
                     aes(x,y,z)
                     )+
      
      geom_tile(aes(fill= asinh(z)))+
      
      geom_contour(bins=5, 
                   aes(z= asinh(z)) , 
                   color="black", 
                   size=0.1)+
      
      geom_point(data = dfVar.list[[i]], 
                 aes(x,y),
                 shape = 21,
                 colour = "red", 
                 size = 0.5, stroke = 0.2, 
                 alpha=0.5)+ 
      
      scale_fill_gradientn(
        colours = terrain.colors(10),
        labels = function(x) formatC(sinh(x), 
                                     format = "e", 
                                     digits = 2),
      )+  
      
      scale_x_continuous(
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2),
        expand = c(0, 0)
        )+
      
      scale_y_continuous(
        labels = function(x) formatC(x, 
                                     format = "e", 
                                     digits = 2),
        expand = c(0,0)
      )+
      
      theme_light()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            text = element_text(size=6),
            plot.title = element_text(size = 7),
            legend.title = element_blank(), 
            legend.text  = element_text(size = 3),
            legend.key.size = unit(0.25, "lines"))+
      guides(shape = guide_legend(override.aes = list(size = 0.75)))+
      guides(color = guide_legend(override.aes = list(size = 0.75)))+
      ggtitle(PlotModelLabel.list[[i]][1]) + 
      xlab(paste0(PlotModelLabel.list[[i]][2])) + 
      ylab(paste0(PlotModelLabel.list[[i]][3]))
  }
  
  for(i in 1:length(p)){
    if(i == 1){visModel2D.res <- p[[i]]}
    else{visModel2D.res <- visModel2D.res / p[[i]]}
  }
  
  if(i == 1){visModel2D.res <- visModel2D.res+
    plot_layout(ncol=1)+ 
    plot_annotation(tag_levels = 'A')}
  
  else {visModel2D.res <- visModel2D.res + plot_layout(ncol=2)+
    plot_annotation(tag_levels = 'A')}
  
  assign('visContour2D', visModel2D.res)
}

visSpot <- function(resAgg){ 
  
  if (resAgg[1] == ''){
    
    index <- 1:length(resAgg)
    
    for(i in head(index, -1)){
      for(j in 1:(ncol(resAgg[[paste0('res', 
                                      toString(i),
                                      sep='')
      ]][[3]])-1)){
        if(i == 1 && j == 1) {c <- 1}
        else {c <- c + 1
        }
      }
    }
    
    if(c==1){
      p2.height <- 3
      p3.height <- 5
      p4.height <- 2.5
      p5.height <- 2.5
    }
    else if(c == 2){
      p2.height <- 3
      p3.height <- 6
      p4.height <- 3
      p5.height <- 3
    }
    else{      
      p2.height <- 6.36
      p3.height <- 6.36
      p4.height <- 6.36
      p5.height <- 6.36
    }
    
    tikz(file = "Overview.tex")
    
    print('')
    
    dev.off()

    tikz(file = "visProgEval.tex", 
         width = 4.5, 
         height = 2.25)
    
    print(visProgEval(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: Progress evaluation printed as visProgEval.tex")
    
    tikz(file = "visPath2D.tex", 
         width = 4.5, 
         height = p2.height)
    
    print(visPath2D(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: 2D pathmap printed as visPath2D.tex")
    
    tikz(file = "visDensMat.tex", 
         width = 4.5, 
         height = p3.height)
    
    print(visDensMat(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: Density Matrix printed as visDensMat.tex")
    
    tikz(file = "visMeanSd2D.tex", 
         width = 4.5, 
         height = p4.height)
    
    print(visMeanSd2D(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: Mean Standard deviation 2D map printed as visMeanSd2D")
    
    tikz(file = "visContour2D.tex", 
         width = 4.5, 
         height = p5.height)
    
    print(visContour2D(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: 2D Contour Model printed as visContour2D.tex ")
    
    
  }
  
  else{
    index <- 1:length(resAgg)
    
    for(i in head(index, -1)){
        for(j in 1:(ncol(resAgg[[paste0('res', 
                                        toString(i),
                                        sep='')
                                 ]][[3]])-1)){
          if(i == 1 && j == 1) {c <- 1}
          else {c <- c + 1
          }
        }
    }
    
    if(c==1){
      p2.height <- 3
      p3.height <- 5
      p4.height <- 2.5
      p5.height <- 2.5
      }
    else if(c == 2){
      p2.height <- 3
      p3.height <- 6
      p4.height <- 3
      p5.height <- 3
      }
    else{      
      p2.height <- 6.36
      p3.height <- 6.36
      p4.height <- 6.36
      p5.height <- 6.36
      }
    
    print(xtable(visInputTable(resAgg)), 
          scalebox='0.5' ,
          floating=FALSE,
          latex.environments=NULL,
          booktabs=TRUE, 
          file = "Overview.tex")
    print("Overview printed as Overview.tex")
   
    tikz(file = "visProgEval.tex", 
         width = 4.5, 
         height = 2.25)
    
    print(visProgEval(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: Progress evaluation printed as visProgEval.tex")
    
    tikz(file = "visPath2D.tex", 
         width = 4.5, 
         height = p2.height)
    
    print(visPath2D(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: 2D pathmap printed as visPath2D.tex")
    
    tikz(file = "visDensMat.tex", 
         width = 4.5, 
         height = p3.height)
    
    print(visDensMat(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: Density Matrix printed as visDensMat.tex")
    
    tikz(file = "visMeanSd2D.tex", 
         width = 4.5, 
         height = p4.height)
    
    print(visMeanSd2D(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: Mean Standard deviation 2D map printed as visMeanSd2D")
    
    tikz(file = "visContour2D.tex", 
         width = 4.5, 
         height = p5.height)
    
    print(visContour2D(resAgg), 
          bottomHeightProportion = 1, 
          leftWidthProportion = 1)
    
    dev.off()
    print("Visualization: 2D Model printed as visContour2D.tex ")
  }
}