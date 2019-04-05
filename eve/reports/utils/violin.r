round.signif <- function(x, p)
{
  ifelse(abs(x)>=1, round(x, p), signif(x, p))
}



add.level <- function( ds1, f){
    f <- f[f != '.']
    stopifnot(all(f %in% colnames(ds1)) && is.data.frame(ds1))

    ds2 <- ds1

    for (i in f){
        if(is.factor(ds1[, i])) ds2[, i] <- as.character(ds2[, i])
        else stopifnot(is.character(ds2[, i]))
    }

    fL1 <- lapply(f, function(i) unique(ds2[, i]) )

    fL2 <- expand.grid(fL1)
    colnames(fL2) <- f

    fL3 <- fL2[!apply(fL2, 1, paste, collapse=' ') %in% apply(ds2[, f], 1, paste, collapse=' '), ,drop=F]
    if(nrow(fL3) == 0) return (ds1)
    else{
        require(plyr)
        ds2 <- rbind.fill(ds2, fL3)
        for (i in f){
            if(is.factor(ds1[, i])) ds2[, i] <- factor(ds2[, i], levels=levels(ds1[, i]))
        }
        return(ds2)
    }

}

violinPlot <- function(
    dsin,    #A data frame object
    x,	#Name of the variable to be plotted on the x axis
    y,	#Name of the variable to be plotted on the y axis
    pch=1,	#Plotting symbol
    col=NA,	#Color for the box, automatically decided if col=NA, otherwise, it has to be a color vector for each level of x. It can be a named vector (e.g: c('MALE'='blue','FEMALE'='green')), where the names will decide factor levels of x. So a level that is not listed will be skipped in plotting.
    dot.col= NA, # color for dot. If NA or missing, no dots are shown 
    y.transf=NA,	#Transformation function for the y variable when plotting.  No transformation applied if =NA. Otherwise the possible values are log2 or log10. Do not use quotation mark
    y.scaling = 'identity',
    xlab='', #Label for the x axis
    ylab='',    #Label for the y axis
    ylim,	#The y limits for the plot. if the y value need be transformed, please provide the ylim in the original scale
    main='',	#A main title for the plot
    sub='',	#A subtitle for the plot
    h=NA,	#The y values for additional horizontal lines. if the y value need be transformed, please provide h in the original scale
    
    h.lty=1,	#The line type for additional horizontal lines
    h.col='black',	#The color for the additional horizontal lines
    h.lwd=1,    #Line width for additional horizontal lines. Can be a vector.
    
    h.facets='.',	#Conditioning variable for horizontal facets. If specified, records with NA h.facets will be ignored
    v.facets='.',	#Conditioning variable for vertical facets. If specified, records with NA v.facets will be ignored
    add.n = T,	#whether to show n in each group
    pMethod='',	#'nonparametric' (Kruskal-Wallis or wilcoxon test),'anova' (for permutation based anova test),''(no p-value shown). Analysis was performed on the origina data before transformation
    pTestSubset=NA,    #run the p value method in a subset of observations
    yRangeAdjust = list(n=0.05, p=0.05, default=0.05), # n: proportion of y range added to show number of observations per group, p: portion of y range added to show p-value, default: added to y by default
    main.size=1,    #Magnification of the main title (main) relative to base font size. Default = 1.
    xlab.size=1,    #Magnification of the x axis label (xlab) relative to base font size. Default = 1.
    ylab.size=1,	#Magnification of the y axis label (ylab) relative to base font size. Default = 1.
    x.jitter=0.25,  #random noise applied to the x axis
    x.anns.size=1,    #Magnification of x axis annotations (x.anns) relative to base font size. Default = 1
    x.anns.angle=0,	#Angle of x axis annotations (x.anns) (degrees to be rotated counter-clockwise). Default = 0 (i.e., parallel to the x axis). If = 90, then annotations are perpendicular to the x axis.
    x.anns.align='center',	#Alignment of x axis annotations (x.anns). Possible values are "left", "center" (default), and "right".
    y.anns.size=1,   #Magnification of y axis annotations (y.anns) relative to base font size. Default = 1
    y.anns.angle=0,	#Angle of y axis annotations (y.anns) (degrees to be rotated counter-clockwise). Default = 0 (i.e., perpendicular to the y axis). If = 90, then annotations are parallel to the y axis.
    y.anns.align='right',	#Alignment of y axis annotations (y.anns). Possible values are "left", "center" (default), and "right"
    
    force.missing.category=F,    # whether to show all x axis levels in subgroups where certain levels are missing
    scaleMethod='count',
    medianColor='white', #color for the plus signs that indicate median,
    medianSymbolSize = 10 #reduce this number will reduce the size of the median symbol
){


    # check input variable
    stopifnot(pMethod %in% c('', 'anova','nonparametric'))
    if(is.na(xlab) || xlab=='') xlab=x
    if(is.na(ylab) || ylab=='') ylab=y

    if(!missing(y.transf)) stopifnot(is.function(y.transf))

    # check input data and columns
    stopifnot(is.data.frame(dsin) && all(c(x, y) %in% colnames(dsin)))

    if(!missing(h.facets) && !is.na(h.facets) && h.facets != '.') {
        stopifnot(h.facets %in% colnames(dsin))
        dsin <- subset(dsin, !is.na(dsin[, h.facets]))
    }

    if(!missing(v.facets) && !is.na(v.facets) && v.facets != '.') {
        stopifnot(v.facets %in% colnames(dsin))
        dsin <- subset(dsin, !is.na(dsin[, v.facets]))
    }

    if(!is.numeric(dsin[, y])) dsin[, y] <- as.numeric(dsin[, y])

    if(!is.factor(dsin[, x])) dsin[, x] <- factor(as.character(dsin[, x]))

    if(!missing(col) && !is.na(col) && nchar(col)> 0){
        stopifnot(length(col) == nlevels(dsin[, x]))
        if(!is.null(names(col))) dsin[, x] <- factor(as.character(dsin[, x]), levels=names(col))
    }else{
        col <- rainbow(nlevels(dsin[, x]))
    }

    ds1 <- droplevels(dsin[!(is.na(dsin[, x]) | is.na(dsin[, y]) ), ,drop=F ] )

    if(nrow(ds1) == 0){
        return('no non-missing data')
    }

    if( force.missing.category && (v.facets != '.' || h.facets != '.')) ds1 <- add.level( ds1, c(v.facets, h.facets, x))

    pv <- NA
    if(pMethod != '' && v.facets == '.' && h.facets == '.') {
        if(is.na(pTestSubset)) {
            selTest <- rep(T, nrow(ds1))
        }else{
            stopifnot(is.vector(pTestSubset) && all(names(pTestSubset) %in% c(x, y)) )
            fl <- list()
            for (i in names(pTestSubset))  fl[[i]] <- eval(parse(text=paste('ds1[, "', i,'"] %in%', pTestSubset[i], sep="")))
            fm <- do.call(cbind, fl)
            selTest <- apply(fm, 1, all)
        }

        if(pMethod == 'nonparametric') pv <- ifelse(nlevels(ds1[selTest, x]) > 1,  kruskal.test(ds1[selTest, y] ~ds1[selTest, x])$p.value, NA)

        if(pMethod == 'anova'){
            require(coin)
            pv <- oneway_test(ds1[selTest, y] ~ds1[selTest, x], distribution=approximate(5000) )
            if (grepl('IndependenceTest', class(pv))) {
                pv <- pvalue(pv)
            }else{
                pv <- NA
            }
        }
    }

     if(is.function(y.transf)){
        y.transf.txt <- deparse(substitute(y.transf))
        if(!grepl(y.transf.txt, ylab, ignore.case = T)) ylab = paste(y.transf.txt, ylab)
        ds1[, y] <- sapply(ds1[, y], y.transf)
        if(any(!is.na(h)) ) h <- sapply(h, y.transf)
    }

  if( missing(ylim) || is.na(ylim) || length(ylim) != 2){
        ylim <- range(ds1[, y],na.rm=T)
        ylim[1] <- ylim[1] - diff(ylim) * yRangeAdjust$default
        ylim[2] <- ylim[2] + diff(ylim) * yRangeAdjust$default
    }else{
        if(is.function(y.transf)) ylim <- sapply(ylim, y.transf)
    }

  #        ylim[1] <- ylim[1] - diff(ylim) * 0.05
    if(!is.na(pv))     ylim[2] <- ylim[2] + diff(ylim) * yRangeAdjust$p # extra space for p-value

   if(add.n) {
        ylim[1] <- ylim[1] - diff(ylim) * yRangeAdjust$n
        xlevels= levels(ds1[, x])

        n_fun <- function(x){

            datan =data.frame(label = paste(length(x)), y=ylim[1], stringsAsFactors=F )
            datan[1,'label'] <- paste('n =',datan[1,'label'])
            return(datan)
        }
    }


    require(ggplot2)
    require(grid)
    set.seed(27519)

    col <- unname(col)

    p1 <- ggplot(ds1, aes_string(x, y))+ geom_violin(aes_string(fill=x), scale = scaleMethod) + scale_fill_manual(values=col) + theme_bw() + labs( x=xlab, y=ylab ) + ggtitle(bquote(atop(.(main), atop(italic(.(sub)), "")))) + theme(legend.position = "none") + scale_y_continuous(limits=ylim, trans=y.scaling) 
    
    if(!is.na(dot.col)) p1 <- p1 + geom_jitter(position = position_jitter(height = 0, width=x.jitter), pch=pch, col=dot.col)

    if(add.n) p1 <- p1 + stat_summary(fun.data = n_fun, geom = "text"  )
  
    #    if(srt != 0) p1 <- p1 + theme(axis.text.x = element_text(angle = srt, hjust = 1))
    p1 <- p1 + stat_summary(fun.y="median", geom="point", pch='+', size=rel(xlab.size* medianSymbolSize), col=medianColor) +  theme(axis.title.x=element_text(size=rel(xlab.size)), axis.text.x=element_text(angle=x.anns.angle, hjust=switch(x.anns.align,left=0,center=0.5,right=1), vjust=0.5, size=rel(x.anns.size)),
                     axis.title.y=element_text(size=rel(ylab.size)), axis.text.y=element_text(angle=y.anns.angle, hjust=switch(y.anns.align,left=0,center=0.5,right=1), vjust=0.5, size=rel(y.anns.size)), plot.title = element_text(size = rel(main.size)))

    if(!is.na(pv)) p1 <- p1 + annotate('text', x=1, y=ifelse(ylim[2] > 0, ylim[2] * 0.95, ylim[2] * 1.05), label=paste('p=', round.signif(pv, 2)), size=rel(x.anns.size*5))

    if(any(!is.na(h))){
        if( any (h < ylim[1] | h > ylim[2]) ) warning(paste('some or all horizontal reference lines that are out of the plotting region will be skipped'))
        p1 <- p1 + geom_hline(yintercept = h, linetype=h.lty, color=h.col, size=h.lwd )
    }

    if(v.facets != '.' || h.facets != '.') p1 <- p1 + facet_grid( as.formula(paste( h.facets, v.facets ,sep='~')), scales = ifelse(add.n,'free', "fixed"))

    gt <- ggplot_gtable(ggplot_build(p1))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    grid.draw(gt)

    invisible(p1)

}

if (F){
dsin <- read.csv('c:/temp/VignetteData.csv', stringsAsFactors=F)
x <- 'BESRSPI'
y <- 'BAGE'
}

# example
#yyboxplot(dsin=dsin, x = 'BESRSPI', y='BAGE', main='baseline', sub='resp', y.transf=log2, h=c(30, 40), pMethod='nonparametric', h.facet ='TRT', x.anns.size=3,  x.anns.angle=45)

#yyboxplot(dsin=dsin, x = 'BESRSPI', y='BAGE', main='baseline', sub='resp',  h=c(30, 40), pMethod='nonparametric', h.facet ='TRT', xlab.size=2,  x.anns.angle=10)

#yyboxplot(dsin=dsin, x = 'BESRSPI', y='BAGE', main='baseline', sub='resp',  h=c(30, 40), pMethod='nonparametric', h.facet ='TRT', xlab.size=2,  x.anns.angle=45, ylim=c(30, 90))

#yyboxplot(dsin=dsin, x = 'BESRSPI', y='RNA', main='baseline', sub='resp',   pMethod='nonparametric', h.facet ='TRT', v.facet='RACE', xlab.size=2,  x.anns.angle=45, y.transf=log2,force.missing.category=F)

