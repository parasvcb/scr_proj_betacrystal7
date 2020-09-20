args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements systemDir and output_append_withoutext", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}
library(ggplot2)
library(reshape)

 
directories=list.dirs(path = systemFile, full.names = TRUE, recursive=F)

cropDF <- function(dfpre,polymer) {
    dffile="processed/dataframe_vnew.tsv"
    dfadd=file.path(dfpre,polymer,dffile)
    dataforce=read.csv(dfadd,sep = "\t",check.names = FALSE)
    colnames(dataforce)[which(names(dataforce) == "frames-raw")] <- "frames"
    dataforce <- subset(dataforce, select = -c(frames))
    df=melt(data = dataforce, id.vars = "displacement")
    newColNames <- c("category", "color")
    newCols <- colsplit(df$variable, "-", newColNames)
    afterpre <- cbind(df, newCols)
    afterpre$polymer <- polymer
    return (afterpre)
}
#pullinginterctaion
resname=paste0(out,"resdf.csv")
if (file.exists(resname)){
  dfgen=read.csv(resname,check.names = FALSE)
} else {
dfgen=cropDF(systemFile,"polyglycine")
dfgen = dfgen[FALSE,]
for (i in directories) {
    #print (i)
    #print (systemFile)
    if (grepl("poly[a-z].*",i)) {
        polymer=basename(i)
        #print (polymer)
        retdf=cropDF(systemFile,polymer)
        dfgen=rbind(dfgen,retdf)
        #this function will get the dataframes appended together
 }
}
write.csv(dfgen,resname,sep = "\t")
}

#dfgen<-dfgen[!(dfgen$color=="ra20" | dfgen$color=="ra250"),]
dfgen=dfgen[!(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]

# dfgen$aatype <- ifelse(dfgen$polymer == 'polyglycine' | dfgen$polymer == 'polyalanine_default' | dfgen$polymer == 'polyala-gly', 'small_group',
#                   ifelse(dfgen$polymer == 'polythreonine' | dfgen$polymer == 'polyasparagine', 'polar',
#                          ifelse(dfgen$polymer == 'polyisoleucine' | dfgen$polymer == 'polyvaline', 'hydrophobic', 'amendment')))

dfgen$aatype <- ifelse(dfgen$polymer == 'polyglycine' | dfgen$polymer == 'polyalanine_default' | dfgen$polymer == 'polyala-gly', 'SAA',
                  ifelse(dfgen$polymer == 'polythreonine' | dfgen$polymer == 'polyasparagine', 'PAA',
                         ifelse(dfgen$polymer == 'polyisoleucine' | dfgen$polymer == 'polyvaline', 'HAA', 'amendment')))

dfgen$homopolymer <- ifelse(dfgen$polymer == 'polyglycine', 'p[G]',
                  ifelse(dfgen$polymer == 'polyalanine_default', 'p[A]',
                  ifelse(dfgen$polymer == 'polyala-gly', 'p[AG]',
                  ifelse(dfgen$polymer == 'polyasparagine', 'p[N]',
                  ifelse(dfgen$polymer == 'polythreonine', 'p[T]',
                  ifelse(dfgen$polymer == 'polyvaline', 'p[V]',
                  ifelse(dfgen$polymer == 'polyisoleucine', 'p[I]','unknown')))))))


#dfstart=df[(df$category=='start'),]
#dfend= df[(df$category=='end'),]
dfgen=dfgen[!(dfgen$category=='start') & !(dfgen$category=='end'),]
dfgen=dfgen[!(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyala_constrained') & !(dfgen$polymer=='polygly_constrained') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]


plotDF <- function(df,extractele,outapp) {
    #three different types of datasets iam looking over here,
    #df having color as only raw, so that line fitting can be done, 
    #df with only ra20, as one of the panle but with line plot, thin line
    #df with only dispav, so that one of set will have it,
    #ablines are needed and hence pending

    outfile=file.path(outapp,extractele)
    dfnew=df[df$category==extractele,]
    dfraw=dfnew[(dfnew$color=='raw'),]
    dfra20=dfnew[(dfnew$color=='ra20'),]
    dfdispav=dfnew[(dfnew$color=='dispav'),]
    
    y_label = switch(extractele, "for" = "Force (pN)", "vel" = "Velocity (Ang/ps)")
    y_label = ifelse(identical(NULL,y_label) & substr(extractele, 1, 3)=="pse","Angle (degree)",
            ifelse(identical(NULL,y_label) & substr(extractele, 1, 2)=="hb","H-Bond count",
            ifelse(identical(NULL,y_label) & substr(extractele, 1, 3)=="com","Distance (Ang)",y_label)))
      
    
    gg <- ggplot(data=dfra20[!is.na(dfra20$value),], aes(x=displacement,y=value,color=homopolymer))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab(y_label)
    ggsave(filename = paste0(outfile,"ra20_line.pdf"),gg, height = 5, width= 15)

    gg <- ggplot()
    gg <- gg + geom_smooth(data=dfraw, aes(x=displacement,y=value,color=homopolymer), size=0.7, method = "auto", fullrange=TRUE) 
    gg <- gg + facet_wrap(~aatype,nrow=1)  + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab(y_label)
    ggsave(filename = paste0(outfile,"fit_auto.pdf"),gg, height = 5, width= 15)

    # gg <- ggplot()
    # #gg <- gg + geom_point(data=dfdispav, aes(x=displacement,y=value,color=homopolymer), size=0.7)
    # gg <- gg + geom_line(data=dfdispav, aes(x=displacement,y=value,color=homopolymer), size=0.5)
    gg <- ggplot(data=dfdispav[!is.na(dfdispav$value),], aes(x=displacement,y=value,color=homopolymer)) 
    gg <- gg + geom_line(size=0.6) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab(y_label)

    ggsave(filename = paste0(outfile,"dispav_line.pdf"),gg, height = 5, width= 15)
}

plotframes <- function(df,extractele,outapp) {
    outfile=file.path(outapp,extractele)
    dfnew=df[df$category==extractele,]
    dfraw=dfnew[(dfnew$color=='raw'),]
 
    gg <- ggplot(data=dfraw[!is.na(dfraw$value),], aes(x=displacement,y=value*0.5,color=homopolymer))
    gg <- gg + geom_line(size=0.6) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + ylab("Time (ps)")
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"raw_line.pdf"),gg, height = 5, width= 15)

    gg <- ggplot()
    gg <- gg + geom_smooth(data=dfraw[!is.na(dfraw$value),], aes(x=displacement,y=value*0.5,color=homopolymer), size=0.7, method = "auto", fullrange=TRUE) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + ylab("Time (ps)")
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"fit_auto.pdf"),gg, height = 5, width= 15)
    
    gg <- ggplot(data=dfraw[!is.na(dfraw$value),], aes(x=displacement,y=value*0.5,color=homopolymer))
    gg <- gg + geom_point(size=0.6) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab("Time (ps)")

    ggsave(filename = paste0(outfile,"raw_point.pdf"),gg, height = 5, width= 15)
}



plotforce <- function(df,extractele,outapp) {
    outfile=file.path(outapp,extractele)
    dfnew=df[df$category==extractele,]
    dfraw=dfnew[(dfnew$color=='raw'),]

    gg <- ggplot(data=dfraw[!is.na(dfraw$value),], aes(x=displacement,y=value,color=homopolymer))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~homopolymer,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab("Force (pN)")
    ggsave(filename = paste0(outfile,"force_raw_line_facet.pdf"),gg, height = 5, width= 25)

    gg <- ggplot(data=dfraw[!is.na(dfraw$value),], aes(x=displacement,y=value,color=homopolymer))
    gg <- gg + geom_point(size=0.6, ) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~homopolymer,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab("Force (pN)")
    ggsave(filename = paste0(outfile,"force_raw_point_facet.pdf"),gg, height = 5, width= 25)

    uni=unique(dfraw$homopolymer)
    for (i in uni) {
    col <- ifelse(i == 'p[G]','darkgreen',
                  ifelse(i ==  'p[A]','brown',
                  ifelse(i ==  'p[AG]','black',
                  ifelse(i == 'p[N]','blue',
                  ifelse(i == 'p[T]','darkorange',
                  ifelse(i == 'p[V]','azure4',
                  ifelse(i == 'p[I]','red','unknown')))))))
    dfcustom=dfraw[(dfraw$homopolymer==i),]
    
    gg <- ggplot(data=dfcustom[!is.na(dfcustom$value),], aes(x=displacement,y=value,color=col))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c(col))
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + scale_y_continuous(name="Force (pN)",limits=c(-100,3500))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    ggsave(filename = paste0(outfile,"force_line_",i,".pdf"),gg,height=5,width=8)

    gg <- ggplot(data=dfcustom[!is.na(dfcustom$value),], aes(x=displacement,y=value,color=col))
    gg <- gg + geom_point(size=0.6, ) + scale_color_manual(values=c(col))
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + scale_y_continuous(name="Force (pN)",limits=c(-100,3500))
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    ggsave(filename = paste0(outfile,"force_point_",i,".pdf"),gg,height=5,width=8)

    }
    }

#for       vel       frametrav pse6      pse4      com1      com2     
# com3      hball     hbadj     hballbbbb hballscsc hballscbb hbadjbbbb
# hbadjscsc hbadjscbb   

plot_com <- function(df,extractele,outapp) {

    outfile=file.path(outapp,extractele)
    dfnew=df[df$category==extractele,]
    dfraw=dfnew[(dfnew$color=='raw'),]
    dfra20=dfnew[(dfnew$color=='ra20'),]
    dfdispav=dfnew[(dfnew$color=='dispav'),]
    
    y_label ="Distance (Ang)"
      
    
    gg <- ggplot(data=dfra20[!is.na(dfra20$value),], aes(x=displacement,y=value,color=homopolymer))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4)) + ylim(-1.5,2.5)
    #gg <- gg + scale_y_continuous(name=y_label,breaks=seq(-1.5,2.5,0.5))
    gg <- gg + ylab(y_label)
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"ra20_line.pdf"),gg, height = 5, width= 15)

    gg <- ggplot()
    gg <- gg + geom_smooth(data=dfraw, aes(x=displacement,y=value,color=homopolymer), size=0.7, method = "auto", fullrange=TRUE) 
    gg <- gg + facet_wrap(~aatype,nrow=1)  + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4)) + ylim(-1.5,2.5)
    #gg <- gg + scale_y_continuous(name=y_label,breaks=seq(-1.5,2.5,0.5))
    gg <- gg + ylab(y_label)
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"fit_auto.pdf"),gg, height = 5, width= 15)

    gg <- ggplot(data=dfdispav[!is.na(dfdispav$value),], aes(x=displacement,y=value,color=homopolymer)) 
    gg <- gg + geom_line(size=0.6) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4)) + ylim(-1.5,2.5)
    #gg <- gg + scale_y_continuous(name=y_label,breaks=seq(-1.5,2.5,0.5))
    gg <- gg + ylab(y_label)
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    
    ggsave(filename = paste0(outfile,"dispav_line.pdf"),gg, height = 5, width= 15)
}

    

uni=unique(dfgen$category)
#uni %in% 'frametrav'
uni = uni[!(uni %in% c('start','end',"com1","com2","com3"))]
print (uni)
#uni = c("frametrav")
#uni=c("for")
#print (head(dfgen))
#print (uni)

# plotforce(dfgen,'for',out)
# for (i in uni) {
#     print (i)
#     if (i=='frametrav') {
#       plotframes(dfgen,i,out)  
#     } else {
#     plotDF(dfgen,i,out)
# }
# }

plot_com(dfgen,'com1',out)
plot_com(dfgen,'com2',out)
plot_com(dfgen,'com3',out)
                                                                          

