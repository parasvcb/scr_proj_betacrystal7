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
dfgen=dfgen[!(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]

#dfgen<-dfgen[!(dfgen$color=="ra20" | dfgen$color=="ra250"),]

dfgen$aatype <- ifelse(dfgen$polymer == 'polyglycine' | dfgen$polymer == 'polyalanine_default' | dfgen$polymer == 'polyala-gly', 'small_group',
                  ifelse(dfgen$polymer == 'polythreonine' | dfgen$polymer == 'polyasparagine', 'polar',
                         ifelse(dfgen$polymer == 'polyisoleucine' | dfgen$polymer == 'polyvaline', 'hydrophobic', 'amendment')))


plotDF <- function(df,extractele,outapp) {
    #three different types of datasets iam looking over here,
    #df having color as only raw, so that line fitting can be done, 
    #df with only ra20, as one of the panle but with line plot, thin line
    #df with only dispav, so that one of set will have it,
    #ablines are needed and hence pending

    dfstart=df[(df$category=='start'),]
    dfend= df[(df$category=='end'),]
    df=df[!(df$category=='start') | (df$category=='end'),]

    outfile=file.path(outapp,extractele)
    dfnew=df[df$category==extractele,]
    dfraw=dfnew[(dfnew$color=='raw'),]
    dfra20=dfnew[(dfnew$color=='ra20'),]
    dfdispav=dfnew[(dfnew$color=='dispav'),]
    #dfnewgly <- dfnewgly[order(dfnewgly$distance),] 


    gg <- ggplot(data=dfra20[!is.na(dfra20$value),], aes(x=displacement,y=value,color=polymer))
    gg <- gg + geom_point(size=0.5)
    #gg <- gg + geom_point(data=dfra20, aes(x=displacement,y=value,color=polymer), size=0.7)
    #gg <- gg + geom_line(data=dfra20, aes(x=displacement,y=value,color=polymer), size=0.5)

    #gg <- gg + geom_vline(xintercept=dfstart$value,color=dfstart$polymer,linetype="soild",size=1)
    #gg <- gg + geom_vline(xintercept=dfend$value,color=dfend$polymer,linetype="dotted",size=1)
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Frame/1ps",breaks=seq(-2,25,4))
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"ra20_point.jpeg"),gg, height = 5, width= 15)

    gg <- ggplot()
    gg <- gg + geom_smooth(data=dfraw, aes(x=displacement,y=value,color=polymer), size=0.7, method = "auto", fullrange=TRUE)
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Frame/1ps",breaks=seq(-2,25,4))
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"fit_auto.jpeg"),gg, height = 5, width= 15)

    # gg <- ggplot()
    # #gg <- gg + geom_point(data=dfdispav, aes(x=displacement,y=value,color=polymer), size=0.7)
    # gg <- gg + geom_line(data=dfdispav, aes(x=displacement,y=value,color=polymer), size=0.5)
    gg <- ggplot(data=dfdispav[!is.na(dfdispav$value),], aes(x=displacement,y=value,color=polymer))
    gg <- gg + geom_point(size=0.5)
    gg <- gg + facet_wrap(~aatype,nrow=1) 
    gg <- gg + scale_x_continuous(name="Frame/1ps",breaks=seq(-2,25,4))
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,"dispav_point.jpeg"),gg, height = 5, width= 17)
}

uni=unique(dfgen$category)
#uni %in% 'frametrav'
uni = uni[!(uni %in% c('frametrav','start','end'))]
print (uni)

#uni = c("for")
#uni=c("force2")
#print (head(dfgen))
#print (uni)
for (i in uni) {
    print (i)
    plotDF(dfgen,i,out)
}
                                                                          

