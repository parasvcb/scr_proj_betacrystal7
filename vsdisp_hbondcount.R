'dataframe_vnew2'
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
    dffile="processed/dataframe_vnew2.tsv"
    dfadd=file.path(dfpre,polymer,dffile)
    dataforce=read.csv(dfadd,sep = "\t",check.names = FALSE)
    colnames(dataforce)[which(names(dataforce) == "frames-raw")] <- "frames"
    dataforce <- subset(dataforce, select = -c(frames))
    df=melt(data = dataforce, id.vars = "displacement")
    newColNames <- c("category", "avtype")
    newCols <- colsplit(df$variable, "-", newColNames)
    afterpre <- cbind(df, newCols)
    afterpre$polymer <- polymer
    return (afterpre)
}
#pullinginterctaion
resname=paste0(out,"resdf.csv")
if (file.exists(resname)){
    print('reading')
  dfgen=read.csv(resname,check.names = FALSE)
} else {
dfgen=cropDF(systemFile,"polyglycine")
dfgen = dfgen[FALSE,]
for (i in directories) {
    print (i)
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
print ('read')

plotallad <- function(df,extractele,outapp) {
    print ('0')
    #print (unique(df$variable))
    df$keeps <- ifelse (grepl(extractele,df$variable),'keep',
    ifelse (grepl("for-dispav",df$variable),'keep','discard'))
    print ('1')
    #print (unique(df$keeps))
    #print(head(df))
    df$hbcolor <- ifelse(grepl('mcmc',df$variable), 'mcmc',
                   ifelse(grepl('scsc',df$category), 'scsc',
                ifelse(grepl('mcsc',df$category), 'mcsc',   'force')))
    dfset=df[(df$keeps=='keep'),]
    dfset=dfset[!(dfset$variable=='for-dispav'),]
    #the above can be removed to add thin layout to force in bacground
    
    print ("here")
    #dfset=dfset[(dfset$avtype=='ra20'),]
    dfset=dfset[(dfset$avtype=='dispav'),]
    print ('2')
    #print (head(dfset))
    #print (unique(dfset$variable))
    dfset$hbnature<- ifelse(grepl('p1',dfset$variable), 'partial',
                   ifelse(grepl('p2',dfset$variable), 'partial',   'unbiased'))
    print ('3')
    dfset2=dfset[(dfset$hbnature=='partial'),]
    dfset1=dfset[(dfset$hbnature=='unbiased'),]
    
    # print (head(dfset1))
    # print (unique(dfset1$value))
    # write.csv(dfset1,'check.tsv',sep = "\t")

    print (unique(dfset1$color))
    print ('4')
    #gg <- ggplot(data=dfraw[!is.na(dfraw$value),], aes(x=displacement,y=value,color=homopolymer))
    gg <- ggplot(data=dfset1[!is.na(dfset1$value),], aes(x=displacement,y=value,color=hbcolor))
    #gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("blue","red","green"))

    gg <- gg + facet_wrap(~polymer,ncol=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 7),axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=5)))
    gg <- gg + ylab("hbcount")
    ggsave(filename = paste0(outapp,'unbiased.pdf'),gg, height = 8, width= 3)

    d1='dashed'
    d2='solid'
    c1='blue'
    c11='lightblue'
    c2='black'
    c21='darkgray'
    c3='red'
    c31='maroon'
    c4='orange'
    c41='yellow'
    c5='green'
    c51='lightgreen'
    c6='purple'
    c61='cyan'
    print (unique(dfset2$category))
    gg <- ggplot(data=dfset2[!is.na(dfset2$value),], aes(x=displacement,y=value,color=category))
    #gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("brown","black","darkgreen","red","blue","darkorange","azure4"))
    gg <- gg + geom_line(size=0.4,) 
    #gg <- gg + scale_linetype_manual(values=c(d1,d2,d1,d2,d1,d2,d1,d2,d1,d2,d1,d2))
    gg <- gg + scale_color_manual(values=c(c1,c11,c2,c21,c3,c31,c1,c11,c2,c21,c3,c31))
    gg <- gg + facet_wrap(~polymer,ncol=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 7),axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=5)))
    gg <- gg + ylab("hbcount")
    ggsave(filename = paste0(outapp,'biased.pdf'),gg, height = 8, width= 3)

    }
    print ('5')
plotcombined <- function(df,extractele,outapp) {
    print ('0')
    #print (unique(df$variable))
    df$keeps <- ifelse (grepl(extractele,df$variable),'keep',
    ifelse (grepl("for-dispav",df$variable),'keep','discard'))
    df$hbcolor <- ifelse(grepl('adjmcmc',df$variable), 'adjmcmc',
                   ifelse(grepl('adjscsc',df$category), 'adjscsc',
                ifelse(grepl('adjmcsc',df$category), 'adjmcsc',
                ifelse(grepl('allmcmc',df$variable), 'allmcmc',
                   ifelse(grepl('allscsc',df$category), 'allscsc',
                ifelse(grepl('allmcsc',df$category), 'allmcsc', 'unknown'
                ))))))
    dfset=df[(df$keeps=='keep'),]
    dfset=dfset[!(dfset$variable=='for-dispav'),]
    #the above can be removed to add thin layout to force in bacground

    dfset=dfset[(dfset$avtype=='dispav'),]
    print ('2')
    dfset$hbnature<- ifelse(grepl('p1',dfset$variable), 'partial',
                   ifelse(grepl('p2',dfset$variable), 'partial',   'unbiased'))
    print ('3')
    dfset2=dfset[(dfset$hbnature=='partial'),]
    dfset1=dfset[(dfset$hbnature=='unbiased'),]
    print ('4')
    gg <- ggplot(data=dfset1[!is.na(dfset1$value),], aes(x=displacement,y=value,color=hbcolor))
    gg <- gg + geom_line(size=0.6,alpha=0.5 ) + scale_color_manual(values=c("black","darkgreen","blue","gray","green","lightblue"))
    gg <- gg + facet_wrap(~polymer,ncol=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 7),axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=5)))
    gg <- gg + ylab("hbcount")
    ggsave(filename = paste0(outapp,'unbiased.pdf'),gg, height = 8, width= 3)
    }
    print ('5')
#dfgen<-dfgen[!(dfgen$color=="ra20" | dfgen$color=="ra250"),]
dfgen=dfgen[!(dfgen$polymer=='polyala_constrained') & !(dfgen$polymer=='polygly_constrained') & !(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]
dfgen=dfgen[!(dfgen$category=='start') & !(dfgen$category=='end'),]
print ('6')
plotcombined(dfgen,'hb',paste0(out,'hbcombined_facet'))
plotallad(dfgen,'hbadj',paste0(out,'hbadjacent_facet'))
plotallad(dfgen,'hball',paste0(out,'hball_facet'))