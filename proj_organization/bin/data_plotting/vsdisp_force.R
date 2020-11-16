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

reptype <- ifelse(grepl('menton_set',systemFile), 'menton_set',
                ifelse(grepl('tyroneP',systemFile), 'tyroneP_set','tyrone_set3'))

directories=list.dirs(path = systemFile, full.names = TRUE, recursive=F)

cropDF <- function(dfpre,polymer) {
    dffile="processed/dataframe_vnew3.tsv"
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
resname=paste0(out,reptype,"resdf.csv")

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
dfgen$keeps <- ifelse (grepl('for',dfgen$variable),'keep','discard')
dfgen=dfgen[(dfgen$keeps=='keep'),]
dfgen$aatype <- ifelse(dfgen$polymer == 'polyglycine' | dfgen$polymer == 'polyalanine_default' | dfgen$polymer == 'polyala-gly', '1_SAA',
                  ifelse(dfgen$polymer == 'polythreonine' | dfgen$polymer == 'polyasparagine', '2_PAA',
                         ifelse(dfgen$polymer == 'polyisoleucine' | dfgen$polymer == 'polyvaline', '3_HAA', 'amendment')))

write.csv(dfgen,resname,sep = "\t",row.names=FALSE)
}

print ('read')

plotallad <- function(df,extractele,outapp) {
    dfset1=df[(df$avtype=='dispav') & (df$displacement<=18) ,]
   
    dfset2=df[(df$avtype=='raw'),]
    #print(head(df))
    gg <- ggplot()
    gg <- gg + geom_line(data=dfset2[!is.na(dfset2$value),], aes(x=displacement,y=value/1000),color='gray', size=0.6)
    gg <- gg + geom_line(data=dfset1[!is.na(dfset1$value),], aes(x=displacement,y=value/1000),color='black', size=0.4)
    gg <- gg + facet_wrap(~polymer,ncol=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 6),axis.text.x = element_text( hjust = 1, size = 6, angle = 45),axis.text.y = element_text( hjust = 1, size = 6), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab("Force (nN)")
    ggsave(filename = paste0(outapp,reptype,'Fig2_FvsD_vertical.pdf'),gg,height = 8, width= 3)

    datadf=dfset1[!is.na(dfset1$value),]
    #print (head(datadf))

    #gg1 <- ggplot(datadf, x=displacement,y=value/1000,color=polymer) 
    #above doesnt work
    gg1 <- ggplot(data=datadf, aes(x=displacement,y=value/1000,color=polymer)) 
    gg1 <- gg1 + geom_line(size=0.4)
    gg1 <- gg1 + facet_wrap(~aatype,nrow=1) 
    gg1 <- gg1 + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg1 <- gg1 + theme (strip.text.x = element_text(size = 6),axis.text.x = element_text( hjust = 1, size = 6, angle = 45),axis.text.y = element_text( hjust = 1, size = 6), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg1 <- gg1 + guides(colour = guide_legend(override.aes = list(size=10)))
    gg1 <- gg1 + ylab("Force (nN)")
    ggsave(filename = paste0(outapp,reptype,'Fig2_FvsD_dispav_horizontal.pdf'),gg1, height = 5, width= 14)


}
dfgen=dfgen[!(dfgen$polymer=='polyala_constrained') & !(dfgen$polymer=='polygly_constrained') & !(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]
dfgen=dfgen[!(dfgen$category=='start') & !(dfgen$category=='end'),]
plotallad(dfgen,'for',paste0(out,'force7'))
