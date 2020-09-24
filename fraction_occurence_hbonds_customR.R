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

cropDF <- function(dfpre,polymer,dffile) {
    dfadd=file.path(dfpre,polymer,dffile)
    dataforce=read.csv(dfadd,sep = "\t",check.names = FALSE)
    #colnames(dataforce)[which(names(dataforce) == "frames-raw")] <- "frames"
    dataforce$polymer <- polymer
    #print (head(dataforce))
    return (dataforce)
}

dfgen_adjacent=cropDF(systemFile,"polyglycine","processed/hbonds_adjacent_occurence_definedbins.tsv")
dfgen_adjacent = dfgen_adjacent[FALSE,]

dfgen_nonadj=cropDF(systemFile,"polyglycine","processed/hbonds_nonadj_occurence_definedbins.tsv")
dfgen_nonadj = dfgen_nonadj[FALSE,]
for (i in directories) {
    print (i)
    #print (systemFile)
    if (grepl("poly[a-z].*",i)) {
        polymer=basename(i)
        #print (polymer)
        retdf_adjacent=cropDF(systemFile,polymer,"processed/hbonds_adjacent_occurence_definedbins.tsv")
        retdf_nonadj=cropDF(systemFile,polymer,"processed/hbonds_nonadj_occurence_definedbins.tsv")
        #print (head(dfgen_nonadj))
        #print (head(dfgen_adjacent))
        dfgen_adjacent=rbind(dfgen_adjacent,retdf_adjacent)
        dfgen_nonadj=rbind(dfgen_nonadj,retdf_nonadj)
        #this function will get the dataframes appended together
 }
}
plotfunc <- function (df,outfile) {
    print (head(df))
    gg <- ggplot(data=df, aes(x=Bins,y=FrequencyClass,fill=hbtype,alpha=0.5))
    gg <- gg + geom_bar(stat='identity',position=position_dodge())
    gg <- gg + facet_wrap(~polymer) 
    gg <- gg + theme (strip.text.x = element_text(size = 6),axis.text.x = element_text( hjust = 1, size = 6, angle = 45),axis.text.y = element_text( hjust = 1, size = 6), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,'_subtypeFrequency.pdf'),gg,height = 4, width= 4)

    gg <- ggplot(data=df, aes(x=Bins,y=FrequencyTotal,fill=hbtype,alpha=0.5))
    gg <- gg + geom_bar(stat='identity',position=position_dodge())
    gg <- gg + facet_wrap(~polymer) 
    gg <- gg + theme (strip.text.x = element_text(size = 6),axis.text.x = element_text( hjust = 1, size = 6, angle = 45),axis.text.y = element_text( hjust = 1, size = 6), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = paste0(outfile,'_wholeFrequency.pdf'),gg,height = 4, width= 4)
    
}

dfgen_adjacent=dfgen_adjacent[!(dfgen_adjacent$polymer=='polyala_constrained') & !(dfgen_adjacent$polymer=='polygly_constrained') & !(dfgen_adjacent$polymer=='polyisoleucine_server191') & !(dfgen_adjacent$polymer=='polyisoleucine_78') & !(dfgen_adjacent$polymer=='polyisoleucine_turing'),]
dfgen_nonadj=dfgen_nonadj[!(dfgen_nonadj$polymer=='polyala_constrained') & !(dfgen_nonadj$polymer=='polygly_constrained') & !(dfgen_nonadj$polymer=='polyisoleucine_server191') & !(dfgen_nonadj$polymer=='polyisoleucine_78') & !(dfgen_nonadj$polymer=='polyisoleucine_turing'),]


plotfunc(dfgen_adjacent,paste0(out,'hbond_frequency_custom_adjacent'))
plotfunc(dfgen_nonadj,paste0(out,'hbond_frequency_custom_nonadj'))
