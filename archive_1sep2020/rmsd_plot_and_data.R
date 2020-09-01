args = commandArgs(trailingOnly=TRUE)
if (length(args)<3) {
  stop("Please supply two arguements systemDir frompdbormin (type wither frame0 or min) and output_append_withoutext", call.=FALSE)
} else if (length(args)==3) {
  # default output file
  parentsystem=args[1]
  rmsdtype=args[2]
  out=args[3]
}
library(ggplot2)


filename= ifelse(rmsdtype=='min','rmsd_frompdb.txt','rmsd_firstframe.txt')

cropDF <- function(dfpre,parenttype,polymer) {
    dffile=paste0("processed/",filename)
    dfadd=file.path(dfpre,parenttype,polymer,dffile)
    datarmsd=read.csv(dfadd,sep = "\t",check.names = FALSE)
    datarmsd <- subset(datarmsd, select = -c(grp,protgrp))
    datarmsd$polymer <- polymer
    datarmsd$replica <- parenttype
    return (datarmsd)
}
#pullinginterctaion
resname=paste0(out,"resdf.csv")

dfgen=cropDF(parentsystem,"./menton_set","polyglycine")
dfgen = dfgen[FALSE,]

directories_parent=list.dirs(path = parentsystem, full.names = F, recursive=F)

for (parent in directories_parent) {
    #print (basename(parent))
    #print ('basename')
    if (basename(parent) %in% c('menton_set','tyroneP_set','tyrone_set3')){
        #print (basename(parent))
        directories_child= list.dirs(path = file.path(parentsystem,parent), full.names = F, recursive=F)
        #print (file.path(parentsystem,parent))
        #print (directories_child)
        for (child in directories_child) {
            if (grepl("poly[a-z].*",child) & !child %in% c('polyala_constrained','polyisoleucine_server191','polygly_constrained','polyisoleucine_78','polyisoleucine_turing')) {
                polymer=basename(child)
                system=basename(parent) 
                retdf=cropDF(parentsystem,system,polymer)
                dfgen=rbind(dfgen,retdf)
        #this function will get the dataframes appended together
}
}
}
}

dfgen$aatype <- ifelse(dfgen$polymer == 'polyglycine' | dfgen$polymer == 'polyalanine_default' | dfgen$polymer == 'polyala-gly', 'SAA',
                  ifelse(dfgen$polymer == 'polythreonine' | dfgen$polymer == 'polyasparagine', 'PAA',
                         ifelse(dfgen$polymer == 'polyisoleucine' | dfgen$polymer == 'polyvaline', 'HAA', 'amendment')))

print (unique(dfgen$replica))
dfgen$replicates <- ifelse(dfgen$replica == 'menton_set', 'R1',
                  ifelse(dfgen$replica == 'tyroneP_set', 'R2',
                  ifelse(dfgen$replica == 'tyrone_set3', 'R3',"NULL")))

dfgen$homopolymer <- ifelse(dfgen$polymer == 'polyglycine', 'p[G]',
                  ifelse(dfgen$polymer == 'polyalanine_default', 'p[A]',
                  ifelse(dfgen$polymer == 'polyala-gly', 'p[AG]',
                  ifelse(dfgen$polymer == 'polyasparagine', 'p[N]',
                  ifelse(dfgen$polymer == 'polythreonine', 'p[T]',
                  ifelse(dfgen$polymer == 'polyvaline', 'p[V]',
                  ifelse(dfgen$polymer == 'polyisoleucine', 'p[I]','unknown')))))))
    
    gg <- ggplot(data=dfgen[!is.na(dfgen$rmsd),], aes(x=(frames),y=rmsd,color=replicates))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("blue","red","green"))
    gg <- gg + facet_wrap(~homopolymer) 
    gg <- gg + scale_x_continuous(name="Time (ps)")
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab('RMSD (ang)')
    ggsave(filename = paste0(out,rmsdtype,"_line_RMSD.pdf"),gg, height = 5, width= 8)

    gg <- ggplot(data=dfgen[!is.na(dfgen$rmsd),], aes(x=(frames),y=rmsd,color=replicates))
    gg <- gg + geom_line(size=0.6, ) + scale_color_manual(values=c("blue","red","green"))
    gg <- gg + facet_wrap(~homopolymer) 
    gg <- gg + scale_x_continuous(name="Time (ps)")
    gg <- gg + theme (strip.text.x = element_text(size = 10),axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    gg <- gg + ylab('RMSD (ang)')
    ggsave(filename = paste0(out,rmsdtype,"_point_RMSD.pdf"),gg, height = 5, width= 8)