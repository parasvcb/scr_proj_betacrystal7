# rscript to plot
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements inputDir and outDir", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  inputdir=args[1]
  outdir=args[2]
}

library(ggplot2)

plothbdata_occfreq <- function(inputdataframe,outputfile) {
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)

# reptypes <- unique(dfgen1$reptype)
# for (i in reptypes) {
#     dfgen <- dfgen1[(dfgen1$reptype==i),]

    gg <- ggplot(data=dfgen,aes(x=hbtype,y=occurence,fill=poltype))
    gg <- gg + geom_bar(stat="identity",position=position_dodge())
    gg <- gg + geom_text(aes(label=occurence),vjust=1, color="black",position = position_dodge(0.9), size=3)
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    #gg <- gg + facet_wrap(~polymer,ncol=1,scale='free_x') 
    ggsave(filename = outputfile, height=5, width=15)
#}
}


plothbdata_framewise <- function(inputdataframe,forceFile,outputfile,val) {
    print (inputdataframe)
    print (forceFile)
    print (outputfile)
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
dfforce=read.csv(forceFile,sep = "\t",check.names = FALSE)
dfgen[dfgen < 0] <- NA
# print (head(dfgen))
# print (head(dfforce))
gg <- ggplot(data=dfgen)
gg <- gg + geom_point(aes(x=frame,y=hbval,color=hbname,shape=hstability),size=0.8,alpha=0.8)
#gg <- gg + scale_shape_manual(values=Valshape)
gg <- gg + geom_line(data=dfforce,aes(x=frame, y=(forceavg / 1000), color='Force'), color='gray', size=0.4, alpha=0.8)
gg <- gg + scale_y_continuous(name='Hb_values',sec.axis = sec_axis(~./1, name="Force (nN)"))
#gg <- gg + scale_color_manual(values=val)
gg <- gg + facet_wrap(~poltype,ncol=1)
#gg <- gg + theme(plot.margin=unit(c(2,2,8,1),"cm"))
#gg <- gg + theme(legend.position=c(1,1))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position='right', panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
gg <- gg + theme (legend.text = element_text(size=8), legend.title = element_text(size=10,face="bold"))
gg <- gg + guides (colour = guide_legend(override.aes = list(size=10)))
gg <- gg + guides (shape = guide_legend(override.aes = list(size=10)))
#  legend.position="top",
# gg <- gg + guides(fill = guide_legend(nrow = 3))
# gg <- gg + facet_wrap(~polymer,ncol=1,scale='free_x') 
ggsave(filename = outputfile, height=15, width=10)
}
val8=c("red","blue","orange","black","cyan","green","brown","black")
val9=c("red","blue","orange","black","cyan","green","brown","gray","black")
val10=c("red","blue","orange","black","cyan","green","brown","gray","pink","black")

rep1='menton_set'
rep2='tyroneP_set'
rep3='tyrone_set3'
plothbdata_framewise(file.path(inputdir,paste0(rep1,'_frameWise_hbpresenceBC.tsv')),file.path(inputdir,paste0(rep1,'_forceread.tsv')), file.path(outdir,paste0(rep1,'_frameWise_hbpresenceBC.pdf')),val9)

plothbdata_framewise(file.path(inputdir,paste0(rep1,'_frameWise_hbpresenceDC.tsv')),file.path(inputdir,paste0(rep1,'_forceread.tsv')),file.path(outdir,paste0(rep1,'_frameWise_hbpresenceDC.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep1,'_frameWise_SimpleRing.tsv')),file.path(inputdir,paste0(rep1,'_forceread.tsv')),file.path(outdir,paste0(rep1,'_frameWise_SimpleRing.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep1,'_frameWise_ComplexRing.tsv')),file.path(inputdir,paste0(rep1,'_forceread.tsv')),file.path(outdir,paste0(rep1,'_frameWise_ComplexRing.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep1,'_frameWise_SuperRing.tsv')),file.path(inputdir,paste0(rep1,'_forceread.tsv')),file.path(outdir,paste0(rep1,'_frameWise_SuperRing.pdf')),val10)

plothbdata_occfreq(file.path(inputdir,paste0(rep1,'_occFreq_hbpresenceBC.tsv')),file.path(outdir,paste0(rep1,'_occFreq_hbpresenceBC.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep1,'_occFreq_hbpresenceDC.tsv')),file.path(outdir,paste0(rep1,'_occFreq_hbpresenceDC.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep1,'_occFreq_SimpleRing.tsv')),file.path(outdir,paste0(rep1,'_occFreq_SimpleRing.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep1,'_occFreq_ComplexRing.tsv')),file.path(outdir,paste0(rep1,'_occFreq_ComplexRing.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep1,'_occFreq_SuperRing.tsv')),file.path(outdir,paste0(rep1,'_occFreq_SuperRing.pdf')))


plothbdata_framewise(file.path(inputdir,paste0(rep2,'_frameWise_hbpresenceBC.tsv')) ,file.path(inputdir,paste0(rep2,'_forceread.tsv')),file.path(outdir,paste0(rep2,'_frameWise_hbpresenceBC.pdf')),val9)
plothbdata_framewise(file.path(inputdir,paste0(rep2,'_frameWise_hbpresenceDC.tsv')) ,file.path(inputdir,paste0(rep2,'_forceread.tsv')),file.path(outdir,paste0(rep2,'_frameWise_hbpresenceDC.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep2,'_frameWise_SimpleRing.tsv')) ,file.path(inputdir,paste0(rep2,'_forceread.tsv')),file.path(outdir,paste0(rep2,'_frameWise_SimpleRing.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep2,'_frameWise_ComplexRing.tsv')) ,file.path(inputdir,paste0(rep2,'_forceread.tsv')),file.path(outdir,paste0(rep2,'_frameWise_ComplexRing.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep2,'_frameWise_SuperRing.tsv')) ,file.path(inputdir,paste0(rep2,'_forceread.tsv')),file.path(outdir,paste0(rep2,'_frameWise_SuperRing.pdf')),val10)

plothbdata_occfreq(file.path(inputdir,paste0(rep2,'_occFreq_hbpresenceBC.tsv')),file.path(outdir,paste0(rep2,'_occFreq_hbpresenceBC.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep2,'_occFreq_hbpresenceDC.tsv')),file.path(outdir,paste0(rep2,'_occFreq_hbpresenceDC.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep2,'_occFreq_SimpleRing.tsv')),file.path(outdir,paste0(rep2,'_occFreq_SimpleRing.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep2,'_occFreq_ComplexRing.tsv')),file.path(outdir,paste0(rep2,'_occFreq_ComplexRing.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep2,'_occFreq_SuperRing.tsv')),file.path(outdir,paste0(rep2,'_occFreq_SuperRing.pdf')))


plothbdata_framewise(file.path(inputdir,paste0(rep3,'_frameWise_hbpresenceBC.tsv')),file.path(inputdir,paste0(rep3,'_forceread.tsv')),file.path(outdir,paste0(rep3,'_frameWise_hbpresenceBC.pdf')),val9)
plothbdata_framewise(file.path(inputdir,paste0(rep3,'_frameWise_hbpresenceDC.tsv')),file.path(inputdir,paste0(rep3,'_forceread.tsv')),file.path(outdir,paste0(rep3,'_frameWise_hbpresenceDC.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep3,'_frameWise_SimpleRing.tsv')),file.path(inputdir,paste0(rep3,'_forceread.tsv')),file.path(outdir,paste0(rep3,'_frameWise_SimpleRing.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep3,'_frameWise_ComplexRing.tsv')),file.path(inputdir,paste0(rep3,'_forceread.tsv')),file.path(outdir,paste0(rep3,'_frameWise_ComplexRing.pdf')),val10)
plothbdata_framewise(file.path(inputdir,paste0(rep3,'_frameWise_SuperRing.tsv')),file.path(inputdir,paste0(rep3,'_forceread.tsv')),file.path(outdir,paste0(rep3,'_frameWise_SuperRing.pdf')),val10)

plothbdata_occfreq(file.path(inputdir,paste0(rep3,'_occFreq_hbpresenceBC.tsv')),file.path(outdir,paste0(rep3,'_occFreq_hbpresenceBC.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep3,'_occFreq_hbpresenceDC.tsv')),file.path(outdir,paste0(rep3,'_occFreq_hbpresenceDC.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep3,'_occFreq_SimpleRing.tsv')),file.path(outdir,paste0(rep3,'_occFreq_SimpleRing.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep3,'_occFreq_ComplexRing.tsv')),file.path(outdir,paste0(rep3,'_occFreq_ComplexRing.pdf')))
plothbdata_occfreq(file.path(inputdir,paste0(rep3,'_occFreq_SuperRing.tsv')),file.path(outdir,paste0(rep3,'_occFreq_SuperRing.pdf')))

print ('DONE')

'''
frame   poltype hbname  hstability      hbval
0       3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1
10      3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1
20      3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1
30      3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1

poltype hbtype  occurence
3VAL    B_7_OT2__C_2_N  0.0
2THR    B_5_N__C_4_O    0.64
1ALAGLY B_7_OT2__C_2_N  0.0

'''

# menton_set_forceread.tsv                menton_set_occFreq_SimpleRing.tsv       tyroneP_set_occFreq_hbpresenceBC.tsv    tyrone_set3_frameWise_SuperRing.tsv
# menton_set_frameWise_ComplexRing.tsv    menton_set_occFreq_SuperRing.tsv        tyroneP_set_occFreq_hbpresenceDC.tsv    tyrone_set3_occFreq_ComplexRing.tsv
# menton_set_frameWise_hbpresenceBC.tsv   tyroneP_set_forceread.tsv               tyroneP_set_occFreq_SimpleRing.tsv      tyrone_set3_occFreq_hbpresenceBC.tsv
# menton_set_frameWise_hbpresenceDC.tsv   tyroneP_set_frameWise_ComplexRing.tsv   tyroneP_set_occFreq_SuperRing.tsv       tyrone_set3_occFreq_hbpresenceDC.tsv
# menton_set_frameWise_SimpleRing.tsv     tyroneP_set_frameWise_hbpresenceBC.tsv  tyrone_set3_forceread.tsv               tyrone_set3_occFreq_SimpleRing.tsv
# menton_set_frameWise_SuperRing.tsv      tyroneP_set_frameWise_hbpresenceDC.tsv  tyrone_set3_frameWise_ComplexRing.tsv   tyrone_set3_occFreq_SuperRing.tsv
# menton_set_occFreq_ComplexRing.tsv      tyroneP_set_frameWise_SimpleRing.tsv    tyrone_set3_frameWise_hbpresenceBC.tsv  
# menton_set_occFreq_hbpresenceBC.tsv     tyroneP_set_frameWise_SuperRing.tsv     tyrone_set3_frameWise_hbpresenceDC.tsv  
# menton_set_occFreq_hbpresenceDC.tsv     tyroneP_set_occFreq_ComplexRing.tsv     tyrone_set3_frameWise_SimpleRing.tsv  