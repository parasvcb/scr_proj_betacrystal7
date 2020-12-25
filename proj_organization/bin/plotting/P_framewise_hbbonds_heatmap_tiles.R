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

plothbdata_framewise_singleHb <- function(inputdataframe,forceFile,outputfile) {
    print ('Single_Hb')
    print (inputdataframe)
    print (forceFile)
    print (outputfile)
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
dfforce=read.csv(forceFile,sep = "\t",check.names = FALSE)

#dfgen[dfgen < 0] <- NA
dfgen$newhbname <- ifelse(dfgen$hbname == 'B_7_OT2__C_2_N' | dfgen$hbname == 'C_2_N__D_7_OT2' , 8,
                  ifelse(dfgen$hbname == 'B_7_OT1__C_2_N' | dfgen$hbname == 'C_2_N__D_7_OT1', 7,
                  ifelse(dfgen$hbname == 'B_7_N__C_2_O' | dfgen$hbname == 'C_3_N__D_6_O', 6,
                  ifelse(dfgen$hbname == 'B_5_O__C_4_N' | dfgen$hbname == 'C_3_O__D_6_N', 5, 
                  ifelse(dfgen$hbname == 'B_5_N__C_4_O' | dfgen$hbname == 'C_5_N__D_4_O', 4,
                  ifelse(dfgen$hbname == 'B_3_O__C_6_N' | dfgen$hbname == 'C_5_O__D_4_N', 3,
                  ifelse(dfgen$hbname == 'B_3_N__C_6_O' | dfgen$hbname == 'C_7_N__D_2_O', 2,
                  ifelse(dfgen$hbname == 'B_2_N__C_7_OT1' | dfgen$hbname == 'C_7_OT1__D_2_N', 1,
                  ifelse(dfgen$hbname == 'B_2_N__C_7_OT2' | dfgen$hbname == 'C_7_OT2__D_2_N', 0 ,'NULL')))))))))
dfgen$presence <- ifelse(dfgen$hbval > 0, 1,0)
dfgen$newhbname= sapply(dfgen$newhbname, as.numeric)
gg <- ggplot(dfgen , aes(x = frame, y = newhbname))
gg <- gg + geom_tile(aes(fill = presence, height=0.9))
gg <- gg + scale_fill_gradient('presence', low = "white",high = "#666666")
gg <- gg + geom_line(data=dfforce,aes(x=frame, y=(forceavg / 500), color='Force'), color='red', size=0.2, alpha=0.8)
gg <- gg + scale_y_continuous(name='Hb_index', breaks = seq(0,8,1) ,sec.axis = sec_axis(~./2, name="Force (nN)"))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position='top', panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
gg <- gg + facet_wrap(~poltype,ncol=1)
pdf(file = outputfile, height=15, width=10, colormodel = 'cmyk') 
print (gg)
dev.off()
#ggsave(filename = outputfile, height=15, width=10)
}

plothbdata_framewise_singleRing <- function(inputdataframe,forceFile,outputfile) {
    print ('Single_ring')
    print (inputdataframe)
    print (forceFile)
    print (outputfile)
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
dfforce=read.csv(forceFile,sep = "\t",check.names = FALSE)

#dfgen[dfgen < 0] <- NA
dfgen$newhbname <- ifelse(dfgen$hbname == 'ring_1_BC',1,
                  ifelse(dfgen$hbname == 'ring_2_BC',2,
                  ifelse(dfgen$hbname == 'ring_3_BC',3,
                  ifelse(dfgen$hbname == 'ring_1_DC',1.5,
                  ifelse(dfgen$hbname == 'ring_2_DC',2.5,
                  ifelse(dfgen$hbname == 'ring_3_DC',3.5,NA))))))

dfgen$presence <- ifelse(dfgen$hbval > 0, 1,0)
#dfgen$presence <- ifelse(dfgen$presence ==1,dfgen$presence,NA)


gg <- ggplot(dfgen , aes(x = frame, y = newhbname))
gg <- gg + geom_tile(aes(fill = presence, height=0.4))
gg <- gg + scale_fill_gradient('presence', low = "white",high = "#666666")
gg <- gg + geom_line(data=dfforce,aes(x=frame, y=(forceavg / 1000), color='Force'), color='red', size=0.2, alpha=0.8)
gg <- gg + scale_y_continuous(name='Hb_index', breaks = seq(0,4,0.5) ,sec.axis = sec_axis(~./1, name="Force (nN)"))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position='top', panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
gg <- gg + facet_wrap(~poltype,ncol=1)
pdf(file = outputfile, height=15, width=10, colormodel = 'cmyk') 
print (gg)
dev.off()
#ggsave(filename = outputfile, height=15, width=10)
}

plothbdata_framewise_superRing <- function(inputdataframe,forceFile,outputfile) {
    print ('Complex_Ring')
    print (inputdataframe)
    print (forceFile)
    print (outputfile)
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
dfforce=read.csv(forceFile,sep = "\t",check.names = FALSE)

dfgen$newhbname <- ifelse(dfgen$hbname == 'SuperRing_1', 1,
                  ifelse(dfgen$hbname == 'SuperRing_2', 2,3))

dfgen$presence <- ifelse(dfgen$hbval > 0, 1,0)
dfgen$newhbname= sapply(dfgen$newhbname, as.numeric)

gg <- ggplot(dfgen , aes(x = frame, y = newhbname))
gg <- gg + geom_tile(aes(fill = presence, height=0.9))
gg <- gg + scale_fill_gradient('presence', low = "white",high = "#666666")
gg <- gg + geom_line(data=dfforce,aes(x=frame, y=(forceavg / 1000), color='Force'), color='red', size=0.2, alpha=0.8)
gg <- gg + scale_y_continuous(name='Hb_index', breaks = seq(0,3.5,1) ,sec.axis = sec_axis(~./1, name="Force (nN)"))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position='top', panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
gg <- gg + facet_wrap(~poltype,ncol=1)
#ggsave(filename = outputfile, height=10, width=10)
pdf(file = outputfile, height=10, width=10, colormodel = 'cmyk') 
print (gg)
dev.off()
}

rep1='menton_set'
rep2='tyroneP_set'
rep3='tyrone_set3'

for (rep in c(rep1,rep2,rep3)){
plothbdata_framewise_singleHb(file.path(inputdir,paste0(rep,'_frameWise_hbpresenceBC.tsv')),file.path(inputdir,paste0(rep,'_forceread.tsv')), file.path(outdir,paste0(rep,'_frameWise_hbpresenceBC_HeatM.pdf')))
plothbdata_framewise_singleHb(file.path(inputdir,paste0(rep,'_frameWise_hbpresenceDC.tsv')),file.path(inputdir,paste0(rep,'_forceread.tsv')),file.path(outdir,paste0(rep,'_frameWise_hbpresenceDC_HeatM.pdf')))
plothbdata_framewise_singleRing(file.path(inputdir,paste0(rep,'_frameWise_SimpleRing.tsv')),file.path(inputdir,paste0(rep,'_forceread.tsv')),file.path(outdir,paste0(rep,'_frameWise_SimpleRing_HeatM.pdf')))
plothbdata_framewise_singleRing(file.path(inputdir,paste0(rep,'_frameWise_ComplexRing.tsv')),file.path(inputdir,paste0(rep,'_forceread.tsv')),file.path(outdir,paste0(rep,'_frameWise_ComplexRing_HeatM.pdf')))
plothbdata_framewise_superRing(file.path(inputdir,paste0(rep,'_frameWise_SuperRing.tsv')),file.path(inputdir,paste0(rep,'_forceread.tsv')),file.path(outdir,paste0(rep,'_frameWise_SuperRing_HeatM.pdf')))
#break
}

# '''
# frame   poltype hbname  hstability      hbval
# 0       3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1
# 10      3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1
# 20      3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1
# 30      3VAL    B_7_OT2__C_2_N  -0.1-0.2        -1

# poltype hbtype  occurence
# 3VAL    B_7_OT2__C_2_N  0.0
# 2THR    B_5_N__C_4_O    0.64
# 1ALAGLY B_7_OT2__C_2_N  0.0

# '''

# menton_set_forceread.tsv                menton_set_occFreq_SimpleRing.tsv       tyroneP_set_occFreq_hbpresenceBC.tsv    tyrone_set3_frameWise_SuperRing.tsv
# menton_set_frameWise_ComplexRing.tsv    menton_set_occFreq_SuperRing.tsv        tyroneP_set_occFreq_hbpresenceDC.tsv    tyrone_set3_occFreq_ComplexRing.tsv
# menton_set_frameWise_hbpresenceBC.tsv   tyroneP_set_forceread.tsv               tyroneP_set_occFreq_SimpleRing.tsv      tyrone_set3_occFreq_hbpresenceBC.tsv
# menton_set_frameWise_hbpresenceDC.tsv   tyroneP_set_frameWise_ComplexRing.tsv   tyroneP_set_occFreq_SuperRing.tsv       tyrone_set3_occFreq_hbpresenceDC.tsv
# menton_set_frameWise_SimpleRing.tsv     tyroneP_set_frameWise_hbpresenceBC.tsv  tyrone_set3_forceread.tsv               tyrone_set3_occFreq_SimpleRing.tsv
# menton_set_frameWise_SuperRing.tsv      tyroneP_set_frameWise_hbpresenceDC.tsv  tyrone_set3_frameWise_ComplexRing.tsv   tyrone_set3_occFreq_SuperRing.tsv
# menton_set_occFreq_ComplexRing.tsv      tyroneP_set_frameWise_SimpleRing.tsv    tyrone_set3_frameWise_hbpresenceBC.tsv  
# menton_set_occFreq_hbpresenceBC.tsv     tyroneP_set_frameWise_SuperRing.tsv     tyrone_set3_frameWise_hbpresenceDC.tsv  
# menton_set_occFreq_hbpresenceDC.tsv     tyroneP_set_occFreq_ComplexRing.tsv     tyrone_set3_frameWise_SimpleRing.tsv  