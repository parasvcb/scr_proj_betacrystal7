#will be renamed to Fig5_barplotfractions_For_Tough.R (primary)
#this should get splitted into two components
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements hbprogram and output_append_withoutext", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemDir=args[1]
  outsys=args[2]
}
library(ggplot2)

plotBar <- function(df,outputfile,legendY,divfac) {
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', '1_SAA',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', '2_PAA',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', '3_HAA', 'amendment')))
df$homopolymer <- ifelse(df$polymer == 'polyglycine', 'poly-gly',
                  ifelse(df$polymer == 'polyalanine_default', 'poly-ala',
                  ifelse(df$polymer == 'polyala-gly', 'poly-alagly',
                  ifelse(df$polymer == 'polyasparagine', 'poly-asp',
                  ifelse(df$polymer == 'polythreonine', 'poly-thr',
                  ifelse(df$polymer == 'polyvaline', 'poly-val',
                  ifelse(df$polymer == 'polyisoleucine', 'poly-ile','unknown')))))))

df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))

gg <- ggplot(data=df, aes(x=homopolymer,y=mean/divfac))
gg <- gg + geom_bar(stat="identity", width=0.5)#, color='green',fill='white')
#gg <- gg + geom_text(aes(label=round(mean/divfac,3)),vjust=1,stat="identity", color="black", size=3)
gg <- gg + geom_errorbar(data= df, aes(x=homopolymer,y=(mean/divfac),ymin=(mean/divfac)-(std/divfac), ymax=(mean/divfac)+(std/divfac)), width=.1)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
#gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
gg <- gg + ylab(legendY)
ggsave(filename = outputfile, height=4, width=7)
}

dfForce=read.csv(file.path(systemDir,'out_forcecum.tsv'),sep = "\t",check.names = FALSE)
dfForce=dfForce[!(dfForce$category=='peak2') & !(dfForce$category=='down2') &!(dfForce$category=='down1'),]

dfTough=read.csv(file.path(systemDir,'out_toughnesscum.tsv'),sep = "\t",check.names = FALSE)
dfTough=dfTough[!(dfTough$category=='1stpeak') & !(dfTough$category=='2aadist'),]

plotBar(dfForce,file.path(outsys,"Fig4_Peak1_ForceReading.pdf"),'Force (nN)',1000)
plotBar(dfTough,file.path(outsys,"Fig4_Complete_ToughnessReading.pdf"),'Toughness (MJ/m3)',1)

