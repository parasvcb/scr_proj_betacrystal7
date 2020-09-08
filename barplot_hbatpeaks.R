# allhb_mcmc      allhb_mcsc      allhb_scsc      adjhb_mcmc      adjhb_mcsc      adjhb_scsc      allp1stochhbmcmc        allp1stochhbmcsc        allp1stochhbscsc        allp1stabhbmcmc allp1stabhbmcsc allp1stabhbscsc adjp1stochhbmcmc        adjp1stochhbmcsc        adjp1stochhbscsc        adjp1stabhbmcmc adjp1stabhbmcsc adjp1stabhbscsc allp2stochhbmcmc        allp2stochhbmcsc        allp2stochhbscsc        allp2stabhbmcmc allp2stabhbmcsc allp2stabhbscsc adjp2stochhbmcmc        adjp2stochhbmcsc        adjp2stochhbscsc        adjp2stabhbmcmc adjp2stabhbmcsc adjp2stabhbscsc displacement
# 14.159  5.683   3.016   14.159  1.778   0.333   0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0                                                                                                     2.1
#                                                                                                                                                 0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     

library(ggplot2)
library(reshape)

plothbdata <- function(inputdataframe,outputfile) {
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
dfgen=dfgen[!(dfgen$polymer=='polyala_constrained') & !(dfgen$polymer=='polygly_constrained') & !(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]

newColNames <- c("category", "hbtype")
newCols <- colsplit(dfgen$datatype, "_", newColNames)
dfgen <- cbind(dfgen, newCols)
drops <- c("datatype")
dfgen<-dfgen[ , !(names(dfgen) %in% drops)]

#print (head(dfgen))
dfgen$hbfrom <- ifelse(grepl('all',dfgen$category), '1_all', '2_adj')
dfgen$hbnature <- ifelse(grepl('stoch',dfgen$category), 'stochastic',
                  ifelse(grepl('stab',dfgen$category), 'stable',   'raw'))
dfgen$hbx <- paste(dfgen$hbnature,'_',dfgen$hbfrom)
#print (head(dfgen))

gg <- ggplot(data=dfgen,aes(x=hbx,y=mean,fill=hbtype))
gg <- gg + geom_bar(stat="identity",position=position_dodge())
gg <- gg + geom_text(aes(label=mean),vjust=1, color="black",position = position_dodge(0.9), size=3)
gg <- gg + geom_errorbar(data= dfgen, aes(x=hbx,y=mean,ymin=mean-std, ymax=mean+std),position = position_dodge(0.9),size=0.3,width=0.5)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
gg <- gg + facet_wrap(~polymer,ncol=1,scale='free_x') 
ggsave(filename = outputfile, height=15, width=5)

}
plothbdata("../data_tsv/out_hbonds_peak_1.tsv","../plots/new4sep/peak1_hbdata.png")
plothbdata("../data_tsv/out_hbonds_peak_2.tsv","../plots/new4sep/peak2_hbdata.png")
