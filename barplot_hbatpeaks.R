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
dfgen$hbfrom <- ifelse(grepl('all',dfgen$category), 'nonadj', 'adj')
dfgen$hbnature <- ifelse(grepl('stoch',dfgen$category), 'stochastic',
                  ifelse(grepl('stab',dfgen$category), 'stable',   'raw'))
# this should be segregated such that, raw information will be available in one plot (df)
# stab and stoch will go together for adj and nondadj so in one df (plot) 

dfgen_primary = dfgen[(dfgen$hbnature=='raw'),]
dfgen_primary$hbx <- paste(dfgen_primary$polymer,'_',dfgen_primary$hbfrom)

#dfgen_secondary_adj = dfgen[(dfgen$hbfrom=='adj') & !(dfgen$hbnature=='raw'),]
#dfgen_secondary_nonadj = dfgen[(dfgen$hbfrom=='nonadj') & !(dfgen$hbnature=='raw'),]
dfgen_secondary_comp = dfgen[!(dfgen$hbnature=='raw'),]
dfgen_secondary_comp$facetvar <- paste(dfgen_secondary_comp$hbfrom,'_',dfgen_secondary_comp$hbnature)
dfgen_secondary_comp$hbx <- dfgen_secondary_comp$polymer


#print (head(dfgen))

#Because of vertical stacking, need to add the means such that SCSC=SCSC, MCSC=SCSC+MCSC, MCMC=SCSC+MCSC+MCMC
# now Upp_lim will be different SCSC=SCSC, MCSC=mean_MCSC-std, MCMC=mean_MCMC-std
# and Low_lim SCSC=SCSC, MCSC=mean_MCSC-std, MCMC=mean_MCMC-std
dfgen_primary$vertmean = dfgen_primary$mean

# dfgen_primary$vertmean[dfgen_primary$hbtype == "mcmc"] <- with (dfgen_primary, mean[dfgen_primary$hbtype == "mcmc"] + mean[dfgen_primary$hbtype == "mcsc"] + mean[dfgen_primary$hbtype == "scsc"])
# dfgen_primary$vertmean[dfgen_primary$hbtype == "mcsc"] <- with (dfgen_primary, mean[dfgen_primary$hbtype == "mcsc"] + mean[dfgen_primary$hbtype == "scsc"])
# dfgen_primary$vertmean[dfgen_primary$hbtype == "scsc"] <- with (dfgen_primary, mean[dfgen_primary$hbtype == "scsc"])

# the above and below values are nothing but tweaks to get the ordering wise std dev bars to appear, cant add text in betwen as of now.
dfgen_primary$vertmean[dfgen_primary$hbtype == "scsc"] <- with (dfgen_primary, mean[dfgen_primary$hbtype == "mcmc"] + mean[dfgen_primary$hbtype == "mcsc"] + mean[dfgen_primary$hbtype == "scsc"])
dfgen_primary$vertmean[dfgen_primary$hbtype == "mcsc"] <- with (dfgen_primary, mean[dfgen_primary$hbtype == "mcsc"] + mean[dfgen_primary$hbtype == "mcmc"])
dfgen_primary$vertmean[dfgen_primary$hbtype == "mcmc"] <- with (dfgen_primary, mean[dfgen_primary$hbtype == "mcmc"])


#print (dfgen_primary)
gg <- ggplot(data=dfgen_primary,aes(x=hbx,y=mean,fill=hbtype))
gg <- gg + geom_bar(position="stack", stat="identity") + geom_col(position = position_stack(reverse = TRUE))
gg <- gg + scale_fill_manual(values=c("gray","green","red"))
#gg <- gg + geom_text(aes(label=mean), position = position_stack(vjust = 0.5)) #vjust=1, color="black", size=1)
gg <- gg + geom_errorbar(aes(x=hbx,y=vertmean,ymin=vertmean-std, ymax=vertmean+std,color = hbtype, width=0.3),position="identity") #+ geom_col(position = position_stack(reverse = TRUE))
gg <- gg + scale_color_manual(values=c("brown","black","blue"))
gg <- gg + ylab("Hb-count")
gg <- gg + xlab("Polymer")
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = paste(outputfile,'_','raw_plot.pdf'))

gg <- ggplot(data=dfgen_primary,aes(x=polymer,y=mean,fill=hbtype))
gg <- gg + geom_bar(position="stack", stat="identity") + geom_col(position = position_stack(reverse = TRUE))
gg <- gg + scale_fill_manual(values=c("gray","green","red"))
#gg <- gg + geom_text(aes(label=mean), position = position_stack(vjust = 0.5)) #vjust=1, color="black", size=1)
gg <- gg + geom_errorbar(aes(x=polymer,y=vertmean,ymin=vertmean-std, ymax=vertmean+std,color = hbtype, width=0.3),position="identity") #+ geom_col(position = position_stack(reverse = TRUE))
gg <- gg + scale_color_manual(values=c("brown","black","blue"))
gg <- gg + ylab("Hb-count")
gg <- gg + xlab("Polymer")
gg <- gg + facet_wrap(~hbfrom,nrow=1)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = paste(outputfile,'_','raw_plot_facet.pdf'))


dfgen_secondary_comp$vertmean[dfgen_secondary_comp$hbtype == "scsc"] <- with (dfgen_secondary_comp, mean[dfgen_secondary_comp$hbtype == "mcmc"] + mean[dfgen_secondary_comp$hbtype == "mcsc"] + mean[dfgen_secondary_comp$hbtype == "scsc"])
dfgen_secondary_comp$vertmean[dfgen_secondary_comp$hbtype == "mcsc"] <- with (dfgen_secondary_comp, mean[dfgen_secondary_comp$hbtype == "mcsc"] + mean[dfgen_secondary_comp$hbtype == "mcmc"])
dfgen_secondary_comp$vertmean[dfgen_secondary_comp$hbtype == "mcmc"] <- with (dfgen_secondary_comp, mean[dfgen_secondary_comp$hbtype == "mcmc"])

#print (dfgen_secondary_comp)

gg <- ggplot(data=dfgen_secondary_comp,aes(x=hbx,y=mean,fill=hbtype))
gg <- gg + geom_bar(position="stack", stat="identity") + geom_col(position = position_stack(reverse = TRUE))
gg <- gg + scale_fill_manual(values=c("gray","green","red"))
#gg <- gg + geom_text(aes(label=mean), position = position_stack(vjust = 0.5)) #vjust=1, color="black", size=1)
gg <- gg + geom_errorbar(data= dfgen_secondary_comp, aes(x=hbx,y=vertmean,ymin=vertmean-std, ymax=vertmean+std,color = hbtype, width=0.3),position="identity") #+ geom_col(position = position_stack(reverse = TRUE))
gg <- gg + scale_color_manual(values=c("brown","black","blue"))
gg <- gg + facet_wrap(~facetvar,ncol=2,nrow=2)
gg <- gg + ylab("Hb-count")
gg <- gg + xlab("Polymer")
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = paste(outputfile,'_','segregated_plot.pdf'))
}
plothbdata("../data_tsv/out_hbonds_peak_1.tsv","../plots/new4sep/peak1_hbdata")
plothbdata("../data_tsv/out_hbonds_peak_2.tsv","../plots/new4sep/peak2_hbdata")


# going to reformat the plotting style
# adjacent and non-adjacent will go as one in barplot, with stacking of different scsc and other types
# For each adj and non adj, there will be supplimnetarty stochastic and stable plots, [lay together], and facet them 
# top and bottom,
# Each column will have color bar as mean and std dev together 