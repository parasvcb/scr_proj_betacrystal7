args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements input tsv and outputfile", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}
library(ggplot2)

df=read.csv(systemFile,sep = "\t",check.names = FALSE)
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', 'small_group',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', 'polar',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', 'hydrophobic', 'amendment')))
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))
gg <- ggplot()
    gg <- gg + geom_bar(data=df, aes(x=polymer,y=mean,fill=category), stat="identity", position=position_dodge(), width=1)
    gg <- gg + geom_errorbar(data= df, aes(x=polymer,y=mean,ymin=mean-std, ymax=mean+std, fill=category), width=.1,position=position_dodge(.9))
    gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
     gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    #gg <- gg + guides(colour = guide_legend(override.aes = list(size=10)))
    ggsave(filename = out, height=5, width=10)