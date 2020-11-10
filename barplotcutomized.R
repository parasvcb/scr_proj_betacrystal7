library(ggplot2)
systemFile="../data_tsv/out_forcecum.tsv"
df=read.csv(systemFile,sep = "\t",check.names = FALSE)
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

df=df[!(df$category=='peak2') & !(df$category=='down2') &!(df$category=='down1'),]
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))

gg <- ggplot(data=df, aes(x=homopolymer,y=mean/1000))
gg <- gg + geom_bar(stat="identity", width=0.5)#, color='green',fill='white')
#gg <- gg + geom_text(aes(label=round(mean/1000,3)),vjust=1,stat="identity", color="black", size=3)
gg <- gg + geom_errorbar(data= df, aes(x=homopolymer,y=(mean/1000),ymin=(mean/1000)-(std/1000), ymax=(mean/1000)+(std/1000)), width=.1)

gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
#gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    gg <- gg + ylab("Force (nN)")
ggsave(filename = "../plots/new4sep/Fig2_forcemod.pdf", height=4, width=7)

systemFile="../data_tsv/out_toughnesscum.tsv"
df=read.csv(systemFile,sep = "\t",check.names = FALSE)
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', '1_SAA',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', '2_PAA',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', '3_HAA', 'amendment')))

df=df[!(df$category=='1stpeak') & !(df$category=='2aadist'),]
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))

df$homopolymer <- ifelse(df$polymer == 'polyglycine', 'poly-gly',
                  ifelse(df$polymer == 'polyalanine_default', 'poly-ala',
                  ifelse(df$polymer == 'polyala-gly', 'poly-alagly',
                  ifelse(df$polymer == 'polyasparagine', 'poly-asp',
                  ifelse(df$polymer == 'polythreonine', 'poly-thr',
                  ifelse(df$polymer == 'polyvaline', 'poly-val',
                  ifelse(df$polymer == 'polyisoleucine', 'poly-ile','unknown')))))))
gg <- ggplot(data=df,aes(x=homopolymer,y=mean))
gg <- gg + geom_bar(stat="identity", width=0.5)#,color='green',fill='white')
#gg <- gg + geom_text(aes(label=mean),vjust=1,stat="identity", color="black", size=3)

gg <- gg + geom_errorbar(aes(x=homopolymer,y=mean,ymin=mean-std, ymax=mean+std), width=.1)

gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    gg <- gg + ylab("Toughness (MJ/m3)")
ggsave(filename = "../plots/new4sep/Fig2_toughnessmod.pdf", height=4, width=7)
