systemFile="../out_forcecum.tsv"
df=read.csv(systemFile,sep = "\t",check.names = FALSE)
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', 'SAA',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', 'PAA',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', 'HAA', 'amendment')))
df$homopolymer <- ifelse(df$polymer == 'polyglycine', 'p[G]',
                  ifelse(df$polymer == 'polyalanine_default', 'p[A]',
                  ifelse(df$polymer == 'polyala-gly', 'p[AG]',
                  ifelse(df$polymer == 'polyasparagine', 'p[N]',
                  ifelse(df$polymer == 'polythreonine', 'p[T]',
                  ifelse(df$polymer == 'polyvaline', 'p[V]',
                  ifelse(df$polymer == 'polyisoleucine', 'p[I]','unknown')))))))

df=df[!(df$category=='peak2') & !(df$category=='down2') &!(df$category=='down1'),]
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))

gg <- ggplot()
gg <- gg + geom_bar(data=df, aes(x=homopolymer,y=mean), stat="identity", width=0.5)
gg <- gg + geom_errorbar(data= df, aes(x=homopolymer,y=mean,ymin=mean-std, ymax=mean+std), width=.1)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    gg <- gg + ylab("Force (pN)")
ggsave(filename = "forcemod.pdf", height=4, width=7)


systemFile="../out_toughnesscum.tsv"
df=read.csv(systemFile,sep = "\t",check.names = FALSE)
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', 'SAA',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', 'PAA',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', 'HAA', 'amendment')))

df=df[!(df$category=='Complete') & !(df$category=='2aadist'),]
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))
df$homopolymer <- ifelse(df$polymer == 'polyglycine', 'p[G]',
                  ifelse(df$polymer == 'polyalanine_default', 'p[A]',
                  ifelse(df$polymer == 'polyala-gly', 'p[AG]',
                  ifelse(df$polymer == 'polyasparagine', 'p[N]',
                  ifelse(df$polymer == 'polythreonine', 'p[T]',
                  ifelse(df$polymer == 'polyvaline', 'p[V]',
                  ifelse(df$polymer == 'polyisoleucine', 'p[I]','unknown')))))))
gg <- ggplot()
gg <- gg + geom_bar(data=df, aes(x=homopolymer,y=mean), stat="identity", width=0.5)
gg <- gg + geom_errorbar(data= df, aes(x=homopolymer,y=mean,ymin=mean-std, ymax=mean+std), width=.1)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    gg <- gg + ylab("Toughness")
ggsave(filename = "toughnessmod.pdf", height=4, width=7)



systemFile="../out_hbondsadjacent_peak_1.tsv"
df=read.csv(systemFile,sep = "\t",check.names = FALSE)
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', 'SAA',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', 'PAA',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', 'HAA', 'amendment')))

df=df[!(df$category=='has_hb_adjacent_scbb'),]
df$category_fill <- ifelse(df$category == 'has_hb_adjacent_bbbb', 'MC_MC_adj',
                  ifelse(df$category == 'has_hb_adjacent_scsc', 'SC_SC_adj', 'extra'))
df=df[!(df$category_fill=='extra'),]
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))
df$homopolymer <- ifelse(df$polymer == 'polyglycine', 'p[G]',
                  ifelse(df$polymer == 'polyalanine_default', 'p[A]',
                  ifelse(df$polymer == 'polyala-gly', 'p[AG]',
                  ifelse(df$polymer == 'polyasparagine', 'p[N]',
                  ifelse(df$polymer == 'polythreonine', 'p[T]',
                  ifelse(df$polymer == 'polyvaline', 'p[V]',
                  ifelse(df$polymer == 'polyisoleucine', 'p[I]','unknown')))))))
gg <- ggplot()
gg <- gg + geom_bar(data=df, aes(x=homopolymer,y=mean,fill=category_fill), stat="identity", width=1,position=position_dodge())
gg <- gg + geom_errorbar(data= df, aes(x=homopolymer,y=mean,ymin=mean-std, ymax=mean+std,fill=category_fill), width=.1,position=position_dodge(.9))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    gg <- gg + ylab("H-bonds")
ggsave(filename = "hadp1.pdf", height=3, width=6)

    
    
    
systemFile="../out_hbondsall_peak_1.tsv"
df=read.csv(systemFile,sep = "\t",check.names = FALSE)
df$aatype <- ifelse(df$polymer == 'polyglycine' | df$polymer == 'polyalanine_default' | df$polymer == 'polyala-gly', 'SAA',
                  ifelse(df$polymer == 'polythreonine' | df$polymer == 'polyasparagine', 'PAA',
                         ifelse(df$polymer == 'polyisoleucine' | df$polymer == 'polyvaline', 'HAA', 'amendment')))

df=df[!(df$category=='has_hb_all_scbb'),]
df$category_fill <- ifelse(df$category == 'has_hb_all_bbbb', 'MC_MC_all',
                  ifelse(df$category == 'has_hb_all_scsc', 'SC_SC_all', 'extra'))
df=df[!(df$category_fill=='extra'),]
df=df[!(df$polymer=='polyisoleucine_server191') & !(df$polymer=='polyala_constrained') & !(df$polymer=='polygly_constrained') & !(df$polymer=='polyisoleucine_78') & !(df$polymer=='polyisoleucine_turing'),]
print (head(df))
df$homopolymer <- ifelse(df$polymer == 'polyglycine', 'p[G]',
                  ifelse(df$polymer == 'polyalanine_default', 'p[A]',
                  ifelse(df$polymer == 'polyala-gly', 'p[AG]',
                  ifelse(df$polymer == 'polyasparagine', 'p[N]',
                  ifelse(df$polymer == 'polythreonine', 'p[T]',
                  ifelse(df$polymer == 'polyvaline', 'p[V]',
                  ifelse(df$polymer == 'polyisoleucine', 'p[I]','unknown')))))))
gg <- ggplot()
gg <- gg + geom_bar(data=df, aes(x=homopolymer,y=mean,fill=category_fill), stat="identity", width=1,position=position_dodge())
gg <- gg + geom_errorbar(data= df, aes(x=homopolymer,y=mean,ymin=mean-std, ymax=mean+std,fill=category_fill), width=.1,position=position_dodge(.9))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + facet_wrap(~aatype,nrow=1,scales='free_x') 
    gg <- gg + ylab("H-bonds")
ggsave(filename = "hallp1.pdf", height=3, width=6)