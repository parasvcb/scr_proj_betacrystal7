df=read.csv('../menton_set/polyalanine_default/processed/dataframe_vnew2.tsv',sep='\t',check.names=FALSE)
colnames(df)[which(names(df) == "for-dispav")] <- "ford"
df1=df <- subset(df, select = c(ford))


#tricky R dataframe selection and antiselection, use -c in row3 to select everything other than row that column
