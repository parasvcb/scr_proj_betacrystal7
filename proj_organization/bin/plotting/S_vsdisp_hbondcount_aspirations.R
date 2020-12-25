'dataframe_vnew2'
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements systemDir and output_append_withoutext", call.=FALSE)
} else if (length(args)==2) {
  systemFile=args[1]
  out=args[2]
}
library(ggplot2)
library(reshape)
directories=list.dirs(path = systemFile, full.names = TRUE, recursive=F)

cropDF <- function(dfpre,polymer) {
    dffile="processed/dataframe_vnew3.tsv"
    dfadd=file.path(dfpre,polymer,dffile)
    dataforce=read.csv(dfadd,sep = "\t",check.names = FALSE)
    colnames(dataforce)[which(names(dataforce) == "frames-raw")] <- "frames"
    dataforce <- subset(dataforce, select = -c(frames))
    df=melt(data = dataforce, id.vars = "displacement")
    newColNames <- c("category", "avtype")
    newCols <- colsplit(df$variable, "-", newColNames)
    afterpre <- cbind(df, newCols)
    afterpre$polymer <- polymer
    return (afterpre)
}

#rename the '-' containing variable to other so that it can be dropped later,
#melt using the displacement, and categorise variable types in different columns 

resname=paste0(out,"resdf.csv")
if (file.exists(resname)){
    print('reading')
  dfgen=read.csv(resname,check.names = FALSE)
} else {
dfgen=cropDF(systemFile,"polyglycine")
dfgen = dfgen[FALSE,]
#above two statements will get me a empty df with only column names, but bit costly
for (i in directories) {
    print (i)
    #print (systemFile)
    if (grepl("poly[a-z].*",i)) {
        polymer=basename(i)
        #print (polymer)
        retdf=cropDF(systemFile,polymer)
        dfgen=rbind(dfgen,retdf)
        #this function will get the dataframes appended together
 }
}
#The above block will iterate all the polymers in one replicate and if of desired type, will concatenate them
#write so that this shouldnt be repeated, or remove wirtten file if raw file data has changed
write.csv(dfgen,resname,sep = "\t")
}
print ('read')

plotallad <- function(df,extractele,outapp) {
    df$keeps <- ifelse (grepl(extractele,df$variable),'keep',
    ifelse (grepl("for-dispav",df$variable),'keep','discard'))
    df$hbcolor <- ifelse(grepl('mcmc',df$variable), 'mcmc',
                   ifelse(grepl('scsc',df$category), 'scsc',
                ifelse(grepl('mcsc',df$category), 'mcsc',   'force')))
    dfset=df[(df$keeps=='keep'),]
    dfForce=dfset[(dfset$variable=='for-dispav'),]
    dfset=dfset[!(dfset$variable=='for-dispav'),]
    #the above can be removed to add thin layout to force in bacground
    
    print ("here")
    #dfset=dfset[(dfset$avtype=='ra20'),]
    dfset=dfset[(dfset$avtype=='dispav'),]
    dfset$hbnature<- ifelse(grepl('p1',dfset$variable), 'partial',
                   ifelse(grepl('p2',dfset$variable), 'partial',   'unbiased'))
    dfset2=dfset[(dfset$hbnature=='partial'),]
    dfset1=dfset[(dfset$hbnature=='unbiased'),]
    
    # print (head(dfset1))
    # print (unique(dfset1$value))
    # write.csv(dfset1,'check.tsv',sep = "\t")

    print (unique(dfset1$polymer))
    print ('4')
    maxval=max(dfset1$value,na.rm = TRUE)
    coeff=round(maxval / 3,digits=2)
    gg <- ggplot(data=dfset1[!is.na(dfset1$value),], aes(x=displacement,y=value,color=hbcolor))
    gg <- gg + geom_line(size=0.4, alpha=0.8 )
    gg <- gg + geom_line(data=dfForce[!is.na(dfForce$value),], aes(x=displacement, y=(value / 1000)*coeff), size=0.1, alpha=0.8)
    gg <- gg + scale_color_manual(values=c("red","blue","black","green"))
    gg <- gg + scale_y_continuous(name='Hb_count',sec.axis = sec_axis(~./coeff, name="Force (nN)")) 
    gg <- gg + facet_wrap(~polymer,ncol=1) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + theme (strip.text.x = element_text(size = 7),axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=1)))
    ggsave(filename = paste0(outapp,'unbiased.pdf'),gg, height = 10, width= 5)

    d1='dashed'
    d2='solid'
    c1='blue'
    c11='lightblue'
    c2='black'
    c21='darkgray'
    c3='pink'
    c31='maroon'
    c4='orange'
    c41='yellow'
    c5='green'
    c51='lightgreen'
    c6='purple'
    c61='cyan'
    print (unique(dfset2$category))
    maxval=max(dfset2$value,na.rm = TRUE)
    coeff=round(maxval / 3,digits=2)
    #data$PROFILE <- sub("^([a-z]+[A-B]).*", "\\1", data$PROFILE)
    #dfset2$category <- sub("all", "nonadj", dfset2$category)
    
    dfset2$category <- sub("hball", "", dfset2$category)
    dfset2$category <- sub("hbadj", "", dfset2$category)
    dfset2$category <- sub("hbnad", "", dfset2$category)
    print ('informations')

    print (unique(dfset2$category))
    gg <- ggplot()
    gg <- gg + geom_line(data=dfset2[!is.na(dfset2$value),], aes(x=displacement,y=value,color=category),size=0.4, alpha=0.8) 
    gg <- gg + geom_line(data=dfForce[!is.na(dfForce$value),], aes(x=displacement, y=(value / 1000)*coeff), size=0.1, alpha=0.8)
    gg <- gg + scale_color_manual(values=c(c1,c11,c2,c21,c3,c31,c1,c11,c2,c21,c3,c31,'red'))
    gg <- gg + scale_y_continuous(name='Hb_count',sec.axis = sec_axis(~./coeff, name="Force (nN)")) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + facet_wrap(~polymer,ncol=1) 
    gg <- gg + theme (strip.text.x = element_text(size = 7),axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=0.2)))
    ggsave(filename = paste0(outapp,'biased1.pdf'),gg, height = 10, width= 5)
    }
    print ('5')
plotcombined <- function(df,extractele,outapp) {
    print ('0')
    df$keeps <- ifelse (grepl(extractele,df$variable),'keep',
    ifelse (grepl("for-dispav",df$variable),'keep','discard'))
    df$hbcolor <- ifelse(grepl('adjmcmc',df$variable), 'adj_mcmc',
                   ifelse(grepl('adjscsc',df$category), 'adj_scsc',
                ifelse(grepl('adjmcsc',df$category), 'adj_mcsc',
                ifelse(grepl('allmcmc',df$variable), 'all_mcmc',
                   ifelse(grepl('allscsc',df$category), 'all_scsc',
                   ifelse(grepl('allmcsc',df$category), 'all_mcsc',
                  ifelse(grepl('nadmcmc',df$variable), 'nonad_mcmc',
                   ifelse(grepl('nadscsc',df$category), 'nonad_scsc',
                   ifelse(grepl('nadmcsc',df$category), 'nonad_mcsc', 'force')))))))))
            
    dfset=df[(df$keeps=='keep'),]
    dfForce=dfset[(dfset$variable=='for-dispav'),]
    #adding this new line here, to get a smoother control in different dataframe
    dfset=dfset[!(dfset$variable=='for-dispav'),]
    dfset=dfset[(dfset$avtype=='dispav'),]
    dfset$hbnature<- ifelse(grepl('p1',dfset$variable), 'partial',
                   ifelse(grepl('p2',dfset$variable), 'partial',   'unbiased'))
    dfset2=dfset[(dfset$hbnature=='partial'),]
    dfset1=dfset[(dfset$hbnature=='unbiased'),]
    
    #now here i do have finer control, upper and lowwer limits starts with 0 and 12,
    #F in nN ranges from the 0 to 3, so divide original scale by the three and divide pN/1000 to get nN
    #Moreover, now lets start specifying new axis, 
    maxval=max(dfset1$value,na.rm = TRUE)
    coeff=round(maxval / 3,digits=2)
    gg <- ggplot(data=dfset1[!is.na(dfset1$value),], aes(x=displacement,y=value,color=hbcolor))
    gg <- gg + geom_line(size=0.4,alpha=0.8 )
    gg <- gg + geom_line(data=dfForce[!is.na(dfForce$value),], aes(x=displacement, y=(value / 1000)*coeff), size=0.1, alpha=1)
    gg <- gg + scale_color_manual(values=c("black","darkgreen","blue","cyan","orange","pink","gray","red","lightblue","green"))
    gg <- gg + scale_y_continuous(name='Hb_count',sec.axis = sec_axis(~./coeff, name="Force (nN)")) 
    gg <- gg + scale_x_continuous(name="Displacement (Ang)",breaks=seq(-2,25,4))
    gg <- gg + facet_wrap(~polymer,ncol=1)
    gg <- gg + theme (strip.text.x = element_text(size = 7),axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
    gg <- gg + guides(colour = guide_legend(override.aes = list(size=1)))
    #gg <- gg + ylab("hbcount")
    ggsave(filename = paste0(outapp,'unbiased.pdf'),gg, height = 10, width= 5)
    }

dfgen=dfgen[!(dfgen$polymer=='polyala_constrained') & !(dfgen$polymer=='polygly_constrained') & !(dfgen$polymer=='polyisoleucine_server191') & !(dfgen$polymer=='polyisoleucine_78') & !(dfgen$polymer=='polyisoleucine_turing'),]
dfgen=dfgen[!(dfgen$category=='start') & !(dfgen$category=='end'),]
#removing the extra polymers 
#and start end remarks
#print ('6')
plotcombined(dfgen,'hb',paste0(out,'hbcombined_vericalStack'))
#above function will only deal with the hbond categories (mc sc and combination) other than partial views,
#it will plot both the ajdcanet and non adjacnet and all with the running average in form of diapav an donly 1 plot
plotallad(dfgen,'hbadj',paste0(out,'hbAdj_verticalStack'))
plotallad(dfgen,'hball',paste0(out,'hbAll_verticalStack'))
plotallad(dfgen,'hbnad',paste0(out,'hbNad_verticalStack'))
#above three statements will account for the biased view having stochastic stable fractions with segregations into MCSC tags but only for adjacent Hbonds
# along with biased voew, it will alos presnet the unbiased/unsegregated fractions and will term them as MCMC, SCSC ,MCSC for whole plot 
# above three will have fraction with dsplacemnt not frames

#Rscript vsdisp_hbondcount.R ../menton_set/ ../plots/new4sep/
