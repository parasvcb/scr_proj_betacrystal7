
# This program does the automatic binning of the input data of default ranges (occ of all hb's)
# and name them occurence (without defined bins tag), 
# input filetyep other than custom, does have default segregation and hence will involve, whole and subtype frequency modules,
args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Please supply two arguements systemDir and output_append_withoutext", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  systemFile=args[1]
  out=args[2]
}
library(ggplot2)
library(reshape)


PRIM_func_dealbins <- function (filenameatend) {
  directories=list.dirs(path = systemFile, full.names = TRUE, recursive=F)

  SECOND_func_cropDF <- function(dfpre,polymer,dffile) {
      dfadd=file.path(dfpre,polymer,dffile)
      dataforce=read.csv(dfadd,sep = "\t",check.names = FALSE)
      #colnames(dataforce)[which(names(dataforce) == "frames-raw")] <- "frames"
      dataforce$polymer <- polymer
      #print (head(dataforce))
      return (dataforce)
  }
  
  SECOND_func_plot <- function (df,outfile,tag) {
     
      if (grepl("defined",tag)){
        #add another tag
        gg_default <- ggplot(data=df, aes(x=Bins,y=FrequencyClass,fill=hbtype,alpha=0.5))
        gg_default <- gg_default + geom_bar(stat='identity',position=position_dodge())
        ggtemp <- ggplot(data=df, aes(x=Bins,y=FrequencyTotal,fill=hbtype,alpha=0.5))
        ggtemp <- ggtemp + geom_bar(stat='identity',position=position_dodge())
        ggtemp <- ggtemp + facet_wrap(~polymer) 
        ggtemp <- ggtemp + theme (strip.text.x = element_text(size = 6),axis.text.x = element_text( hjust = 1, size = 6, angle = 45),axis.text.y = element_text( hjust = 1, size = 6), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
        #ggtemp <- ggtemp + guides(colour = guide_legend(override.aes = list(size=10)))
        ggsave(filename = paste0(outfile,'_wholeFrequency.pdf'),ggtemp,height = 4, width= 4)
      } else {
        gg_default <- ggplot(data=df, aes(x=occurence))
        gg_default <- gg_default + stat_bin(data=df,aes(y=..count../sum(..count..)), breaks = seq(-0.1,1.1, by = .1), alpha=0.5)
      }
     
      gg_default <- gg_default + facet_wrap(~polymer)  
      gg_default <- gg_default + theme (strip.text.x = element_text(size = 6),axis.text.x = element_text( hjust = 1, size = 6, angle = 45),axis.text.y = element_text( hjust = 1, size = 6), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
      # gg_default <- gg_default + guides(colour = guide_legend(override.aes = list(size=10)))
      ggsave(filename = paste0(outfile,'_subtypeFrequency.pdf'),gg_default,height = 4, width= 4)

      
  }

  dfgen_adjacent=SECOND_func_cropDF(systemFile,"polyglycine",paste0("processed/hbonds_adjacent_",filenameatend,".tsv"))
  dfgen_adjacent = dfgen_adjacent[FALSE,]

  dfgen_nonadj=SECOND_func_cropDF(systemFile,"polyglycine",paste0("processed/hbonds_nonadj_",filenameatend,".tsv"))
  dfgen_nonadj = dfgen_nonadj[FALSE,]

  dfgen_all=SECOND_func_cropDF(systemFile,"polyglycine",paste0("processed/hbonds_all_",filenameatend,".tsv"))
  dfgen_all = dfgen_all[FALSE,]

  for (i in directories) {
      print (i)
      #print (systemFile)
      if (grepl("poly[a-z].*",i)) {
          polymer=basename(i)
          #print (polymer)
          retdf_adjacent=SECOND_func_cropDF(systemFile,polymer,paste0("processed/hbonds_adjacent_",filenameatend,".tsv"))
          retdf_nonadj=SECOND_func_cropDF(systemFile,polymer,paste0("processed/hbonds_nonadj_",filenameatend,".tsv"))
          retdf_all=SECOND_func_cropDF(systemFile,polymer,paste0("processed/hbonds_all_",filenameatend,".tsv"))
          #print (head(dfgen_nonadj))
          #print (head(dfgen_adjacent))
          dfgen_adjacent=rbind(dfgen_adjacent,retdf_adjacent)
          dfgen_nonadj=rbind(dfgen_nonadj,retdf_nonadj)
          dfgen_all=rbind(dfgen_all,retdf_all)
          #this function will get the dataframes appended together
          print (i)
      }
      # print (head(dfgen_adjacent))
      # print (head(dfgen_nonadj))
      # print (head(dfgen_all))
  }
dfgen_adjacent=dfgen_adjacent[!(dfgen_adjacent$polymer=='polyala_constrained') & !(dfgen_adjacent$polymer=='polygly_constrained') & !(dfgen_adjacent$polymer=='polyisoleucine_server191') & !(dfgen_adjacent$polymer=='polyisoleucine_78') & !(dfgen_adjacent$polymer=='polyisoleucine_turing'),]
dfgen_nonadj=dfgen_nonadj[!(dfgen_nonadj$polymer=='polyala_constrained') & !(dfgen_nonadj$polymer=='polygly_constrained') & !(dfgen_nonadj$polymer=='polyisoleucine_server191') & !(dfgen_nonadj$polymer=='polyisoleucine_78') & !(dfgen_nonadj$polymer=='polyisoleucine_turing'),]
dfgen_all=dfgen_all[!(dfgen_all$polymer=='polyala_constrained') & !(dfgen_all$polymer=='polygly_constrained') & !(dfgen_all$polymer=='polyisoleucine_server191') & !(dfgen_all$polymer=='polyisoleucine_78') & !(dfgen_all$polymer=='polyisoleucine_turing'),]
SECOND_func_plot(dfgen_adjacent,paste0(out,'hbond_frequency_adjacent_',filenameatend),filenameatend)
SECOND_func_plot(dfgen_nonadj,paste0(out,'hbond_frequency_nonadj_',filenameatend),filenameatend)
SECOND_func_plot(dfgen_all,paste0(out,'hbond_frequency_all_',filenameatend),filenameatend)
}
PRIM_func_dealbins('occurence_definedbins')
PRIM_func_dealbins('occurence')



