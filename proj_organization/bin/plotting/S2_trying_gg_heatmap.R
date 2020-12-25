short 
   penetration scc          pi0
1            0   0  0.002545268
2            5   0 -0.408621176
3           10   0 -0.929432006
4           15   0 -1.121309680
5           20   0 -1.587298317
6           25   0 -2.957853131
7           30   0 -5.123329738
8            0  50  1.199748327
9            5  50  0.788581883
10          10  50  0.267771053
11          15  50  0.075893379
12          20  50 -0.390095258
13          25  50 -1.760650073
14          30  50 -3.926126679
15           0 100  2.396951386
16           5 100  1.985784941
17          10 100  1.464974112
18          15 100  1.273096438
19          20 100  0.807107801
20          25 100 -0.563447014
21          30 100 -2.728923621

frame   poltype hbname  hstability      hbval
0       3VAL    SuperRing_1     0.8-1.0 1
10      3VAL    SuperRing_1     0.8-1.0 1
20      3VAL    SuperRing_1     0.8-1.0 1
30      3VAL    SuperRing_1     0.8-1.0 1
40      3VAL    SuperRing_1     0.8-1.0 1


mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")

ggplot(data = df, aes(x = frame, y = hbname)) +
  geom_tile(aes(fill = pi0)) +
  scale_fill_gradientn(colours = mycol)


library(ggplot2)
outputfile= ('../../plots')
inputdataframe = ('../../results/data_tsv_fast/hbdata/menton_set_frameWise_SuperRing.tsv')
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
dfgen$ring <- ifelse(dfgen$hbname == 'SuperRing_1', '1',
                  ifelse(dfgen$hbname == 'SuperRing_2', '2','3'))
dfgen$presence <- ifelse(dfgen$hbval > 0, 1,0)

head(dfgen)
gg <- ggplot(dfgen , aes(x = frame, y = hbname))
#testheatmap1
 #gg <- gg + geom_raster(aes(fill = presence), interpolate=F) 
 #gg <- gg + scale_fill_gradient2(low="white", mid="gray", high="black", midpoint=0.5, limits=range(dfgen$hbval))
#testheatmap2
gg <- gg + geom_tile(aes(fill = presence, height=0.9))
##gg <- gg + scale_fill_gradient1('hbval', low = "white", mid = "gray", high = "black", midpoint = 0.5)
gg <- gg + scale_fill_gradient('presence', low = "white",high = "black")
gg <- gg + facet_wrap(~poltype,ncol=1)
#gg <- gg + facet_grid (~ring)
ggsave(filename = file.path(outputfile,'testheatmap1.pdf'))

#gg <- gg + geom_tile(aes(fill = hbval), interpolate=T)

