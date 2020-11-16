import subprocess
import re
import sys
import os
import numpy as np

hasframes = {
    'menton_set':
    {
        'polyala-gly': 2625,
        'polyalanine_default': 2757,
        'polyasparagine': 2417,
        'polyglycine': 2533,
        'polyisoleucine': 2604,
        'polythreonine': 1977,
        'polyvaline': 2573
    },
    'tyroneP_set':
    {
        'polyala-gly': 2486,
        'polyalanine_default': 2756,
        'polyasparagine': 2539,
        'polyglycine': 2327,
        'polyisoleucine': 2491,
        'polythreonine': 1801,
        'polyvaline': 2116
    },
    'tyrone_set3':
    {
        'polyala-gly': 2800,
        'polyalanine_default': 2782,
        'polyasparagine': 2213,
        'polyglycine': 1897,
        'polyisoleucine': 3104,
        'polythreonine': 1923,
        'polyvaline': 2336
    }
}


def readthehbonds(has, chainsHead, reptype, poltype):
    dirloc = os.path.join(parent, reptype, poltype,
                          'processed', 'hbonds_adjacent')
    # print(dirloc)
    dictonary = {'polyalanine_default': '1ALA', 'polyglycine': '1GLY', 'polyasparagine': '2ASN',
                 'polythreonine': '2THR', 'polyvaline': '3VAL', 'polyisoleucine': '3ILE', 'polyala-gly': '1ALAGLY'}
    val = dictonary[poltype]
    if not val:
        print(val)
        sys.exit()
    temphas = {}
    for i in os.listdir(dirloc):
        # print(i)
        namef = int(i.split(".")[0])
        if namef < hasframes[reptype][poltype]:
            with open(os.path.join(dirloc, i)) as fin:
                dat = [k for k in fin.read().split("\n")[2:] if len(k) > 0]
            for line in dat:
                # print(line)
                line = re.sub(r'7-Side-OT1', '7-Main-OT1', line)
                line = re.sub(r'7-Side-OT2', '7-Main-OT2', line)
                ele = line.split()
                chains = [ele[0][3], ele[1][3]]
                atomtype = [ele[0].split('-')[3], ele[1].split('-')[3]]
                hbtype = [ele[0].split('-')[2], ele[1].split('-')[2]]

                chainssort = chains
                chainssort.sort()
                resids = [ele[0].split('-')[1][3], ele[1].split('-')[1][3]]
                #print(hbtype[0], chainssort, chainsHead, resids)
                #print(hbtype == ['Main', 'Main'])
                if hbtype == ['Main', 'Main']:
                    # print('in')
                    if chainssort == chainsHead:
                        key = [resids[0]+'_'+chains[0]+'_'+atomtype[0],
                               resids[1]+'_'+chains[1]+'_'+atomtype[1]]
                        key.sort()
                        key = tuple(key)
                        # print(key)
                        if key not in temphas:
                            temphas[key] = 0
                        temphas[key] += 1
    if val not in has:
        has[val] = {}
    for i in temphas:
        #print(i, temphas[i])
        if i not in has[val]:
            has[val][i] = []
        #has[val][i] += [round(temphas[i]/hasframes[reptype][poltype], 2)]
        has[val][i] += [round(temphas[i]*0.5, 2)]
    return has


hasmegaBC = {}
hasmegaCD = {}
deleted = ['polyala_constrained', 'polygly_constrained',
           'polyisoleucine_server191', 'polyisoleucine_78', 'polyisoleucine_turing']
parent = "../"
parent = "../"
for reptype in os.listdir(parent):
    if reptype in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        for poltype in os.listdir(os.path.join(parent, reptype)):
            if re.match(r'poly[a-z]+', poltype) and poltype not in deleted:
                hasmegaBC = readthehbonds(
                    hasmegaBC, ['B', 'C'], reptype, poltype)
                hasmegaCD = readthehbonds(
                    hasmegaCD, ['C', 'D'], reptype, poltype)


def writetofile(has, outfile):
    commonkeys = {}
    for aa in has:
        if aa not in commonkeys:
            commonkeys[aa] = []

        for hb in has[aa]:
            if len(has[aa][hb]) == 3:
                commonkeys[aa] += [hb]
    setcommon = set(commonkeys['1ALA'])
    for i in commonkeys:
        setcommon &= set(commonkeys[i])
    with open(outfile, 'w') as fout:
        fout.write('poltype\thbtype\tmean\tstddev\n')
        for aa in has:
            for hb in has[aa]:
                if hb in setcommon:
                    fout.write('%s\t%s\t%s\t%s\n' % (aa, hb, round(
                        np.mean(has[aa][hb]), 2), round(np.std(has[aa][hb]), 2)))


writetofile(hasmegaBC, '../data_tsv/timehbnondstillpeak1BC.tsv')
writetofile(hasmegaCD, '../data_tsv/timehbnondstillpeak1CD.tsv')


'''
#rscript to plot
library(ggplot2)

plothbdata <- function(inputdataframe,outputfile) {
dfgen=read.csv(inputdataframe,sep = "\t",check.names = FALSE)
gg <- ggplot(data=dfgen,aes(x=hbtype,y=mean,fill=poltype))
gg <- gg + geom_bar(stat="identity",position=position_dodge())
gg <- gg + geom_text(aes(label=mean),vjust=1, color="black",position = position_dodge(0.9), size=3)
gg <- gg + geom_errorbar(data= dfgen, aes(x=hbtype,y=mean,ymin=mean-stddev, ymax=mean+stddev),position = position_dodge(0.9),size=0.3,width=0.5)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
#gg <- gg + facet_wrap(~polymer,ncol=1,scale='free_x') 
ggsave(filename = outputfile, height=5, width=15)
}
plothbdata('../data_tsv/timehbnondstillpeak1BC.tsv','../plots/hbondtimechainBC.pdf')
plothbdata('../data_tsv/timehbnondstillpeak1CD.tsv','../plots/hbondtimechainCD.pdf')

'''
