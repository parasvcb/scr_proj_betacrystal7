import subprocess
import re
import sys
import os
import numpy as np
'''
Here i would like to change the way we had analysed H bonds, 
go for the maximal representation of the frames possible per replicate,
1 plot per replicate and hence 3 plots overall,
for each such plot, there will be two polymorphs showing bonds with either B or D chain,
each such polymer window will later have follwing hbonds colors drawn with frame,
stability 0-20, and so on till 100 for frames till peak,  (category 1, with 5 classes)
for each consecutive 3 frames if H bond is presnet in any three frames then 1, else 0
(segregate above such that Hh between 7 and 2 will have different value than later, 0 
'''

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

def getframes(has):
    a=has.values()
    return max(a)+500

def hbbin(val):
    #print (val)
    has={(-0.10,0.20):0,(0.20,0.40):0,(0.40,0.60):0,(0.60,0.80):0,(0.80,1.00):0}
    for i in has:
        if i[0]<val<=i[1]:
            return str(i[0])+'-'+str(i[1])
def hbvalfunc(framerange,totalframes,hbpos,hbneg):
    has={}
    ma=max(framerange)
    mi=min(framerange)
    '''
    for i in totalframes:
        #if i%3==0:
            if i in framerange:
                    has[i]=hbpos
            else:
                #print ('isneg')
                has[i]=hbneg
    '''
    for i in framerange:
        if i%3==0:
            if i<=ma-2 and i >=mi+2:
                if i-1 in framerange and i+1 in framerange:
                    has[i]=hbpos
                else:
                    has[i]=hbneg
            else:
                has[i]=hbpos
    return has,hbpos+0.3,hbneg-0.3

def getforce(confile,logfile,frames):
    with open (confile) as fin:
        dirx,diry,dirz=re.search(r'SMDDir.*',fin.read()).group().split()[-3:]
    with open (logfile) as fin:
        data=[i.split()[2:] for i in fin.read().split("\n") if len(i)>0 and i[:3]=="SMD" and len(i.split())==8 and int(i.split()[1])%500==0]

    distance=[]
    force=[]
    for i in data:
            #print (i)
            distance+=[float(i[0])*float(dirx)+float(i[1])*float(diry)+float(i[2])*float(dirz)]
            force+=[float(i[3])*float(dirx)+float(i[4])*float(diry)+float(i[5])*float(dirz)]
    distane=[i-distance[0] for i in distance]
    forceav=[]
    for ind,val in enumerate(force):
        if ind>=2 and ind<len(force)-3:
            forceav+=[sum(force[ind-2:ind+3])/5]
        else:
            forceav+=[0]
    return {ind:val for ind,val in enumerate(forceav)}

def readthehbonds(has, chainsHead, reptype, poltype, framesupp, parentallowed):
    reverseiter=True if 'D' in chainsHead else False
    dirloc = os.path.join(parent, reptype, poltype,
                          'processed', 'hbonds_adjacent')
    # print(dirloc)
    dictonary = {'polyalanine_default': '1ALA', 'polyglycine': '1GLY', 'polyasparagine': '2ASN',
                 'polythreonine': '2THR', 'polyvaline': '3VAL', 'polyisoleucine': '3ILE', 'polyala-gly': '1ALAGLY'}
    val = dictonary[poltype]
    temphas = {}
    
    for i in os.listdir(dirloc):
        #print(i)
        namef = int(i.split(".")[0])
        if namef <= framesupp:
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

                chainssort = chains.copy()
                chainssort.sort(reverse=reverseiter)
                resids = [ele[0].split('-')[1][3], ele[1].split('-')[1][3]]
                sortresid=list(map(int,resids.copy()))
                sortresid.sort()
                #print(hbtype[0], chainssort, chainsHead, resids)
                #print(hbtype == ['Main', 'Main'])
                if hbtype == ['Main', 'Main'] and not ('C' in atomtype[0] or 'C' in atomtype[1]):
                    #print('in')
                    if chainssort == chainsHead:
                        if tuple(sortresid) in parentallowed:

                            #ABOVE CAN BE REMOVED FOR GREATER h BONDS COVERAGE
                            key = [chains[0]+'_'+resids[0]+'_'+atomtype[0],
                                chains[1]+'_'+ resids[1]+'_'+atomtype[1]]
                            key.sort()
                            key = tuple(key)
                            # print(key)
                            if key not in temphas:
                                temphas[key] = []
                            temphas[key] += [namef]
    specificframes=hasframes[reptype][poltype]
    if val not in has:
        has[val] = {}
    for i in temphas:
        #print(i, temphas[i])
        #print (val,i)
        if i not in has[val]:
            has[val][i] = []
        occurence= round(len([j for j in temphas[i] if j <=specificframes])/specificframes,2)
        has[val][i] = [temphas[i],occurence]

    return has



def writetofile(has, outfile,hbpos,hbneg,framecons,reptype):
    parent = "../"
    totalframes=list(range(0,framecons+1))
    hbhasval={}
    commonkeys = {}
    for aa in has:
        for hb in has[aa]:
            if hb not in commonkeys:
                commonkeys[hb]=0
            commonkeys[hb] += 1
    commonlis=[[commonkeys[i],i] for i in commonkeys]
    commonlis.sort(reverse=True)
    #print (commonlis)
    colorschema={}
    count=1
    for i in commonlis:
        colorschema[i[1]]=str(count)+'_'+i[1][0]+'__'+i[1][1]
        count+=1
    with open(outfile+'_hbdata.tsv', 'w') as fout:
        fout.write('frame\tpoltype\thbname\thstability\thbval\n')
        for aa in has:
            for hb in has[aa]:
                hbname=colorschema[hb]
                poltype=aa
                hstability=hbbin(has[aa][hb][1])
                framerange=has[aa][hb][0]
                #thisframerange only accounts for frames over which this HB was present
                framerange.sort()
                if hb not in hbhasval:
                    hbhasval[hb]=(hbpos,hbneg)
                framerangeref,hbpos,hbneg=hbvalfunc(framerange,totalframes,hbhasval[hb][0],hbhasval[hb][1])
                #print (hb,hbhasval[hb])
                for i in [j for j in framerange if j%6==0]:
                #for i in framerange:
                    hbval=framerangeref[i]
                    #force=forceref[i]
                    fout.write('%s\t%s\t%s\t%s\t%s\n' % (i,poltype, hbname[2:],hstability, hbval))
    dictionary = {'polyalanine_default': '1ALA', 'polyglycine': '1GLY', 'polyasparagine': '2ASN',
                 'polythreonine': '2THR', 'polyvaline': '3VAL', 'polyisoleucine': '3ILE', 'polyala-gly': '1ALAGLY'}
    dicnew={dictionary[i]:i for i in dictionary}
    with open (outfile+'_forceread.tsv','w') as fout:
        fout.write('frame\tpoltype\tforceavg\n')
        for aa in has:
            #print (reptype,poltype)
            forceavgref= getforce (os.path.join(parent, reptype, dicnew[aa],
                          'configurations/force.conf'), os.path.join(parent, reptype, dicnew[aa],
                          'logs/force.log'),framecons)
            for i in totalframes:
                fout.write('%s\t%s\t%s\n'%(i,aa,forceavgref[i]))

        

deleted = ['polyala_constrained', 'polygly_constrained',
           'polyisoleucine_server191', 'polyisoleucine_78', 'polyisoleucine_turing']
parent = "../"

parentresids=[(2,7),(3,6),(4,5)]

for reptype in os.listdir(parent):
    if reptype in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        hasmegaBC = {}
        hasmegaCD = {}
        theframescount=getframes(hasframes[reptype])

        for poltype in os.listdir(os.path.join(parent, reptype)):
            if re.match(r'poly[a-z]+', poltype) and poltype not in deleted:
                hasmegaBC = readthehbonds(hasmegaBC, ['B', 'C'], reptype, poltype, theframescount,parentresids)
                hasmegaCD = readthehbonds(hasmegaCD, ['D', 'C'], reptype, poltype, theframescount,parentresids)
                #break

        hbpos=1
        hbneg=-1
        writetofile(hasmegaBC, '../data_tsv/%s_frameandatbilitywiseBC'%(reptype),hbpos,hbneg,theframescount,reptype)
        hbpos=1
        hbneg=-1
        writetofile(hasmegaCD, '../data_tsv/%s_frameandatbilitywiseDC'%(reptype),hbpos,hbneg,theframescount,reptype)
        #break


'''
#rscript to plot
library(ggplot2)

plothbdata <- function(inputdataframe,outputfile,val) {
dfgen=read.csv(paste0(inputdataframe,'_hbdata.tsv'),sep = "\t",check.names = FALSE)
dfforce=read.csv(paste0(inputdataframe,'_forceread.tsv'),sep = "\t",check.names = FALSE)

gg <- ggplot(data=dfgen)
gg <- gg + geom_point(aes(x=frame,y=hbval,color=hbname,shape=hstability),size=0.8,alpha=0.8)
gg <- gg + scale_shape_manual(values=c(0,1,2,3,8))
gg <- gg + geom_line(data=dfforce,aes(x=frame, y=(forceavg / 1000), color='Force'), size=0.1, alpha=0.8)
gg <- gg + scale_y_continuous(name='Hb_values',sec.axis = sec_axis(~./1, name="Force (nN)"))
gg <- gg + scale_color_manual(values=val)
gg <- gg + facet_wrap(~poltype,ncol=1)
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
#gg <- gg + facet_wrap(~polymer,ncol=1,scale='free_x') 
ggsave(filename = outputfile, height=20, width=10)
}
val8=c("red","blue","orange","black","cyan","green","brown","black")
val9=c("red","blue","orange","black","cyan","green","brown","gray","black")
val10=c("red","blue","orange","black","cyan","green","brown","gray","pink","black")

plothbdata('../data_tsv/menton_set_frameandatbilitywiseBC','../plots/menton_set_frameandatbilitywiseBC.pdf',val9)
plothbdata('../data_tsv/tyroneP_set_frameandatbilitywiseBC','../plots/tyroneP_set_frameandatbilitywiseBC.pdf',val9)
plothbdata('../data_tsv/tyrone_set3_frameandatbilitywiseBC','../plots/tyrone_set3_frameandatbilitywiseBC.pdf',val9)
plothbdata('../data_tsv/menton_set_frameandatbilitywiseDC','../plots/menton_set_frameandatbilitywiseDC.pdf',val10)
plothbdata('../data_tsv/tyroneP_set_frameandatbilitywiseDC','../plots/tyroneP_set_frameandatbilitywiseDC.pdf',val9)
plothbdata('../data_tsv/tyrone_set3_frameandatbilitywiseDC','../plots/tyrone_set3_frameandatbilitywiseDC.pdf',val9)
'''
