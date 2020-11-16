import os
import re
import sys
import pandas as pd
import numpy as np

if len(sys.argv) != 4:
    print("Please enter correct cmd options 1forcefilehavingpeaks 2directroy having three sets 3outputfileappend")
    sys.exit()
prog, sourcefile, parent, outfile = sys.argv

def get_frames(pandasdf,p1,d1):
    df=pd.read_csv(pandasdf,sep='\t')
    disp=df['displacement']
    #print ('d1',d1)
    #print ('p1',p1)
    datastream = [i for i in disp if p1 <= i <= d1]
    return len(datastream)

def frame_reads(has_forces, has, polymer,reptype):
    f1=has_forces[polymer]['peak1'][reptype]
    d1=has_forces[polymer]['down1'][reptype]
    framescount=get_frames(os.path.join(parent,reptype,polymer,'processed','dataframe_vnew.tsv'),f1,d1)
    if polymer not in has[reptype]:
        has[reptype][polymer]=framescount
    return has

has_forces_vals={}
with open (sourcefile) as fin:
    dat=[i for i in fin.read().split("\n")[1:] if len(i)>0] 
    for i in dat:
        ele=i.split('\t')
        key = ele[0]
        if key == 'peak1' or key =='down1':
            polymer=ele[1]
            if polymer not in ['polyisoleucine_server191','polyisoleucine_78','polyisoleucine_turing']:
                #print (ele)
                menton,tyrone3,tyroneP=list(map(float,ele[5:8]))
                #print (menton)
                if polymer not in has_forces_vals:
                    has_forces_vals[polymer]={'down1':{},'peak1':{}}
                has_forces_vals[polymer][key]={'menton_set':menton,'tyroneP_set':tyroneP,'tyrone_set3':tyrone3}

has_frames={'menton_set':{},'tyroneP_set':{},'tyrone_set3':{}}
for reptype in os.listdir(parent):
    if reptype in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        for polymer in os.listdir(os.path.join(parent, reptype)):
            if re.match(r'poly[a-z]+', polymer):
                if polymer not in ['polyisoleucine_server191','polyisoleucine_78','polyisoleucine_turing']:
                    #print (polymer)
                    has_frames = frame_reads(has_forces_vals, has_frames, polymer,reptype)

with open(outfile+"_framesfrom_P1toD1.tsv", 'w') as fin:
    fin.write("polymer\tmenton_set\ttyroneP_set\ttyrone_set3\n")
    categories = ['raw', 'minimized', 'equilibrated']
    polymer = ['polyala-gly', 'polyalanine_default', 'polyglycine',
               'polythreonine', 'polyasparagine', 'polyisoleucine', 'polyvaline']
    for pol in polymer:
        val1 = np.mean(has_frames['menton_set'][pol])
        val2 = np.mean(has_frames['tyroneP_set'][pol])
        val3 = np.mean(has_frames['tyrone_set3'][pol])
        fin.write("%s\t%s\t%s\t%s\n" % (pol, val1,val2,val3))
