import sys
import os
import re
import subprocess
import module_processing as mdpro
import modules_computing as mdcom
import pandas as pd
import time
# this program taken in the arguements from user where the polymeric simulation trajectory has to be referred twice (they can be different also), this directory if have processed folder, then all the raw files will be added to that folder and thereafter will be used as file directory where the presence of output files will matter and raw data generation or computing from those will be continued later.

inst1=time.time()
if len(sys.argv) != 4:
    print("Please enter correct cmd arguements \
        1:simdir 2:outdir 3:analyze the production run (0 no, 1 yes)")
    sys.exit()

def timeoutput(start,end,message,exitflag=False):
    print (message,round(end-start,3))
    if exitflag:
        sys.exit()

program, dirsim, outdir, analyse_production = sys.argv
del (program)
analyse_production = int(analyse_production)

# general operations

processingdir = os.path.join(dirsim, "processed")
outputfinalfile = os.path.join(processingdir, "dataframe_vnew3.tsv")
outputfinalfileOnlyForce = os.path.join(processingdir, "dataframe_Force.tsv")

if os.path.isfile(outputfinalfile) and 0:
    sys.exit()

outdir = processingdir
mdpro.makedir(processingdir)

for i in os.listdir(dirsim):
    source = os.path.join(dirsim, i)
    destination = os.path.join(processingdir, i)
    if os.path.isfile(source) and i != 'par_all36_prot.prm':
        mdpro.movefile(source, destination)
        print("moved", source)
# general cleanup of main sim dir


# common file variable declarataions:
catdcd = "/home/paras/bin/catdcd"
centerofmasstclscript = "/home/paras/bin/centerofmassgeom.tcl"

waterpsf = os.path.join(dirsim, 'before_mini', "ionized.psf")
waterpdb = os.path.join(dirsim, 'before_mini', "ionized.pdb")

waterminipdb = os.path.join(dirsim, 'before_mini', "coordmini.pdb")
minimizedpdb = os.path.join(processingdir, "minimized_protein.pdb")

rawpdb = os.path.join(processingdir, "raw_protein.pdb")
psf = os.path.join(processingdir, "wowater.psf")

waterdcdpull = os.path.join(dirsim, "dcd_outputs", "pull", "force_pull.dcd")
waterdcdprod = os.path.join(dirsim, "dcd_outputs", "press_concatenate.dcd")

dcdpull = os.path.join(processingdir, "wowaterpulling.dcd")
dcdprod = os.path.join(processingdir, "wowaterprod.dcd")

# waterpsf,waterminipdb,wasterpdb,waterdcdpull are must needed files
# -> check their presemce and exit the program if not there with a warning,

timeoutput(inst1,time.time(),'filecreation')

#sys.exit()
listFileEssential = [catdcd, centerofmasstclscript,
                     waterpsf, waterpdb, waterminipdb, waterdcdpull]

if not mdcom.checkFiles(listFileEssential):
    print("make sure follwing basic files were present and run again")
    for i in listFileEssential:
        print(i)
    sys.exit()

toughnessFile = os.path.join(processingdir, "toughness.tsv")
dimensionsFile = os.path.join(processingdir, "dimensions.tsv")
peakdistancesfile = os.path.join(processingdir, "peaks_distances.tsv")

timeoutput(inst1,time.time(),'Basics')

if analyse_production:
    # -> so many tweaks are pending
    # print (how many frames are usually taken for delving deep into the equilibrium run)
    if not os.path.isfile(dcdprod):
        if not os.path.isfile(waterdcdprod):
            waterdcdprod = mdcom.concatenate_dcd(os.path.join(dirsim, "dcd_outputs"), [
                "press_equil1/equil_p.dcd", "press_equil2/equil_p.dcd", "press_equil3/equil_p.dcd", "press_equil4/equil_p.dcd", "press_equil5/equil_p.dcd"], waterdcdprod, catdcd)
            # -> this over here needs the changes, for the slower runs, this needs to be handled differently segregate it from non equilibrium simulations
        if waterdcdprod:
            mdcom.removewaterfromdcd(waterpsf, waterpdb, waterdcdprod,
                                     dcdprod, catdcd, processingdir)
        else:
            dcdprod = False
    if not (os.path.isfile(os.path.join(processingdir, "rmsd_firstframe.txt")) and os.path.isfile(os.path.join(processingdir, "rmsd_frompdb.txt"))):
        mdcom.rmsd(dcdprod, minimizedpdb, processingdir)

# -> arrow below needs to be fulfilled and otther dimesnions should be recorded but only for ra20 force peaks
#
has_fast={
    'polyalanine_default':
            {
                'menten': {'up': [(1.5, 3), (8, 10)],'down': [(6, 8), (10, 14)]},
                'tyroneP': {'up': [(1.5, 3), (8, 10)], 'down': [(6, 8), (10, 13)]},
                'tyroneR': {'up': [(1.5, 3), (8, 10)], 'down': [(6, 8), (10, 14)]}
            },
    'polyala-gly': 
            {
                'menten': {'up': [(1.5, 3), (8, 10)],'down': [(5, 8), (11, 14)]},
                'tyroneP': {'up': [(1.5, 3), (8, 10)], 'down': [(5, 8), (11, 14)]},
                'tyroneR': {'up': [(1.5, 3), (8, 10)], 'down': [(5, 8), (11, 14)]}
            },
    'polyglycine':
            {
                'menten': {'up': [(1.5, 3), (5, 7)],'down': [(2, 6), (6.5, 9)]},
                'tyroneP': {'up': [(1.5, 3), (8, 10)], 'down': [(4, 8), (10, 12)]},
                'tyroneR': {'up': [(1.5, 3), (5, 6.5)], 'down': [(2.5, 5), (6, 8)]}
            },
    'polyasparagine': 
            {
                'menten': {'up': [(1.5, 4), (8, 12)],'down': [(6, 10), (10, 13)]}, #second peak is bit confusing
                'tyroneP': {'up': [(1.5, 4), (8, 12)], 'down': [(6, 10), (10, 13)]},
                'tyroneR': {'up': [(1.5, 4), (8, 12)], 'down': [(6, 10), (10, 13)]}
            },
    'polythreonine': 
            {
                'menten': {'up': [(1.5, 3), (8, 12)],'down': [(4, 8), (12, 16)]},
                'tyroneP': {'up': [(1.5, 3), (7, 12)], 'down': [(4, 8), (12, 16)]},
                'tyroneR': {'up': [(1.5, 3), (6, 10)], 'down': [(4, 8), (10, 12)]}
            },
    'polyisoleucine': 
            {   'menten': {'up': [(2, 6), (8, 10)],'down': [(6, 10), (10, 14)]},
                'tyroneP': {'up': [(2, 6), (8, 12)], 'down': [(6, 10), (10, 14)]},
                'tyroneR': {'up': [(2, 8), (10, 14)], 'down': [(8, 10), (12, 16)]}
            },
    'polyvaline': 
        {
            'menten': {'up': [(2, 6), (8, 12)],'down': [(6, 8), (12, 14)]},
            'tyroneP': {'up': [(2, 6), (8, 12)], 'down': [(6, 8), (12, 14)]},
            'tyroneR': {'up': [(2, 6), (8, 10)], 'down': [(6, 8), (8, 11)]}
        }
}
has_slow={
    'polyalanine_default':
            {
                'menten': {'up': [(0, 4), (7, 10)],'down': [(5, 8), (12, 14)]},
                'tyroneP': {'up': [(0, 4), (7, 10)],'down': [(5, 8), (12, 14)]},
                'tyroneR':  {'up': [(0, 4), (7, 10)],'down': [(5, 8), (10, 12)]}
            },
    'polyala-gly': 
            {
                'menten': {'up': [(0, 4), (8, 10)],'down': [(5, 7), (12, 14)]},  #pending
                'tyroneP': {'up': [(0, 4), (8, 10)],'down': [(5, 8), (12, 14)]},
                'tyroneR': {'up': [(0, 4), (8, 10)],'down': [(5, 8), (11, 14)]}
            },
    'polyglycine':
            {
                'menten': {'up': [(0, 4), (4, 6)],'down': [(3, 5), (6,9)]}, #pending
                'tyroneP':{'up': [(0, 4), (4, 6)],'down': [(3, 5), (6, 9)]},
                'tyroneR': {'up': [(0, 4), (4, 6)],'down': [(3, 5), (6, 9)]}
            },
    'polythreonine': 
            {
                'menten': {'up': [(0, 2), (8, 10)],'down': [(4, 6), (10, 11)]},
                'tyroneP': {'up': [(0, 2), (7, 10)],'down': [(5, 7), (11, 13)]},
                'tyroneR': {'up': [(0, 4), (6, 10)], 'down': [(4, 8), (10, 12)]}
            },
  
    'polyasparagine': 
            {
                'menten': {'up': [(2, 4), (6, 8)],'down': [(4, 7), (8, 10)]}, #done
                'tyroneP': {'up': [(2, 5), (6, 8)], 'down': [(4, 6), (7, 10)]},
                'tyroneR': {'up': [(0, 4), (8, 12)], 'down': [(6, 9), (10, 12)]}
            },
    'polyisoleucine': 
            {   'menten': {'up': [(2, 6), (8, 12)],'down': [(6, 10), (14, 15)]},
                'tyroneP': {'up': [(2, 6), (8, 12)],'down': [(5, 8), (11, 13)]},
                'tyroneR': {'up': [(2, 6), (8, 12)],'down': [(6, 10), (11, 13)]}# peak 1 and 2 are confusing
            },
    'polyvaline': 
            {
                'menten': {'up': [(2, 6), (8, 12)],'down': [(6, 10), (12, 12.5)]},
                'tyroneP': {'up': [(2, 6), (8, 12)],'down': [(5, 8), (11, 13)]},
                'tyroneR': {'up': [(2, 6), (8, 10)],'down': [(4, 8), (10, 12)]}
            }
}
has_coordinate_range= has_fast if 'fast' in os.path.abspath(dirsim) else has_slow
# up and downs dict labels have tuple of two values as list, \
# where they represent the bin range from to where search the peak/descent
# calc. for only first two peaks, (only 1 is considered later)\
#  [Based on visual depictions]
poltype = ''
for i in has_coordinate_range:
    if i in dirsim:
        poltype = i
reptype = "menten" if "menton" in os.path.abspath(dirsim) else "tyroneP"\
    if "tyroneP" in os.path.abspath(dirsim) else "tyroneR"
print (reptype,poltype)

if not mdcom.checkFiles([rawpdb]) or not mdcom.checkFiles([minimizedpdb]):
    # both are without water versions
    # if 0:
    mdcom.removewaterfrompdb(
        waterpsf, waterpdb, os.path.join(processingdir, "raw_protein"))
    mdcom.removewaterfrompdb(waterpsf, waterminipdb, os.path.join(
        processingdir, "minimized_protein"))

if not mdcom.checkFiles([psf]):
    psfsource = os.path.join(dirsim, 'before_mini/no_clash_z_autopsf.psf')
    subprocess.check_output(['cp', psfsource, psf])

if not mdcom.checkFiles([dcdpull]):
    mdcom.removewaterfromdcd(
        waterpsf, waterpdb, waterdcdpull, dcdpull, catdcd, processingdir)

if not mdcom.checkFiles([dimensionsFile]):
    rawdim = mdcom.nanocrystal_dimensions(
        waterpsf, waterpdb, centerofmasstclscript)
    mindim = mdcom.nanocrystal_dimensions(waterpsf, os.path.join(
        dirsim, "before_mini/coordmini.pdb"), centerofmasstclscript)
    equildim = mdcom.nanocrystal_dimensions(waterpsf, os.path.join(
        dirsim, "before_mini/ready_p2.pdb"), centerofmasstclscript)

    df = pd.DataFrame([rawdim, mindim, equildim], index=[
        'raw', 'minimized', 'equilibrated']).transpose()
    df['dim_categories'] = df.index
    df.to_csv(dimensionsFile, index=False, sep='\t')

timeoutput(inst1,time.time(),'timebeforeComputing')

# main variables computation (1st instance)
dispint = 0.1
# displacement interval over which displacement wise averaging of property X will be computed
distC, forceC = mdpro.forcedistance(dirsim)
# distance and force readings were taken, and commented in module
forceC, force_ra20, forcedispav = mdpro.compute_averages(
    forceC, distC, dispint)
# made a crucial refinement in ra20 averaging, instead of adding 0 in variable keys
# anyhow, even in previous calculations, no such bias was added as emphasis was put on the
# displacement wise averaged plot

timeoutput(inst1,time.time(),'Force_calculation()')

#sys.exit()
keystemp = list(forcedispav.keys())
keystemp.sort()
newarrtemp = [[i, forcedispav[i]] for i in keystemp]
#print(newarrtemp)
#print(has_coordinate_range[poltype][reptype]['up'][1])
#print(has_coordinate_range[poltype][reptype]['down'][1])
# each element is a sublist having displacement bins averaged variable (force) values.
peak1 = mdcom.peakdistances(
    newarrtemp, has_coordinate_range[poltype][reptype]['up'][0], "max")
down1 = mdcom.peakdistances(
    newarrtemp, has_coordinate_range[poltype][reptype]['down'][0], "min")
peak2 = mdcom.peakdistances(
    newarrtemp, has_coordinate_range[poltype][reptype]['up'][1], "max")
down2 = mdcom.peakdistances(
    newarrtemp, has_coordinate_range[poltype][reptype]['down'][1], "min")
print (peak1,down1,peak2,down2)
# above four values have 2 values packed, displacemnet and force


#if not mdcom.checkFiles([peakdistancesfile]):
if 1:
    with open(peakdistancesfile, "w") as fin:
        fin.write("peaktype\tdistance\tforce\tpoltype\treptype\n")
        fin.write("peak1\t%s\t%s\n" % (peak1[0], peak1[1]))
        fin.write("down1\t%s\t%s\n" % (down1[0], down1[1]))
        fin.write("peak2\t%s\t%s\n" % (peak2[0], peak2[1]))
        fin.write("down2\t%s\t%s\n" % (down2[0], down2[1]))

peakhas = {peak1[0]: peak1[1], peak2[0]: peak2[1]}
downhas = {down1[0]: down1[1], down2[0]: down2[1]}
# so here we stored, displacement peak values as keys and force as their values

timeoutput(inst1,time.time(),'peakwriting')

#if not (mdcom.checkFiles([toughnessFile])):
if 1:
    volume = mdcom.nanocrystal_volume(waterpsf, waterpdb)
    # print (volume,poltype)
    toughness_data_all, toughness_data_1stpeak, \
        toughness_data_2aa = mdpro.toughness(forcedispav, down1[0], volume)
    with open(toughnessFile, 'w') as fout:
        fout.write("Peaktype\tArea\n")
        fout.write("Complete\t%s\n" % (toughness_data_all))
        fout.write("1stpeak\t%s\n" % (toughness_data_1stpeak))
        fout.write("2aadist\t%s\n" % (toughness_data_2aa))

timeoutput(inst1,time.time(),'toughness_data')

# refined
# check first of data exists, if yes ,, call for processing, else for computing
# tdirhball = os.path.join(processingdir, 'hbonds_all')

cuttoff_hbstable = 0.3
if poltype == "polyglycine" and 0:
    floatdigit = len(str(peak1[0]).split('.')[-1])
    # print(peak1[0], floatdigit)
    peak1[0] = round(peak1[0] - 0.5, floatdigit)
    # print(peak1[0])
# if poltype == 'polyglycine':
#     cuttoff_hbstable = 0.5
#     print('in', 2)
# print(peak1[0], 'after')
fullformhbs = {'hbnad': 'hbonds_nonadj',
               'hbadj': 'hbonds_adjacent', 'hball': 'hbonds_all'}

hashbdf = {}
extractout = {}

# main variables computation (2nd instance)
countHb = 1


temp_peak1 = peak1[0]
temp_down1 = down1[0]
temp_peak2 = peak2[0]

framesrange_first_ascent = [i for i in distC if distC[i] <= temp_peak1+0.1]
framesrange_second_ascent = [i for i in distC if temp_down1-0.1 <= distC[i] <= temp_peak2+0.1]
framesrange_slip1 =  [i for i in distC if temp_peak1-0.1 <= distC[i] <= temp_down1+0.1]
if 1:
    #writing basic frames information once
    fileframerange = os.path.join(processingdir, 'frame_ranges.txt')
    with open (fileframerange,'w') as fout:
        fout.write("Till first peak, Totalframes:%s (%s to %s)"%(len(framesrange_first_ascent),min(framesrange_first_ascent),max(framesrange_first_ascent)))
        fout.write("From 1st to onset of 2nd, Totalframes:%s (%s to %s)"%(len(framesrange_slip1),min(framesrange_slip1),max(framesrange_slip1)))
        fout.write("from onset of 2nd to its peak, Totalframes:%s (%s to %s)"%(len(framesrange_second_ascent),min(framesrange_second_ascent),max(framesrange_second_ascent)))


for hbtype in ['hbadj', 'hbnad', 'hball']:
    inst1=time.time()
    dirtemp = os.path.join(processingdir, fullformhbs[hbtype])
    mdpro.makedir(dirtemp)
    if not len(os.listdir(dirtemp)) > 10:
        # arbitrary cuttoff_hbstable
        mdcom.hbonds_calculator3layer(
            psf=psf, dcd=dcdpull, outfile=dirtemp, mode=hbtype)
        # if directory files are not present
        # compute them else merge them into dataframe
    print("***",poltype, reptype)
    hbondData = mdpro.hbonds_calculator3layer(
        dirtemp, hbtype, distC, peak1[0], down1[0], peak2[0], dispint, framesrange_first_ascent,framesrange_second_ascent, cuttofffrommain=cuttoff_hbstable
)
    temphas = {i: hbondData[i] for i in hbondData if 'dispav' in i}
    hbondData = {i: hbondData[i] for i in hbondData if 'dispav' not in i}
    extractout.update(temphas)
    hashbdf[hbtype] = pd.DataFrame(hbondData)
    timeoutput ( inst1, time.time(), 'Hbdata_computing %s'%hbtype)

inst1=time.time()

hblisdf = list(hashbdf.values())
hbcompdf = hblisdf[0]
for i in hblisdf[1:]:
    hbcompdf = pd.merge(hbcompdf, i, how='outer',
                        left_index=True, right_index=True)

timeoutput ( inst1, time.time(), 'pdmerging hbtype')
inst1=time.time()

print(poltype, reptype)


extractout['for-dispav'] = forcedispav
dfdispav = pd.DataFrame(extractout)
dfdispav['displacement'] = dfdispav.index
# hbcompdf needs frame column,
hbcompdf['frames-raw'] = hbcompdf.index

# dfdispav needs to be matched with displacemenet values

dfmain = pd.DataFrame([forceC, force_ra20, distC],
                      index=['for-raw', 'for-ra20', 'displacement']).transpose()


dfmain['frames-raw'] = dfmain.index
dfmain.to_csv(outputfinalfileOnlyForce,
              index=False, sep='\t')


dfmain = pd.merge(dfmain, hbcompdf, on=['frames-raw'], how='outer')
dfmain = pd.merge(dfmain, dfdispav, on=['displacement'], how='outer')

dfmain = dfmain.sort_values(['frames-raw'])

dfmain.to_csv(outputfinalfile,
              index=False, sep='\t')
timeoutput ( inst1, time.time(), 'writing file to dffinal')
inst1=time.time()

peak1has = {}
peak2has = {}

for keys in extractout:
    if 'hb' in keys:
        # print(keys)
        newkey = keys.split('-')[0]
        # newkey = keys
        if 'p2' in keys:
            peak2has[newkey] = {}
            peak2has[newkey][peak2[0]] = extractout[keys][peak2[0]]
            continue
        if 'p1' in keys:
            peak1has[newkey] = {}
            peak1has[newkey][peak1[0]] = extractout[keys][peak1[0]]
            continue
        peak2has[newkey] = {}
        peak2has[newkey][peak2[0]] = extractout[keys][peak2[0]]
        peak1has[newkey] = {}
        peak1has[newkey][peak1[0]] = extractout[keys][peak1[0]]
# print(peak1has)

dfhbp1 = pd.DataFrame(peak1has)
dfhbp1['displacement'] = dfhbp1.index
dfhbp1.to_csv(os.path.join(outdir, "hbonds_peaks1.tsv"),
              index=False, sep='\t')
timeoutput ( inst1, time.time(), 'hbpeaks and dome')
inst1=time.time()