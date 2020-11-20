from progress import bar as Bar
import sys
import os
import re
import subprocess
import module_processing as mdpro
import modules_computing as mdcom
import pandas as pd
import numpy as np
# this program taken in the arguements from user where the polymeric simulation trajectory has to be referred twice (they can be different also), this directory if have processed folder, then all the raw files will be added to that folder and thereafter will be used as file directory where the presence of output files will matter and raw data generation or computing from those will be continued later.
if len(sys.argv) != 3:
    print("Please enter correct cmd arguements 1:simdir 2:outdir")
    sys.exit()
program, dirsim, outdir = sys.argv
processingdir = os.path.join(dirsim, "processed")
outputfinalfile = os.path.join(outdir, "dataframe_vnew3.tsv")

outdir = processingdir
mdpro.makedir(processingdir)
del program
prottype = re.search(r'poly\w+\-?\w*', dirsim).group()
reptype = "menten" if "menton" in dirsim else "tyroneP" if "tyroneP" in dirsim else "tyroneR"
if os.path.isfile(outputfinalfile):
    sys.exit()
    #
# -> arrow below needs to be fulfilled and otther dimesnions should be recorded but only for ra20 force peaks
#
has_testing = {'menten': {'polyalanine_default': {'up': [(1.5, 2.5), (8, 10)], 'down': [(6, 8), (10, 14)]},
                          'polyala-gly': {'up': [(1.5, 3), (8, 10)], 'down': [(4, 8), (12, 14)]},
                          'polyvaline': {'up': [(2, 6), (8, 12)], 'down': [(6, 8), (12, 14)]},
                          'polygly_constrained': {'up': [(1.5, 2.5), (5, 6.5)], 'down': [(2.5, 4.5), (6.5, 9)]},
                          'polythreonine': {'up': [(1.5, 2.5), (8, 12)], 'down': [(4, 8), (12, 16)]},
                          'polyala_constrained': {'up': [(1.5, 2.5), (8, 10)], 'down': [(6, 8), (10, 14)]},
                          'polyglycine': {'up': [(1.5, 2.5), (5, 6.5)], 'down': [(2.5, 4.5), (6.5, 9)]},
                          'polyasparagine': {'up': [(1.5, 4), (8, 12)], 'down': [(6, 10), (10, 12)]},
                          'polyisoleucine': {'up': [(2, 6), (8, 10)], 'down': [(6, 10), (10, 14)]}},
               'tyroneP': {'polyalanine_default': {'up': [(1.5, 2.5), (8, 10)], 'down': [(6, 8), (10, 12)]},
                           'polyala-gly': {'up': [(1.5, 3), (8, 10)], 'down': [(4, 8), (12, 14)]},
                           'polyvaline': {'up': [(2, 6), (8, 12)], 'down': [(6, 8), (12, 14)]},
                           'polygly_constrained': {'up': [(1.5, 2.5), (5, 6.5)], 'down': [(2.5, 4.5), (6.5, 9)]},
                           'polythreonine': {'up': [(1.5, 2.5), (8, 12)], 'down': [(4, 8), (12, 16)]},
                           'polyala_constrained': {'up': [(1.5, 2.5), (8, 10)], 'down': [(6, 8), (10, 14)]},
                           'polyglycine': {'up': [(1.5, 2.5), (5, 6.5)], 'down': [(2.5, 4.5), (6.5, 9)]},
                           'polyasparagine': {'up': [(1.5, 4), (8, 12)], 'down': [(6, 10), (10, 13)]},
                           'polyisoleucine': {'up': [(2, 6), (8, 12)], 'down': [(6, 10), (10, 14)]}},
               'tyroneR': {'polyalanine_default': {'up': [(1.5, 2.5), (8, 10)], 'down': [(6, 8), (10, 14)]},
                           'polyala-gly': {'up': [(1.5, 3), (8, 10)], 'down': [(4, 8), (12, 14)]},
                           'polyvaline': {'up': [(2, 6), (8, 10)], 'down': [(6, 8), (10, 11)]},
                           'polygly_constrained': {'up': [(1.5, 2.5), (5, 6.5)], 'down': [(2.5, 4.5), (6.5, 9)]},
                           'polythreonine': {'up': [(1.5, 2.5), (8, 12)], 'down': [(6, 8), (8, 10)]},
                           'polyala_constrained': {'up': [(1.5, 2.5), (8, 10)], 'down': [(6, 8), (10, 14)]},
                           'polyglycine': {'up': [(1.5, 2.5), (5, 6.5)], 'down': [(2.5, 4.5), (6.5, 10)]},
                           'polyasparagine': {'up': [(1.5, 4), (8, 12)], 'down': [(6, 10), (10, 12)]},
                           'polyisoleucine': {'up': [(2, 8), (10, 14)], 'down': [(8, 10), (12, 14)]}}
               }
# up and downs dict labels have tuple of two values as list, where they represent the bin range from to where search the peak/descent
# calculations for only first two peaks, (should be only 1)


poltype = ''
for i in has_testing['menten']:
    if i in dirsim:
        poltype = i


for i in os.listdir(dirsim):
    # print(i)
    if os.path.isfile(os.path.join(dirsim, i)) and i != 'par_all36_prot.prm':
        source = os.path.join(dirsim, i)
        destination = os.path.join(processingdir, i)
        mdpro.movefile(source, destination)
        print("moved", source)
# general cleanup of main sim dir

centerofmasstclscript = "/home/paras/bin/centerofmassgeom.tcl"
catdcd = "$HOME/bin/catdcd"
waterpsf = os.path.join(dirsim, 'before_mini', "ionized.psf")
waterpdb = os.path.join(dirsim, 'before_mini', "ionized.pdb")
waterminipdb = os.path.join(dirsim, 'before_mini', "coordmini.pdb")

minimizedpdb = os.path.join(processingdir, "minimized_protein.pdb")
rawpdb = os.path.join(processingdir, "raw_protein.pdb")
waterdcdpull = os.path.join(dirsim, "dcd_outputs", "pull", "force_pull.dcd")
waterdcdprod = os.path.join(dirsim, "dcd_outputs", "press_concatenate.dcd")
psf = os.path.join(processingdir, "wowater.psf")
dcdpull = os.path.join(processingdir, "wowaterpulling.dcd")
dcdprod = os.path.join(processingdir, "wowaterprod.dcd")

if not os.path.isfile(dcdprod):
    # if 0:
    if not os.path.isfile(waterdcdprod):
        waterdcdprod = mdcom.concatenate_dcd(os.path.join(dirsim, "dcd_outputs"), [
                                             "press_equil1/equil_p.dcd", "press_equil2/equil_p.dcd", "press_equil3/equil_p.dcd", "press_equil4/equil_p.dcd", "press_equil5/equil_p.dcd"], waterdcdprod, catdcd)
        # -> this over here needs the changes, for the slower runs, this needs to be handled differently segregate it from non equilibrium simulations
    if waterdcdprod:
        mdcom.removewaterfromdcd(waterpsf, waterpdb, waterdcdprod,
                                 dcdprod, catdcd, processingdir)
    else:
        dcdprod = False

if not os.path.isfile(rawpdb) or not os.path.isfile(minimizedpdb):
    # if 0:
    mdcom.removewaterfrompdb(
        waterpsf, waterpdb, os.path.join(processingdir, "raw_protein"))
    mdcom.removewaterfrompdb(waterpsf, waterminipdb, os.path.join(
        processingdir, "minimized_protein"))

if not (os.path.isfile(os.path.join(processingdir, "rmsd_firstframe.txt")) and os.path.isfile(os.path.join(processingdir, "rmsd_frompdb.txt"))):
    # if 0:
    mdcom.rmsd(dcdprod, minimizedpdb, processingdir)
# sys.exit()
# usually placed for fast group calculations

if not os.path.isfile(psf):
    psfsource = os.path.join(dirsim, 'before_mini/no_clash_z_autopsf.psf')
    subprocess.check_output(['cp', psfsource, psf])

if not os.path.isfile(dcdpull):
    mdcom.removewaterfromdcd(
        waterpsf, waterpdb, waterdcdpull, dcdpull, catdcd, processingdir)


# sys.exit()
# if not os.path.isfile(os.path.join(processingdir, "dimensions.tsv")):
if 0:
    rawdim = mdcom.nanocrystal_dimensions(
        waterpsf, waterpdb, centerofmasstclscript)
    mindim = mdcom.nanocrystal_dimensions(waterpsf, os.path.join(
        dirsim, "before_mini/coordmini.pdb"), centerofmasstclscript)
    equildim = mdcom.nanocrystal_dimensions(waterpsf, os.path.join(
        dirsim, "before_mini/ready_p2.pdb"), centerofmasstclscript)

    df = pd.DataFrame([rawdim, mindim, equildim], index=[
        'raw', 'minimized', 'equilibrated']).transpose()
    df['dim_categories'] = df.index
    df.to_csv(os.path.join(processingdir, "dimensions.tsv"),
              index=False, sep='\t')

dispint = 0.1
# displacement interval over which displacemnet wise averaging of property X will be computed
distC, forceC = mdpro.forcedistance(dirsim)
forceC, force_ra20, forcedispav = mdpro.compute_averages(
    forceC, distC, dispint)
# print(forcedispav)
# sys.exit()
keystemp = list(forcedispav.keys())
keystemp.sort()
newarrtemp = [[i, forcedispav[i]] for i in keystemp]
peak1 = mdcom.peakdistances(
    newarrtemp, has_testing[reptype][poltype]['up'][0], "max")
down1 = mdcom.peakdistances(
    newarrtemp, has_testing[reptype][poltype]['down'][0], "min")
peak2 = mdcom.peakdistances(
    newarrtemp, has_testing[reptype][poltype]['up'][1], "max")
down2 = mdcom.peakdistances(
    newarrtemp, has_testing[reptype][poltype]['down'][1], "min")
# above four values have 2 values packed, displacemnet and force

with open(os.path.join(processingdir, "peaks_distances.tsv"), "w") as fin:
    fin.write("peaktype\tdistance\tforce\tpoltype\treptype\n")
    fin.write("peak1\t%s\t%s\n" % (peak1[0], peak1[1]))
    fin.write("down1\t%s\t%s\n" % (down1[0], down1[1]))
    fin.write("peak2\t%s\t%s\n" % (peak2[0], peak2[1]))
    fin.write("down2\t%s\t%s\n" % (down2[0], down2[1]))

peakhas = {peak1[0]: peak1[1], peak2[0]: peak2[1]}
downhas = {down1[0]: down1[1], down2[0]: down2[1]}
# so here we stored, displacement peak values as keys and force as their values

# if not (os.path.isfile(os.path.join(processingdir, "toughness.tsv"))):
if 0:
    volume = mdcom.nanocrystal_volume(waterpsf, waterpdb)
    # print (volume,poltype)
    toughness_data_all, toughness_data_1stpeak, toughness_data_2aa = mdpro.toughness(
        forcedispav, down1[0], volume)
    with open(os.path.join(processingdir, "toughness.tsv"), 'w') as fout:
        fout.write("Peaktype\tArea\n")
        fout.write("Complete\t%s\n" % (toughness_data_all))
        fout.write("1stpeak\t%s\n" % (toughness_data_1stpeak))
        fout.write("2aadist\t%s\n" % (toughness_data_2aa))

# sys.exit()
framestravelled = mdpro.framestravelled(distC, dispint=0.5)


# refined
# check first of data exists, if yes ,, call for processing, else for computing
# tdirhball = os.path.join(processingdir, 'hbonds_all')
cuttoff = 0.3
if poltype == "polyglycine" and 0:
    floatdigit = len(str(peak1[0]).split('.')[-1])
    # print(peak1[0], floatdigit)
    peak1[0] = round(peak1[0] - 0.5, floatdigit)
    # print(peak1[0])
# if poltype == 'polyglycine':
#     cuttoff = 0.5
#     print('in', 2)
# print(peak1[0], 'after')
fullformhbs = {'hbnad': 'hbonds_nonadj',
               'hbadj': 'hbonds_adjacent', 'hball': 'hbonds_all'}
hashbdf = {}
extractout = {}
for hbtype in ['hbadj', 'hbnad', 'hball']:
    dirtemp = os.path.join(processingdir, fullformhbs[hbtype])
    mdpro.makedir(dirtemp)
    if not len(os.listdir(dirtemp)) > 10:
        # arbitrary cuttoff
        print('inside')
        mdcom.hbonds_calculator3layer(
            psf=psf, dcd=dcdpull, outfile=dirtemp, mode=hbtype)
        # if directory files are not present
        # compute them else merge them into dataframe
    print(poltype, reptype)
    hbondData = mdpro.hbonds_calculator3layer(
        dirtemp, hbtype, distC, peak1[0], down1[0], peak2[0], dispint=dispint, cuttofffrommain=cuttoff)
    sys.exit()
    temphas = {i: hbondData[i] for i in hbondData if 'dispav' in i}
    hbondData = {i: hbondData[i] for i in hbondData if 'dispav' not in i}
    extractout.update(temphas)
    hashbdf[hbtype] = pd.DataFrame(hbondData)
hblisdf = list(hashbdf.values())
hbcompdf = hblisdf[0]
for i in hblisdf[1:]:
    hbcompdf = pd.merge(hbcompdf, i, how='outer',
                        left_index=True, right_index=True)

if not mdcom.checkFiles(['com_gcm_pull.log'], directory=processingdir):
    mdcom.com_calc(psf, dcdpull, centerofmasstclscript, processingdir)

print(poltype, reptype)
# print(peak1[0], 'before')


com1st, com1_ra20, com1_dispav, com2nd, com2_ra20, com2_dispav, com3rd, com3_ra20, com3_dispav = mdpro.centreofmasscalc(
    os.path.join(processingdir, "com_gcm_pull.log"), distC, dispint)

extractout['com3-dispav'] = com3_dispav
extractout['com2-dispav'] = com2_dispav
extractout['com1-dispav'] = com1_dispav
extractout['for-dispav'] = forcedispav

dfdispav = pd.DataFrame(extractout)
dfdispav['displacement'] = dfdispav.index

# hbcompdf needs frame column,
hbcompdf['frames-raw'] = hbcompdf.index

# dfdispav needs to be matched with displacemenet values

dfmain = pd.DataFrame([forceC, force_ra20, distC,   framestravelled,
                       com1st, com1_ra20,  com2nd, com2_ra20,  com3rd, com3_ra20],
                      index=['for-raw', 'for-ra20', 'displacement', 'frametrav-raw',
                             'com1-raw', 'com1-ra20', 'com2-raw', 'com2-ra20', 'com3-raw', 'com3-ra20', ]).transpose()


dfmain['frames-raw'] = dfmain.index

dfmain = pd.merge(dfmain, hbcompdf, on=['frames-raw'], how='outer')
dfmain = pd.merge(dfmain, dfdispav, on=['displacement'], how='outer')

dfmain = dfmain.sort_values(['frames-raw'])

dfmain.to_csv(outputfinalfile,
              index=False, sep='\t')

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
dfhbp2 = pd.DataFrame(peak2has)

dfhbp1['displacement'] = dfhbp1.index
dfhbp2['displacement'] = dfhbp2.index

dfhbp1.to_csv(os.path.join(outdir, "hbonds_peaks1.tsv"),
              index=False, sep='\t')
dfhbp2.to_csv(os.path.join(outdir, "hbonds_peaks2.tsv"),
              index=False, sep='\t')
