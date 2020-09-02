from progress import bar as Bar
import sys
import os
import re
import subprocess
import module_processing as mdpro
import modules_computing as mdcom
import pandas as pd
if len(sys.argv) != 3:
    print("Please enter correct cmd arguements 1:simdir 2:outfile")
    sys.exit()
program, dirsim, outfile = sys.argv
processingdir = os.path.join(dirsim, "processed")
outfile = processingdir
mdpro.makedir(processingdir)
del program
prottype = re.search(r'poly\w+\-?\w*', dirsim).group()
reptype = "menten" if "menton" in dirsim else "tyroneP" if "tyroneP" in dirsim else "tyroneR"
if os.path.isfile(os.path.join(outfile, "dataframe_vnew.tsv")) and 0:
    sys.exit()
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
# up will have bins here, from whree to whree and same for down
# calculations for only first two peaks


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

print(1)
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
    if not os.path.isfile(waterdcdprod):
        waterdcdprod = mdcom.concatenate_dcd(os.path.join(dirsim, "dcd_outputs"), [
                                             "press_equil1/equil_p.dcd", "press_equil2/equil_p.dcd", "press_equil3/equil_p.dcd", "press_equil4/equil_p.dcd", "press_equil5/equil_p.dcd"], waterdcdprod, catdcd)
    if waterdcdprod:
        mdcom.removewaterfromdcd(waterpsf, waterpdb, waterdcdprod,
                                 dcdprod, catdcd, processingdir)
    else:
        dcdprod = False

if not os.path.isfile(rawpdb) or not os.path.isfile(minimizedpdb):
    mdcom.removewaterfrompdb(
        waterpsf, waterpdb, os.path.join(processingdir, "raw_protein"))
    mdcom.removewaterfrompdb(waterpsf, waterminipdb, os.path.join(
        processingdir, "minimized_protein"))

if not (os.path.isfile(os.path.join(processingdir, "rmsd_firstframe.txt")) and os.path.isfile(os.path.join(processingdir, "rmsd_frompdb.txt"))):
    mdcom.rmsd(dcdprod, minimizedpdb, processingdir)
# sys.exit()
# usually placed for fast group calculations

if not os.path.isfile(psf):
    psfsource = os.path.join(dirsim, 'before_mini/no_clash_z_autopsf.psf')
    subprocess.check_output(['cp', psfsource, psf])

if not os.path.isfile(dcdpull):
    mdcom.removewaterfromdcd(
        waterpsf, waterpdb, waterdcdpull, dcdpull, catdcd, processingdir)


print(2)
# sys.exit()
if not os.path.isfile(os.path.join(processingdir, "dimensions.tsv")):
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

if not os.path.isfile(os.path.join(processingdir, "toughness.tsv")):
    volume = mdcom.nanocrystal_volume(waterpsf, waterpdb)
    toughness_data_all, toughness_data_1stpeak, toughness_data_2aa = mdpro.toughness(
        forcedispav, down1[0], volume)
    with open(os.path.join(processingdir, "toughness.tsv"), 'w') as fout:
        fout.write("Peaktype\tArea\n")
        fout.write("Complete\t%s\n" % (toughness_data_all))
        fout.write("1stpeak\t%s\n" % (toughness_data_1stpeak))
        fout.write("2aadist\t%s\n" % (toughness_data_2aa))

framestravelled = mdpro.framestravelled(distC, dispint=0.5)

print(4)


# refined
# check first of data exists, if yes ,, call for processing, else for computing
tdirhball = os.path.join(processingdir, 'hbonds_all')
tdirhbadj = os.path.join(processingdir, 'hbonds_adjacent')
mdpro.makedir(tdirhball)
mdpro.makedir(tdirhbadj)
if not (len(os.listdir(tdirhball)) > 1000 and len(os.listdir(tdirhbadj)) > 1000):
    print('inside')
    mdcom.hbonds_calculator3layer(psf=psf, dcd=dcdpull, outfile=processingdir)

if not mdcom.checkFiles(['com_gcm_pull.log'], directory=processingdir):
    mdcom.com_calc(psf, dcdpull, centerofmasstclscript, processingdir)

hbondData = mdpro.hbonds_calculator3layer(
    processingdir, distC, peak1[0], down1[0], peak2[0], dispint=dispint)
print(hbondData.keys())
# sys.exit()


['all_mcmc_raw', 'all_mcsc_raw', 'all_scsc_raw',
 'all_mcmc_ra20', 'all_mcsc_ra20', 'all_scsc_ra20',
 'all_mcmc_dispav', 'all_mcsc_dispav', 'all_scsc_dispav',
 'adj_mcmc_raw', 'adj_mcsc_raw', 'adj_scsc_raw',
 'adj_mcmc_ra20', 'adj_mcsc_ra20', 'adj_scsc_ra20',
 'adj_mcmc_dispav', 'adj_mcsc_dispav', 'adj_scsc_dispav',
 'all_p1_mcmc_raw', 'all_p1_mcsc_raw', 'all_p1_scsc_raw',
 'all_p1_mcmc_ra20', 'all_p1_mcsc_ra20', 'all_p1_scsc_ra20',
 'all_p1_mcmc_dispav', 'all_p1_mcsc_dispav', 'all_p1_scsc_dispav',
 'adj_p1_mcmc_raw', 'adj_p1_mcsc_raw', 'adj_p1_scsc_raw',
 'adj_p1_mcmc_ra20', 'adj_p1_mcsc_ra20', 'adj_p1_scsc_ra20',
 'adj_p1_mcmc_dispav', 'adj_p1_mcsc_dispav', 'adj_p1_scsc_dispav',
 'all_p2_mcmc_raw', 'all_p2_mcsc_raw', 'all_p2_scsc_raw',
 'all_p2_mcmc_ra20', 'all_p2_mcsc_ra20', 'all_p2_scsc_ra20',
 'all_p2_mcmc_dispav', 'all_p2_mcsc_dispav', 'all_p2_scsc_dispav',
 'adj_p2_mcmc_raw', 'adj_p2_mcsc_raw', 'adj_p2_scsc_raw',
 'adj_p2_mcmc_ra20', 'adj_p2_mcsc_ra20', 'adj_p2_scsc_ra20',
 'adj_p2_mcmc_dispav', 'adj_p2_mcsc_dispav', 'adj_p2_scsc_dispav']

com1st, com1_ra20, com1_dispav, com2nd, com2_ra20, com2_dispav, com3rd, com3_ra20, com3_dispav = mdpro.centreofmasscalc(
    os.path.join(processingdir, "com_gcm_pull.log"), distC, dispint)

df = pd.DataFrame([forceC, force_ra20, distC,   framestravelled,
                   com1st, com1_ra20,  com2nd, com2_ra20,  com3rd, com3_ra20,
                   hbondData['all_mcmc_raw'], hbondData['all_mcmc_ra20'], hbondData['all_mcsc_raw'], hbondData[
                       'all_mcsc_ra20'], hbondData['all_scsc_raw'], hbondData['all_scsc_ra20'],
                   hbondData['adj_mcmc_raw'], hbondData['adj_mcmc_ra20'], hbondData['adj_mcsc_raw'], hbondData[
                       'adj_mcsc_ra20'], hbondData['adj_scsc_raw'], hbondData['adj_scsc_ra20'],
                   hbondData['all_p1_mcmc_raw'], hbondData['all_p1_mcmc_ra20'], hbondData['all_p1_mcsc_raw'], hbondData[
                       'all_p1_mcsc_ra20'], hbondData['all_p1_scsc_raw'], hbondData['all_p1_scsc_ra20'],
                   hbondData['adj_p1_mcmc_raw'], hbondData['adj_p1_mcmc_ra20'], hbondData['adj_p1_mcsc_raw'], hbondData[
                       'adj_p1_mcsc_ra20'], hbondData['adj_p1_scsc_raw'], hbondData['adj_p1_scsc_ra20'],
                   hbondData['all_p2_mcmc_raw'], hbondData['all_p2_mcmc_ra20'], hbondData['all_p2_mcsc_raw'], hbondData[
                       'all_p2_mcsc_ra20'], hbondData['all_p2_scsc_raw'], hbondData['all_p2_scsc_ra20'],
                   hbondData['adj_p2_mcmc_raw'], hbondData['adj_p2_mcmc_ra20'], hbondData['adj_p2_mcsc_raw'], hbondData[
                       'adj_p2_mcsc_ra20'], hbondData['adj_p2_scsc_raw'], hbondData['adj_p2_scsc_ra20']],
                  index=['for-raw', 'for-ra20', 'displacement', 'frametrav-raw',
                         'com1-raw', 'com1-ra20', 'com2-raw', 'com2-ra20', 'com3-raw', 'com3-ra20',
                         'all_mcmc-raw', 'all_mcmc-ra20', 'all_mcsc-raw', 'all_mcsc-ra20', 'all_scsc-raw', 'all_scsc-ra20',
                         'adj_mcmc-raw', 'adj_mcmc-ra20', 'adj_mcsc-raw', 'adj_mcsc-ra20', 'adj_scsc-raw', 'adj_scsc-ra20',
                         'all_p1_mcmc-raw', 'all_p1_mcmc-ra20', 'all_p1_mcsc-raw', 'all_p1_mcsc-ra20', 'all_p1_scsc-raw', 'all_p1_scsc-ra20',
                         'adj_p1_mcmc-raw', 'adj_p1_mcmc-ra20', 'adj_p1_mcsc-raw', 'adj_p1_mcsc-ra20', 'adj_p1_scsc-raw', 'adj_p1_scsc-ra20',
                         'all_p2_mcmc-raw', 'all_p2_mcmc-ra20', 'all_p2_mcsc-raw', 'all_p2_mcsc-ra20', 'all_p2_scsc-raw', 'all_p2_scsc-ra20',
                         'adj_p2_mcmc-raw', 'adj_p2_mcmc-ra20', 'adj_p2_mcsc-raw', 'adj_p2_mcsc-ra20', 'adj_p2_scsc-raw', 'adj_p2_scsc-ra20']).transpose()

df['frames-raw'] = df.index
df['all_mcmc-dispav'] = df['displacement'].map(hbondData['all_mcmc_dispav'])
df['all_mcsc-dispav'] = df['displacement'].map(hbondData['all_mcsc_dispav'])
df['all_scsc-dispav'] = df['displacement'].map(hbondData['all_scsc_dispav'])

df['adj_mcmc-dispav'] = df['displacement'].map(hbondData['adj_mcmc_dispav'])
df['adj_mcsc-dispav'] = df['displacement'].map(hbondData['adj_mcsc_dispav'])
df['adj_scsc-dispav'] = df['displacement'].map(hbondData['adj_scsc_dispav'])
print(hbondData['all_p1_mcmc_dispav'])
df['all_p1_mcmc-dispav'] = df['displacement'].map(
    hbondData['all_p1_mcmc_dispav'])
df['all_p1_mcsc-dispav'] = df['displacement'].map(
    hbondData['all_p1_mcsc_dispav'])
df['all_p1_scsc-dispav'] = df['displacement'].map(
    hbondData['all_p1_scsc_dispav'])

df['adj_p1_mcmc-dispav'] = df['displacement'].map(
    hbondData['adj_p1_mcmc_dispav'])
df['adj_p1_mcsc-dispav'] = df['displacement'].map(
    hbondData['adj_p1_mcsc_dispav'])
df['adj_p1_scsc-dispav'] = df['displacement'].map(
    hbondData['all_p1_scsc_dispav'])

df['all_p2_mcmc-dispav'] = df['displacement'].map(
    hbondData['all_p2_mcmc_dispav'])
df['all_p2_mcsc-dispav'] = df['displacement'].map(
    hbondData['all_p2_mcsc_dispav'])
df['all_p2_scsc-dispav'] = df['displacement'].map(
    hbondData['all_p2_scsc_dispav'])

df['adj_p2_mcmc-dispav'] = df['displacement'].map(
    hbondData['adj_p2_mcmc_dispav'])
df['adj_p2_mcsc-dispav'] = df['displacement'].map(
    hbondData['adj_p2_mcsc_dispav'])
df['adj_p2_scsc-dispav'] = df['displacement'].map(
    hbondData['all_p2_scsc_dispav'])


df['com1-dispav'] = df['displacement'].map(com1_dispav)
df['com2-dispav'] = df['displacement'].map(com2_dispav)
df['com3-dispav'] = df['displacement'].map(com1_dispav)
df['for-dispav'] = df['displacement'].map(forcedispav)

df.to_csv(os.path.join(outfile, "dataframe_vnew2.tsv"), index=False, sep='\t')
sys.exit()
has_hb_all = {peak1[0]: hb_all_dispav[peak1[0]],
              peak2[0]: hb_all_dispav[peak2[0]]}
has_hb_all_bbbb = {peak1[0]: hb_allbbbb_dispav[peak1[0]],
                   peak2[0]: hb_allbbbb_dispav[peak2[0]]}
has_hb_all_scsc = {peak1[0]: hb_allscsc_dispav[peak1[0]],
                   peak2[0]: hb_allscsc_dispav[peak2[0]]}
has_hb_all_scbb = {peak1[0]: hb_allscbb_dispav[peak1[0]],
                   peak2[0]: hb_allscbb_dispav[peak2[0]]}
has_hb_adjacent = {peak1[0]: hb_adjacent_dispav[peak1[0]],
                   peak2[0]: hb_adjacent_dispav[peak2[0]]}
has_hb_adjacent_bbbb = {peak1[0]: hb_adjacentbbbb_dispav[peak1[0]],
                        peak2[0]: hb_adjacentbbbb_dispav[peak2[0]]}
has_hb_adjacent_scsc = {peak1[0]: hb_adjacentscsc_dispav[peak1[0]],
                        peak2[0]: hb_adjacentscsc_dispav[peak2[0]]}
has_hb_adjacent_scbb = {peak1[0]: hb_adjacentscbb_dispav[peak1[0]],
                        peak2[0]: hb_adjacentscbb_dispav[peak2[0]]}
dfhb = pd.DataFrame([has_hb_all, has_hb_all_bbbb, has_hb_all_scsc, has_hb_all_scbb, has_hb_adjacent,
                     has_hb_adjacent_bbbb, has_hb_adjacent_scsc, has_hb_adjacent_scbb],
                    index=['has_hb_all', 'has_hb_all_bbbb', 'has_hb_all_scsc', 'has_hb_all_scbb', 'has_hb_adjacent',
                           'has_hb_adjacent_bbbb', 'has_hb_adjacent_scsc', 'has_hb_adjacent_scbb']).transpose()
# print(dfhb)
dfhb['displacement'] = dfhb.index
dfhb.to_csv(os.path.join(outfile, "hbonds_peaks.tsv"), index=False, sep='\t')
