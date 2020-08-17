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
               }  # up will have bins here, from whree to whree and same for down
# above needs to be reverified and change in values


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

print(1)
centerofmasstclscript = "/home/paras/bin/centerofmassgeom.tcl"
catdcd = "$HOME/bin/catdcd"
waterpsf = os.path.join(dirsim, 'before_mini', "ionized.psf")
waterpdb = os.path.join(dirsim, 'before_mini', "ionized.pdb")
waterminipdb = os.path.join(dirsim, 'before_mini', "coordmini.pdb")

minimizedpdb = os.path.join(processingdir, "minimized_protein.pdb")
rawpdb= os.path.join(processingdir, "raw_protein.pdb")
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
    mdcom.removewaterfrompdb(waterpsf,waterpdb, os.path.join(processingdir,"raw_protein"))
    mdcom.removewaterfrompdb(waterpsf,waterminipdb, os.path.join(processingdir,"minimized_protein"))

if not (os.path.isfile(os.path.join(processingdir,"rmsd_firstframe.txt")) and os.path.isfile(os.path.join(processingdir,"rmsd_frompdb.txt"))):
    mdcom.rmsd(dcdprod, minimizedpdb, processingdir)
sys.exit()

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
distC, forceC = mdpro.forcedistance(dirsim)
forcedispav = mdpro.dispAvg(distC, forceC, dispint)
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

peakhas = {peak1[0]: peak1[1], peak2[0]: peak2[1]}
downhas = {down1[0]: down1[1], down2[0]: down2[1]}

if not os.path.isfile(os.path.join(processingdir, "toughness.tsv")):
    toughness_data_all, toughness_data_1stpeak, toughness_data_2aa = mdpro.toughness(
        forcedispav, down1[0])
    with open(os.path.join(processingdir, "toughness.tsv"), 'w') as fout:
        fout.write("Peaktype\tArea\n")
        fout.write("Complete\t%s\n" % (toughness_data_all))
        fout.write("1stpeak\t%s\n" % (toughness_data_1stpeak))
        fout.write("2aadist\t%s\n" % (toughness_data_2aa))

velocityC = mdpro.velocity(distC)
velocitydispav = mdpro.dispAvg(distC, velocityC, dispint)

framestravelled = mdpro.framestravelled(distC, dispint=0.5)

print(down1)
with open(os.path.join(processingdir, "peaks_distances.tsv"), "w") as fin:
    fin.write("peaktype\tdistance\tforce\tpoltype\treptype\n")
    fin.write("peak1\t%s\t%s\n" % (peak1[0], peak1[1]))
    fin.write("down1\t%s\t%s\n" % (down1[0], down1[1]))
    fin.write("peak2\t%s\t%s\n" % (peak2[0], peak2[1]))
    fin.write("down2\t%s\t%s\n" % (down2[0], down2[1]))

print(4)
# ->forcepeaks are pending and corresponding calls

forceC, force_ra20, force_ra250 = mdpro.compute_averages(
    forceC, distC)

velocityC, velocity_ra20, velocity_ra250 = mdpro.compute_averages(
    velocityC, distC)

# refined
# check first of data exists, if yes ,, call for processing, else for computing
if not mdcom.checkFiles(['hbonds_C_protein.dat', 'hbonds_C_BD.dat', 'hbonds_C_protein_bbbb.dat',
                         'hbonds_C_BD_bbbb.dat', 'hbonds_C_protein_scsc.dat', 'hbonds_C_BD_scsc.dat',
                         'hbonds_C_protein_scbb.dat', 'hbonds_C_BD_scbb.dat'], directory=processingdir):
    mdcom.hbonds_calculator3layer(psf=psf, dcd=dcdpull, outfile=processingdir)
if not mdcom.checkFiles(['angle_resid4', 'angle_resid2'], directory=processingdir):
    print("yes")
    mdcom.angle_cal(psf, dcdpull, processingdir)
if not mdcom.checkFiles(['com_gcm_pull.log'], directory=processingdir):
    mdcom.com_calc(psf, dcdpull, centerofmasstclscript, processingdir)

hbondData = mdpro.hbonds_calculator3layer(processingdir, distC)
hb_all_raw, hb_all_ra20, hb_all_ra250 = hbondData['all']
hb_adjacent_raw, hb_adjacent_ra20, hb_adjacent_ra250 = hbondData['adjacent']
hb_allbbbb_raw, hb_allbbbb_ra20, hb_allbbbb_ra250 = hbondData['allbbbb']
hb_adjacentbbbb_raw, hb_adjacentbbbb_ra20, hb_adjacentbbbb_ra250 = hbondData[
    'adjacentbbbb']
hb_allscsc_raw, hb_allscsc_ra20, hb_allscsc_ra250 = hbondData['allscsc']
hb_adjacentscsc_raw, hb_adjacentscsc_ra20, hb_adjacentscsc_ra250 = hbondData[
    'adjacentscsc']
hb_allscbb_raw, hb_allscbb_ra20, hb_allscbb_ra250 = hbondData['allscbb']
hb_adjacentscbb_raw, hb_adjacentscbb_ra20, hb_adjacentscbb_ra250 = hbondData[
    'adjacentscbb']


pse6, pse6ra20, pse6ra250 = mdpro.angleaverages(
    os.path.join(processingdir, 'angle_resid6'), distC)
pse4, pse4ra20, pse4ra250 = mdpro.angleaverages(
    os.path.join(processingdir, 'angle_resid4'), distC)


com1st, com1_ra20, com1_ra250, com2nd, com2_ra20, com2_ra250, com3rd,\
    com3_ra20, com3_ra250 = mdpro.centreofmasscalc(
        os.path.join(processingdir, "com_gcm_pull.log"), distC)

#print("com1st", list(com1st.keys())[:10])
hb_all_dispav = mdpro.dispAvg(distC, hb_all_raw, dispint)
hb_adjacent_dispav = mdpro.dispAvg(distC, hb_adjacent_raw, dispint)
hb_allbbbb_dispav = mdpro.dispAvg(distC, hb_allbbbb_raw, dispint)
hb_allscsc_dispav = mdpro.dispAvg(distC, hb_allscsc_raw, dispint)
hb_allscbb_dispav = mdpro.dispAvg(distC, hb_allscbb_raw, dispint)
hb_adjacentbbbb_dispav = mdpro.dispAvg(
    distC, hb_adjacentbbbb_raw, dispint)
hb_adjacentscsc_dispav = mdpro.dispAvg(
    distC, hb_adjacentscsc_raw, dispint)
hb_adjacentscbb_dispav = mdpro.dispAvg(
    distC, hb_adjacentscbb_raw, dispint)
pse6_dispav = mdpro.dispAvg(distC, pse6, dispint)
pse4_dispav = mdpro.dispAvg(distC, pse4, dispint)
com1_dispav = mdpro.dispAvg(distC, com1st, dispint)
vel_dispav = mdpro.dispAvg(distC, velocityC, dispint)
# print(com1_dispav)
print("done1")
com2_dispav = mdpro.dispAvg(distC, com2nd, dispint)
print("done2")
com3_dispav = mdpro.dispAvg(distC, com3rd, dispint)
print("done3")
framecount = 8001
for i in hb_allbbbb_dispav.keys():
    if i not in distC.values():
        distC[framecount] = i
        framecount += 1

df = pd.DataFrame([forceC, force_ra20, distC, velocityC, velocity_ra20,
                   framestravelled, pse6, pse6ra20, pse4, pse4ra20, com1st, com1_ra20,
                   com2nd, com2_ra20, com3rd, com3_ra20, hb_all_raw, hb_all_ra20,
                   hb_adjacent_raw, hb_adjacent_ra20, hb_allbbbb_raw, hb_allbbbb_ra20,
                   hb_allscsc_raw, hb_allscsc_ra20, hb_allscbb_raw, hb_allscbb_ra20,
                   hb_adjacentbbbb_raw, hb_adjacentbbbb_ra20, hb_adjacentscsc_raw,
                   hb_adjacentscsc_ra20, hb_adjacentscbb_raw, hb_adjacentscbb_ra20],
                  index=['for-raw', 'for-ra20', 'displacement', 'vel-raw', 'vel-ra20',
                         'frametrav-raw', 'pse6-raw', 'pse6-ra20', 'pse4-raw', 'pse4-ra20', 'com1-raw', 'com1-ra20',
                         'com2-raw', 'com2-ra20', 'com3-raw', 'com3-ra20', 'hball-raw', 'hball-ra20',
                         'hbadj-raw', 'hbadj-ra20', 'hballbbbb-raw', 'hballbbbb-ra20',
                         'hballscsc-raw', 'hballscsc-ra20', 'hballscbb-raw', 'hballscbb-ra20',
                         'hbadjbbbb-raw', 'hbadjbbbb-ra20', 'hbadjscsc-raw',
                         'hbadjscsc-ra20', 'hbadjscbb-raw', 'hbadjscbb-ra20']).transpose()
df['frames-raw'] = df.index
df['hball-dispav'] = df['displacement'].map(hb_all_dispav)
df['hbadj-dispav'] = df['displacement'].map(hb_adjacent_dispav)
df['hballbbbb-dispav'] = df['displacement'].map(hb_allbbbb_dispav)
df['hballscsc-dispav'] = df['displacement'].map(hb_allscsc_dispav)
df['hballscbb-dispav'] = df['displacement'].map(hb_allscbb_dispav)
df['hbadjbbbb-dispav'] = df['displacement'].map(hb_adjacentbbbb_dispav)
df['hbadjscsc-dispav'] = df['displacement'].map(hb_adjacentscsc_dispav)
df['hbadjscbb-dispav'] = df['displacement'].map(hb_adjacentscbb_dispav)
df['pse6-dispav'] = df['displacement'].map(pse6_dispav)
df['pse4-dispav'] = df['displacement'].map(pse4_dispav)
df['com1-dispav'] = df['displacement'].map(com1_dispav)
df['com2-dispav'] = df['displacement'].map(com2_dispav)
df['com3-dispav'] = df['displacement'].map(com1_dispav)
df['for-dispav'] = df['displacement'].map(forcedispav)
df['vel-dispav'] = df['displacement'].map(vel_dispav)
df['start-raw'] = df['displacement'].map(peakhas)
df['end-raw'] = df['displacement'].map(downhas)
df['vel-dispav']
# print(pse6_dispav)
# this will be added
df.to_csv(os.path.join(outfile, "dataframe_vnew.tsv"), index=False, sep='\t')

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
