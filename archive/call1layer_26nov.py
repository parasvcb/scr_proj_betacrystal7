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
mdpro.makedir(processingdir)
del program

prottype = re.search(r'poly\w+\-?\w*', dirsim).group()
reptype = "menten" if "menton" in dirsim else "tyroneP" if "tyroneP" in dirsim else "tyroneR"

# -> arrow below needs to be fulfilled and otther dimesnions should be recorded but only for ra20 force peaks

has_testing = {'menten': {'polyalanine': {'up': [1, 2.5], 'down': [4, 6]},
                          'polyala-gly': {'up': [1, 2], 'down': [2, 3]},
                          'polyvaline': {'up': [1, 4], 'down': [4, 6]},
                          'polythreonine': {'up': [1, 6], 'down': [4, 8]},
                          'polyglycine': {'up': [1, 2], 'down': [2, 3]},
                          'polyasparagine': {'up': [2, 6], 'down': [4, 8]},
                          'polyisoleucine': {'up': [1, 4], 'down': [4, 6]},
               'tyroneP': {'polyalanine': {'up': [1, 4], 'down': [4, 6]},
                          'polyala-gly': {'up': [1, 2], 'down': [2, 4]},
                          'polyvaline': {'up': [1, 4], 'down': [4, 6]},
                          'polythreonine': {'up': [1, 6], 'down': [4, 8]},
                          'polyglycine': {'up': [1, 2], 'down': [2, 3]},
                          'polyasparagine': {'up': [2, 6], 'down': [4, 8]},
                          'polyisoleucine': {'up': [2, 4], 'down': [4, 6]},
               'tyroneR': {'polyalanine': {'up': [1, 4], 'down': [4, 6]},
                          'polyala-gly': {'up': [1, 2], 'down': [2, 3]},
                          'polyvaline': {'up': [1, 4], 'down': [4, 6]},
                          'polythreonine': {'up': [1, 6], 'down': [4, 8]},
                          'polyglycine': {'up': [1, 2], 'down': [2, 3]},
                          'polyasparagine': {'up': [2, 6], 'down': [4, 8]},
                          'polyisoleucine': {'up': [2, 6], 'down': [4, 8]},
               }  # up here means, upper peak range, where that hsould be found
               # constrast from 3 layer system, wher its split in two
               # this depicts bin, should move in bins of updown values
# above needs to be reverified and change in values
poltype = ''
for i in has_testing['menton']:
    if i in dirsim:
        poltype = i

for i in os.listdir(dirsim):
    if os.path.isfile(os.path.join(dirsim, i)) and i != 'par_all36_prot.prm':
        source = os.path.join(dirsim, i)
        destination = os.path.join(processingdir, i)
        mdpro.movefile(source, destination)

centerofmasstclscript = "$HOME/bin/centerofmassgeom.tcl"
catdcd = "$HOME/bin/catdcd"
waterpsf = os.path.join(dirsim, 'before_mini', "ionized.psf")
waterpdb = os.path.join(dirsim, 'before_mini', "ionized.pdb")
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


if not os.path.isfile(psf):
    psfsource = os.path.join(dirsim, 'before_mini/no_clash_z_autopsf.psf')
    subprocess.check_output(['cp', psfsource, psf])

if not os.path.isfile(dcdpull):
    mdcom.removewaterfromdcd(
        waterpsf, waterpdb, waterdcdpull, dcdpull, catdcd, processingdir)

dispint = 0.1
distC, forceC = mdpro.forcedistance(dirsim)
forcedispav = mdpro.dispAvg(distC, forceC, dispint)
keystemp = list(distC.keys())
keystemp.sort()
newarrtemp = [[i, forcedispav[i]] for i in keystemp]
peak1 = mdcom.peakdistances(
    newarrtemp, has_testing[reptype][poltype]['up'], "max")
down1 = mdcom.peakdistances(
    newarrtemp, has_testing[reptype][poltype]['down'], "min")


if not os.path.isfile(os.path.join(processingdir, "_toughness.tsv")):
    toughness_data_all, toughness_data_1stpeak, toughness_data_2aa = mdpro.toughness(
        forcedispav, down1[0])
    with open(os.path.join(processingdir, "_toughness.tsv"), 'w') as fout:
        fout.write("reptype\tpoltype\tPeaktype\tArea\n")
        fout.write("%s\t%s\tComplete\t%s\n" %
                   (reptype, poltype, toughness_data_all))
        fout.write("%s\t%s\t1stpeak\t%s\n" %
                   (reptype, poltype, toughness_data_1stpeak))
        fout.write("%s\t%s\t2aadist\t%s\n" %
                   (reptype, poltype, toughness_data_2aa))

velocityC = mdpro.velocity(distC)
velocitydispav = mdpro.dispAvg(distC, velocityC, dispint)

framestravelled = mdpro.framestravelled(distC)

with open(os.path.join(processingdir, "peaks_distances.tsv"), "w") as fin:
    fin.write("peaktype\tdistance\tforce\tpoltype\treptype\n")
    fin.write("peak1\t%s\t%s\t%s\t%s\n" %
              (peak1[0], peak1[1], poltype, reptype))
    fin.write("down1\t%s\t%s\t%s\t%s\n" %
              (down1[0], down1[1], poltype, reptype))

# ->forcepeaks are pending and corresponding calls

forceC, force_ra20, force_ra250, toughness_data = mdpro.compute_averages(
    forceC, distC)
velocityC, velocity_ra20, velocity_ra250 = mdpro.compute_averages(
    velocityC, distC)

# refined
# check first of data exists, if yes ,, call for processing, else for computing
if not mdcom.checkFiles(['hbonds_C_BD.dat', 'hbonds_C_BD_bbbb.dat', 'hbonds_C_BD_scsc.dat',
                         'hbonds_C_BD_scbb.dat'], directory=processingdir):
    mdcom.hbonds_calculator3layer(psf=psf, dcd=dcdpull, outfile=processingdir)

if not mdcom.checkFiles(['angle_resid4', 'angle_resid2'], directory=processingdir):
    print("yes")
    mdcom.angle_cal(psf, dcdpull, processingdir)

hbondData = mdpro.hbonds_calculator3layer(processingdir, distC)
hb_adjacent_raw, hb_adjacent_ra20, hb_adjacent_ra250 = hbondData['adjacent']
hb_adjacentbbbb_raw, hb_adjacentbbbb_ra20, hb_adjacentbbbb_ra250 = hbondData[
    'adjacentbbbb']
hb_adjacentscsc_raw, hb_adjacentscsc_ra20, hb_adjacentscsc_ra250 = hbondData[
    'adjacentscsc']
hb_adjacentscbb_raw, hb_adjacentscbb_ra20, hb_adjacentscbb_ra250 = hbondData[
    'adjacentscbb']
pse2, pse2ra20, pse2ra250 = mdpro.angleaverages(
    os.path.join(processingdir, 'angle_resid4'), distC)
pse4, pse4ra20, pse4ra250 = mdpro.angleaverages(
    os.path.join(processingdir, 'angle_resid2'), distC)


hb_adjacent_dispav = mdpro.dispAvg(distC, hb_adjacent_raw, dispint)
hb_adjacentbbbb_dispav = mdpro.dispAvg(
    distC, hb_adjacentbbbb_raw, dispint)
hb_adjacentscsc_dispav = mdpro.dispAvg(
    distC, hb_adjacentscsc_raw, dispint)
hb_adjacentscbb_dispav = mdpro.dispAvg(
    distC, hb_adjacentscbb_raw, dispint)
pse2_dispav = mdpro.dispAvg(distC, pse2, dispint)
pse4_dispav = mdpro.dispAvg(distC, pse4, dispint)

framecount = 8001
for i in hb_adjacentbbbb_dispav.keys():
    if i not in distC.values():
        distC[framecount] = i
        framecount += 1

df = pd.DataFrame([forceC, force_ra20, distC, velocityC, velocity_ra20,
                   framestravelled, pse2, pse2ra20, pse4, pse4ra20,
                   hb_adjacent_raw, hb_adjacent_ra20,
                   hb_adjacentbbbb_raw, hb_adjacentbbbb_ra20, hb_adjacentscsc_raw,
                   hb_adjacentscsc_ra20, hb_adjacentscbb_raw, hb_adjacentscbb_ra20],
                  index=['for-raw', 'for-ra20', 'displacement', 'vel-raw', 'vel-ra20',
                         'frametrav-raw', 'pse2-raw', 'pse2-ra20', 'pse4-raw', 'pse4-ra20',
                         'hbadjacent-raw', 'hbadjacent-ra20', 'hbadjbbbb-raw', 'hbadjbbbb-ra20', 'hbadjscsc-raw',
                         'hbadjscsc-ra20', 'hbadjscbb-raw', 'hbadjscbb-ra20']).transpose()
df['frames'] = df.index
df['hbadjdispav'] = df['displacement'].map(hb_adjacent_dispav)
df['hbadjbbbbdispav'] = df['displacement'].map(hb_adjacentbbbb_dispav)
df['hbadjscscdispav'] = df['displacement'].map(hb_adjacentscsc_dispav)
df['hbadjscbbdispav'] = df['displacement'].map(hb_adjacentscbb_dispav)
df['pse2dispav'] = df['displacement'].map(pse2_dispav)
df['pse4dispav'] = df['displacement'].map(pse4_dispav)
df['forcedispav'] = df['displacement'].map(forcedispav)
# print(pse2_dispav)
# this will be added
df.to_csv(os.path.join(outfile, "dataframe_vnew.tsv"),
          index=False, sep='\t')
'''
old and should be removed
keys_1 = list(distC.keys())
keys_1.sort()
with open(outfile+".tsv", "w") as fin:
    hbondBDstring = "1_HbBD-raw\t1_HbBD-ra20\t1_HbBD-ra250"
    velocitystring = "2_velocity-raw\t2_velocity-ra20\t2_velocity-ra250"

    pse2string = "5_pseudo2-raw\t5_pseudo2-ra20\t5_pseudo2-ra250"
    pse4string = "6_pseudo4-raw\t6_pseudo4-ra20\t6_pseudo4-ra250"

    phi2string = "7_phi2-raw\t7_phi2-ra20\t7_phi2-ra250"
    phi4string = "8_phi4-raw\t8_phi4-ra20\t8_phi4-ra250"

    psi2string = "9_psi2-raw\t9_psi2-ra20\t9_psi2-ra250"
    psi4string = "10_psi4-raw\t10_psi4-ra20\t10_psi4-ra250"
    framestravelledstring = "12_framestravelled-raw"
    framestravelledstring2 = "13_framestravelled2-raw"
    force1string = "14_force1-raw\t13_force1-ra20\t13_force1-ra250"

    force2string = "15_force2-raw\t14_force2-ra20\t14_force2-ra250"
    # not writing frames here that was in variable i

    fin.write("distance\tframes\t%s\n" % ("\t".join(
        [hbondBDstring, velocitystring,
         pse2string, pse4string, phi2string, phi4string, psi2string,
         psi4string, framestravelledstring, framestravelledstring2,
         force1string, force2string])))

    bar = Bar('Processing', max=len(keys_1))
    for i in keys_1:
        hbondBDdata = "%s\t%s\t%s" % (
            hbondBD[i], avgframeBD_ra20[i], avgframeBD_ra250[i])
        pse2data = "%s\t%s\t%s" % (pse2[i], pse2ra20[i], pse2ra250[i])
        pse4data = "%s\t%s\t%s" % (pse4[i], pse4ra20[i], pse2ra250[i])
        phi2data = "%s\t%s\t%s" % (phi2[i], phi2ra20[i], phi2ra250[i])
        phi4data = "%s\t%s\t%s" % (phi4[i], phi4ra20[i], phi2ra250[i])
        psi2data = "%s\t%s\t%s" % (psi2[i], psi2ra20[i], psi2ra250[i])
        psi4data = "%s\t%s\t%s" % (psi4[i], psi4ra20[i], psi2ra250[i])

        velocitydata = "%s\t%s\t%s" % (
            velocityC[i], velocity_ra20[i], velocity_ra250[i])
        framestravelleddata = "%s" % framestravelled[i]
        force1data = "%s\t%s\t%s" % (
            forceC[i], force_ra20[i], force_ra250[i])

        force2data = "%s\t%s\t%s" % (
            forceC[i], force_ra20[i], force_ra250[i])
        fin.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            distC[i], i, hbondBDdata, velocitydata,
            pse2data, pse4data, phi2data, phi4data, psi2data, psi4data,
            framestravelleddata, framestravelleddata, force1data, force2data))
        bar.next()
    bar.finish()
topkey = list(hbondBD.keys())
topkey.sort()
frames = {i: i for i in range(0, topkey[-1])}

print(toughness_data)
fout2 = open(os.path.join(outfile, "_toughness.tsv"), 'w')
fout2.write("color\tPeaktype\tArea\n")
for i in toughness_data:
    fout2.write("%s\tComplete\t%s\n" % (i, toughness_data[i][0]))
    fout2.write("%s\t1stpeak\t%s\n" % (i, toughness_data[i][1]))
    fout2.write("%s\t2aadist\t%s\n" % (i, toughness_data[i][2]))
fout2.close()
'''
