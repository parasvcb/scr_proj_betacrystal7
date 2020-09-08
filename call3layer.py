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

if not (os.path.isfile(os.path.join(processingdir, "toughness.tsv")) and 0):
    volume = mdcom.nanocrystal_volume(waterpsf, waterpdb)
    #print (volume,poltype)
    toughness_data_all, toughness_data_1stpeak, toughness_data_2aa = mdpro.toughness(
        forcedispav, down1[0], volume)
    with open(os.path.join(processingdir, "toughness.tsv"), 'w') as fout:
        fout.write("Peaktype\tArea\n")
        fout.write("Complete\t%s\n" % (toughness_data_all))
        fout.write("1stpeak\t%s\n" % (toughness_data_1stpeak))
        fout.write("2aadist\t%s\n" % (toughness_data_2aa))
sys.exit()
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
#print(hbondData.keys())
# sys.exit()
{'all_mcmc_raw' 'all_mcsc_raw' 'all_scsc_raw'
'all_mcmc_ra20' 'all_mcsc_ra20' 'all_scsc_ra20'
'all_mcmc_dispav' 'all_mcsc_dispav' 'all_scsc_dispav'
'adj_mcmc_raw' 'adj_mcsc_raw' 'adj_scsc_raw'
'adj_mcmc_ra20' 'adj_mcsc_ra20' 'adj_scsc_ra20'
'adj_mcmc_dispav' 'adj_mcsc_dispav' 'adj_scsc_dispav'

'allp1mcmcstochraw' 'allp1mcscstochraw' 'allp1scscstochraw'
'allp1mcmcstochra20' 'allp1mcscstochra20' 'allp1scscstochra20'
'allp1mcmcstochdispav' 'allp1mcscstochdispav' 'allp1scscstochdispav'
'adjp1mcmcstochraw' 'adjp1mcscstochraw' 'adjp1scscstochraw'
'adjp1mcmcstochra20' 'adjp1mcscstochra20' 'adjp1scscstochra20'
'adjp1mcmcstochdispav' 'adjp1mcscstochdispav' 'adjp1scscstochdispav'

'allp1mcmcstabraw' 'allp1mcscstabraw' 'allp1scscstabraw'
'allp1mcmcstabra20' 'allp1mcscstabra20' 'allp1scscstabra20'
'allp1mcmcstabdispav' 'allp1mcscstabdispav' 'allp1scscstabdispav'
'adjp1mcmcstabraw' 'adjp1mcscstabraw' 'adjp1scscstabraw'
'adjp1mcmcstabra20' 'adjp1mcscstabra20' 'adjp1scscstabra20'
'adjp1mcmcstabdispav' 'adjp1mcscstabdispav' 'adjp1scscstabdispav'

'allp2mcmcstochraw' 'allp2mcscstochraw' 'allp2scscstochraw'
'allp2mcmcstochra20' 'allp2mcscstochra20' 'allp2scscstochra20'
'allp2mcmcstochdispav' 'allp2mcscstochdispav' 'allp2scscstochdispav'
'adjp2mcmcstochraw' 'adjp2mcscstochraw' 'adjp2scscstochraw'
'adjp2mcmcstochra20' 'adjp2mcscstochra20' 'adjp2scscstochra20'
'adjp2mcmcstochdispav' 'adjp2mcscstochdispav' 'adjp2scscstochdispav'

'allp2mcmcstabraw' 'allp2mcscstabraw' 'allp2scscstabraw'
'allp2mcmcstabra20' 'allp2mcscstabra20' 'allp2scscstabra20'
'allp2mcmcstabdispav' 'allp2mcscstabdispav' 'allp2scscstabdispav'
'adjp2mcmcstabraw' 'adjp2mcscstabraw' 'adjp2scscstabraw'
'adjp2mcmcstabra20' 'adjp2mcscstabra20' 'adjp2scscstabra20'
'adjp2mcmcstabdispav' 'adjp2mcscstabdispav' 'adjp2scscstabdispav'
}



com1st, com1_ra20, com1_dispav, com2nd, com2_ra20, com2_dispav, com3rd, com3_ra20, com3_dispav = mdpro.centreofmasscalc(
    os.path.join(processingdir, "com_gcm_pull.log"), distC, dispint)

df = pd.DataFrame([forceC, force_ra20, distC,   framestravelled,
                   com1st, com1_ra20,  com2nd, com2_ra20,  com3rd, com3_ra20,
                   hbondData['all_mcmc_raw'], hbondData['all_mcmc_ra20'], hbondData['all_mcsc_raw'], hbondData[
                       'all_mcsc_ra20'], hbondData['all_scsc_raw'], hbondData['all_scsc_ra20'],
                   hbondData['adj_mcmc_raw'], hbondData['adj_mcmc_ra20'], hbondData['adj_mcsc_raw'], hbondData[
                       'adj_mcsc_ra20'], hbondData['adj_scsc_raw'], hbondData['adj_scsc_ra20'],
                   
                   hbondData['allp1mcmcstochraw'], hbondData['allp1mcmcstochra20'], hbondData['allp1mcmcstabraw'], hbondData['allp1mcmcstabra20'], 
                   hbondData['allp1mcscstochraw'], hbondData['allp1mcscstochra20'], hbondData['allp1mcscstabraw'], hbondData['allp1mcscstabra20'], 
                   hbondData['allp1scscstochraw'], hbondData['allp1scscstochra20'], hbondData['allp1scscstabraw'], hbondData['allp1scscstabra20'], 

                   hbondData['adjp1mcmcstochraw'], hbondData['adjp1mcmcstochra20'], hbondData['adjp1mcmcstabraw'], hbondData['adjp1mcmcstabra20'], 
                   hbondData['adjp1mcscstochraw'], hbondData['adjp1mcscstochra20'], hbondData['adjp1mcscstabraw'], hbondData['adjp1mcscstabra20'], 
                   hbondData['adjp1scscstochraw'], hbondData['adjp1scscstochra20'], hbondData['adjp1scscstabraw'], hbondData['adjp1scscstabra20'], 

                   hbondData['allp2mcmcstochraw'], hbondData['allp2mcmcstochra20'], hbondData['allp2mcmcstabraw'], hbondData['allp2mcmcstabra20'], 
                   hbondData['allp2mcscstochraw'], hbondData['allp2mcscstochra20'], hbondData['allp2mcscstabraw'], hbondData['allp2mcscstabra20'], 
                   hbondData['allp2scscstochraw'], hbondData['allp2scscstochra20'], hbondData['allp2scscstabraw'], hbondData['allp2scscstabra20'], 

                   hbondData['adjp2mcmcstochraw'], hbondData['adjp2mcmcstochra20'], hbondData['adjp2mcmcstabraw'], hbondData['adjp2mcmcstabra20'], 
                   hbondData['adjp2mcscstochraw'], hbondData['adjp2mcscstochra20'], hbondData['adjp2mcscstabraw'], hbondData['adjp2mcscstabra20'], 
                   hbondData['adjp2scscstochraw'], hbondData['adjp2scscstochra20'], hbondData['adjp2scscstabraw'], hbondData['adjp2scscstabra20'], 
                   ],
                  index=['for-raw', 'for-ra20', 'displacement', 'frametrav-raw',
                         'com1-raw', 'com1-ra20', 'com2-raw', 'com2-ra20', 'com3-raw', 'com3-ra20',
                         'hballmcmc-raw', 'hballmcmc-ra20', 'hballmcsc-raw', 'hballmcsc-ra20', 'hballscsc-raw', 'hballscsc-ra20',
                         'hbadjmcmc-raw', 'hbadjmcmc-ra20', 'hbadjmcsc-raw', 'hbadjmcsc-ra20', 'hbadjscsc-raw', 'hbadjscsc-ra20',
                        
                   'hballp1mcmcstoch-raw', 'hballp1mcmcstoch-ra20', 'hballp1mcmcstab-raw', 'hballp1mcmcstab-ra20', 
                   'hballp1mcscstoch-raw', 'hballp1mcscstoch-ra20', 'hballp1mcscstab-raw', 'hballp1mcscstab-ra20', 
                   'hballp1scscstoch-raw', 'hballp1scscstoch-ra20', 'hballp1scscstab-raw', 'hballp1scscstab-ra20', 

                   'hbadjp1mcmcstoch-raw', 'hbadjp1mcmcstoch-ra20', 'hbadjp1mcmcstab-raw', 'hbadjp1mcmcstab-ra20', 
                   'hbadjp1mcscstoch-raw', 'hbadjp1mcscstoch-ra20', 'hbadjp1mcscstab-raw', 'hbadjp1mcscstab-ra20', 
                   'hbadjp1scscstoch-raw', 'hbadjp1scscstoch-ra20', 'hbadjp1scscstab-raw', 'hbadjp1scscstab-ra20', 

                   'hballp2mcmcstoch-raw', 'hballp2mcmcstoch-ra20', 'hballp2mcmcstab-raw', 'hballp2mcmcstab-ra20', 
                   'hballp2mcscstoch-raw', 'hballp2mcscstoch-ra20', 'hballp2mcscstab-raw', 'hballp2mcscstab-ra20', 
                   'hballp2scscstoch-raw', 'hballp2scscstoch-ra20', 'hballp2scscstab-raw', 'hballp2scscstab-ra20', 

                   'hbadjp2mcmcstoch-raw', 'hbadjp2mcmcstoch-ra20', 'hbadjp2mcmcstab-raw', 'hbadjp2mcmcstab-ra20', 
                   'hbadjp2mcscstoch-raw', 'hbadjp2mcscstoch-ra20', 'hbadjp2mcscstab-raw', 'hbadjp2mcscstab-ra20', 
                   'hbadjp2scscstoch-raw', 'hbadjp2scscstoch-ra20', 'hbadjp2scscstab-raw', 'hbadjp2scscstab-ra20']).transpose()

df['frames-raw'] = df.index
df['hballmcmc-dispav'] = df['displacement'].map(hbondData['all_mcmc_dispav'])
df['hballmcsc-dispav'] = df['displacement'].map(hbondData['all_mcsc_dispav'])
df['hballscsc-dispav'] = df['displacement'].map(hbondData['all_scsc_dispav'])

df['hbadjmcmc-dispav'] = df['displacement'].map(hbondData['adj_mcmc_dispav'])
df['hbadjmcsc-dispav'] = df['displacement'].map(hbondData['adj_mcsc_dispav'])
df['hbadjscsc-dispav'] = df['displacement'].map(hbondData['adj_scsc_dispav'])
#print(hbondData['allp1mcmc_dispav'])
df['hballp1mcmcstoch-dispav'] = df['displacement'].map(hbondData['allp1mcmcstochdispav'])
df['hballp1scscstoch-dispav'] = df['displacement'].map(hbondData['allp1scscstochdispav'])
df['hballp1mcscstoch-dispav'] = df['displacement'].map(hbondData['allp1mcscstochdispav'])

df['hballp1mcmcstab-dispav'] = df['displacement'].map(hbondData['allp1mcmcstabdispav'])
df['hballp1scscstab-dispav'] = df['displacement'].map(hbondData['allp1scscstabdispav'])
df['hballp1mcscstab-dispav'] = df['displacement'].map(hbondData['allp1mcscstabdispav'])

df['hbadjp1mcmcstoch-dispav'] = df['displacement'].map(hbondData['adjp1mcmcstochdispav'])
df['hbadjp1scscstoch-dispav'] = df['displacement'].map(hbondData['adjp1scscstochdispav'])
df['hbadjp1mcscstoch-dispav'] = df['displacement'].map(hbondData['adjp1mcscstochdispav'])

df['hbadjp1mcmcstab-dispav'] = df['displacement'].map(hbondData['adjp1mcmcstabdispav'])
df['hbadjp1scscstab-dispav'] = df['displacement'].map(hbondData['adjp1scscstabdispav'])
df['hbadjp1mcscstab-dispav'] = df['displacement'].map(hbondData['adjp1mcscstabdispav'])

df['hballp2mcmcstoch-dispav'] = df['displacement'].map(hbondData['allp2mcmcstochdispav'])
df['hballp2scscstoch-dispav'] = df['displacement'].map(hbondData['allp2scscstochdispav'])
df['hballp2mcscstoch-dispav'] = df['displacement'].map(hbondData['allp2mcscstochdispav'])

df['hballp2mcmcstab-dispav'] = df['displacement'].map(hbondData['allp2mcmcstabdispav'])
df['hballp2scscstab-dispav'] = df['displacement'].map(hbondData['allp2scscstabdispav'])
df['hballp2mcscstab-dispav'] = df['displacement'].map(hbondData['allp2mcscstabdispav'])

df['hbadjp2mcmcstoch-dispav'] = df['displacement'].map(hbondData['adjp2mcmcstochdispav'])
df['hbadjp2scscstoch-dispav'] = df['displacement'].map(hbondData['adjp2scscstochdispav'])
df['hbadjp2mcscstoch-dispav'] = df['displacement'].map(hbondData['adjp2mcscstochdispav'])

df['hbadjp2mcmcstab-dispav'] = df['displacement'].map(hbondData['adjp2mcmcstabdispav'])
df['hbadjp2scscstab-dispav'] = df['displacement'].map(hbondData['adjp2scscstabdispav'])
df['hbadjp2mcscstab-dispav'] = df['displacement'].map(hbondData['adjp2mcscstabdispav'])

df['com1-dispav'] = df['displacement'].map(com1_dispav)
df['com2-dispav'] = df['displacement'].map(com2_dispav)
df['com3-dispav'] = df['displacement'].map(com1_dispav)
df['for-dispav'] = df['displacement'].map(forcedispav)

df.to_csv(os.path.join(outfile, "dataframe_vnew2.tsv"), index=False, sep='\t')
#sys.exit()


allp1hb_mcmc = {peak1[0]: hbondData['all_mcmc_dispav'][peak1[0]]}
allp1hb_mcsc = {peak1[0]: hbondData['all_mcsc_dispav'][peak1[0]]}
allp1hb_scsc = {peak1[0]: hbondData['all_scsc_dispav'][peak1[0]]}
adjp1hb_mcmc = {peak1[0]: hbondData['adj_mcmc_dispav'][peak1[0]]}
adjp1hb_mcsc = {peak1[0]: hbondData['adj_mcsc_dispav'][peak1[0]]}
adjp1hb_scsc = {peak1[0]: hbondData['adj_scsc_dispav'][peak1[0]]}

allp2hb_mcmc = {peak2[0]: hbondData['all_mcmc_dispav'][peak2[0]]}
allp2hb_mcsc = {peak2[0]: hbondData['all_mcsc_dispav'][peak2[0]]}
allp2hb_scsc = {peak2[0]: hbondData['all_scsc_dispav'][peak2[0]]}
adjp2hb_mcmc = {peak2[0]: hbondData['adj_mcmc_dispav'][peak2[0]]}
adjp2hb_mcsc = {peak2[0]: hbondData['adj_mcsc_dispav'][peak2[0]]}
adjp2hb_scsc = {peak2[0]: hbondData['adj_scsc_dispav'][peak2[0]]}

print (hbondData['allp1mcmcstochdispav'])
allp1stochhbmcmc = {peak1[0]: hbondData['allp1mcmcstochdispav'][peak1[0]]}
allp1stochhbmcsc = {peak1[0]: hbondData['allp1mcscstochdispav'][peak1[0]]}
allp1stochhbscsc = {peak1[0]: hbondData['allp1scscstochdispav'][peak1[0]]}
allp1stabhbmcmc = {peak1[0]: hbondData['allp1mcmcstabdispav'][peak1[0]]}
allp1stabhbmcsc = {peak1[0]: hbondData['allp1mcscstabdispav'][peak1[0]]}
allp1stabhbscsc = {peak1[0]: hbondData['allp1scscstabdispav'][peak1[0]]}

adjp1stochhbmcmc = {peak1[0]: hbondData['adjp1mcmcstochdispav'][peak1[0]]}
adjp1stochhbmcsc = {peak1[0]: hbondData['adjp1mcscstochdispav'][peak1[0]]}
adjp1stochhbscsc = {peak1[0]: hbondData['adjp1scscstochdispav'][peak1[0]]}
adjp1stabhbmcmc = {peak1[0]: hbondData['adjp1mcmcstabdispav'][peak1[0]]}
adjp1stabhbmcsc = {peak1[0]: hbondData['adjp1mcscstabdispav'][peak1[0]]}
adjp1stabhbscsc = {peak1[0]: hbondData['adjp1scscstabdispav'][peak1[0]]}

allp2stochhbmcmc = {peak2[0]: hbondData['allp2mcmcstochdispav'][peak2[0]]}
allp2stochhbmcsc = {peak2[0]: hbondData['allp2mcscstochdispav'][peak2[0]]}
allp2stochhbscsc = {peak2[0]: hbondData['allp2scscstochdispav'][peak2[0]]}
allp2stabhbmcmc = {peak2[0]: hbondData['allp2mcmcstabdispav'][peak2[0]]}
allp2stabhbmcsc = {peak2[0]: hbondData['allp2mcscstabdispav'][peak2[0]]}
allp2stabhbscsc = {peak2[0]: hbondData['allp2scscstabdispav'][peak2[0]]}

adjp2stochhbmcmc = {peak2[0]: hbondData['adjp2mcmcstochdispav'][peak2[0]]}
adjp2stochhbmcsc = {peak2[0]: hbondData['adjp2mcscstochdispav'][peak2[0]]}
adjp2stochhbscsc = {peak2[0]: hbondData['adjp2scscstochdispav'][peak2[0]]}
adjp2stabhbmcmc = {peak2[0]: hbondData['adjp2mcmcstabdispav'][peak2[0]]}
adjp2stabhbmcsc = {peak2[0]: hbondData['adjp2mcscstabdispav'][peak2[0]]}
adjp2stabhbscsc = {peak2[0]: hbondData['adjp2scscstabdispav'][peak2[0]]}

dfhbp1 = pd.DataFrame([allp1hb_mcmc,allp1hb_mcsc,allp1hb_scsc,\
adjp1hb_mcmc,adjp1hb_mcsc,adjp1hb_scsc,\
allp1stochhbmcmc,allp1stochhbmcsc,allp1stochhbscsc,\
allp1stabhbmcmc,allp1stabhbmcsc,allp1stabhbscsc,\
adjp1stochhbmcmc,adjp1stochhbmcsc,adjp1stochhbscsc,\
adjp1stabhbmcmc,adjp1stabhbmcsc,adjp1stabhbscsc],
index=['allhb_mcmc','allhb_mcsc','allhb_scsc',\
'adjhb_mcmc','adjhb_mcsc','adjhb_scsc',\
'allp1stochhb_mcmc','allp1stochhb_mcsc','allp1stochhb_scsc',\
'allp1stabhb_mcmc','allp1stabhb_mcsc','allp1stabhb_scsc',\
'adjp1stochhb_mcmc','adjp1stochhb_mcsc','adjp1stochhb_scsc',\
'adjp1stabhb_mcmc','adjp1stabhb_mcsc','adjp1stabhb_scsc']).transpose()

dfhbp2 = pd.DataFrame([allp2hb_mcmc,allp2hb_mcsc,allp2hb_scsc,\
adjp2hb_mcmc,adjp2hb_mcsc,adjp2hb_scsc,\
allp2stochhbmcmc,allp2stochhbmcsc,allp2stochhbscsc,\
allp2stabhbmcmc,allp2stabhbmcsc,allp2stabhbscsc,\
adjp2stochhbmcmc,adjp2stochhbmcsc,adjp2stochhbscsc,\
adjp2stabhbmcmc,adjp2stabhbmcsc,adjp2stabhbscsc],
index=['allhb_mcmc','allhb_mcsc','allhb_scsc',\
'adjhb_mcmc','adjhb_mcsc','adjhb_scsc',\
'allp2stochhb_mcmc','allp2stochhb_mcsc','allp2stochhb_scsc',\
'allp2stabhb_mcmc','allp2stabhb_mcsc','allp2stabhb_scsc',\
'adjp2stochhb_mcmc','adjp2stochhb_mcsc','adjp2stochhb_scsc',\
'adjp2stabhb_mcmc','adjp2stabhb_mcsc','adjp2stabhb_scsc']).transpose()

dfhbp1['displacement'] = dfhbp1.index
dfhbp2['displacement'] = dfhbp2.index

dfhbp1.to_csv(os.path.join(outfile, "hbonds_peaks1.tsv"), index=False, sep='\t')
dfhbp2.to_csv(os.path.join(outfile, "hbonds_peaks2.tsv"), index=False, sep='\t')
