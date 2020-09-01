import os
import re
import sys
import pandas as pd
import numpy as np
hbonds = '''
reptype poltype Peaktype        Area
menten  polyalanine_default     Complete        17557.26125
menten  polyalanine_default     1stpeak 11400.791

has_hb_all      has_hb_all_bbbb has_hb_all_scsc has_hb_adjacent has_hb_adjacent_bbbb    has_hb_adjacent_scsc    has_hb_adjacent_scbb    displacement    reptype poltype
18.762  14.127  1.905   14.429  14.127  0.206   0.0     2.1     menten  polyalanine_default
'''
if len(sys.argv) != 3:
    print("Please enter correct cmd options 1 directroy having three sets 2 outputfileappend")
    sys.exit()
prog, parent, outfile = sys.argv

has_toughness = {'Complete': {}, '1stpeak': {}, '2aadist': {}}
has_force = {'peak1': {}, 'down1': {}, 'peak2': {}, 'down2': {}}

has_hbonds = {'has_hb_all': {}, 'has_hb_all_bbbb': {}, 'has_hb_all_scsc': {}, 'has_hb_all_scbb': {},
              'has_hb_adjacent': {}, 'has_hb_adjacent_bbbb': {}, 'has_hb_adjacent_scsc': {}, 'has_hb_adjacent_scbb': {}}


def toughness_reads(filename, has, key):
    dftough = pd.read_csv(filename, sep='\t')
    for seckeys, area in zip(dftough['Peaktype'], dftough['Area']):
        if key not in has[seckeys]:
            has[seckeys][key] = []
        has[seckeys][key] += [area]
    return has


def hbonds_reads(filename, has, key):
    dfHbond = pd.read_csv(filename, sep='\t')
    for columns in dfHbond.columns.values:
        if columns != 'displacement':
            for ind, val in enumerate(dfHbond[columns]):
                string = 'peak_%s' % (ind+1)
                if key not in has[columns]:
                    has[columns][key] = {}
                if string not in has[columns][key]:
                    has[columns][key][string] = []
                has[columns][key][string] += [val]
    return has


def force_reads(filename, has, key):
    dfforce = pd.read_csv(filename, sep='\t')
    #print (dfforce.columns.values)
    #sys.exit()
    for seckeys, force, disp in zip(dfforce['peaktype'], dfforce['force'], dfforce['distance']):
        if key not in has[seckeys]:
            has[seckeys][key] = []
        has[seckeys][key] += [(force,disp)]
    return has


for reptype in os.listdir(parent):
    if reptype in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        for polymer in os.listdir(os.path.join(parent, reptype)):
            if re.match(r'poly[a-z]+', polymer):
                has_toughness = toughness_reads(os.path.join(
                    parent, reptype, polymer, 'processed', 'toughness.tsv'), has_toughness, polymer)
                has_hbonds = hbonds_reads(os.path.join(
                    parent, reptype, polymer, 'processed', 'hbonds_peaks.tsv'), has_hbonds, polymer)
                has_force = force_reads(os.path.join(
                    parent, reptype, polymer, 'processed', 'peaks_distances.tsv'), has_force, polymer)
                # break
# print(has_hbonds)
with open(outfile+"_forcecum.tsv", 'w') as fin:
    fin.write("category\tpolymer\tforce-menton_set\tforce-tyroneP_set\tforce-tyrone_set3\tdistance-menton_set\tdistance-tyroneP_set\tdistance-tyrone_set3\tmean\tstd\n")
    for peak in has_force:
        for polymer in has_force[peak]:
            values_force="\t".join(map(str,[i[0] for i in has_force[peak][polymer]]))
            values_disp="\t".join(map(str,[i[1] for i in has_force[peak][polymer]]))
            fin.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (peak, polymer,values_force,values_disp, np.mean([i[0] for i in 
                has_force[peak][polymer]]), np.std([i[0] for i in has_force[peak][polymer]])))

with open(outfile+"_toughnesscum.tsv", 'w') as fin:
    fin.write("category\tpolymer\tmean\tstd\n")
    for typeof in has_toughness:
        for polymer in has_toughness[typeof]:
            fin.write("%s\t%s\t%s\t%s\n" % (typeof, polymer, np.mean(
                has_toughness[typeof][polymer]), np.std(has_toughness[typeof][polymer])))

with open(outfile+"_hbondsall_peak_1.tsv", 'w') as fin:
    fin.write("category\tpolymer\tmean\tstd\n")
    for typeof in has_hbonds:
        if typeof in ['has_hb_all', 'has_hb_all_bbbb', 'has_hb_all_scsc', 'has_hb_all_scbb']:
            for polymer in has_hbonds[typeof]:
                fin.write("%s\t%s\t%s\t%s\n" % (typeof, polymer, np.mean(
                    has_hbonds[typeof][polymer]['peak_1']), np.std(has_hbonds[typeof][polymer]['peak_1'])))

with open(outfile+"_hbondsall_peak_2.tsv", 'w') as fin:
    fin.write("category\tpolymer\tmean\tstd\n")
    for typeof in has_hbonds:
        if typeof in ['has_hb_all', 'has_hb_all_bbbb', 'has_hb_all_scsc', 'has_hb_all_scbb']:
            for polymer in has_hbonds[typeof]:
                fin.write("%s\t%s\t%s\t%s\n" % (typeof, polymer, np.mean(
                    has_hbonds[typeof][polymer]['peak_2']), np.std(has_hbonds[typeof][polymer]['peak_2'])))

with open(outfile+"_hbondsadjacent_peak_1.tsv", 'w') as fin:
    fin.write("category\tpolymer\tmean\tstd\n")
    for typeof in has_hbonds:
        if typeof in ['has_hb_adjacent', 'has_hb_adjacent_bbbb', 'has_hb_adjacent_scsc', 'has_hb_adjacent_scbb']:
            for polymer in has_hbonds[typeof]:
                fin.write("%s\t%s\t%s\t%s\n" % (typeof, polymer, np.mean(
                    has_hbonds[typeof][polymer]['peak_1']), np.std(has_hbonds[typeof][polymer]['peak_1'])))

with open(outfile+"_hbondsadjacent_peak_2.tsv", 'w') as fin:
    fin.write("category\tpolymer\tmean\tstd\n")
    for typeof in has_hbonds:
        if typeof in ['has_hb_adjacent', 'has_hb_adjacent_bbbb', 'has_hb_adjacent_scsc', 'has_hb_adjacent_scbb']:
            for polymer in has_hbonds[typeof]:
                fin.write("%s\t%s\t%s\t%s\n" % (typeof, polymer, np.mean(
                    has_hbonds[typeof][polymer]['peak_2']), np.std(has_hbonds[typeof][polymer]['peak_2'])))
