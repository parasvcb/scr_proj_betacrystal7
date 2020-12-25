import os
import re
import sys
import pandas as pd
import numpy as np

if len(sys.argv) != 3:
    print("Please enter correct cmd options 1 directroy having three sets 2 outputfileappend")
    sys.exit()
prog, parent, outfile = sys.argv

has_dimesnions = {'raw': {}, 'minimized': {}, 'equilibrated': {}}


def force_reads(filename, has, key):
    dfforce = pd.read_csv(filename, sep='\t')
    for cat in ['raw', 'minimized', 'equilibrated']:
        for valcat, dim in zip(dfforce[cat], dfforce['dim_categories']):
            if dim == 'ca':
                if key not in has[cat]:
                    has[cat][key] = []
                has[cat][key] += [valcat]
    return has


for reptype in os.listdir(parent):
    if reptype in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        for polymer in os.listdir(os.path.join(parent, reptype)):
            if re.match(r'poly[a-z]+', polymer):
                has_dimesnions = force_reads(os.path.join(
                    parent, reptype, polymer, 'processed', 'dimensions.tsv'), has_dimesnions, polymer)

with open(outfile+"_dimensions_table.tsv", 'w') as fin:
    fin.write("polymer\traw\tminimized\tequilibrated\n")
    categories = ['raw', 'minimized', 'equilibrated']
    polymer = ['polyala-gly', 'polyalanine_default', 'polyglycine',
               'polythreonine', 'polyasparagine', 'polyisoleucine', 'polyvaline']
    for pol in polymer:
        val1 = np.mean(has_dimesnions[categories[0]][pol])
        val2 = np.mean(has_dimesnions[categories[1]][pol])
        val3 = np.mean(has_dimesnions[categories[2]][pol])
        valstd = np.std(has_dimesnions[categories[2]][pol])
        fin.write("%s\t%s\t%s\t%s+-%s\n" % (pol, np.round(val1, 3),
                                            np.round(val2, 3), np.round(val3, 3), np.round(valstd, 3)))
