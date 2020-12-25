import re
import sys
import os
from collections import Counter

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
'''
Pending:

running and updating all layers
3 reads and improvemnets
'''
if len(sys.argv) != 3:
    print("Please enter correct cmd arguements \
        1:parentDir(having three replicates) 2:outdir")
    sys.exit()

program, dirsim, outdir = sys.argv
del (program)

# follwing will be group of two hashes, has fast and has slow,
hasframes_fast = {
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

hasframes_Slow = {
    'menton_set':
    {
        'polyala-gly': 23097,
        'polyalanine_default': 23483,
        'polyasparagine': 21563,
        'polyglycine': 20118,
        'polyisoleucine': 30127,
        'polythreonine': 13887,
        'polyvaline': 25901
    },
    'tyroneP_set':
    {
        'polyala-gly': 23097,
        'polyalanine_default': 27032,
        'polyasparagine': 21641,
        'polyglycine': 20118,
        'polyisoleucine': 22533,
        'polythreonine': 13386,
        'polyvaline': 23696
    },
    'tyrone_set3':
    {
        'polyala-gly': 25892,
        'polyalanine_default': 25414,
        'polyasparagine': 16029,
        'polyglycine': 17459,
        'polyisoleucine': 29205,
        'polythreonine': 18394,
        'polyvaline': 20268
    }
}

hasFrames = hasframes_fast if 'fast' in os.path.abspath(
    dirsim) else hasframes_Slow


def getframes(has):
    a = has.values()
    return max(a) + 500


def hbbin(val):
    has = {(-0.10, 0.20): 0, (0.20, 0.40): 0, (0.40, 0.60)
            : 0, (0.60, 0.80): 0, (0.80, 1.00): 0}
    for i in has:
        if i[0] < val <= i[1]:
            return str(i[0]) + '-' + str(i[1])


def hbvalfunc(framerange, totalframes, hbpos, hbneg):
    # framerange is list of frames over the hb is presnet and total frames i sthe list of consensun frames
    '''
    framerangeref, hbpos, hbneg = hbvalfunc(
        framerange, totalframes, hbhasval[hb][0], hbhasval[hb][1])
    iterates the consensusFramerange (specofc for hb in which it was found present),
        for every 10th frame of it
            will check if prev and next 5 frames are present in specific framerange,
            give pos value and neg value based on if qualified

    '''
    totalframes.sort()
    ma = totalframes[-1]
    mi = totalframes[0]
    setFramerange = set(framerange)
    framequalifier = {
        iteration: hbpos if iteration in setFramerange else hbneg for iteration in range(mi, ma+1)}
    '''
    below code was used before to record and represent value at every 10th of frame
    total = 0
    pos = 0
    vals = []
    checkEvery = 10  # 1  # should be even interval
    beforeAfter = int(checkEvery / 2)
    valPres = 0.5
    for iteration in range(mi, ma + 1):
        # iterating the consensus range, and choosing the every 10th frame divisioble by 10
        if iteration % checkEvery == 0:
            total += 1
            setInterest = set(range(iteration - beforeAfter,
                                    iteration + 1 + (beforeAfter)))
            presentCount = hbpos if len(
                setInterest & setFramerange) / len(setInterest) >= valPres else hbneg
            vals += [len(setInterest & setFramerange)]
            pos = pos + 1 if presentCount == hbpos else pos + 0
            framequalifier[iteration] = presentCount
    '''

    if framerange and 0:
        print('totalframes:%s (%s to %s)' % (len(totalframes), min(
            totalframes), max(totalframes)))
        print('framerange:%s (%s to %s)' % (len(framerange), min(
            framerange), max(framerange)))
        print(total, pos, Counter(vals), len(
            setFramerange & set(totalframes)), '\n')

    return framequalifier, hbpos + 0.3, hbneg - 0.3


def ring_structure(temphashb):
    # temphas_hb[chainKey][key] += [frameNo]
    # print(temphashb)
    has = {ch: {("__".join(("_".join(hbid[0]), "_".join(hbid[1])))):
                temphashb[ch][hbid] for hbid in temphashb[ch]} for ch in temphashb}
    has_simplerings = {(('B_7_N__C_2_O'), ('B_7_OT1__C_2_N')): 'ring_3_BC',
                       (('B_7_N__C_2_O'), ('B_7_OT2__C_2_N')): 'ring_3_BC',
                       (('B_5_N__C_4_O'), ('B_5_O__C_4_N')): 'ring_2_BC',
                       (('B_3_O__C_6_N'), ('B_5_N__C_4_O')): 'ring_1_BC',
                       (('C_3_N__D_6_O'), ('C_3_O__D_6_N')): 'ring_3_DC',
                       (('C_5_N__D_4_O'), ('C_5_O__D_4_N')): 'ring_2_DC',
                       (('C_7_N__D_2_O'), ('C_7_OT1__D_2_N')): 'ring_1_DC',
                       (('C_7_N__D_2_O'), ('C_7_OT2__D_2_N')): 'ring_1_DC',
                       }
    has_combined_rings = {
        ((('B_7_N__C_2_O'), ('B_7_OT1__C_2_N')), (('C_3_N__D_6_O'), ('C_3_O__D_6_N'))): 'SuperRing_3',
        ((('B_7_N__C_2_O'), ('B_7_OT2__C_2_N')), (('C_3_N__D_6_O'), ('C_3_O__D_6_N'))): 'SuperRing_3',
        ((('B_5_N__C_4_O'), ('B_5_O__C_4_N')), (('C_5_N__D_4_O'), ('C_5_O__D_4_N'))): 'SuperRing_2',
        ((('B_3_O__C_6_N'), ('B_5_N__C_4_O')), (('C_7_N__D_2_O'), ('C_7_OT1__D_2_N'))): 'SuperRing_1',
        ((('B_3_O__C_6_N'), ('B_5_N__C_4_O')), (('C_7_N__D_2_O'), ('C_7_OT2__D_2_N'))): 'SuperRing_1',
    }
    keys = has.keys()
    # print(has.values())
    # print(list(keys)[:5])
    unionset = set()
    for key1 in keys:
        for key2 in has[key1].keys():
            unionset |= set(has[key1][key2])

    hasReverse = {i: [] for i in unionset}
    for hb in has['BC']:
        for frame in has['BC'][hb]:
            hasReverse[frame] += [hb]
    for hb in has['DC']:
        for frame in has['DC'][hb]:
            hasReverse[frame] += [hb]

    has_outSimple_only1 = {'ring_1_BC': [], 'ring_2_BC': [], 'ring_3_BC': [
    ], 'ring_1_DC': [], 'ring_2_DC': [], 'ring_3_DC': []}
    has_outSimple_both = {'ring_1_BC': [], 'ring_2_BC': [], 'ring_3_BC': [
    ], 'ring_1_DC': [], 'ring_2_DC': [], 'ring_3_DC': []}
    has_outCombine = {'SuperRing_1': [], 'SuperRing_2': [], 'SuperRing_3': []}

    # print(len(unionset), 'setunion')
    for frame in unionset:
        # 1,2,3,4 and so on
        for ring in has_simplerings:
            # for 1, ietarte all, hb combs, got the problem, set is needed, termini combinations are being added twice
            if ring[0] in hasReverse[frame] or ring[1] in hasReverse[frame]:
                has_outSimple_only1[has_simplerings[ring]] += [frame]
            if ring[0] in hasReverse[frame] and ring[1] in hasReverse[frame]:
                has_outSimple_both[has_simplerings[ring]] += [frame]
        for ring in has_combined_rings:
            p1 = ring[0]
            p2 = ring[1]
            if (p1[0] in hasReverse[frame] or p1[1] in hasReverse[frame]) and (p2[0] in hasReverse[frame] or p2[1] in hasReverse[frame]):
                has_outCombine[has_combined_rings[ring]] += [frame]
                # print('yes')
    return has_outSimple_only1, has_outSimple_both, has_outCombine, has


def getforce(confile, logfile, frames):
    with open(confile) as fin:
        dirx, diry, dirz = re.search(
            r'SMDDir.*', fin.read()).group().split()[-3:]
    with open(logfile) as fin:
        data = [i.split()[2:] for i in fin.read().split("\n") if len(
            i) > 0 and i[:3] == "SMD" and len(i.split()) == 8 and int(i.split()[1]) % 500 == 0]

    distance = []
    force = []
    for i in data:
        # print (i)
        distance += [float(i[0]) * float(dirx) + float(i[1])
                     * float(diry) + float(i[2]) * float(dirz)]
        force += [float(i[3]) * float(dirx) + float(i[4]) *
                  float(diry) + float(i[5]) * float(dirz)]
    distane = [i-distance[0] for i in distance]
    forceav = []
    for ind, val in enumerate(force):
        if ind >= 2 and ind < len(force)-3:
            forceav += [sum(force[ind-2:ind+3])/5]
        else:
            forceav += [0]
    return {ind: val for ind, val in enumerate(forceav)}


def readthehbonds(has, chainsHead, reptype, poltype, framesupp, parentallowed):
    chainsHead = set(['B', 'C', 'D'])
    # readthehbonds(hasmegaBC, ['B', 'C'], reptype, poltype, theframescount,parentresids)
    # parentresids=[(2,7),(3,6),(4,5)]
    '''
    analysis has focussed on only 1 directory, that is the adjacent hbonded category,
    has, single HB, and then chains BC and CD,
    has, rings, cooperativity, if 1 is presnet and both is present
    has, super_rings, atleast 1 of two different rings
    '''
    dirloc = os.path.join(dirsim, reptype, poltype,
                          'processed', 'hbonds_adjacent')
    # print(dirloc)
    dictionary = {'polyalanine_default': '1ALA', 'polyglycine': '1GLY', 'polyasparagine': '2ASN',
                  'polythreonine': '2THR', 'polyvaline': '3VAL', 'polyisoleucine': '3ILE', 'polyala-gly': '1ALAGLY'}
    polymorph = dictionary[poltype]
    temphas_hb = {'DC': {}, 'BC': {}}
    # temphas_hb will have key as hb idenituty tuple and list as frame counts
    for i in os.listdir(dirloc):
        # print(i)
        frameNo = int(i.split(".")[0])
        if frameNo <= framesupp:
            # framesupp refers to the upper limit [from all polymorphs from peak 1] of frames till peak 1 is reached, processed from the has frame defined above
            with open(os.path.join(dirloc, i)) as fin:
                '''
                Found 11 hbonds.
                donor                    acceptor                occupancy
                SegCP1-ALA2-Main-N       SegBP1-ALA7-Side-OT2    100.00%
                SegCP1-ALA3-Main-N       SegDP1-ALA6-Main-O      100.00%
                SegCP1-ALA4-Main-N       SegBP1-ALA5-Main-O      100.00%
                '''
                raw = fin.read()
                raw = re.sub(r'7-Side-OT1', '7-Main-OT1', raw)
                raw = re.sub(r'7-Side-OT2', '7-Main-OT2', raw)
                dat = [k for k in raw.split("\n")[2:] if len(k) > 0]
            for line in dat:
                ele = line.split()
                hbtype = [ele[0].split('-')[2], ele[1].split('-')[2]]
                atomtype = [ele[0].split('-')[3], ele[1].split('-')[3]]
                chains = [ele[0][3], ele[1][3]]
                chainssort = chains.copy()
                reverseiter = True if 'D' in chainssort else False
                # print(chainssort, reverseiter)
                chainssort.sort(reverse=reverseiter)
                # print(chainssort)
                # ReverseIter = True if D in list, means CD,DC will have CD
                # else False if BC,CB will result BC
                resids = [ele[0].split('-')[1][3], ele[1].split('-')[1][3]]
                sortresid = list(map(int, resids.copy()))
                sortresid.sort()
                cond1 = hbtype == ['Main', 'Main']
                cond2 = not ('C' in atomtype[0] or 'C' in atomtype[1])
                cond3 = set(chainssort) < chainsHead
                cond4 = tuple(sortresid) in parentallowed
                cond5 = 'C' in chainssort
                chainKey = "".join(chainssort)
                if cond1 and cond2 and cond3 and cond4 and cond5:
                    # ABOVE CAN BE REMOVED FOR GREATER h BONDS COVERAGE
                    key = [(chains[0], resids[0], atomtype[0]),
                           (chains[1], resids[1], atomtype[1])]
                    key.sort()
                    key = tuple(key)
                    # print(key)
                    if key not in temphas_hb[chainKey]:
                        temphas_hb[chainKey][key] = []
                    temphas_hb[chainKey][key] += [frameNo]
    # print(temphas_hb['BC'].keys(), temphas_hb['DC'].keys())
    temphas_outSimple_only1, temphas_outSimple_both, temphas_outCombine, temphas_hb = ring_structure(
        temphas_hb)
    # print(temphas_hb['BC'].keys(), temphas_hb['DC'].keys())
    # sys.exit()

    # temphas_hb has been read with care,
    # temphas_ring, temphas_supering will be restored from it later

    # above ran till consensusn max highest frame + 500 i think

    if polymorph not in has['hb']['BC']:
        has['hb']['BC'][polymorph] = {}
    if polymorph not in has['hb']['DC']:
        has['hb']['DC'][polymorph] = {}
    if polymorph not in has['ringOr']:
        has['ringOr'][polymorph] = {}
    if polymorph not in has['ringAnd']:
        has['ringAnd'][polymorph] = {}
    if polymorph not in has['super_ring']:
        has['super_ring'][polymorph] = {}

    def updateHas(temphas, has, polymorph):
        specificframes = hasFrames[reptype][poltype]
        for hb in temphas:
            if hb not in has[polymorph]:
                # this has must be has with subtypes
                # same for the temphas version
                has[polymorph][hb] = []
            occurence = round(
                len([j for j in set(temphas[hb]) if j <= specificframes]) / specificframes, 2)
            if occurence > 1:
                print('greater than 1 freq detected, be careful')
                print(reptype, poltype, hb, occurence,
                      specificframes, len(temphas[hb]), len(set(temphas[hb])))
            has[polymorph][hb] = [list(set(temphas[hb])), occurence]
        return has

    has['hb']['BC'] = updateHas(temphas_hb['BC'], has['hb']['BC'], polymorph)
    has['hb']['DC'] = updateHas(temphas_hb['DC'], has['hb']['DC'], polymorph)
    has['ringOr'] = updateHas(temphas_outSimple_only1,
                              has['ringOr'], polymorph)
    has['ringAnd'] = updateHas(
        temphas_outSimple_both, has['ringAnd'], polymorph)
    has['super_ring'] = updateHas(
        temphas_outCombine, has['super_ring'], polymorph)
    # so this has now will have both the frames list and its occurence value listed
    return has


def writetofile(has, outfile, hbpos, hbneg, consensusFrames, reptype):
    '''
    writetofile(hasmegaBC, '../data_tsv/%s_frameandatbilitywiseBC' %
                    (reptype), hbpos, hbneg, theframescount, reptype)
    '''
    # transforming to previous version
    totalframes = list(range(0, consensusFrames + 1))
    hbhasval = {}
    commonkeys = {}
    for aa in has:
        for hb in has[aa]:
            if hb not in commonkeys:
                commonkeys[hb] = 0
            commonkeys[hb] += 1
    commonlis = [[commonkeys[i], i] for i in commonkeys]
    # common keys will have dict of union of hbkeys with their overall occurence as values from 3 replicates
    # commonlis will have their count and identity in list
    commonlis.sort(reverse=True)
    # descending order
    # print (commonlis)
    colorschema = {}
    count = 1
    for hbinfo in commonlis:
        # hbinfo = [23,(idneityhb)]
        colorschema[hbinfo[1]] = str(count) + '_' + hbinfo[1]
        count += 1
    with open(outfile, 'w') as fout:
        fout.write('frame\tpoltype\thbname\thstability\thbval\n')
        for aa in has:
            for hb in has[aa]:
                hbname = colorschema[hb]
                # hbname will be 1_idenityhbpartner1_idenityhbpartner2
                hstability = hbbin(has[aa][hb][1])
                # print(has[aa][hb][1], hstability)
                # hstability is occurence in given framerange
                framerange = has[aa][hb][0]
                # thisframerange only accounts for list of frames over which this HB was present
                framerange.sort()
                if hb not in hbhasval:
                    hbhasval[hb] = (hbpos, hbneg)
                # print(hb, aa)
                framerangeref, hbpos, hbneg = hbvalfunc(
                    framerange, totalframes, hbhasval[hb][0], hbhasval[hb][1])
                # totalFrames is overallaconsensus range
                # print (hb,hbhasval[hb])
                for tenthFrame in framerangeref:
                    #
                    # will have pos and neg values
                    hbval = framerangeref[tenthFrame]
                    fout.write('%s\t%s\t%s\t%s\t%s\n' %
                               (tenthFrame, aa, hbname[2:], hstability, hbval))
    # above is a master file


def forceWriting(has, outfile, consensusFrames, dirsim, reptype):
    totalframes = list(range(0, consensusFrames + 1))
    dictionary = {'polyalanine_default': '1ALA', 'polyglycine': '1GLY', 'polyasparagine': '2ASN',
                  'polythreonine': '2THR', 'polyvaline': '3VAL', 'polyisoleucine': '3ILE', 'polyala-gly': '1ALAGLY'}
    dicnew = {dictionary[i]: i for i in dictionary}
    with open(outfile, 'w') as fout:
        fout.write('frame\tpoltype\tforceavg\n')
        for aa in has:
            # print (reptype,poltype)
            forceavgref = getforce(os.path.join(dirsim, reptype, dicnew[aa],
                                                'configurations/force.conf'), os.path.join(dirsim, reptype, dicnew[aa],
                                                                                           'logs/force.log'), consensusFrames)
            for i in totalframes:
                fout.write('%s\t%s\t%s\n' % (i, aa, forceavgref[i]))


def getTheOccurenecFrequency(has, outfile):
    # print(has)
    with open(outfile, 'w') as fout:
        fout.write('poltype\thbtype\toccurence\n')
        for aa in has:
            for hb in has[aa]:
                fout.write('%s\t%s\t%s\n' %
                           (aa, hb, has[aa][hb][1]))


deleted = ['polyala_constrained', 'polygly_constrained',
           'polyisoleucine_server191', 'polyisoleucine_78', 'polyisoleucine_turing']

parentresids = [(2, 7), (3, 6), (4, 5)]

print(os.listdir(dirsim))
for reptype in os.listdir(dirsim):
    if reptype in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        hasmega = {
            'hb': {'DC': {}, 'BC': {}},
            'ringOr': {},
            'ringAnd': {},
            'super_ring': {},
        }
        theframescount = getframes(hasFrames[reptype])
        # will give the maximal frame count from all polymorphs in system +500, such that overall force value can be measured.

        for poltype in os.listdir(os.path.join(dirsim, reptype)):
            if re.match(r'poly[a-z]+', poltype) and poltype not in deleted:
                # poltype = 'polyisoleucine'
                hasmega = readthehbonds(
                    hasmega, ['B', 'C'], reptype, poltype, theframescount, parentresids)

                # this has now will have both the frames list and its occurence value listed as list of two keys, first reptype and then poltype
                # break

        hbpos = 1
        hbneg = -1
        writetofile(hasmega['super_ring'], os.path.join(
            outdir, '%s_frameWise_SuperRing.tsv' % (reptype)), hbpos, hbneg, theframescount, reptype)
        # sys.exit()
        hbpos = 1
        hbneg = -1
        writetofile(hasmega['hb']['BC'], os.path.join(
            outdir, '%s_frameWise_hbpresenceBC.tsv' % (reptype)), hbpos, hbneg, theframescount, reptype)
        hbpos = 1
        hbneg = -1
        writetofile(hasmega['hb']['DC'], os.path.join(outdir, '%s_frameWise_hbpresenceDC.tsv' % (
            reptype)), hbpos, hbneg, theframescount, reptype)
        hbpos = 1
        hbneg = -1
        writetofile(hasmega['ringOr'], os.path.join(outdir, '%s_frameWise_SimpleRing.tsv' % (
            reptype)), hbpos, hbneg, theframescount, reptype)
        hbpos = 1
        hbneg = -1
        writetofile(hasmega['ringAnd'], os.path.join(outdir, '%s_frameWise_ComplexRing.tsv' % (
            reptype)), hbpos, hbneg, theframescount, reptype)
        hbpos = 1
        hbneg = -1

        getTheOccurenecFrequency(hasmega['hb']['BC'], os.path.join(
            outdir, '%s_occFreq_hbpresenceBC.tsv' % (reptype)))
        getTheOccurenecFrequency(hasmega['hb']['DC'], os.path.join(outdir, '%s_occFreq_hbpresenceDC.tsv' % (
            reptype)))
        getTheOccurenecFrequency(hasmega['ringOr'], os.path.join(outdir, '%s_occFreq_SimpleRing.tsv' % (
            reptype)))
        getTheOccurenecFrequency(hasmega['ringAnd'], os.path.join(outdir, '%s_occFreq_ComplexRing.tsv' % (
            reptype)))
        getTheOccurenecFrequency(hasmega['super_ring'], os.path.join(
            outdir, '%s_occFreq_SuperRing.tsv' % (reptype)))

        forceWriting(hasmega['hb']['BC'], os.path.join(
            outdir, '%s_forceread.tsv') % (reptype), theframescount, dirsim, reptype)
