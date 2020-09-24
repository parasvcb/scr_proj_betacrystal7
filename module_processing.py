# will process all the dummyfiles to help create dataframe needed for both the 3 and as well as for 3 layer system
import sys
import os
import re
import numpy as np
from progress.bar import Bar
from scipy.spatial import distance as scdist
import subprocess
from collections import Counter


def movefile(source, destination):
    subprocess.check_output(["mv", source, destination])


def makedir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)


def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def vec_sub(vec1lis, vec2lis):
    vec1lis = list(vec1lis)
    vec2lis = list(vec2lis)
    return (scdist.euclidean(vec1lis, vec2lis))


def running_mean_refiner(lis, avgcount, val):
    if val:
        lisav = running_mean(lis, avgcount)
    else:
        lisav = lis
    # print (len(lisav))
    lisavg = ['NA']*(int(avgcount/2 - 1)) + \
        list(lisav) + ['NA']*(int(avgcount/2))
    return lisavg


def toughness(has, toughnessDistance, crysvol):
    # units here are pn, ang and ang^3
    '''
    old and can be removed
    lisra20 = [[lis[i][0], lisra20[i]] for i in range(0, len(lis))]
    # added distance as first element and averaged component as second
    areacompra20 = np.trapz(x=[i[0] for i in lisra20 if i[1] != 'NA'], y=[
        i[1] for i in lisra20 if i[1] != 'NA'])
    area1stra20 = np.trapz(x=[i[0] for i in lisra20 if i[1] != 'NA' and i[0] <= toughnessDistance], y=[
        i[1] for i in lisra20 if i[1] != 'NA' and i[0] <= toughnessDistance])

    area2aagenra20 = np.trapz(x=[i[0] for i in lisra20 if i[1] != 'NA' and i[0] <= 7.6], y=[
        i[1] for i in lisra20 if i[1] != 'NA' and i[0] <= 7.6])

    return raw, avgframe_ra20, avgframe_ra250, {'ra20': (areacompra20, area1stra20, area2aagenra20)}
    '''
    def return_jm3(pnm_forceDist, vol_ang3):
        # print ('force:%s,vol:%s'%(pnm_forceDist,vol_ang3))
        # will return in Mega joules
        return (((pnm_forceDist*1e-22) / (vol_ang3*1e-30))*1e-6)
    tkeys = list(has.keys())
    tkeys.sort()
    x = []
    y = []
    xy = []
    for i in tkeys:
        x += [i]
        y += [has[i]]
        xy += [(i, has[i])]
    areall = np.trapz(x=x, y=y)
    area1stpeak = np.trapz(x=[i[0] for i in xy if i[0] <= toughnessDistance], y=[
        i[1] for i in xy if i[0] <= toughnessDistance])
    area2aagen = np.trapz(x=[i[0] for i in xy if i[0] <= 7.6], y=[
        i[1] for i in xy if i[0] <= 7.6])

    return (return_jm3(areall, crysvol), return_jm3(area1stpeak, crysvol), return_jm3(area2aagen, crysvol))


def angleaverages(filename, dis):
    # this will feed in the H bond data from file to hash
    with open(filename) as fin:
        lis = []
        for i in fin.read().split('\n')[1:]:
            if len(i) > 0:
                ele = i.split()
                frame = int(ele[0])
                angle = float(ele[2])
                lis += [(frame, angle)]
        raw = {i[0]: i[1] for i in lis[:8000]}
    raw, avgframe_ra20, avgframe_ra250 = compute_averages(raw, dis)
    return raw, avgframe_ra20, avgframe_ra250


def dispAvg(distC, varible, dispint=1):
    distval = list(distC.values())
    distval.sort()
    bins = np.arange(-4, int(distval[-1])+6, dispint)
    binshas = {(bins[i], bins[i+1]): [] for i in range(0, len(bins)-1)}
    for i in varible:
        if i in distC:
            dispval = distC[i]
            for j in binshas:
                if dispval >= j[0] and dispval < j[1]:
                    binshas[j] += [varible[i]]
                    break
        else:
            # print(i)
            pass
    # remove the unfilled keys
    binshas = {i: np.mean(binshas[i]) for i in binshas if binshas[i]}
    hastemp = {round(i[0], 3): round(binshas[i], 3) for i in binshas}

    return hastemp


def compute_averages(raw, dis, dispint, flag=True):
    # flag is for distance wise averaging, else for production runs
    tdis = list(dis.keys())
    traw = list(raw.keys())
    tdis.sort()
    traw.sort()
    if flag:
        if dis.keys() != raw.keys():
            # print ("true")
            for i in dis.keys():
                if i not in raw.keys():
                    raw[i] = 0

        lis = [[dis[i], raw[i], i] for i in dis]
        # lis= [distance of frame, property value at that frame and then the frame]
        lis.sort()
        listempn = [i[1] for i in lis]
        # got the properties as per sorted list

        lisra20 = running_mean_refiner(listempn, 20, 1)
        # property is averaged, here for instance force after sorting accoridng to distnace
        avgframe_ra20 = {lis[i][2]: lisra20[i] for i in range(0, len(lisra20))}
        # then they are arranged as per frames
        dispavg = dispAvg(dis, raw, dispint)
        return raw, avgframe_ra20, dispavg
    else:
        # production run
        lis = [[i, raw[i]] for i in raw]
        # lis= [distance of frame, property value at that frame and then the frame]
        lis.sort()
        listempn = [i[1] for i in lis]
        # got the properties as per sorted list
        lisra20 = running_mean_refiner(listempn, 20, 1)
        # property is averaged, here for instance force after sorting accoridng to distnace
        avgframe_ra20 = {lis[i][0]: lisra20[i] for i in range(0, len(lisra20))}
        # then they are arranged as per frames
        return raw, avgframe_ra20


def framestravelled(distC, dispint):
    distval = list(distC.values())
    distval.sort()
    bins = [i for i in np.arange(-4, int(distval[-1])+6, dispint)]
    binshas = {(bins[i], bins[i+1]): 0 for i in range(0, len(bins)-1)}
    for i in distval:
        for j in binshas:
            if i >= j[0] and i < j[1]:
                binshas[j] += 1
                break
    distframes = {}
    ins = 0
    out = 0
    # print(binshas)
    for i in distC:
        ins += 1
        for j in binshas:
            if distC[i] >= j[0] and distC[i] < j[1]:
                out += 1
                distframes[i] = binshas[j]
                break
    return distframes


def forcedistance(dirsim):
    # will read smd file
    with open(os.path.join(dirsim, "configurations/force.conf")) as fin:
        dirx, diry, dirz = re.search(
            r'SMDDir.*', fin.read()).group().split()[-3:]
    data = []
    with open(os.path.join(dirsim, "logs/force.log")) as fin:
        for i in fin.read().split("\n"):
            if len(i) > 0 and i[:3] == "SMD" and len(i.split()) == 8 and int(i.split()[1]) % 500 == 0:
                data += [i.split()[2:]]
    distance = []
    force = []
    for i in data[:8000]:
        distance += [float(i[0]) * float(dirx) + float(i[1]) *
                     float(diry) + float(i[2]) * float(dirz)]
        force += [float(i[3]) * float(dirx) + float(i[4]) *
                  float(diry) + float(i[5]) * float(dirz)]

    distC = {i: val - distance[0] for i, val in enumerate(distance)}
    forceC = {i: val for i, val in enumerate(force)}
    return distC, forceC


def parseEnergyFiles(filename):
    # Frame\tTime\tElec\tVdW\tNonbond\tTotal\tVdWForce\tElecForce\tTotalForce
    with open(filename) as fin:
        dat = [i for i in fin.read().split("\n")[1:] if len(i) > 0]
    hasElec = {}
    hasVdw = {}
    hasTotal = {}
    hasElecForce = {}
    hasVdwForce = {}
    hasTotalForce = {}
    for i in dat:
        ele = i.split()
        frame = int(ele[0])
        hasElec[frame] = float(ele[2])
        hasVdw[frame] = float(ele[3])
        hasTotal[frame] = float(ele[5])
        hasElecForce[frame] = float(ele[6])
        hasVdwForce[frame] = float(ele[7])
        hasTotalForce[frame] = float(ele[8])
    return hasElec, hasVdw, hasTotal, hasElecForce, hasVdwForce, hasTotalForce


def hbondaverages(filename, dis, dispint):
    # this will feed in the H bond data from file to hash
    # print(filename)
    with open(filename) as fin:
        lis = []
        for i in fin.read().split('\n'):
            if len(i) > 0 and re.match(r'\d+\s+[-]?\d+.?\d*', i):
                ele = i.split()
                at = int(ele[0])
                bt = float(ele[1])
                lis += [(at, bt)]
        raw = {i[0]: i[1] for i in lis[:8000]}
    raw, avgframe_ra20, avgframe_ra250 = compute_averages(raw, dis)
    return raw, avgframe_ra20, avgframe_ra250


def writehbhastocsv(has, filename):
    binshas = {(0, 0.1): [], (0.1, 0.2): [], (0.2, 0.3): [], (0.3, 0.4): [], (0.4, 0.5): [
    ], (0.5, 0.6): [], (0.6, 0.7): [], (0.7, 0.8): [], (0.9, 1.0): [], (1.0, 1.1): []}
    hascustom = {'Main_Main': {}, 'Main_Side': {}, 'Side_Side': {}}
    hascustom['Main_Main'] = {i: 0 for i in binshas}
    hascustom['Main_Side'] = {i: 0 for i in binshas}
    hascustom['Side_Side'] = {i: 0 for i in binshas}
    with open(filename + '_occurence.tsv', 'w') as fout:
        fout.write('hbtype\tatoms\toccurence\n')
        if has and len(has) > 0:
            for i in has:
                identity = i[0]
                value = has[i]
                accdon = "_".join(i[2])
                hbtype = "_".join(i[1])
                fout.write('%s\t%s\t%s\n' % (hbtype, accdon, value))
                for j in binshas:
                    if j[0] <= value < j[1]:
                        binshas[j] += [(accdon, i[1], i[0])]
                        hascustom[hbtype][j] += 1
                        break
        else:
            fout.write('Main_Main\tNULL_NULL\0.0\n')

    keys = list(binshas.keys())
    keys.sort()
    with open(filename + '_occurence_definedbins.tsv', 'w') as fout:
        fout.write('Bins\thbtype\tvalues\tFrequencyClass\tFrequencyTotal\n')
        def getdiv(numerator, denominator): return round(
            numerator / denominator, 2) if numerator > 0 and denominator > 0 else 0
        totalSum = sum([sum(hascustom[i].values()) for i in hascustom])
        for hbtype in hascustom:
            classSum = sum(hascustom[hbtype].values())
            for bins in keys:
                binval = hascustom[hbtype][bins]
                fout.write("%s\t%s\t%s\t%s\t%s\n" % ("_".join(map(str, bins)), hbtype, binval, getdiv(
                    binval, classSum), getdiv(binval, totalSum)))

    with open(filename + '_atomdetailed.log', 'w') as fout:
        fout.write('Hbtype\ttotalBonds\n')
        for hbtype in hascustom:
            classSum = sum(hascustom[hbtype].values())
            fout.write("%s\t%s\n" % (hbtype, classSum))

        fout.write('\nOcc_range\tatomshbtype\thbtype\thb_identity\n')
        for i in keys:
            sorteditems = binshas[i]
            # print(sorteditems)
            sorteditems = sorted(sorteditems, key=lambda x: x[1])
            # print(sorteditems)
            for j in sorteditems:
                hbtype = "_".join(j[1])
                fout.write('%s\t%s\t%s\t%s\n' %
                           ("_".join(map(str, i)), j[0], hbtype, " : ".join(j[-1])))
            fout.write('\n')


def hbondaverages_new(directory, dis, dispint, cuttoff=0.6, framerange=False, plot=True, mind=-1, maxd=22):
    '''
    Lets see how i have deisgned this progrAM
    '''
    def refine_atomtype(key_hbatom, typebond):
        # this function will get called only if, theres going to be OT in one of the atom
        if 'OT' in key_hbatom[0] and 'OT' in key_hbatom[1]:
            typebond = ['Main', 'Main']
        elif 'OT' in key_hbatom[0]:
            typebond[0] = 'Main'
        else:
            typebond[1] = 'Main'
        return typebond

    def firstsetocalculations(filerange, directory):
        hbhas = {}
        framewise = {}
        raw_mcmc = {}
        raw_mcsc = {}
        raw_scsc = {}
        # SegCP1-ASN2-Main-N 	 SegBP1-ASN7-Side-OT1 	 100.00%
        # print(filerange)
        for fil in filerange:
            frame = int(fil.split('.')[0])
            framewise[frame] = []
            if frame <= 8000:
                hbtemp = {('Main', 'Side'): 0, ('Main', 'Main')                          : 0, ('Side', 'Side'): 0}
                with open(os.path.join(directory, fil)) as fin:
                    for line in [i for i in fin.read().split('\n')[2:] if len(i) > 0]:
                        acc, don, occ = line.split()
                        # e.g it should be 1000.hbdata
                        # add condition, if both C,C continue statemnet
                        # of OT1 or OT2, then add main
                        accele = acc.split("-")
                        donele = don.split("-")
                        typebond = [accele[2], donele[2]]
                        key_hbatom = [accele[3], donele[3]]
                        if 'C' in key_hbatom[0] and 'C' in key_hbatom[1]:
                            continue
                        if 'OT' in key_hbatom[0] or 'OT' in key_hbatom[1]:
                            typebond = refine_atomtype(key_hbatom, typebond)
                        typebond.sort()
                        typebond = tuple(typebond)
                        identity = [acc, don]
                        identity.sort()
                        identity = tuple(identity)
                        hbhaskey = (identity, typebond, tuple(key_hbatom))
                        if hbhaskey not in hbhas:
                            hbhas[hbhaskey] = [frame]
                        hbhas[hbhaskey] += [frame]
                        framewise[frame] += [(identity, typebond)]
                        hbtemp[typebond] += 1

                raw_mcmc[frame] = hbtemp[('Main', 'Main')]
                raw_mcsc[frame] = hbtemp[('Main', 'Side')]
                raw_scsc[frame] = hbtemp[('Side', 'Side')]
        # print(Counter([len(i) for i in framewise.values()]))
        # print('framewise above')
        # sys.exit()
        # print(Counter([len(i) for i in hbhas.values()]))
        # print('hbhas above')
        # # sys.exit()
        return hbhas, framewise, raw_mcmc, raw_mcsc, raw_scsc

    def stoch_stable_hunter(hbhas, framewise, framerange, cuttoff):
        mcmcstoch_raw = {}
        mcmcstab_raw = {}
        scscstoch_raw = {}
        scscstab_raw = {}
        mcscstoch_raw = {}
        mcscstab_raw = {}
        # print (len(set(framerange)&set(framewise.keys())), list(set(framerange)&set(framewise.keys()))[:5])
        # creating checkpoint
        if set(framerange) < set(framewise.keys()):
            print('passed CK framerange subset of framewise')
        for frame in framerange:
            templishb = []
            templishb = [(list(hbdata[0]), hbdata[1])
                         for hbdata in framewise[frame]]
            # where 0th element is identity tuple and 1st element in tuple of hbtype

            mcmc = []
            scsc = []
            mcsc = []
            for hb in templishb:
                hb_id = hb[0]
                hb_id.sort()
                if hb[1] == ('Main', 'Main'):
                    mcmc += [tuple(hb_id)]
                elif hb[1] == ('Side', 'Side'):
                    scsc += [tuple(hb_id)]
                else:
                    mcsc += [tuple(hb_id)]
            mcmcstoch_raw[frame] = len([i for i in mcmc if (
                i not in hbhas or (i in hbhas and hbhas[i] < cuttoff))])
            mcmcstab_raw[frame] = len(
                [i for i in mcmc if i in hbhas and hbhas[i] >= cuttoff])
            scscstoch_raw[frame] = len([i for i in scsc if (
                i not in hbhas or (i in hbhas and hbhas[i] < cuttoff))])
            scscstab_raw[frame] = len(
                [i for i in scsc if i in hbhas and hbhas[i] >= cuttoff])
            mcscstoch_raw[frame] = len([i for i in mcsc if (
                i not in hbhas or (i in hbhas and hbhas[i] < cuttoff))])
            mcscstab_raw[frame] = len(
                [i for i in mcsc if i in hbhas and hbhas[i] >= cuttoff])
            # the above defined should be modified hence to include data
            # for frames not captured till peak
            # The above is modified

        return mcmcstoch_raw, mcmcstab_raw, scscstoch_raw, scscstab_raw, mcscstoch_raw, mcscstab_raw

    # framerangewillbe a lilst of eligible frames
    filerange = os.listdir(directory) if not framerange \
        else [str(i) + '.hbdata' for i in framerange]

    # removed chunk of caluculations to sub () called firstsetocalculations
    hbhas, framewise, raw_mcmc, raw_mcsc, raw_scsc = firstsetocalculations(
        filerange, directory)
    # the above for the p1 and p2 data is limited to the framelist defined

    if framerange:
        # print(framewise)
        # sys.exit()
        # get in form of stochastic and stable values
        # print(cuttoff, 'cuttoff')
        hbhasoccurence = {hbond: round(len(set(hbhas[hbond]) & set(
            framerange)) / len(framerange), 2) for hbond in hbhas}
        # print(np.histogram(list(hbhasoccurence.values()), bins = np.arange(0, 1.1, 0.1)))
        # for hbond in hbhas:
        #     # here i need to think, about gettingthi data to csv
        #     presence =
        #     # instead of creating new record, add this to hbhas as new val
        #     hbhas[hbond] = presence
        #     # # presence above shoudl be dealth with two frame conditions
        #     # # one going default to the all frames from start to p1 and second going
        #     # # so that it can be written to csv file and
        #     # # second to cropped framerange if any such that default to
        # # # this process above of calculating presence should be
        # # dealt such that set framerange should calculate distance less than some cuttoff
        if plot:
            filetocsvname = os.path.join(os.path.dirname(
                os.path.normpath(directory)), os.path.basename(os.path.normpath(directory)))
            # print(hbhasoccurence)
            writehbhastocsv(hbhasoccurence, filetocsvname)
        mcmcstoch_raw, mcmcstab_raw, scscstoch_raw,\
            scscstab_raw, mcscstoch_raw, mcscstab_raw = stoch_stable_hunter(
                hbhas, framewise, framerange, cuttoff)
        mcmcstoch_raw, mcmcstoch_ra20, mcmcstochdispav = compute_averages(
            mcmcstoch_raw, dis, dispint)
        # print (mcmcstoch_ra20)
        # print (mcmcstochdispav)
        # sys.exit()

        mcmcstab_raw, mcmcstab_ra20, mcmcstabdispav = compute_averages(
            mcmcstab_raw, dis, dispint)
        scscstoch_raw, scscstoch_ra20, scscstochdispav = compute_averages(
            scscstoch_raw, dis, dispint)
        scscstab_raw, scscstab_ra20, scscstabdispav = compute_averages(
            scscstab_raw, dis, dispint)
        mcscstoch_raw, mcscstoch_ra20, mcscstochdispav = compute_averages(
            mcscstoch_raw, dis, dispint)
        mcscstab_raw, mcscstab_ra20, mcscstabdispav = compute_averages(
            mcscstab_raw, dis, dispint)

        # distrange = [dis[i] for i in framerange]
        # maxd = max(distrange)
        # mind = min(distrange)
        # print (maxd,mind)
        # sys.exit()
        # print(raw_mcmcdispav)
        mcmcstochdispav = {i: mcmcstochdispav[i]
                           for i in mcmcstochdispav if mind <= i <= maxd}
        mcmcstabdispav = {i: mcmcstabdispav[i]
                          for i in mcmcstabdispav if mind <= i <= maxd}
        scscstochdispav = {i: scscstochdispav[i]
                           for i in scscstochdispav if mind <= i <= maxd}
        scscstabdispav = {i: scscstabdispav[i]
                          for i in scscstabdispav if mind <= i <= maxd}
        mcscstochdispav = {i: mcscstochdispav[i]
                           for i in mcscstochdispav if mind <= i <= maxd}
        mcscstabdispav = {i: mcscstabdispav[i]
                          for i in mcscstabdispav if mind <= i <= maxd}

        return mcmcstoch_raw, mcmcstoch_ra20, mcmcstochdispav, \
            mcmcstab_raw, mcmcstab_ra20, mcmcstabdispav, \
            scscstoch_raw, scscstoch_ra20, scscstochdispav, \
            scscstab_raw, scscstab_ra20, scscstabdispav, \
            mcscstoch_raw, mcscstoch_ra20, mcscstochdispav, \
            mcscstab_raw, mcscstab_ra20, mcscstabdispav
    else:
        raw_mcmc, raw_mcmc_ra20, raw_mcmcdispav = compute_averages(
            raw_mcmc, dis, dispint)
        raw_mcsc, raw_mcsc_ra20, raw_mcscdispav = compute_averages(
            raw_mcsc, dis, dispint)
        raw_scsc, raw_scsc_ra20, raw_scscdispav = compute_averages(
            raw_scsc, dis, dispint)
        return raw_mcmc, raw_mcmc_ra20, raw_mcmcdispav, raw_mcsc, raw_mcsc_ra20, raw_mcscdispav, raw_scsc, raw_scsc_ra20, raw_scscdispav


def hbonds_calculator3layer(dirsim, dis, p1, d1, p2, dispint, cuttofffrommain):
    # print(p1, d1, p2, "peaks")
    # folder_all = os.path.join(dirsim, 'hbonds_all')
    folder_all = os.path.join(dirsim, 'hbonds_nonadj')
    folder_adj = os.path.join(dirsim, 'hbonds_adjacent')

    # p1 is definately final major chunk of the data we are expecting to discuss,
    # send framerange to find stoch and stab as different range and extrapltae stable informatyiion to the peak

    # above two folders will have the data for H bond types amd per frame files
    allmcmc, allmcmcra20, allmcmcdispav, allmcsc, allmcscra20, allmcscdispav, allscsc, allscscra20, allscscdispav = hbondaverages_new(
        folder_all, dis, dispint, framerange=False)
    adjmcmc, adjmcmcra20, adjmcmcdispav, adjmcsc, adjmcscra20, adjmcscdispav, adjscsc, adjscscra20, adjscscdispav = hbondaverages_new(
        folder_adj, dis, dispint, framerange=False)

    framesrange_first_ascent = [i for i in dis if dis[i] <= p1+0.1]
    framesrange_second_ascent = [i for i in dis if d1-0.1 <= dis[i] <= p2+0.1]
    allp1mcmcstoch_raw, allp1mcmcstoch_ra20, allp1mcmcstochdispav, \
        allp1mcmcstab_raw, allp1mcmcstab_ra20, allp1mcmcstabdispav, \
        allp1scscstoch_raw, allp1scscstoch_ra20, allp1scscstochdispav, \
        allp1scscstab_raw, allp1scscstab_ra20, allp1scscstabdispav, \
        allp1mcscstoch_raw, allp1mcscstoch_ra20, allp1mcscstochdispav, \
        allp1mcscstab_raw, allp1mcscstab_ra20, allp1mcscstabdispav = hbondaverages_new(folder_all, dis, dispint, cuttofffrommain,
                                                                                       framerange=framesrange_first_ascent, maxd=p1)
    # sys.exit()
    adjp1mcmcstoch_raw, adjp1mcmcstoch_ra20, adjp1mcmcstochdispav, \
        adjp1mcmcstab_raw, adjp1mcmcstab_ra20, adjp1mcmcstabdispav, \
        adjp1scscstoch_raw, adjp1scscstoch_ra20, adjp1scscstochdispav, \
        adjp1scscstab_raw, adjp1scscstab_ra20, adjp1scscstabdispav, \
        adjp1mcscstoch_raw, adjp1mcscstoch_ra20, adjp1mcscstochdispav, \
        adjp1mcscstab_raw, adjp1mcscstab_ra20, adjp1mcscstabdispav = hbondaverages_new(folder_adj, dis, dispint, cuttofffrommain,
                                                                                       framerange=framesrange_first_ascent, maxd=p1)
    allp2mcmcstoch_raw, allp2mcmcstoch_ra20, allp2mcmcstochdispav, \
        allp2mcmcstab_raw, allp2mcmcstab_ra20, allp2mcmcstabdispav, \
        allp2scscstoch_raw, allp2scscstoch_ra20, allp2scscstochdispav, \
        allp2scscstab_raw, allp2scscstab_ra20, allp2scscstabdispav, \
        allp2mcscstoch_raw, allp2mcscstoch_ra20, allp2mcscstochdispav, \
        allp2mcscstab_raw, allp2mcscstab_ra20, allp2mcscstabdispav = hbondaverages_new(folder_all, dis, dispint, cuttofffrommain,
                                                                                       framerange=framesrange_second_ascent, plot=False, mind=d1, maxd=p2)
    adjp2mcmcstoch_raw, adjp2mcmcstoch_ra20, adjp2mcmcstochdispav, \
        adjp2mcmcstab_raw, adjp2mcmcstab_ra20, adjp2mcmcstabdispav, \
        adjp2scscstoch_raw, adjp2scscstoch_ra20, adjp2scscstochdispav, \
        adjp2scscstab_raw, adjp2scscstab_ra20, adjp2scscstabdispav, \
        adjp2mcscstoch_raw, adjp2mcscstoch_ra20, adjp2mcscstochdispav, \
        adjp2mcscstab_raw, adjp2mcscstab_ra20, adjp2mcscstabdispav = hbondaverages_new(folder_adj, dis, dispint, cuttofffrommain,
                                                                                       framerange=framesrange_second_ascent, plot=False, mind=d1, maxd=p2)

    has = {'all_mcmc_raw': allmcmc, 'all_mcsc_raw': allmcsc, 'all_scsc_raw': allscsc,
           'all_mcmc_ra20': allmcmcra20, 'all_mcsc_ra20': allmcscra20, 'all_scsc_ra20': allscscra20,
           'all_mcmc_dispav': allmcmcdispav, 'all_mcsc_dispav': allmcscdispav, 'all_scsc_dispav': allscscdispav,
           'adj_mcmc_raw': adjmcmc, 'adj_mcsc_raw': adjmcsc, 'adj_scsc_raw': adjscsc,
           'adj_mcmc_ra20': adjmcmcra20, 'adj_mcsc_ra20': adjmcscra20, 'adj_scsc_ra20': adjscscra20,
           'adj_mcmc_dispav': adjmcmcdispav, 'adj_mcsc_dispav': adjmcscdispav, 'adj_scsc_dispav': adjscscdispav,

           'allp1mcmcstochraw': allp1mcmcstoch_raw, 'allp1mcscstochraw': allp1mcscstoch_raw, 'allp1scscstochraw': allp1scscstoch_raw,
           'allp1mcmcstochra20': allp1mcmcstoch_ra20, 'allp1mcscstochra20': allp1mcscstoch_ra20, 'allp1scscstochra20': allp1scscstoch_ra20,
           'allp1mcmcstochdispav': allp1mcmcstochdispav, 'allp1mcscstochdispav': allp1mcscstochdispav, 'allp1scscstochdispav': allp1scscstochdispav,
           'adjp1mcmcstochraw': adjp1mcmcstoch_raw, 'adjp1mcscstochraw': adjp1mcscstoch_raw, 'adjp1scscstochraw': adjp1scscstoch_raw,
           'adjp1mcmcstochra20': adjp1mcmcstoch_ra20, 'adjp1mcscstochra20': adjp1mcscstoch_ra20, 'adjp1scscstochra20': adjp1scscstoch_ra20,
           'adjp1mcmcstochdispav': adjp1mcmcstochdispav, 'adjp1mcscstochdispav': adjp1mcscstochdispav, 'adjp1scscstochdispav': adjp1scscstochdispav,

           'allp1mcmcstabraw': allp1mcmcstab_raw, 'allp1mcscstabraw': allp1mcscstab_raw, 'allp1scscstabraw': allp1scscstab_raw,
           'allp1mcmcstabra20': allp1mcmcstab_ra20, 'allp1mcscstabra20': allp1mcscstab_ra20, 'allp1scscstabra20': allp1scscstab_ra20,
           'allp1mcmcstabdispav': allp1mcmcstabdispav, 'allp1mcscstabdispav': allp1mcscstabdispav, 'allp1scscstabdispav': allp1scscstabdispav,
           'adjp1mcmcstabraw': adjp1mcmcstab_raw, 'adjp1mcscstabraw': adjp1mcscstab_raw, 'adjp1scscstabraw': adjp1scscstab_raw,
           'adjp1mcmcstabra20': adjp1mcmcstab_ra20, 'adjp1mcscstabra20': adjp1mcscstab_ra20, 'adjp1scscstabra20': adjp1scscstab_ra20,
           'adjp1mcmcstabdispav': adjp1mcmcstabdispav, 'adjp1mcscstabdispav': adjp1mcscstabdispav, 'adjp1scscstabdispav': adjp1scscstabdispav,

           'allp2mcmcstochraw': allp2mcmcstoch_raw, 'allp2mcscstochraw': allp2mcscstoch_raw, 'allp2scscstochraw': allp2scscstoch_raw,
           'allp2mcmcstochra20': allp2mcmcstoch_ra20, 'allp2mcscstochra20': allp2mcscstoch_ra20, 'allp2scscstochra20': allp2scscstoch_ra20,
           'allp2mcmcstochdispav': allp2mcmcstochdispav, 'allp2mcscstochdispav': allp2mcscstochdispav, 'allp2scscstochdispav': allp2scscstochdispav,
           'adjp2mcmcstochraw': adjp2mcmcstoch_raw, 'adjp2mcscstochraw': adjp2mcscstoch_raw, 'adjp2scscstochraw': adjp2scscstoch_raw,
           'adjp2mcmcstochra20': adjp2mcmcstoch_ra20, 'adjp2mcscstochra20': adjp2mcscstoch_ra20, 'adjp2scscstochra20': adjp2scscstoch_ra20,
           'adjp2mcmcstochdispav': adjp2mcmcstochdispav, 'adjp2mcscstochdispav': adjp2mcscstochdispav, 'adjp2scscstochdispav': adjp2scscstochdispav,

           'allp2mcmcstabraw': allp2mcmcstab_raw, 'allp2mcscstabraw': allp2mcscstab_raw, 'allp2scscstabraw': allp2scscstab_raw,
           'allp2mcmcstabra20': allp2mcmcstab_ra20, 'allp2mcscstabra20': allp2mcscstab_ra20, 'allp2scscstabra20': allp2scscstab_ra20,
           'allp2mcmcstabdispav': allp2mcmcstabdispav, 'allp2mcscstabdispav': allp2mcscstabdispav, 'allp2scscstabdispav': allp2scscstabdispav,
           'adjp2mcmcstabraw': adjp2mcmcstab_raw, 'adjp2mcscstabraw': adjp2mcscstab_raw, 'adjp2scscstabraw': adjp2scscstab_raw,
           'adjp2mcmcstabra20': adjp2mcmcstab_ra20, 'adjp2mcscstabra20': adjp2mcscstab_ra20, 'adjp2scscstabra20': adjp2scscstab_ra20,
           'adjp2mcmcstabdispav': adjp2mcmcstabdispav, 'adjp2mcscstabdispav': adjp2mcscstabdispav, 'adjp2scscstabdispav': adjp2scscstabdispav,

           }
    return has


def hbonds_calculator1layer(dirsim, dis):
    filelis = ['hbonds_C_BD.dat', 'hbonds_C_BD_bbbb.dat',
               'hbonds_C_BD_scsc.dat', 'hbonds_C_BD_scbb.dat']
    has_map = {'hbonds_C_BD.dat': 'adjacent', 'hbonds_C_BD_bbbb.dat': 'adjacentbbbb',
               'hbonds_C_BD_scsc.dat': 'adjacentscsc', 'hbonds_C_BD_scbb.dat': 'adjacentscbb'}
    has = {}
    for i in filelis:
        raw, ra20, ra250 = hbondaverages(
            os.path.join(dirsim, i), dis)
        has[has_map[i]] = (raw, ra20, ra250)
    return has


def centreofmasscalc(filename, dis, dispint):
    with open(filename) as fin:
        dat = [i for i in fin.read().split('\n')[1:] if len(i) > 0]
    comupsublow1st = []
    comupsublow2nd = []
    comupsublow3rd = []

    for i in dat:
        # print(i.dat.split("\t"))
        frame, comup1st, comlw1st, comup2nd, comlw2nd, comup3rd, comlw3rd, = i.split(
            '\t')
        comup1st = map(float, comup1st.split())
        comlw1st = map(float, comlw1st.split())
        comup2nd = map(float, comup2nd.split())
        comlw2nd = map(float, comlw2nd.split())
        comup3rd = map(float, comup3rd.split())
        comlw3rd = map(float, comlw3rd.split())
        frame = int(frame)
        comupsublow1st += [vec_sub(comup1st, comlw1st)]
        comupsublow2nd += [vec_sub(comup2nd, comlw2nd)]
        comupsublow3rd += [vec_sub(comup3rd, comlw3rd)]

    com1st = {i: val-comupsublow1st[0]
              for i, val in enumerate(comupsublow1st)}
    com1st, com1_ra20, com1_dispav = compute_averages(com1st, dis, dispint)
    com2nd = {i: val-comupsublow2nd[0]
              for i, val in enumerate(comupsublow2nd)}
    com2nd, com2_ra20, com2_dispav = compute_averages(com2nd, dis, dispint)

    com3rd = {i: val-comupsublow3rd[0]
              for i, val in enumerate(comupsublow3rd)}
    com3rd, com3_ra20, com3_dispav = compute_averages(com3rd, dis, dispint)
    # print("com1st", comupsublow1st[0])
    # print("com2nd", comupsublow2nd[0])
    # print("com3rd", comupsublow3rd[0])
    return com1st, com1_ra20, com1_dispav, com2nd, com2_ra20, com2_dispav, com3rd, com3_ra20, com3_dispav
