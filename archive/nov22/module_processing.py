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


def compute_averages(raw, dis, dispint, notProd=True):
    # notProd is for distance wise averaging, else for production runs
    if notProd:
        if dis.keys() != raw.keys():
            notInRaw = set(dis.keys()) - set(raw.keys())
            dis = {i: dis[i] for i in dis if i not in notInRaw}

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
    # will read smd file, to get the pulling direction coordinates
    with open(os.path.join(dirsim, "configurations/force.conf")) as fin:
        dirx, diry, dirz = re.search(
            r'SMDDir.*', fin.read()).group().split()[-3:]
    data = []
    with open(os.path.join(dirsim, "logs/force.log")) as fin:
        for i in fin.read().split("\n"):
            if len(i) > 0 and i[:3] == "SMD" and len(i.split()) == 8 and int(i.split()[1]) % 500 == 0:
                data += [i.split()[2:]]
    # took and recorded the SMD atom pullout readings at every 0.5ps
    distance = []
    force = []
    # for i in data[:8000]:
    for i in data:
        distance += [float(i[0]) * float(dirx) + float(i[1]) *
                     float(diry) + float(i[2]) * float(dirz)]
        force += [float(i[3]) * float(dirx) + float(i[4]) *
                  float(diry) + float(i[5]) * float(dirz)]
    distC = {i: val - distance[0] for i, val in enumerate(distance)}
    forceC = {i: val for i, val in enumerate(force)}
    # distance and force data were before subsetted to take only the first 8000 entries (for faster pull),
    # I am removing that condition (will be interesting to see, if inconsistencies in number of frames will emerge hence)
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


def getdiv(numerator, denominator, prec): return round(
    numerator / denominator, prec) if numerator > 0 and denominator > 0 else 0


def uniqueHbondasFuncofdisplacement(distC, framewise, filename):
    listrelevant = [[distC[frame], frame]
                    for frame in set(distC.keys()) & set(framewise.keys())]
    listrelevant.sort()
    allbin = set()
    # print(list(framewise.values())[1][1])
    totalhbonds = len(set([j[0] for i in framewise.values() for j in i]))
    with open(filename + '_uniqueAsDisplacementFunc.tsv', 'w') as fout:
        fout.write('Frame\tdisplacement\tuniqueBonds\tfracUniqueBonds\n')
        for value in listrelevant:
            frame = value[1]
            hbhas = {hb[0]: hb[1]
                     for hb in framewise[frame]}
            # for unique count i need to get
            uniqueHb = len(set(hbhas.keys())-allbin)
            allbin |= set(hbhas.keys())
            frac = getdiv(uniqueHb, totalhbonds, 4)
            fout.write('%s\t%s\t%s\t%s\n' % (frame, value[0], uniqueHb, frac))


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
    frac3plus = 0
    frac4plus = 0
    with open(filename + '_occurence_definedbins.tsv', 'w') as fout:
        fout.write('Bins\thbtype\tvalues\tFrequencyClass\tFrequencyTotal\n')
        totalSum = sum([sum(hascustom[i].values()) for i in hascustom])
        for hbtype in hascustom:
            classSum = sum(hascustom[hbtype].values())
            for bins in keys:
                binval = hascustom[hbtype][bins]
                if bins[1] > 0.3:
                    frac3plus += getdiv(binval, totalSum, 2)
                if bins[1] > 0.4:
                    frac4plus += getdiv(binval, totalSum, 2)
                fout.write("%s\t%s\t%s\t%s\t%s\n" % ("_".join(map(str, bins)), hbtype, binval, getdiv(
                    binval, classSum, 2), getdiv(binval, totalSum, 2)))

    print(os.path.basename(os.path.normpath(filename)))
    print('frac3plus=%s,frac4plus=%s' %
          (round(frac3plus, 3), round(frac4plus, 3)))

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


def hbondaverages_new(directory, dis, dispint, cuttoff=0.3, framerange=False, plot=True, mind=-1, maxd=22):
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
                        hbhasKeyAndFramewiseValue = (
                            identity, typebond, tuple(key_hbatom))
                        if hbhasKeyAndFramewiseValue not in hbhas:
                            hbhas[hbhasKeyAndFramewiseValue] = [frame]
                        hbhas[hbhasKeyAndFramewiseValue] += [frame]
                        framewise[frame] += [hbhasKeyAndFramewiseValue]
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

    def stoch_stable_hunter(hbhasoccurence, framewise, framerange, cuttoff):
        mcmcstoch_raw = {}
        mcmcstab_raw = {}
        scscstoch_raw = {}
        scscstab_raw = {}
        mcscstoch_raw = {}
        mcscstab_raw = {}
        # print (len(set(framerange)&set(framewise.keys())), list(set(framerange)&set(framewise.keys()))[:5])
        # creating checkpoint
        # if set(framerange) == set(framewise.keys()):
        #     print('passed CK framerange equal to framewise')
        # hbhasoccurence key has (identity, typebond, tuple(key_hbatom)) as key
        # and frameiwse will have them as value elemnet
        for frame in framerange:
            templishb = []
            templishb = list(framewise[frame])
            # they will cross refer to hbhasoccirce keys
            # where 0th element is identity tuple and 1st element in tuple of hbtype

            mcmc = []
            scsc = []
            mcsc = []
            for hb in templishb:
                if hb[1] == ('Main', 'Main'):
                    mcmc += [hb]
                elif hb[1] == ('Side', 'Side'):
                    scsc += [hb]
                else:
                    mcsc += [hb]
            # print(cuttoff)
            # print([hbhasoccurence[i] for i in mcmc if i in hbhasoccurence])
            mcmcstoch_raw[frame] = len([i for i in mcmc if (
                i not in hbhasoccurence or (i in hbhasoccurence and hbhasoccurence[i] < cuttoff))])
            # print(mcmcstoch_raw[frame])
            mcmcstab_raw[frame] = len(
                [i for i in mcmc if i in hbhasoccurence and hbhasoccurence[i] >= cuttoff])
            # print(mcmcstab_raw[frame])
            scscstoch_raw[frame] = len([i for i in scsc if (
                i not in hbhasoccurence or (i in hbhasoccurence and hbhasoccurence[i] < cuttoff))])
            scscstab_raw[frame] = len(
                [i for i in scsc if i in hbhasoccurence and hbhasoccurence[i] >= cuttoff])
            mcscstoch_raw[frame] = len([i for i in mcsc if (
                i not in hbhasoccurence or (i in hbhasoccurence and hbhasoccurence[i] < cuttoff))])
            mcscstab_raw[frame] = len(
                [i for i in mcsc if i in hbhasoccurence and hbhasoccurence[i] >= cuttoff])
            # break
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
            uniqueHbondasFuncofdisplacement(dis, framewise, filetocsvname)
        mcmcstoch_raw, mcmcstab_raw, scscstoch_raw,\
            scscstab_raw, mcscstoch_raw, mcscstab_raw = stoch_stable_hunter(
                hbhasoccurence, framewise, framerange, cuttoff)
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


def hbonds_calculator3layer(dirsim, appendhbtype, dis, p1, d1, p2, dispint, cuttofffrommain):
    # print(p1, d1, p2, "peaks")
    # folder_all = os.path.join(dirsim, 'hbonds_all')

    # p1 is definately final major chunk of the data we are expecting to discuss,
    # send framerange to find stoch and stab as different range and extrapltae stable informatyiion to the peak

    # above two folders will have the data for H bond types amd per frame files
    mcmc, mcmcra20, mcmcdispav, mcsc, mcscra20, mcscdispav, scsc, scscra20, scscdispav = hbondaverages_new(
        dirsim, dis, dispint, framerange=False)

    framesrange_first_ascent = [i for i in dis if dis[i] <= p1+0.1]
    framesrange_second_ascent = [i for i in dis if d1-0.1 <= dis[i] <= p2+0.1]
    print(max(framesrange_first_ascent), min(framesrange_first_ascent))
    sys.exit()
    p1mcmcstoch_raw, p1mcmcstoch_ra20, p1mcmcstochdispav, \
        p1mcmcstab_raw, p1mcmcstab_ra20, p1mcmcstabdispav, \
        p1scscstoch_raw, p1scscstoch_ra20, p1scscstochdispav, \
        p1scscstab_raw, p1scscstab_ra20, p1scscstabdispav, \
        p1mcscstoch_raw, p1mcscstoch_ra20, p1mcscstochdispav, \
        p1mcscstab_raw, p1mcscstab_ra20, p1mcscstabdispav = hbondaverages_new(dirsim, dis, dispint, cuttofffrommain,
                                                                              framerange=framesrange_first_ascent, maxd=p1)
    # sys.exit()
    p2mcmcstoch_raw, p2mcmcstoch_ra20, p2mcmcstochdispav, \
        p2mcmcstab_raw, p2mcmcstab_ra20, p2mcmcstabdispav, \
        p2scscstoch_raw, p2scscstoch_ra20, p2scscstochdispav, \
        p2scscstab_raw, p2scscstab_ra20, p2scscstabdispav, \
        p2mcscstoch_raw, p2mcscstoch_ra20, p2mcscstochdispav, \
        p2mcscstab_raw, p2mcscstab_ra20, p2mcscstabdispav = hbondaverages_new(dirsim, dis, dispint, cuttofffrommain,
                                                                              framerange=framesrange_second_ascent, plot=False, mind=d1, maxd=p2)

    has = {appendhbtype + 'mcmc-raw': mcmc, appendhbtype + 'mcsc-raw': mcsc, appendhbtype + 'scsc-raw': scsc,
           appendhbtype + 'mcmc-ra20': mcmcra20, appendhbtype + 'mcsc-ra20': mcscra20, appendhbtype + 'scsc-ra20': scscra20,
           appendhbtype + 'mcmc-dispav': mcmcdispav, appendhbtype + 'mcsc-dispav': mcscdispav, appendhbtype + 'scsc-dispav': scscdispav,

           appendhbtype + 'p1mcmcstoch-raw': p1mcmcstoch_raw, appendhbtype + 'p1mcscstoch-raw': p1mcscstoch_raw, appendhbtype + 'p1scscstoch-raw': p1scscstoch_raw,
           appendhbtype + 'p1mcmcstoch-ra20': p1mcmcstoch_ra20, appendhbtype + 'p1mcscstoch-ra20': p1mcscstoch_ra20, appendhbtype + 'p1scscstoch-ra20': p1scscstoch_ra20,
           appendhbtype + 'p1mcmcstoch-dispav': p1mcmcstochdispav, appendhbtype + 'p1mcscstoch-dispav': p1mcscstochdispav, appendhbtype + 'p1scscstoch-dispav': p1scscstochdispav,

           appendhbtype + 'p1mcmcstab-raw': p1mcmcstab_raw, appendhbtype + 'p1mcscstab-raw': p1mcscstab_raw, appendhbtype + 'p1scscstab-raw': p1scscstab_raw,
           appendhbtype + 'p1mcmcstab-ra20': p1mcmcstab_ra20, appendhbtype + 'p1mcscstab-ra20': p1mcscstab_ra20, appendhbtype + 'p1scscstab-ra20': p1scscstab_ra20,
           appendhbtype + 'p1mcmcstab-dispav': p1mcmcstabdispav, appendhbtype + 'p1mcscstab-dispav': p1mcscstabdispav, appendhbtype + 'p1scscstab-dispav': p1scscstabdispav,

           appendhbtype + 'p2mcmcstoch-raw': p2mcmcstoch_raw, appendhbtype + 'p2mcscstoch-raw': p2mcscstoch_raw, appendhbtype + 'p2scscstoch-raw': p2scscstoch_raw,
           appendhbtype + 'p2mcmcstoch-ra20': p2mcmcstoch_ra20, appendhbtype + 'p2mcscstoch-ra20': p2mcscstoch_ra20, appendhbtype + 'p2scscstoch-ra20': p2scscstoch_ra20,
           appendhbtype + 'p2mcmcstoch-dispav': p2mcmcstochdispav, appendhbtype + 'p2mcscstoch-dispav': p2mcscstochdispav, appendhbtype + 'p2scscstoch-dispav': p2scscstochdispav,

           appendhbtype + 'p2mcmcstab-raw': p2mcmcstab_raw, appendhbtype + 'p2mcscstab-raw': p2mcscstab_raw, appendhbtype + 'p2scscstab-raw': p2scscstab_raw,
           appendhbtype + 'p2mcmcstab-ra20': p2mcmcstab_ra20, appendhbtype + 'p2mcscstab-ra20': p2mcscstab_ra20, appendhbtype + 'p2scscstab-ra20': p2scscstab_ra20,
           appendhbtype + 'p2mcmcstab-dispav': p2mcmcstabdispav, appendhbtype + 'p2mcscstab-dispav': p2mcscstabdispav, appendhbtype + 'p2scscstab-dispav': p2scscstabdispav,

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
