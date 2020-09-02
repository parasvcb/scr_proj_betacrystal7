# will process all the dummyfiles to help create dataframe needed for both the 3 and as well as for 3 layer system
import sys
import os
import re
import numpy as np
from progress.bar import Bar
from scipy.spatial import distance as scdist
import subprocess


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
    def return_jm3(vol_ang3, pnm_forceDist):
        return ((pnm_forceDist**-22) / (vol_ang3**-30))

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
        distance += [float(i[0])*float(dirx)+float(i[1]) *
                     float(diry)+float(i[2])*float(dirz)]
        force += [float(i[3])*float(dirx)+float(i[4]) *
                  float(diry)+float(i[5])*float(dirz)]

    distC = {i: val-distance[0] for i, val in enumerate(distance)}
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


def hbondaverages_new(directory, dis, dispint, framerange=False):
    # framerangewillbe a lilst of eligible frames
    # this will feed in the H bond data from frame file to hash
    filerange = os.listdir(directory) if not framerange \
        else [str(i)+'.hbdata' for i in framerange]
    # abvove can be costly step
    # SegCP1-ASN2-Main-N 	 SegBP1-ASN7-Side-OT1 	 100.00%
    hbhas = {}
    raw_mcmc = {}
    raw_mcsc = {}
    raw_scsc = {}
    for fil in filerange:
        frame = int(fil.split('.')[0])
        if frame <= 8000:
            with open(os.path.join(directory, fil)) as fin:
                hbtemp = {('Main', 'Side'): 0, ('Main', 'Main')
                           : 0, ('Side', 'Side'): 0}
                for line in [i for i in fin.read().split('\n')[2:] if len(i) > 0]:
                    acc, don, occ = line.split()
                    # e.g it should be 1000.hbdata
                    typebond = [acc.split('-')[2], don.split('-')[2]]
                    typebond.sort()
                    identity = [acc, don]
                    identity.sort()
                    if identity not in hbhas:
                        hbhas[identity] = 0
                    hbhas[identity] += 1
                    hbtemp[tuple(typebond)] += 1

            raw_mcmc[frame] = hbtemp[('Main', 'Main')]
            raw_mcsc[frame] = hbtemp[('Main', 'Side')]
            raw_scsc[frame] = hbtemp[('Side', 'Side')]
    raw_mcmc, raw_mcmc_ra20, raw_mcmcdispav = compute_averages(
        raw_mcmc, dis, dispint)
    raw_mcsc, raw_mcsc_ra20, raw_mcscdispav = compute_averages(
        raw_mcsc, dis, dispint)
    raw_scsc, raw_scsc_ra20, raw_scscdispav = compute_averages(
        raw_scsc, dis, dispint)
    if framerange:
        distrange = [dis[i] for i in framerange]
        maxd = max(distrange)
        mind = min(distrange)
        raw_mcmcdispav = [raw_mcmcdispav[i]
                          for i in raw_mcmcdispav if i not in mind <= i <= maxd]
        raw_mcscdispav = [raw_mcscdispav[i]
                          for i in raw_mcscdispav if i not in mind <= i <= maxd]
        raw_scscdispav = [raw_scscdispav[i]
                          for i in raw_scscdispav if i not in mind <= i <= maxd]

    return raw_mcmc, raw_mcmc_ra20, raw_mcmcdispav, raw_mcsc, raw_mcsc_ra20, raw_mcscdispav, raw_scsc, raw_scsc_ra20, raw_scscdispav


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


def hbonds_calculator3layer(dirsim, dis, p1, d1, p2, dispint):
    folder_all = os.path.join(dirsim, 'hbonds_all')
    folder_adj = os.path.join(dirsim, 'hbonds_adjacent')
    # above two folders will have the data for H bond types amd per frame files
    allmcmc, allmcmcra20, allmcmcdispav, allmcsc, allmcscra20, allmcscdispav, allscsc, allscscra20, allscscdispav = \
        hbondaverages_new(folder_all, dis, dispint, framerange=False)
    adjmcmc, adjmcmcra20, adjmcmcdispav, adjmcsc, adjmcscra20, adjmcscdispav, adjscsc, adjscscra20, adjscscdispav = \
        hbondaverages_new(folder_adj, dis, dispint, framerange=False)

    framesrange_first_ascent = [i for i in dis if dis[i] <= p1]
    framesrange_second_ascent = [i for i in dis if d1 <= dis[i] <= p2]
    all_p1mcmc, all_p1mcmcra20, all_p1mcmcdispav, all_p1mcsc, all_p1mcscra20, all_p1mcscdispav, all_p1scsc, all_p1scscra20, all_p1scscdispav = \
        hbondaverages_new(folder_all, dis, dispint,
                          framerange=framesrange_first_ascent)
    adj_p1mcmc, adj_p1mcmcra20, adj_p1mcmcdispav, adj_p1mcsc, adj_p1mcscra20, adj_p1mcscdispav, adj_p1scsc, adj_p1scscra20, adj_p1scscdispav = \
        hbondaverages_new(folder_adj, dis, dispint,
                          framerange=framesrange_first_ascent)
    all_p2mcmc, all_p2mcmcra20, all_p2mcmcdispav, all_p2mcsc, all_p2mcscra20, all_p2mcscdispav, all_p2scsc, all_p2scscra20, all_p2scscdispav = \
        hbondaverages_new(folder_all, dis, dispint,
                          framerange=framesrange_second_ascent)
    adj_p2mcmc, adj_p2mcmcra20, adj_p2mcmcdispav, adj_p2mcsc, adj_p2mcscra20, adj_p2mcscdispav, adj_p2scsc, adj_p2scscra20, adj_p2scscdispav = \
        hbondaverages_new(folder_adj, dis, dispint,
                          framerange=framesrange_second_ascent)

    has = {'hbonds_C_protein.dat': 'all', 'hbonds_C_BD.dat': 'adjacent',
           'hbonds_C_protein_bbbb.dat': 'allbbbb', 'hbonds_C_BD_bbbb.dat': 'adjacentbbbb',
           'hbonds_C_protein_scsc.dat': 'allscsc', 'hbonds_C_BD_scsc.dat': 'adjacentscsc',
           'hbonds_C_protein_scbb.dat': 'allscbb', 'hbonds_C_BD_scbb.dat': 'adjacentscbb'}
    has = {'all_mcmc_raw': allmcmc, 'all_mcsc_raw': allmcsc, 'all_scsc_raw': allscsc,
           'all_mcmc_ra20': allmcmcra20, 'all_mcsc_ra20': allmcscra20, 'all_scsc_ra20': allscscra20,
           'all_mcmc_dispav': allmcmcdispav, 'all_mcsc_dispav': allmcscdispav, 'all_scsc_dispav': allscscdispav,
           'adj_mcmc_raw': adjmcmc, 'adj_mcsc_raw': adjmcsc, 'adj_scsc_raw': adjscsc,
           'adj_mcmc_ra20': adjmcmcra20, 'adj_mcsc_ra20': adjmcscra20, 'adj_scsc_ra20': adjscscra20,
           'adj_mcmc_dispav': adjmcmcdispav, 'adj_mcsc_dispav': adjmcscdispav, 'adj_scsc_dispav': adjscscdispav,

           'all_p1_mcmc_raw': all_p1mcmc, 'all_p1_mcsc_raw': all_p1mcsc, 'all_p1_scsc_raw': all_p1scsc,
           'all_p1_mcmc_ra20': all_p1mcmcra20, 'all_p1_mcsc_ra20': all_p1mcscra20, 'all_p1_scsc_ra20': all_p1scscra20,
           'all_p1_mcmc_dispav': all_p1mcmcdispav, 'all_p1_mcsc_dispav': all_p1mcscdispav, 'all_p1_scsc_dispav': all_p1scscdispav,
           'adj_p1_mcmc_raw': adj_p1mcmc, 'adj_p1_mcsc_raw': adj_p1mcsc, 'adj_p1_scsc_raw': adj_p1scsc,
           'adj_p1_mcmc_ra20': adj_p1mcmcra20, 'adj_p1_mcsc_ra20': adj_p1mcscra20, 'adj_p1_scsc_ra20': adj_p1scscra20,
           'adj_p1_mcmc_dispav': adj_p1mcmcdispav, 'adj_p1_mcsc_dispav': adj_p1mcscdispav, 'adj_p1_scsc_dispav': adj_p1scscdispav,

           'all_p2_mcmc_raw': all_p2mcmc, 'all_p2_mcsc_raw': all_p2mcsc, 'all_p2_scsc_raw': all_p2scsc,
           'all_p2_mcmc_ra20': all_p2mcmcra20, 'all_p2_mcsc_ra20': all_p2mcscra20, 'all_p2_scsc_ra20': all_p2scscra20,
           'all_p2_mcmc_dispav': all_p2mcmcdispav, 'all_p2_mcsc_dispav': all_p2mcscdispav, 'all_p2_scsc_dispav': all_p2scscdispav,
           'adj_p2_mcmc_raw': adj_p2mcmc, 'adj_p2_mcsc_raw': adj_p2mcsc, 'adj_p2_scsc_raw': adj_p2scsc,
           'adj_p2_mcmc_ra20': adj_p2mcmcra20, 'adj_p2_mcsc_ra20': adj_p2mcscra20, 'adj_p2_scsc_ra20': adj_p2scscra20,
           'adj_p2_mcmc_dispav': adj_p2mcmcdispav, 'adj_p2_mcsc_dispav': adj_p2mcscdispav, 'adj_p2_scsc_dispav': adj_p2scscdispav,
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


def centreofmasscalc(filename, dis):
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
    com1st, com1_ra20, com1_ra250 = compute_averages(com1st, dis)
    com2nd = {i: val-comupsublow2nd[0]
              for i, val in enumerate(comupsublow2nd)}
    com2nd, com2_ra20, com2_ra250 = compute_averages(com2nd, dis)

    com3rd = {i: val-comupsublow3rd[0]
              for i, val in enumerate(comupsublow3rd)}
    com3rd, com3_ra20, com3_ra250 = compute_averages(com3rd, dis)
    print("com1st", comupsublow1st[0])
    print("com2nd", comupsublow2nd[0])
    print("com3rd", comupsublow3rd[0])
    return com1st, com1_ra20, com1_ra250, com2nd, com2_ra20, com2_ra250, com3rd, com3_ra20, com3_ra250
