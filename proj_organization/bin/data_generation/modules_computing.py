import re
import os
import subprocess
import sys
# will create all the dummy files ever needed for rest of processing

def callcatdcd(catdcd, indexfile, inpdcd, outdcd):
    subprocess.check_output(
        catdcd+" -o %s -i %s %s" % (outdcd, indexfile, inpdcd), shell=True)
    return

def callscript(script, mode="vmd"):
    if mode == "vmd":
        with open("tempscript.vmd", 'w') as fin:
            fin.write("%s" % script)
        res = subprocess.check_output(
            "vmd -dispdev text -e tempscript.vmd", shell=True)
    else:
        # it should be R
        with open("tempscript.R", 'w') as fin:
            fin.write("%s" % script)
        res = subprocess.check_output(
            "Rscript tempscript.R", shell=True)
    # os.remove("tempscript.vmd")
    return res

def checkFiles(filenamelis, directory=False):
    if directory:
        flag = True
        for i in filenamelis:
            if not os.path.isfile(os.path.join(directory, i)):
                flag = False
                break
        return flag
    else:
        flag = True
        for i in filenamelis:
            if not os.path.isfile(i):
                flag = False
                break
        return flag

def concatenate_dcd(sourcedir, listaftersource, outputdcd, catdcd):
    newlis = []
    for i in listaftersource:
        compadd = os.path.join(sourcedir, i)
        if os.path.isfile(compadd):
            newlis += [compadd]
        else:
            print("file %s not present, complete the production runs" % compadd)
            return False
    stringdcd = " ".join(newlis)
    subprocess.check_output(
        catdcd+" -o %s %s" % (outputdcd, stringdcd), shell=True)
    return outputdcd

def hbonds_calculator3layer(psf, dcd, outfile, mode):
    hbgroup = {'hball': '[atomselect top "protein and not chain C"]',
               'hbnad': '[atomselect top "protein and not chain C B D"]',
               'hbadj': '[atomselect top "chain B D"]'}
    script = '''
        package require hbonds
        mol load psf %s dcd %s
        set g1 [atomselect top "chain C"]
        set g2 %s

        set mol [molinfo top]
        set nf [molinfo $mol get numframes]
    
        for {set i 0} {$i < $nf} {incr i} {
        $g1 frame $i
        $g2 frame $i
        hbonds -sel1 $g1 -sel2 $g2 -frames $i:$i -dist 3.5 -ang 30 -writefile yes -upsel yes -plot no -type unique -detailout %s/$i.hbdata
        }
        exit
        ''' % (psf, dcd, hbgroup[mode], outfile)
    if os.path.isfile('hbonds.dat'):
        os.remove('hbonds.dat')
    res = callscript(script)
    del res
    return

def nanocrystal_dimensions(psf, pdb, centermassfile):
    def fetchdistance(text):
        text = text.decode('utf-8')
        # print(text)
        ab = float(re.search(r'ab=\-?\d+\.?\d+', text).group().split("=")[-1])
        abx = float(re.search(r'abx=\-?\d+\.?\d+',
                              text).group().split("=")[-1])
        cb = float(re.search(r'cb=\-?\d+\.?\d+', text).group().split("=")[-1])
        cbx = float(re.search(r'cbx=\-?\d+\.?\d+',
                              text).group().split("=")[-1])
        ca = float(re.search(r'ca=\-?\d+\.?\d+', text).group().split("=")[-1])
        cax = float(re.search(r'cax=\-?\d+\.?\d+',
                              text).group().split("=")[-1])
        ca = float(re.search(r'ca=\-?\d+\.?\d+', text).group().split("=")[-1])
        aainx = float(re.search(r'aainx=\-?\d+\.?\d+',
                                text).group().split("=")[-1])
        has = {'ab': ab, 'abx': abx, 'cb': cb, 'cbx': cbx,
               'ca': ca, 'cax': cax, 'aainx': aainx}
        return has

    script = '''
    source %s
    mol load psf %s pdb %s
    set a [atomselect top "(chain K L M N O T U) and (name CA) and protein"]
    set b [atomselect top "(chain A B C D E P Q) and (name CA) and protein"]
    set c [atomselect top "(chain F G H J I R S) and (name CA) and protein"]

    set acom [center_of_mass $a]
    set bcom [center_of_mass $b]
    set ccom [center_of_mass $c]

    set acomx [lindex $acom 0]
    set bcomx [lindex $bcom 0]
    set ccomx [lindex $ccom 0]

    set dim [atomselect top "(resid 4 and chain C)"]
    set dimx [measure minmax $dim]
    set outdim [expr [lindex $dimx 0 0] - [lindex $dimx 1 0]]

    set abx [expr $acomx - $bcomx]
    set cbx [expr $ccomx - $bcomx]
    set cax [expr $ccomx - $acomx]

    puts "aainx=$outdim"
    puts "abx=$abx"
    puts "cbx=$cbx"
    puts "cax=$cax"
    puts "ab=[vecdist $acom $bcom]"
    puts "cb=[vecdist $ccom $bcom]"
    puts "ca=[vecdist $ccom $acom]"
    exit
    ''' % (centermassfile, psf, pdb)
    res = callscript(script)
    hasdist = fetchdistance(res)
    return hasdist


def nanocrystal_volume(psf, pdb):
    def fetchvolume(text):
        text = text.decode('utf-8')
        # print(text)
        vol = float(re.search(r'vol=\-?\d+\.?\d+',
                              text).group().split("=")[-1])
        return abs(vol)

    script = '''
    mol load psf %s pdb %s
    set a [atomselect top "protein and backbone"]
    
    set dim [measure minmax $a]
    set x [expr [lindex $dim 0 0] - [lindex $dim 1 0]]
    set y [expr [lindex $dim 0 1] - [lindex $dim 1 1]]
    set z [expr [lindex $dim 0 2] - [lindex $dim 1 2]]
    
    set vol [expr $x * $y *$z]
    
    puts "vol=$vol"
    exit
    ''' % (psf, pdb)
    res = callscript(script)
    volume = fetchvolume(res)
    return volume


def peakdistances(datastream, lis, needmaxorlow):
    upper, lower = lis
    # datastream will be 2d array, first value displacemnet and second value propertyX (which is force over here [to get peaks])
    datastream = [i for i in datastream if upper <= i[0] <= lower]
    # get values withing displacement range
    maxval = datastream[0]
    minval = datastream[0]
    # getting their disp,force pairs
    for i in datastream:
        if i[1] > maxval[1]:
            maxval = i
        if i[1] < minval[1]:
            minval = i
        # comparing force at particular displacement intervals
    if needmaxorlow == 'max':
        return maxval
    else:
        return minval


def rmsd(dcd, minimizedpdb, output):
    script = '''
    library('bio3d')
    dcd <- read.dcd('%s')
    pdb <- read.pdb('%s')
    ca.inds <- atom.select(pdb, 'calpha', chain=c("I","H","G","E","D","C","B","A","L","M","N"))
    xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
    rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
    rdfrompdb=rmsd(pdb$xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
    dfr = data.frame(rmsd=rd,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
    dfrfrompdb = data.frame(rmsd=rdfrompdb,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
    write.table(dfr,"%s/rmsd_firstframe.txt" ,sep = "\\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    write.table(dfrfrompdb,"%s/rmsd_frompdb.txt" ,sep = "\\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
    ''' % (dcd, minimizedpdb, output, output)
    res = callscript(script, mode="R")
    return


def removewaterfromdcd(waterpsf, waterpdb, inpdcd, outdcd, catdcd, tempadd):
    filelis = [inpdcd]
    if checkFiles(filelis):
        script = '''
        mol load psf %s pdb %s
        set a [atomselect top protein]
        set atoms [$a get index]
        set out [open "%s/indexfile.txt" w]
        puts $out $atoms
        close $out
        exit''' % (waterpsf, waterpdb, tempadd)
        # write that script and run it and remove it
        res = callscript(script)
        indexfile = '%s/indexfile.txt' % (tempadd)
        callcatdcd(catdcd, indexfile, inpdcd, outdcd)
        return
    else:
        print("dcd %s is missing.. exiting" % (inpdcd))
        sys.exit()


def removewaterfrompdb(psf, pdb, outapp):
    script = '''
    mol load psf %s pdb %s
    set a [atomselect top protein]
    $a writepdb %s.pdb
    $a writepdb %s.psf
    exit''' % (psf, pdb, outapp, outapp)
    # write that script and run it and remove it
    res = callscript(script)
    return
