import re
import os
import subprocess
import sys
import math
from Bio.PDB import Vector
from Bio.PDB.vectors import calc_dihedral
# will create all the dummy files ever needed for rest of processing


def angle_cal(psf, dcd, output):
    outfilelis = ['%s/coods_resid4' % (output), '%s/coods_resid6' % (output)]
    if not (os.path.isfile(outfilelis[0]) and os.path.isfile(outfilelis[1])):
        # will be calculating pseudo dihedral angles + actual dihedral angles for 2 and 4
        string2 = '''
        mol load psf %s dcd %s
        set outdatafile6 [open "%s/coods_resid6" w]
        set outdatafile4 [open "%s/coods_resid4" w]
        puts $outdatafile6 "frame\tresid6CAxyz\tresid6CBxyz"
        puts $outdatafile4 "frame\tresid4CAxyz\tresid4CBxyz"

        set g6_1 [atomselect top "resid 6 and name C and chain C"]
        set g6_2 [atomselect top "resid 6 and name O and chain C"]
        set g4_1 [atomselect top "resid 4 and name C and chain C"]
        set g4_2 [atomselect top "resid 4 and name O and chain C"]
        set nf [molinfo top get numframes]

        for {set i 0 } {$i < $nf} {incr i} {
        $g6_1 frame $i
        $g6_1 update
        $g6_2 frame $i
        $g6_2 update
        $g4_1 frame $i
        $g4_1 update
        $g4_2 frame $i
        $g4_2 update
        set a1 [$g6_1 get {x y z}]
        set a2 [$g6_2 get {x y z}]
        set b1 [$g4_1 get {x y z}]
        set b2 [$g4_2 get {x y z}]

        puts $outdatafile6 "$i\t$a1\t$a2"
        puts $outdatafile4 "$i\t$b1\t$b2"
        unset a1
        unset a2
        unset b1
        unset b2
        }
        close $outdatafile6
        close $outdatafile4
        exit
        ''' % (psf, dcd, output, output)

        res = callscript(string2)
        del res
    for i in outfilelis:
        refineangleoutput(i)
        print("call", i)
    return


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


def checkFiles(filenamelis, directory="./"):
    flag = True
    for i in filenamelis:
        if not os.path.isfile(os.path.join(directory, i)):
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


def com_calc(psf, dcd, centermassfile, output):
    # this script needs a refresher, unset and remove
    script = '''
    source %s
    mol load psf %s dcd %s
    set filename "%s/com_gcm_pull.log"
    set nf [molinfo top get numframes]
    set outDataFile [open $filename w]
    puts $outDataFile "frameno\tcomupp1st\tcomlow1st\tcomupp2nd\tcomlow2nd\tcomupp3rd\tcomlow3rd"

    set upper1st [atomselect top "((chain G or chain I) and ((resid 7 or resid 6) and name CA)) or (chain H and ((resid 2 or resid 3) and name CA))"]
    set lower1st [atomselect top "((chain L or chain N) and ((resid 7 or resid 6) and name CA)) or (chain M and ((resid 2 or resid 3) and name CA))"]
    set upper2nd [atomselect top "((chain G or chain I) and ((resid 5 or resid 4) and name CA)) or (chain H and ((resid 5 or resid 4) and name CA))"]
    set lower2nd [atomselect top "((chain L or chain N) and ((resid 5 or resid 4) and name CA)) or (chain M and ((resid 5 or resid 4) and name CA))"]
    set upper3rd [atomselect top "((chain G or chain I) and ((resid 3 or resid 2) and name CA)) or (chain H and ((resid 6 or resid 7) and name CA))"]
    set lower3rd [atomselect top "((chain L or chain N) and ((resid 3 or resid 2) and name CA)) or (chain M and ((resid 6 or resid 7) and name CA))"]

    for {set i 0 } {$i < $nf} {incr i} {
    $upper1st frame $i
    $upper1st update
    $lower1st frame $i
    $lower1st update
    $upper2nd frame $i
    $upper2nd update
    $lower2nd frame $i
    $lower2nd update
    $upper3rd frame $i
    $upper3rd update
    $lower3rd frame $i
    $lower3rd update

    set uppveccom1st [center_of_mass $upper1st]
    set lowveccom1st [center_of_mass $lower1st]

    set uppveccom2nd [center_of_mass $upper2nd]
    set lowveccom2nd [center_of_mass $lower2nd]

    set uppveccom3rd [center_of_mass $upper3rd]
    set lowveccom3rd [center_of_mass $lower3rd] 
    puts $outDataFile "$i\t$uppveccom1st\t$lowveccom1st\t$uppveccom2nd\t$lowveccom2nd\t$uppveccom3rd\t$lowveccom3rd"
    unset uppveccom1st
    unset lowveccom1st
    unset uppveccom2nd
    unset lowveccom2nd
    unset uppveccom3rd
    unset lowveccom3rd
    }
    close $outDataFile
    exit
    ''' % (centermassfile, psf, dcd, output)
    res = callscript(script)
    del res
    return


def energy_calc_3layer(psf, pdb, dcd, output, pulling=True):

    stringPulling = '''
    mol load psf %s dcd %s
    set chC [atomselect top "chain C and protein"]
    set chCbb [atomselect top "chain C and protein and backbone"]

    set chBD [atomselect top "(chain D B) and protein"]
    set chBDbb [atomselect top "(chain D B) and protein and backbone"]

    set up [atomselect top "(chain I H G) and protein"]
    set upbb [atomselect top "(chain I H G) and protein and backbone"]

    set lw [atomselect top "(chain L M N) and protein"]
    set lwbb [atomselect top "(chain L M N) and protein and backbone"]

    package require namdenergy
    set totalE [namdenergy -vdw -elec -nonb -sel $chC $chBD -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandBD.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $chBDbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandBD_backbone.interaction"]

    set totalE [namdenergy -vdw -elec -nonb -sel $chC $up -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandUP.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $upbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandUP_backbone.interaction"]

    set totalE [namdenergy -vdw -elec -nonb -sel $chC $lw -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandLOW.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $lwbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandLOW_backbone.interaction"]
    exit
    ''' % (psf, dcd, output, output, output, output, output, output)

    stringStatic = '''
    mol load psf %s dcd %s

    set midd [atomselect top "(chain D A E C B) and protein"]
    set middbb [atomselect top "(chain D A E C B) and protein and backbone"]
    set up [atomselect top "(chain I H G) and protein"]
    set upbb [atomselect top "(chain I H G) and protein and backbone"]
    set lw [atomselect top "(chain L M N) and protein"]
    set lwbb [atomselect top "(chain L M N) and protein and backbone"]

    package require namdenergy
    set totalE [namdenergy -vdw -elec -nonb -sel $midd $up -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/prodprodMiddUp.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $middbb $upbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/prodMiddUp_backbone.interaction"]

    set totalE [namdenergy -vdw -elec -nonb -sel $midd $lw -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/prodMiddLow.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $middbb $lwbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/prodMiddLow_backbone.interaction"]
    exit
    ''' % (psf, dcd, output, output, output, output)
    string = stringPulling if pulling else stringStatic
    res = callscript(string)
    if 'ERROR' in res.decode('utf-8'):
        print("ERROR temppullInt-->>", psf, pdb, output)


def energy_calc_1layer(psf, pdb, dcd, output, pulling=True):

    stringPulling = '''
    mol load psf %s dcd %s
    set chC [atomselect top "chain C and protein"]
    set chCbb [atomselect top "chain C and protein and backbone"]

    set chBD [atomselect top "(chain D B) and protein"]
    set chBDbb [atomselect top "(chain D B) and protein and backbone"]

    package require namdenergy
    set totalE [namdenergy -vdw -elec -nonb -sel $chC $chBD -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandBD.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $chBDbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "%s/pullChainCandBD_backbone.interaction"]
    exit
    ''' % (psf, dcd, output, output, output, output, output, output)

    string = stringPulling
    res = callscript(string)
    if 'ERROR' in res.decode('utf-8'):
        print("ERROR temppullInt-->>", psf, pdb, output)


def hbonds_calculator3layer(psf, dcd, outfile):
    script = '''
        package require hbonds
        mol load psf %s dcd %s
        set g1 [atomselect top "chain C"]
        set g2 [atomselect top "protein and not chain C"]
        set g3 [atomselect top "chain B D"]
        
        set mol [molinfo top]
        set nf [molinfo $mol get numframes]
    
        for {set i 0} {$i < $nf} {incr i} {
        $g1 frame $i
        $g2 frame $i
        $g3 frame $i

        hbonds -sel1 $g1 -sel2 $g2 -frames $i:$i -dist 3.5 -ang 30 -writefile yes -upsel yes -plot no -type unique -detailout %s/hbonds_all/$i.hbdata
        hbonds -sel1 $g1 -sel2 $g3 -frames $i:$i -dist 3.5 -ang 30 -writefile yes -upsel yes -plot no -type unique -detailout %s/hbonds_adjacent/$i.hbdata
        }
        exit
        ''' % (psf, dcd, outfile, outfile)
    os.remove('hbonds.dat')
    res = callscript(script)
    del res
    return


def hbonds_calculator1layer(psf, dcd, outfile):

    script = '''
        package require hbonds
        mol load psf %s dcd %s
        set g1 [atomselect top "chain C"]
        set g3 [atomselect top "chain B D"]
        hbonds -sel1 $g1 -sel2 $g3 -dist 3.5 -ang 25 -outfile %s/hbonds_C_BD.dat -writefile yes -upsel yes -plot no

        set g1 [atomselect top "chain C and backbone"]
        set g3 [atomselect top "chain B D and backbone"]
        hbonds -sel1 $g1 -sel2 $g3 -dist 3.5 -ang 25 -outfile %s/hbonds_C_BD_bbbb.dat -writefile yes -upsel yes -plot no

        set g1 [atomselect top "chain C  and not backbone"]
        set g3 [atomselect top "chain B D  and not backbone"]
        hbonds -sel1 $g1 -sel2 $g3 -dist 3.5 -ang 25 -outfile %s/hbonds_C_BD_scsc.dat -writefile yes -upsel yes -plot no
        exit
        ''' % (psf, dcd, outfile, outfile, outfile)
    fileparent = "%s/hbonds_C_BD.dat" % (outfile)
    filechild1 = "%s/hbonds_C_BD_bbbb.dat" % (outfile)
    filechild2 = "%s/hbonds_C_BD_scsc.dat" % (outfile)
    filenewchild = "%s/hbonds_C_BD_scbb.dat" % (outfile)
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


def refineangleoutput(filename):
    with open(filename) as fin:
        dat = fin.read()
        dat = re.sub(r'{', '', dat)
        dat = re.sub(r'}', '', dat)
        dat = [i for i in dat.split('\n')[1:] if len(i) > 0]
        print(dat[0], 'dat[0]')
    newname = re.sub(r'coods', 'angle', filename)
    with open(newname, 'w') as fin:
        fin.write("Frame\tDihedralAngleRaw\tDihedralAngle360\n")
        cacoods = []
        cbcoods = []
        for i in dat:
            # print(i)
            frame, ca, cb = i.split("\t")
            ca = Vector(list(map(float, ca.split())))
            cb = Vector(list(map(float, cb.split())))
            cacoods += [ca]
            cbcoods += [cb]
        refca = cacoods[0]
        refcb = cbcoods[0]
        sub360 = False
        subraw = False
        for i in range(0, len(cacoods)):
            angleraw = math.degrees(calc_dihedral(
                refca, refcb, cbcoods[i], cacoods[i]))
            #print(i, angleraw)
            angle360 = angleraw if angleraw >= 0 else angleraw+360
            if i == 0:
                # sub360 = angle360
                # subraw = angleraw
                sub360 = 0
                subraw = 0
            fin.write("%s\t%s\t%s\n" % (i, angleraw-subraw, angle360-sub360))


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
