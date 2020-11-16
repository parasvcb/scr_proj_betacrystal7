import sys
import re,os
import subprocess
print ("Please make sure wowater.psf pullwowater.dcd and wowaterprod.dcd present outside else catdcd should be in path, wo wowater.psf, program will end")
if len(sys.argv)!=2:
    print("Please enter correct 0:args python ....py 1:directoryStructure(this should have appropriate tree structure with dirs starting with poly)")
    sys.exit()
prog,directoryStructure=sys.argv

def getindexfile():
    string='''
    mol load psf before_mini/ionized.psf pdb before_mini/ionized.pdb
    set a [atomselect top "protein"]
    set atoms [$a get index]
    set out [open "indexfile.txt" w]
    $a writepdb protein.pdb
    puts $out $atoms
    close $out
    mol load psf before_mini/ionized.psf pdb before_mini/coordmini.pdb
    set a [atomselect top "protein"]
    $a writepdb minimized_protein.pdb
    exit   
    '''
    with open ("temprotind.tcl",'w') as fout:
        fout.write("%s"%string)
    subprocess.check_output(['vmd','-dispdev','text','-e','temprotind.tcl'])
    os.remove('temprotind.tcl')
    return 'indexfile.txt'

stringRmsd='''
library('bio3d')
dcd <- read.dcd('wowaterprod.dcd')
pdb <- read.pdb('minimized_protein.pdb')
ca.inds <- atom.select(pdb, 'calpha', chain=c("I","H","G","E","D","C","B","A","L","M","N"))
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
rdfrompdb=rmsd(pdb$xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
dfr = data.frame(rmsd=rd,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
dfrfrompdb = data.frame(rmsd=rdfrompdb,frames=c(1: dim(dcd)[1]),grp='whole_protein',protgrp='_all')
write.table(dfr,"rmsd_firstframe.txt" ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
write.table(dfrfrompdb,"rmsd_frompdb.txt" ,sep = "\t" ,row.names = FALSE, col.names = TRUE, dec=".", quote=FALSE)
'''

stringPulling='''
    mol load psf wowater.psf dcd wowater.dcd
    set chC [atomselect top "chain C and protein"]
    set chCbb [atomselect top "chain C and protein and backbone"]

    set chBD [atomselect top "(chain D B) and protein"]
    set chBDbb [atomselect top "(chain D B) and protein and backbone"]

    set up [atomselect top "(chain I H G) and protein"]
    set upbb [atomselect top "(chain I H G) and protein and backbone"]

    set lw [atomselect top "(chain L M N) and protein"]
    set lwbb [atomselect top "(chain L M N) and protein and backbone"]

    package require namdenergy
    set totalE [namdenergy -vdw -elec -nonb -sel $chC $chBD -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "pullChainCandBD.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $chBDbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "pullChainCandBD_backbone.interaction"]

    set totalE [namdenergy -vdw -elec -nonb -sel $chC $up -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "pullChainCandUP.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $upbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "pullChainCandUP_backbone.interaction"]

    set totalE [namdenergy -vdw -elec -nonb -sel $chC $lw -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "pullChainCandLOW.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $chCbb $lwbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "pullChainCandLOW_backbone.interaction"]
    exit
    '''
    
stringStatic='''
    mol load psf wowater.psf dcd wowaterprod.dcd
    
    set midd [atomselect top "(chain D A E C B) and protein"]
    set middbb [atomselect top "(chain D A E C B) and protein and backbone"]
    set up [atomselect top "(chain I H G) and protein"]
    set upbb [atomselect top "(chain I H G) and protein and backbone"]
    set lw [atomselect top "(chain L M N) and protein"]
    set lwbb [atomselect top "(chain L M N) and protein and backbone"]

    package require namdenergy
    set totalE [namdenergy -vdw -elec -nonb -sel $midd $up -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "prodprodMiddUp.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $middbb $upbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "prodMiddUp_backbone.interaction"]

    set totalE [namdenergy -vdw -elec -nonb -sel $midd $lw -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "prodMiddLow.interaction"]
    set totalE [namdenergy -vdw -elec -nonb -sel $middbb $lwbb -par "par_all36_prot.prm" -keepforce -projforce -exe  "~/bin/namd2_m" -ofile "prodMiddLow_backbone.interaction"]
    exit
    '''
    
os.chdir(directoryStructure)
for i in [os.path.abspath(j) for j in os.listdir(".") if j[0:4]=="poly" and not re.search(r'poly[A-Z]G',j)]:
    os.chdir(i)
    print (i)
    if not os.path.isfile('wowater.psf'):
        subprocess.check_output(['cp','before_mini/no_clash_z_autopsf.psf','wowater.psf'])
    
    if not os.path.isfile('wowater.dcd'):
        indexfile='indexfile.txt' if os.path.isfile('indexfile.txt') else getindexfile()
        print ("No wowater.dcd, creating one")
        out=subprocess.call(['/home/paras/bin/catdcd','-o','wowater.dcd','-i','%s'%indexfile,'-last','80001','dcd_outputs/pull/force_pull.dcd'])
        print (out)
    if not os.path.isfile('wowaterprod.dcd'):
        print ("No wowaterprod.dcd, creating one")
        indexfile='indexfile.txt' if os.path.isfile('indexfile.txt') else getindexfile()
        out=subprocess.call(['/home/paras/bin/catdcd','-o','wowaterprod.dcd','-i','%s'%indexfile,'dcd_outputs/press_equil2/equil_p.dcd',\
        'dcd_outputs/press_equil3/equil_p.dcd','dcd_outputs/press_equil4/equil_p.dcd','dcd_outputs/press_equil5/equil_p.dcd'])
        print (out)
    
    print ("HERE")
    
    indexfile=getindexfile()
    with open("temppullInt.tcl","w") as fout:
        fout.write("%s"%stringPulling)

    with open("tempprodInt.tcl","w") as fout:
        fout.write("%s"%stringStatic)

    with open("calcrmsd.R","w") as fout:
        fout.write("%s"%stringRmsd)
    
    res=subprocess.check_output(["vmd","-dispdev","text","-e","temppullInt.tcl"])
    if 'ERROR' in res.decode('utf-8'):
        print ("ERROR temppullInt-->>",i)
    res1=subprocess.check_output(["vmd","-dispdev","text","-e","tempprodInt.tcl"])
    
    if 'ERROR' in res1.decode('utf-8'):
        print ("ERROR tempprodInt-->>",i)    
    res2=subprocess.check_output(["Rscript","calcrmsd.R"])
    
    if 'ERROR' in res2.decode('utf-8'):
        print ("ERROR rscript-->>",i)

    print (i)
    #break

