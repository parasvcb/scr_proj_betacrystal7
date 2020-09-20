import os, subprocess,re
parent = "../"
has_data={}
for i in os.listdir(parent):
    if i in ["menton_set", "tyroneP_set", "tyrone_set3"]:
        has_data[i]={}
        for j in os.listdir(os.path.join(parent,i)):
            if re.match(r'poly[a-z]+',j):
                key=i+'_'+j
                has_data[i][j]={}
                directory=os.path.join(parent,i,j,'processed','hbonds_all')
                for fil in os.listdir(directory):
                    with open (os.path.join(directory,fil)) as fin:
                        dat=[k for k in fin.read().split("\n")[2:] if len(k)>0]
                    for line in dat:
                        acc,don,occ=line.split()
                        accele=acc.split("-")
                        donele=don.split("-")
                        key_hbtype=[accele[2],donele[2]]
                        key_hbatom=[accele[3],donele[3]]
                        key_hbtype.sort()
                        key_hbatom.sort()
                        key_hbtype=tuple(key_hbtype)
                        key_hbatom=tuple(key_hbatom)
                        if key_hbtype!=('Main','Main'):
                            if key_hbtype not in has_data[i][j]:
                                has_data[i][j][key_hbtype]={}
                            if key_hbatom not in has_data[i][j][key_hbtype]:
                                has_data[i][j][key_hbtype][key_hbatom]=[]
                            has_data[i][j][key_hbtype][key_hbatom]+=[int(fil.split(".")[0])]
with open ("data_temp_hjbondtype",'w') as fout:
    fout.write("Rep\tpol\th-type\tatompair\ttotalframeocc\toccupancy\n")
    for rep in has_data:
        for pol in has_data[rep]:
            for hbtype in has_data[rep][pol]:
                for atompair in has_data[rep][pol][hbtype]:
                    val=len(has_data[rep][pol][hbtype][atompair])
                    fout.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(rep,pol,hbtype,atompair,val,(val/8000)*100))
                fout.write("---------\n")
            fout.write("\n")
        fout.write("\n\n")