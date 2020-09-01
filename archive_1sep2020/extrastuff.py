if not os.path.isfile('wowater.psf'):
        subprocess.check_output(
            ['cp', 'before_mini/no_clash_z_autopsf.psf', 'wowater.psf'])

    if not os.path.isfile('wowater.dcd'):
        indexfile = 'indexfile.txt' if os.path.isfile(
            'indexfile.txt') else getindexfile()
        print("No wowater.dcd, creating one")
        out = subprocess.call(['/home/paras/bin/catdcd', '-o', 'wowater.dcd', '-i', '%s' %
                               indexfile, '-last', '80001', 'dcd_outputs/pull/force_pull.dcd'])
        print(out)
    if not os.path.isfile('wowaterprod.dcd'):
        print("No wowaterprod.dcd, creating one")
        indexfile = 'indexfile.txt' if os.path.isfile(
            'indexfile.txt') else getindexfile()
        out = subprocess.call(['/home/paras/bin/catdcd', '-o', 'wowaterprod.dcd', '-i', '%s' % indexfile, 'dcd_outputs/press_equil2/equil_p.dcd',
                               'dcd_outputs/press_equil3/equil_p.dcd', 'dcd_outputs/press_equil4/equil_p.dcd', 'dcd_outputs/press_equil5/equil_p.dcd'])
        print(out)

    print("HERE")

    indexfile = getindexfile()



    with open("calcrmsd.R", "w") as fout:
        fout.write("%s" % stringRmsd)
 
    res2 = subprocess.check_output(["Rscript", "calcrmsd.R"])

    if 'ERROR' in res2.decode('utf-8'):
        print("ERROR rscript-->>", i)

    print(i)
