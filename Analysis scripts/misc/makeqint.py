def makeqint(zval):
    
    """
    edits lammps data files to get whole integer charges (not yet generalised) for solute molecule   
    """
    #address and name of old datafiles 
    adress='/home/kselvon/Amber-amino-acids/aminoacids/newGLU/'
    oldfile = 'data.elba_glu_sca_neutral_z'+str(zval)
    lines = open(adress+oldfile).readlines()
    #output new file name 
    newfile = 'data.elba_glu_sca_neutral_z'+str(zval)+'_custom'    
    #6173-6179 to change for lammps asp datafile  which is lines 6172-178 in here (python counts from zero)
    #write the start of the file     
    open(adress+newfile,'w').writelines(lines[0:6172])
    # These are the correct charges values to give interger charge
    mycharges =[-0.114411, 0.041715, 0.041715, 0.041715, -0.031111, 0.021963, 0.021963, 0.724490, -0.621691, 0.411711, -0.538059]
    #edit the resid 129 charge section    #resid 129 is lines 6173-6182 in data file 
    # below 6172 is used in loop because python counts from zero, 6183 is used because range(n) terminates at n-1  
    for n, i in enumerate(range(6172, 6183)):
        print(lines[i])
        print(n, mycharges[n])
        line = lines[i].split()
        line[7] = str(mycharges[n])   #7th element of line is the charge 
        custline = '  '+"   ".join(line)
        open(adress+newfile, 'a').writelines(custline+'\n')
    #write the end of the file 
    open(adress+newfile,'a').writelines(lines[6183:-1])
    return
