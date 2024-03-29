    
def nudgeXY(zval):
    
    """
    function to nudge solvent in the xy plane by modifying the coordinates, not yet generalised
    need to edit adresses to use. 
    
    """
    #address and name of old datafiles 
    adress='/home/kselvon/peturbdata/'
    oldfile = 'data.elba_glu_sca_neutral_z'+str(zval)+'_custom'
    lines = open(adress+oldfile).readlines()
    
    #output new file name EDIT FOR WHATS NEEDED
    newfile = 'data.elba_glu_sca_neutral_z'+str(zval)+'_custom_XYedit'    
    
    
    #6173-6179 to change for lammps asp datafile  which is lines 6172-178 in here (python counts from zero)

    #write the start of the file     
    open(adress+newfile, 'w').writelines(lines[0:6172])

    
    #edit the resid 129 charge section    #resid 129 is lines 6173-6182 in data file 
    # below 6172 is used in loop because python counts from zero, 6183 is used because range(n) terminates at n-1  
    for n, i in enumerate(range(6172, 6183)):
        line = lines[i].split()
        #edit 2nd and 3rd element of line corrosponding to x y coordinate 
        line[2] = str(float(line[2])+1)
        line[3] = str(float(line[3])+1)  
        custline = '  '+"   ".join(line)
        open(adress+newfile, 'a').writelines(custline+'\n')
    
    
    #write the end of the file 
    open(adress+newfile, 'a').writelines(lines[6183:10082])
    return 

