def Hbonds(molecule,zval,atomid1,atomid2,inout,twotb=False,save=False,plotting=False):
    """
    Get internal H bonds between atomid1 = bynum of polar h, atomid2 = bynum of heavy atom
    """

    #Get the trajectory into universe object
    if twotb == True:
        u = newUin(molecule, inout, zval)
    if twotb == False:
        u = readinU("/home/kselvon/Downloads/directedruns/", molecule, zval, inout)

    freqstore = np.zeros(69)
    count = 0
    pos_array1 = np.zeros((len(u.trajectory), 3))
    pos_array2 = np.zeros((len(u.trajectory), 3))
    r = np.zeros((len(u.trajectory), 1))
    atom1 = u.select_atoms(atomid1)
    atom2 = u.select_atoms(atomid2)
    print('getting atom positions')
    for n, ts in enumerate(u.trajectory):

        pos_array1[n] = atom1.get_positions()
        pos_array2[n] = atom2.get_positions()

    xyzdiff = pos_array1 - pos_array2
    xyzdiff = xyzdiff**2

    for i in range(len(u.trajectory)):
        r[i] = sum(xyzdiff[i])
    r = r**0.5
    freq, bin_edges = np.histogram(r, bins=np.arange(0,30.25,0.25))
    #normalise by traj length 
    freq = freq/np.float(len(u.trajectory))
    #count number of bonds within the H bond length 3.25 
    bondval  = sum(r <= 3.25)/float(len(u.trajectory))
    return bondval

def graph_Hbonds(molecule,id1,id2,inout,zrange,twotb=False,saverdf=False,save=False,plot=False):
    """  plot number of H bonds as function of Z """
    bvalstore = []
    for i in zrange:
        print('on z val ',i)
        bvalstore.append(Hbonds_new(molecule, i, id1, id2, inout, twotb, saverdf, plotting=plot))
    plt.figure()
    plt.plot(bvalstore)
    plt.ylabel('Intramolecular hydrogen bonds')
    plt.xlabel('z Angstrom')
    plt.title('Internal H bonding between '+id1+' and '+id2)
    if save == True:
        plt.savefig('/media/kselvon/D20492AA049290D9/Iridis-runs/E-decomposed/'+molecule+'/InternalHbonds_'+id1+'_'+id2)
    return bvalstore

def get_hvy_indicies(molecule, typelist, inout, twotb, zval=0):
    """Get the heavy atom indicies for Hbonding functions"""    
    if twotb == True:
        u = newUin(molecule, inout, zval)
    if twotb == False:
        u = readinU("/home/kselvon/Downloads/directedruns/", molecule , zval, inout)

    solute = u.select_atoms('resid 129')
    indlist = [i for i, x in enumerate(solute.names) if x in typelist]

    hvylist = []
    namelist = []
    for i in indlist:
        hvylist.append(solute.indices[i]+1) #the plus 1 is needed to get correct indexing sol.indicies retuns a number 1 too small
        namelist.append(solute.names[i])
    return hvylist, namelist

def internal_hydrogen_bonding(molecule, inout, dnr, polarh, custnames=False, indexs=[], twotb=False, savehbonds=False):
    """Wrapper function to call all necesary functions for intH bonding """
    if custnames == False:
        hatomlist, names = get_hvy_indicies(molecule, dnr, inout, twotb, zval=0)
    if custnames == True:
        hatomlist = indexs    
    mydict = {}
    for n, j in enumerate(polarh):
        for m, i in enumerate(hatomlist):
             print('on ', n, ' of ', len(polarh), 'polarh, on ', m,' of ', len(hatomlist),' hatoms' )
             mydict['atoms '+str(i)+' '+str(j)] = graph_Hbonds(molecule, 'bynum '+str(i), 'bynum '+str(j), inout, range(31), twotb, saverdf=False, save=savehbonds ,plot=False)
    return mydict

def intH_inout_compare_plotter(molecule, dictin, dictout, typedict, save=False):
    """ plot intH for a each heavy atom, on seperate graphs comparing in and out"""
    inkeys = dictin.keys()
    outkeys = dictout.keys()
    assert (inkeys == outkeys), 'keys from in and out dictionaries differ'
    for key in inkeys:
        title = 'Int H bonding'+' '+molecule+' for: '+str(typedict[key[:10]]+' h'+key[11:15])
        plt.figure()
        plt.plot(dictin[key],label='in')
        plt.plot(dictout[key],label='out')
        plt.title(title)
        plt.legend(loc='best')
        plt.xlabel('z depth Angstom')
        plt.ylabel('H bonds')
        if save == True:
            savedir = '/home/kselvon/Downloads/directedruns/'+molecule+'/intH'
            plt.savefig(savedir+'/'+title)
