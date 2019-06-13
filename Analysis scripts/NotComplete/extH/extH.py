
def new_watHbonds(molecule,zval,atomid1,inout,save=False,plotting=False,plotfbf=False):
    """
    Get internal H bonds between atomid1 = bynum of polar h, atomid2 = bynum of heavy atom
    """
    #Get the trajectory into universe object
    u = newUin(molecule,inout,zval)





    atom1 = u.select_atoms(atomid1)
    water = u.select_atoms('resid 0')
    pos_array1 = np.zeros((len(u.trajectory),3))
    pos_array2 = np.zeros((len(u.trajectory),len(water),3))
    r = np.zeros((len(u.trajectory),len(water)))
    xyzdiff = np.zeros_like(pos_array2)



    print('getting atom positions')
    for n,ts in enumerate(u.trajectory):

        pos_array1[n] = atom1.get_positions()
        pos_array2[n] = water.get_positions()

    for i in range(len(u.trajectory)):
        xyzdiff[i] = pos_array2[i][:] - pos_array1[i]



    xyzdiff = xyzdiff**2


    for j in range(len(u.trajectory)):
        r[j] = np.sum(xyzdiff[j],axis=1)

    r = r**0.5



    fBfHbonds = framebyframeHbond(r)#,len(u.trajectory))
    #print(fBfHbonds)
#    if molecule == 'piroxicam' and zval ==7:
#        bondval = sum(fBfHbonds[:-2]) / np.float(len(u.trajectory)-1)
    if molecule == 'acyclovir' and zval == 7:
        bondval = sum(fBfHbonds[:-1]) / np.float(len(u.trajectory)-1)
    if molecule == 'pindolol' and zval == 2:
        bondval = sum(fBfHbonds[:-1]) / np.float(len(u.trajectory)-1)
    if molecule == 'pindolol' and zval == 6:
        bondval = sum(fBfHbonds[:-1]) / np.float(len(u.trajectory)-1)
    if molecule == 'pindolol' and zval == 22:
        bondval = sum(fBfHbonds[:-1]) / np.float(len(u.trajectory)-1)
    else:
        bondval = sum(fBfHbonds[:-1]) / np.float(len(u.trajectory)-1)


    #return r, bondval

    #r= r.flatten()
    #print('lenth of r: ',len(r))


    freq,bin_edges = np.histogram(r,bins=np.arange(0,30.25,0.25))

    freq = freq/np.float(len(u.trajectory))
    if plotting == True:
        print('plotting')
        plt.figure()
        plt.plot(bin_edges[0:len(bin_edges)-1],freq)
        plt.xlabel('z Angstrom')
        plt.ylabel('Frequency')
        title='External H bonds '+atomid1+ ' '+str(zval)
        plt.title(title)
        if save == True:
            plt.savefig('/media/kselvon/D20492AA049290D9/Iridis-runs/E-decomposed/'+molecule+'/'+title)

    #bondval = np.trapz(freq[0:13],bin_edges[0:13],dx=0.25)
    #bondval  = sum(r <= 3.25)/float(len(u.trajectory))

    return bondval#,fBfHbonds


def new_watgraph_Hbonds(molecule,id1,zrange,inout,saverdf=False,save=False,plot=False):

    bvalstore = []
    for i in zrange:
        print('on z val ',i)
        bvalstore.append(new_watHbonds(molecule,i,id1,inout,saverdf,plotting=plot))
    plt.figure()
    plt.plot(bvalstore)
    plt.ylabel('External hydrogen bonds')
    plt.xlabel('z Angstrom')
    plt.title('External H bonding '+id1)
    if save == True:
        plt.savefig('/media/kselvon/D20492AA049290D9/Iridis-runs/E-decomposed/'+molecule+'/External H bonding '+id1)
    return bvalstore



def new_watinternal_hydrogen_bonding(molecule,inout,dnr,poalrh,savegraph=False):


    hatomlist,names = new_get_hvy_indicies(molecule,dnr,inout,zval=0)
    mydict ={}
    for i in hatomlist:
             mydict['atoms '+str(i)] = new_watgraph_Hbonds(molecule,'bynum '+str(i),range(31),inout,saverdf=False,save=savegraph,plot=False)


    return mydict

def new_watinternal_hydrogen_bonding2(molecule,inout,dnr,poalrh,savegraph=False):


    hatomlist,names = new_get_hvy_indicies(molecule,dnr,inout,zval=0)
    mydict ={}
    for i in hatomlist:
             mydict['atoms '+str(i)] = new_watgraph_Hbonds(molecule,'bynum '+str(i),range(27,28),'out',saverdf=False,save=savegraph,plot=False)


    return mydict



def new_get_hvy_indicies(molecule,typelist,inout,zval=0):
    u = newUin(molecule,inout,zval)
    solute = u.select_atoms('resid 129')
    indlist = [i for i, x in enumerate(solute.names) if x in typelist]

    hvylist =[]
    namelist =[]
    for i in indlist:
        hvylist.append(solute.indices[i]+1)#the plus 1 is needed to get correct indexing sol.indicies retuns a number 1 too small
        namelist.append(solute.names[i])
    return hvylist,namelist
