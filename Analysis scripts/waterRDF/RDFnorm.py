def watRDF_normalised(molecule, moldnr, moltype, inout, zval=0, save=False):
    """
    edited version of watRDF to get water RDFs + normalised for MD simulations 
    """
    #Get the trajectory into universe object
    u = readinU("/home/kselvon/Downloads/directedruns/", molecule, zval, inout)

    idlist = get_hvy_indicies(molecule, moldnr, 'in')[0]
    plt.figure()
    for hvyatom in idlist:
        atom1 = u.select_atoms('bynum '+str(hvyatom))
        water = u.select_atoms('resid 0')
        pos_array1 = np.zeros((len(u.trajectory), 3))
        pos_array2 = np.zeros((len(u.trajectory), len(water), 3))
        r = np.zeros((len(u.trajectory), len(water)))
        xyzdiff = np.zeros_like(pos_array2)

        print('getting atom positions')
        for n, ts in enumerate(u.trajectory):
            pos_array1[n] = atom1.get_positions()
            pos_array2[n] = water.get_positions()
        for i in range(len(u.trajectory)):
            xyzdiff[i] = pos_array2[i][:] - pos_array1[i]
        print('xyz shape ', np.shape(xyzdiff))
        xyzdiff = xyzdiff**2
        for j in range(len(u.trajectory)):
            r[j] = np.sum(xyzdiff[j],axis=1)
        r = r**0.5
        r = np.ndarray.flatten(r)
        hist, binedges = np.histogram(r, bins=range(31))
        hist = hist/np.float(len(u.trajectory))
        for n, i in enumerate(hist):
            hist[n] = i/(4/3. * np.pi * (n+1)**3)
        print('histshape ', np.shape(hist))
        plt.plot(hist, label=str(moltype['atoms '+str(hvyatom)]))
    plt.title('Water count at zval: '+str(zval)+' '+inout)
    plt.legend(loc='best')
    plt.xlabel('r')
    plt.ylabel('Waters')
    plt.savefig('/home/kselvon/Downloads/directedruns/'+molecule+'/waterRDF/'+molecule+' zval: '+str(zval)+' normalised'+' '+inout)
    return
