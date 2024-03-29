
def dihedral2(p):
    """(2nd version of function with imporved speed) compute dihedral angle from input quartet
    returns dihedral in degrees """
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2( y, x ))


def PBCFIX_dihedralplusangles(molecule,quartetlist,zval): 
    """Function to perform analysis of dihedrals during simulations with PBXFIX applied """
    u_in = readinU('/media/kselvon/7A82B43582B3F42B/iridisruns2019/trajectoriesForMD/'
, molecule, zval, 'in')
    u_out = readinU('/media/kselvon/7A82B43582B3F42B/iridisruns2019/trajectoriesForMD/'
, molecule, zval, 'out')
    
    #uncomment this line to check in/out length consistency 
    #assert (len(u_in.trajectory)==len(u_out.trajectory)), 'In/out trajectories different length'
    framelist_in = []
    framelist_out= []
    globaldhstore_in = []
    globaldhstore_out= []
    #Get frame by frame dimensions
    dims_in = np.zeros((len(u_in.trajectory), 3))
    dims_out = np.zeros((len(u_out.trajectory), 3))
    for ts in u_in.trajectory:
        dims_in[ts.frame] = u_in.dimensions[0:3]
    for ts in u_out.trajectory:
        dims_out[ts.frame] = u_out.dimensions[0:3]
    drug_in = u_in.select_atoms('resid 129')
    drug_out = u_out.select_atoms('resid 129')
    for ts in u_in.trajectory:
        #ref atom xyz
        refatom_in =  drug_in.get_positions()[0]
        #new position array
        newpos_in = np.zeros_like(drug_in.get_positions())
        #loop through all other atoms in drug
        for n,i in enumerate(drug_in.get_positions()):
            vector_in = refatom_in - i
            vector_in = vector_in - dims_in[ts.frame]*np.rint(vector_in/dims_in[ts.frame])
            newpos_in[n] = refatom_in + vector_in
            drug_in.set_positions(newpos_in)
        framelist_in.append(newpos_in)
    for ts in u_out.trajectory:
        refatom_out= drug_out.get_positions()[0]
        newpos_out = np.zeros_like(drug_out.get_positions())
        for n,i in enumerate(drug_out.get_positions()):
            vector_out = refatom_out - i
            vector_out = vector_out - dims_out[ts.frame]*np.rint(vector_out/dims_out[ts.frame])
            newpos_out[n] = refatom_out + vector_out
            drug_out.set_positions(newpos_out)
        framelist_out.append(newpos_out)
    #number of dihedrals
    dhnum = np.shape(quartetlist)[0]
    #loop through each quartet
    for i in range(dhnum):
        dhatom1in = []
        dhatom2in = []
        dhatom3in = []
        dhatom4in = []
        dhatom1out = []
        dhatom2out = []
        dhatom3out = []
        dhatom4out = []
        #loop through the framelist of coordinates to get dihedral atom coordinates
        for frame in range(len(framelist_in)):
            dhatom1in.append(framelist_in[frame][quartetlist[i][0]-1])  
            dhatom2in.append(framelist_in[frame][quartetlist[i][1]-1])   #need the -1 because pdb start at ATOM 1, python array starts at zero
            dhatom3in.append(framelist_in[frame][quartetlist[i][2]-1])
            dhatom4in.append(framelist_in[frame][quartetlist[i][3]-1])
        for frame in range(len(framelist_out)):
            dhatom1out.append(framelist_out[frame][quartetlist[i][0]-1])
            dhatom2out.append(framelist_out[frame][quartetlist[i][1]-1])
            dhatom3out.append(framelist_out[frame][quartetlist[i][2]-1])
            dhatom4out.append(framelist_out[frame][quartetlist[i][3]-1])
        dhstore_in =[]
        dhstore_out =[]
        for j in range(len(dhatom1in)):
            dhstore_in.append(dihedral2(np.vstack((dhatom1in[j], dhatom2in[j], dhatom3in[j], dhatom4in[j]))))
        for j in range(len(dhatom1out)):
            dhstore_out.append(dihedral2(np.vstack((dhatom1out[j], dhatom2out[j], dhatom3out[j], dhatom4out[j]))))

        globaldhstore_in.append(dhstore_in)
        globaldhstore_out.append(dhstore_out)
        
        #to normalise by framelength
        normalised_gdhs_in = []
        normalised_gdhs_out = []
         
        #function does not currently return these values, use for debugging 
        for i in globaldhstore_in:
            normalised_gdhs_in.append(np.array(i)/np.float(len(u_in.trajectory)))
        for i in globaldhstore_out:
            normalised_gdhs_out.append(np.array(i)/np.float(len(u_out.trajectory)))
return globaldhstore_in,globaldhstore_out 

def plotgauch(instore,outstore,dhlist,molecule,save,shift=False):
    """Function to plot dihedral angles as a function of z depth from PBCFIX_dihedralplusangles output"""
    for dh in range(len(dhlist)):
        if shift == True:
            freqin = np.histogram(instore[dh], np.arange(0, 361, 1))[0]
            freqout = np.histogram(outstore[dh], np.arange(0, 361, 1))[0]
        else:
            freqin = np.histogram(instore[dh], np.arange(-180, 181, 1))[0]
            freqout = np.histogram(outstore[dh], np.arange(-180, 181, 1))[0]

        plt.figure()
        if shift == True:
            title = molecule+' SHIFTED'+str(dhlist[dh])
            plt.title(title)
        else:
            title = molecule+' '+str(dhlist[dh])
            plt.title(title)
        if shift == True:
            plt.plot(np.arange(0, 360, 1), freqin, label='in')
            plt.plot(np.arange(0, 360, 1), freqout, label='out')
        else:
            plt.plot(np.arange(-180, 180, 1), freqin, label='in')
            plt.plot(np.arange(-180, 180, 1), freqout, label='out')
        plt.xlabel('angle in degrees')
        plt.ylabel('freqency')
        plt.legend(loc='best')
        if save == True:
            plt.savefig('/home/kselvon/Downloads/directedruns/'+molecule+'/Gauchness/'+title)
        plt.show()

def extractdhacrossz(dhstore, zran):
    """ Function to extract values from PBCFIX_dihedralplusangles output  """ 
    dhnum = np.shape(dhstore)[1]
    dh_overz = []
    for dh in range(dhnum):
        temp=[]
        for z in zran:
            for i in range(len(dhstore[z][dh])-1):
                temp.append(dhstore[z][dh][i])
        dh_overz.append(temp)
    return dh_overz
