def Orientation_spherical2(molecule, bynum, selectedzeds=range(31), saveplots=False, twotb=False):
    """ Wrapper function to perform orientational analysis of molecule in simulation
        2nd version of function with modifications to the order parameter produced 
    """ 

    #set up save folder
    if saveplots is True:
        if not os.path.exists('/home/kselvon/Downloads/directedruns/'+molecule):
            os.mkdir('/home/kselvon/Downloads/directedruns/'+molecule)
        if not os.path.exists('/home/kselvon/Downloads/directedruns/'+molecule+'/OrientationPHI2019'):
            os.mkdir('/home/kselvon/Downloads/directedruns/'+molecule+'/OrientationPHI2019')

    savepath = '/home/kselvon/Downloads/directedruns/'+molecule+'/OrientationPHI/'

    rin, thetain, phiin, rout, thetaout, phiout = Axis_spherical(molecule, bynum, selectedzeds, twotb)

    #plot theta in (av cos2 as in legendre polynomial ) as function of z? done?
    P2phiin_avstore = []
    P2phiout_avstore = []
    Cosphiin_avstore = []
    Cosphiout_avstore = []
    for i in selectedzeds:
        #S order parameter of phi
        P2phiin_avstore.append(np.average(0.5*(3*np.cos(np.deg2rad(phiin[i]))**2 -1)))
        P2phiout_avstore.append(np.average(0.5*(3*np.cos(np.deg2rad(phiout[i]))**2 -1)))
        # cosine of phi
        Cosphiin_avstore.append(np.average(np.cos(np.deg2rad(phiin[i]))))
        Cosphiout_avstore.append(np.average(np.cos(np.deg2rad(phiout[i]))))
    #plot s order for in and out
    plt.figure()
    title1 = molecule+' S order parameter of vector '+str(bynum)
    plt.title(title1)
    plt.plot(P2phiin_avstore, label='in')
    plt.plot(P2phiout_avstore, label='out')
    plt.xlabel('Z val')
    plt.ylabel('Average S order of phi')
    plt.legend(loc='best')
    if saveplots ==True:
        plt.savefig(savepath+title1)
    #plot cos for in and out
    plt.figure()
    title2 = molecule+' Cos of vector '+str(bynum)
    plt.title(title2)
    plt.plot(Cosphiin_avstore, label='in')
    plt.xlabel('Z val')
    plt.ylabel('Average cosine of phi')
    plt.plot(Cosphiout_avstore, label='out')
    plt.legend(loc='best')
    if saveplots is True:
        plt.savefig(savepath+title2)
    return 

def Axis_spherical(molecule, bynum_vec, selectedzeds, twotb=False):
    """Function to get spherical coordinates of molecules as defined by a vector across the molecule """
    theta_instore = []
    phi_instore = []
    r_instore = []
    theta_outstore = []
    phi_outstore = []
    r_outstore = []
    for zval in selectedzeds:
        print('on z ', zval)
        atom1in, atom2in, frlenin = Get_atomcoords(molecule, zval, 'in', bynum_vec, twotb)
        atom1out, atom2out, frlenout = Get_atomcoords(molecule, zval, 'out', bynum_vec, twotb)
        #shift the vector to touch the origin
        shifted_atom2in = translate2origin(atom1in, atom2in)
        shifted_atom2out = translate2origin(atom1out, atom2out)
        r_in, theta_in, phi_in = cart2polar(shifted_atom2in)
        r_out, theta_out, phi_out = cart2polar(shifted_atom2out)
        theta_instore.append(theta_in)
        phi_instore.append(phi_in)
        r_instore.append(r_in)
        theta_outstore.append(theta_out)
        phi_outstore.append(phi_out)
        r_outstore.append(r_out)
    return r_instore, theta_instore, phi_instore, r_outstore, theta_outstore, phi_outstore
    
def Get_atomcoords(molecule, zval, inout, bynum_vec, twotb):
    """ Function to get cartesian coords of molecules during simulation"""
    if twotb is False:
        u = readinU("/home/kselvon/Downloads/directedruns/", molecule, zval, inout)
    if twotb is True:
        u = newUin(molecule, inout, zval)

    #define the two atoms forming the molecule vector
    atom1 = u.select_atoms('bynum '+bynum_vec[0])
    atom2 = u.select_atoms('bynum '+bynum_vec[1])
    #stores for dimensions and atom coords
    dims = np.zeros((len(u.trajectory), 3))
    atom1store = np.zeros_like(dims)
    atom2store = np.zeros_like(dims)


    #get the dimensions, correct the positions and store the coords
    for ts in u.trajectory:
        dims[ts.frame] = u.dimensions[0:3]

    drug = u.select_atoms('resid 129')
    for ts in u.trajectory:
        #ref atom xyz
        refatom=  drug.get_positions()[0]
        #new position array
        newpos = np.zeros_like(drug.get_positions())
        #loop through all other atoms in drug
        for n, i in enumerate(drug.get_positions()):
            vector = refatom - i
            vector = vector - dims[ts.frame]*np.rint(vector/dims[ts.frame])
            newpos[n] = refatom + vector
            drug.set_positions(newpos)
        #Below stores the bynum_vec atom coordinates
        atom1store[ts.frame] = atom1.positions
        atom2store[ts.frame] = atom2.positions
    return atom1store, atom2store, len(u.trajectory)
