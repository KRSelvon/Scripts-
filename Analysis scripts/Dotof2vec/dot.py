def dot(v1, v2):
    """
    compute the angle between two 3d vectors using the dot product rule.
    Returns answer in degrees.     
    """
    magv1 = (v1[0]**2 + v1[1]**2 + v1[2]**2)**0.5
    magv2 = (v2[0]**2 + v2[1]**2 + v2[2]**2)**0.5
    costheta = np.dot(v1,v2) / float((magv1*magv2)) 
    theta_rad = np.arccos(costheta)
    theta_deg = np.rad2deg(theta_rad)
    return theta_deg 

def plotangles(molecule, bynumv1, bynumv2, z, in_store, out_store, zeds, saveplots, savepath, lenin, lenout, specbins, bins):
    """ Function to plot/histogram angles as function of z from dotoftwovects function""" 
     plt.figure()
     title = molecule+' z'+str(zeds[z])+' angle between vects: '+str(bynumv1)+' and '+str(bynumv2)
     
     if specbins == False:
         freqin = np.histogram(in_store[z], np.arange(0,190,10))[0]
         freqout = np.histogram(out_store[z], np.arange(0,190,10))[0]
     if specbins == True:       
         freqin = np.histogram(in_store[z], bins)[0]
         freqout = np.histogram(out_store[z], bins)[0]
         plot_bins = np.arange(bins[1]/2., max(bins), bins[3]-bins[2]) 
     #normalise by frame length
     freqin = freqin/float(lenin[z])
     freqout = freqout /float(lenout[z])
     if specbins == False:
         plt.plot(range(5, 185, 10), freqin, label='in')
         plt.plot(range(5, 185, 10), freqout, label='out')
         plt.xticks(range(5, 185, 10))
     if specbins == True:
         plt.plot(plot_bins, freqin, label='in')
         plt.plot(plot_bins, freqout, label='out')
         plt.xticks(plot_bins)
     plt.legend(loc='best')
     plt.title(title)
     plt.xlabel('Angle in degrees')
     plt.ylabel('Frequency')
     plt.grid()
    
     if saveplots == True:
         if specbins == False:
             plt.savefig(savepath+title)
             return
         if specbins == True:
            plt.savefig(savepath+title+str(bins))
            return plot_bins
    
            
    
def dotoftwovects(molecule, bynum_vec1, bynum_vec2, selectedzeds=[0,5,10,15,20,25,30], twotb=False, saveplots=False, specbins=False, bins=[]):
    """
    get the angle between two vectors (BA and CD) as defined by 4 atoms across a molecule 
    
    A------B--[bits of molecule]----C------D
    
    
    A <------B    
    find angle between these vectors 
    C ------> D   
    
    
    vector BA defined by atoms in bynum_vec1 eg bynum_vec1[0] = index of A, bynum_vec1[1] = index of B
    vector CD defined by atoms in bynum_vec1 eg bynum_vec2[0] = index of c, bynum_vec1[1] = index of D
    
    with this convention in mind an angle of 0 degrees would be total overlap and 180 degrees total 
    unfolding of the molecule 
    
    """   
    if saveplots==True:
        if not os.path.exists('/home/kselvon/Downloads/directedruns/'+molecule):
            os.mkdir('/home/kselvon/Downloads/directedruns/'+molecule)
        if not os.path.exists('/home/kselvon/Downloads/directedruns/'+molecule+'/Vector-folding'):
            os.mkdir('/home/kselvon/Downloads/directedruns/'+molecule+'/Vector-folding')

    savepath = '/home/kselvon/Downloads/directedruns/'+molecule+'/Vector-folding/'

    zstore_in = []
    zstore_out = []      
    
    #store lengths of trajs to normalise frequnecy in histogram plots 
    lenstore_in = []
    lenstore_out = [] 

    for zval in selectedzeds:

        print('on z ', zval)
        #get atom coords for vect1 and vect2 for both in and out 
        atom1in_v1, atom2in_v1, frlenin_v1 = Get_atomcoords2(molecule, zval, 'in', bynum_vec1, twotb)
        atom1in_v2, atom2in_v2, frlenin_v2 = Get_atomcoords2(molecule, zval, 'in', bynum_vec2, twotb)
        atom1out_v1, atom2out_v1, frlenout_v1 = Get_atomcoords2(molecule, zval, 'out', bynum_vec1, twotb)
        atom1out_v2, atom2out_v2, frlenout_v2 = Get_atomcoords2(molecule, zval, 'out', bynum_vec2, twotb)
        
        #check in and out are same length, save lengths to pass to plotting function to norm
        assert (frlenin_v1==frlenin_v2), 'trajectories for vectors are different lengths'        
        assert (frlenout_v1==frlenout_v2), 'trajectories for vectors are different lengths'  
        lenstore_in.append(frlenin_v1)
        lenstore_out.append(frlenout_v1)
        if frlenout_v1 != frlenin_v1:
            print('IN/OUT have different traj lengths')
            
        #get vectors for in
        vector1_in = atom1in_v1 - atom2in_v1  # this is vector BA  <----- 
        vector2_in = atom2in_v2 - atom1in_v2  # this is vector CD  -----> 
        #get vectors for out
        vector1_out = atom1out_v1 - atom2out_v1
        vector2_out = atom2out_v2 - atom1out_v2
 
        instore = [] 
        outstore = [] 
        
        #for in 
        for i in range(frlenin_v1):
            instore.append(dot(vector1_in[i], vector2_in[i]))
        #for out 
        for i in range(frlenout_v1):
            outstore.append(dot(vector1_out[i], vector2_out[i]))
        
        zstore_in.append(instore)
        zstore_out.append(outstore)
        
    
    #graph the results
    for n, z in enumerate(selectedzeds):    
     plotangles(molecule, bynum_vec1, bynum_vec2, n, zstore_in, zstore_out, selectedzeds, saveplots, savepath, lenstore_in, lenstore_out, specbins, bins)
             
    return 
