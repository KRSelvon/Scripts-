def thermocycle(molecule1,molecule2,savename,pka,save=False):
    """molecule2 should be the neutral species"""
    
    #read in acid
    datain = np.genfromtxt('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule1+'/wham.pmf')
    normalisedin =  normalisedata(datain)
    #read in neutral    
    dataout = np.genfromtxt('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule2+'/wham.pmf')
    normalisedout =  normalisedata(dataout)
    
    #pka correction on neutral 
    normalisedout[:,1] = normalisedout[:,1] - pKa_correction(pka)
        
    #interpolate the PMFS     
    acid = interpolate.interp1d(normalisedin[:,0],normalisedin[:,1]) 
    neutral = interpolate.interp1d(normalisedout[:,0],normalisedout[:,1]) 
    
    
    for i in range(31):
        neutral(i) + 
    
    
    
    acid(0)    
    neutral(0)   
   
   
   
    plt.figure()
    title2 = savename+' PMF DOPC'
    plt.title(title2)
    plt.ylabel('kJ/mol')
    plt.xlabel('Z depth in Angstrom from membrane core')
    plt.errorbar(normalisedin[:,0],normalisedin[:,1],normalisedin[:,2],errorevery=5,label=molecule1)
    #plt.errorbar(dataout[:,0],dataout[:,1],dataout[:,2],errorevery=5,label=molecule2)    
    plt.errorbar(normalisedout[:,0],normalisedout[:,1],normalisedout[:,2],errorevery=5,label=molecule2)
    plt.legend(loc='best')
    plt.xlim(-0.2,30)
    if save == True:
        plt.savefig('/media/kselvon/7A82B43582B3F42B/wham2019/'+savename+'PMF')
 
 
 
 
 
def Goffset(pka):
    R = 8.31445e-3
    T=298.0
    ph=7
    G = (pka -ph)/ (2.303*R*T)
    return G
    
 

def plotbothpmf(molecule,save=False):
    datain = np.genfromtxt('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule+'/in/wham.pmf')
    normalisedin =  normalisedata(datain)
    dataout = np.genfromtxt('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule+'/out/wham.pmf')
    normalisedout =  normalisedata(dataout)
   
   
    plt.figure()
    title2 = molecule+' PMF DOPC'
    plt.title(title2)
    plt.ylabel('kJ/mol')
    plt.xlabel('Z depth in Angstrom from membrane core')
    plt.errorbar(normalisedin[:,0],normalisedin[:,1],normalisedin[:,2],errorevery=5,label='in')
    plt.errorbar(normalisedout[:,0],normalisedout[:,1],normalisedout[:,2],errorevery=5,label='out')
    plt.legend(loc='best')
    plt.xlim(-0.2,30)
    if save == True:
        plt.savefig('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule+'/'+molecule+'PMF')

def plotbothpmf_amino(molecule1,molecule2,savename,pka=0,pkacorrect=False,save=False):
    """molecule2 should be the neutral species"""
    
    datain = np.genfromtxt('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule1+'/wham.pmf')
    normalisedin =  normalisedata(datain)
    dataout = np.genfromtxt('/media/kselvon/7A82B43582B3F42B/wham2019/'+molecule2+'/wham.pmf')
    normalisedout =  normalisedata(dataout)
    
    if pkacorrect==True:
        normalisedout[:,1] = normalisedout[:,1] - pKa_correction(pka)
        #dataout[:,1] = dataout[:,1]+ pKa_correction(pka)
        
   
   
    plt.figure()
    title2 = savename+' PMF DOPC'
    plt.title(title2)
    plt.ylabel('kJ/mol')
    plt.xlabel('Z depth in Angstrom from membrane core')
    plt.errorbar(normalisedin[:,0],normalisedin[:,1],normalisedin[:,2],errorevery=5,label=molecule1)
    #plt.errorbar(dataout[:,0],dataout[:,1],dataout[:,2],errorevery=5,label=molecule2)    
    plt.errorbar(normalisedout[:,0],normalisedout[:,1],normalisedout[:,2],errorevery=5,label=molecule2)
    plt.legend(loc='best')
    plt.xlim(-0.2,30)
    if save == True:
        plt.savefig('/media/kselvon/7A82B43582B3F42B/wham2019/'+savename+'PMF')


def pKa_correction(pKa,T=298):
    
    #henderson hasslebalch equation ph = pka +log10(Ka)    
    # Ka =  10^(ph - pKa) 
    ph = 7
    Ka = 10**(ph-np.float(pKa))
    
    #deltaG = -RT lnKeq   R has units of kJ  K−1 mol−1  so G has KJ mol-1 
    R = 8.31445e-3  
    #T= 298
    G = -R * T * np.log(Ka)
    
    return G 
