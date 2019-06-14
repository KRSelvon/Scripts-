def makecols(table, index, y_index, names=0):
    """Function to make columns for catagorical analysis table """
    blankcol = []
    ycol = []
    blankdict = {}
    blankydict = {}
    for i in range(len(table)):
        print('data being read into returned colx')
        if table[i][index] != 'nan':
            if table[i][y_index] != 'nan':
                blankcol.append(round(np.float(table[i][index]), 47))
                ycol.append(table[i][y_index])
                print(round(np.float(table[i][index]), 46), table[i][0])
                blankdict[round(np.float(table[i][index]), 46)] = table[i][names]
                blankydict[round(np.float(table[i][y_index]), 10)] = table[i][names]
        print('rankings of xcol...index:', index)
    for n, k in enumerate(np.sort(blankdict.keys())):
        print('RANK ', n, ' ', blankdict[k])
    print('rankings of ycol....index:', y_index)    
    for n, j in enumerate(np.sort(blankydict.keys())):
        print('RANK ', n, ' ', blankydict[j]) 
    xs = np.zeros(len(blankcol))
    ys = np.zeros(len(blankcol))
    for i in range(len(blankcol)):
        xs[i] = np.float(blankcol[i])
        ys[i] = np.float(ycol[i])
    return xs, ys, blankdict, blankydict
    
def spearmans_rank_correlation(xs, ys):
    """ Function to compute spearmans rank correlation for input array of xs and ys"""
    # Calculate the rank of x's
    xranks = pd.Series(xs).rank()
    # Caclulate the ranking of the y's
    yranks = pd.Series(ys).rank()
    # Calculate Pearson's correlation coefficient on the ranked versions of the data
    return scipy.stats.pearsonr(xranks, yranks)   


def LOBFbasic(xs, ys, table, dic, title='', xlab='', ylab='', save=False):
    """get line of best fit and statistics for a set of xy coordinates
       basic version of LOBFcatagorical function"""
    m, c, r_val, p_val, stder = stats.linregress(xs, ys)
    colordic = {'a': 'r', 'b': 'b', 'n': 'g'}
    plt.figure()    
    for i in range(len(xs)): 
        plt.scatter(xs[i], ys[i], color=colordic[table[i][27]]) # pass acid base col for colour
        #plt.annotate(dic[xs[i]][1:5],(xs[i],ys[i]))
    if title == '':
       title = ylab+' vs '+xlab    ## CURRENTLY OVERRIDING FUNCTION ARG ABOVE
    def f(x, gradient, intercept):
        return gradient*x + intercept
    lobf = f(xs, m, c)
    plt.plot(xs, lobf, label='P :'+str(round(p_val, 3))+'\n'+r'$r^{2} :$'+str(round(r_val**2, 2)))

    plt.legend(loc='best')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)

    return m, c, r_val, p_val, stder
    
def LOBFcatagorical(xs, ys, table, dic, save, title='', xlab='', ylab=''):
    """get line of best fit and statistics for a set of xy coordinates, perform catagorical analysis"""
    m, c, r_val, p_val, stder = stats.linregress(xs, ys)
    colordic = {'a':'r','b':'b','n':'g'}

    plt.figure(figsize=(7.5, 6.5))
    for i in range(len(xs)):
        plt.scatter(xs[i], ys[i], color=colordic[table[i][27]]) #,marker=table[i][30]) # pass acid base col for colour
        #plt.annotate(dic[xs[i]][1:5],(xs[i],ys[i]))
    title = xlab+' vs '+ylab    ## CURRENTLY OVERRIDING FUNCTION ARG ABOVE
    def f(x, gradient, intercept):
        return gradient*x + intercept
    lobf = f(xs, m, c)
    plt.plot(xs, lobf, label='P :'+str(p_val)+'\n'+'r^2 :'+str(r_val**2))
    
    
    if xlab[0:3] == 'Log':
            plt.plot([np.log10(1)]*len(ys), ys, color='k')
            plt.plot([np.log10(5)]*len(ys), ys, color='k')
    else:
        plt.plot([1]*len(ys), ys, color='k')
        plt.plot([5]*len(ys), ys, color='k')
    
    if ylab[0:3] =='Log':
        plt.plot(xs, [m*np.log10(1) +c]*len(xs), color='k')
        plt.plot(xs, [m*np.log10(5) +c]*len(xs), color='k')
    else:
        plt.plot(xs, [m*1 +c]*len(xs), color='k')
        plt.plot(xs, [m*5 +c]*len(xs), color='k')
        
    
    plt.legend(loc='best')
    #plt.legend(bbox_to_anchor=(0.8, -0.13))
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    #plt.ylim(-1.9,-0.5)
    plt.title(title)
    if save == True:
        plt.savefig('/home/kselvon/Desktop/AZcorrelations_groupcolour/'+title, bbox_extra_artists=(0.8, -0.13))
    print('r :',r_val, 'r^2: ',r_val**2)
    print('p :',p_val)
    print('m',m, 'c ',c)
    return m, c, r_val, p_val, stder
    
    
    
def wham_data(molecule, save=False):
    """ Hacky unfinished Function to extract PR1 and PR2 from drug PMFS"""
    if save==True:
       savedir2 = '/home/kselvon/Downloads/directedruns/'+molecule
       if not os.path.exists(savedir2):
          os.mkdir(savedir2)

    if save==True:
       savedir = '/home/kselvon/Downloads/directedruns/'+molecule+'/PMF'
       if not os.path.exists(savedir):
          os.mkdir(savedir)

    #read in the data
    adress_in = '/media/kselvon/D20492AA049290D9/Iridis-runs/directedruns/'+molecule+'/in/'
    adress_out = '/media/kselvon/D20492AA049290D9/Iridis-runs/directedruns/'+molecule+'/out/'
    data_in = np.genfromtxt(adress_in+'wham.pmf')
    data_out = np.genfromtxt(adress_out+'wham.pmf')

    #normalise to be zero at water boundary
    norm_in = normalisedata(data_in)
    norm_out = normalisedata(data_out)

    #get an interpolated average
    xav,yav = interpolatedaverage(norm_in,norm_out)

    #get the in/out average pr1 and pr2 
    PR1av = PR1(xav,yav)
    PR2av = PR2(xav,yav)
    PR2av_full = PR2_KBT(xav,yav)
    print('PR2av: ',PR2av,' PR2av_full: ',PR2av_full)
    
    #get the in/out PR1     
    PR1in = PR1(norm_in[:,0], norm_in[:,1])
    PR1out = PR1(norm_out[:,0], norm_out[:,1])

    #get the in/out PR2
    PR2in = PR2(norm_in[:,0], norm_in[:,1])
    PR2out = PR2(norm_out[:,0], norm_out[:,1])

    PR2in_full = PR2_KBT(norm_in[:,0], norm_in[:,1])
    PR2out_full = PR2_KBT(norm_out[:,0], norm_out[:,1])    
    print('PR2in: ', PR2in,' PR2in_full: ', PR2in_full,' PR2out: ', PR2out, ' PR2out_full: ', PR2out_full)    
    
    #get the TFE at z0 using from scipy import interpolate
    
    fin = interpolate.interp1d(norm_in[:,0], norm_in[:,1]) 
    TFEin = fin(0)
    
    fout = interpolate.interp1d(norm_out[:,0], norm_out[:,1]) 
    TFEout = fout(0)
    
    fav = interpolate.interp1d(xav, yav)
    TFEav = fav(0)
    
    #plot raw data with error bars
    plt.figure()
    title2 = molecule+' PMF DOPC'
    plt.title(title2)
    plt.ylabel('kJ/mol')  #kK jJ
    plt.xlabel('Z depth in Angstrom from membrane core')
    plt.errorbar(norm_in[:,0], norm_in[:,1], norm_in[:,2], errorevery=5, label='in')
    plt.errorbar(norm_out[:,0], norm_out[:,1], norm_out[:,2], errorevery=5, label='out')
    plt.xlim(-0.2, 30)
    plt.legend(loc='best')
    if save == True:
        plt.savefig(savedir+'/'+title2)
        plt.savefig('/home/kselvon/Desktop/betablockers/PMFS-all/'+title2)
    #plot raw data no error bars along with the interpolated average
    plt.figure()
    plt.xlim(0, 30)
    plt.plot(xav, yav, label='average PR1:'+str(PR1av)[0:5]+' PR2:'+str(PR2av)[0:5])
    #plt.plot(norm_in[:,0],data_in[:,1],label='in')
    #plt.plot(norm_out[:,0],data_out[:,1],label='out')
    title = molecule+' interpolated inout average PMF DOPC'
    plt.title(title)
    plt.ylabel('kJ/mol')
    plt.xlabel('Z depth in Angstrom from membrane core')
    plt.legend(loc='best')
    if save==True:
        plt.savefig(savedir+'/'+title)
        plt.savefig('/home/kselvon/Desktop/betablockers/PMFS-all/'+title)
    return PR2out_full

