def diheds_zcompare(molecule, dhstore, zran, qlist, savedir, inout, save, angrange=np.arange(-180, 190, 10)):
    """ Function to take dihedral store and plot each dihedral as function of z"""
    zstore = []
    #loop each dihedral
    for dh in np.arange(np.shape(dhstore)[1]):
        plt.figure()
        title = 'Z compare ' + molecule+' '+inout+' '+' dihedral: '+str(qlist[dh])
        plt.title(title)
        for z in zran:
            freqin = np.histogram(dhstore[z][dh], angrange)[0]
            plt.plot(range(-175, 185, 10), freqin, label=str(z))
        plt.legend(loc='best')
        if save is True:
            plt.savefig(savedir+'/'+title)
