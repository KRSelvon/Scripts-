

def trajlencheck(molecule):
    """ Function to check the trajectory length of simulations"""
    for i in range(31):
        uin = newUin(molecule, 'in', i)    
        uout = newUin(molecule, 'out', i)    
        lenin= len(uin.trajectory)
        lenout =len(uout.trajectory)
        print('for zval: '+str(i))
        print('in is '+str(lenin)+' frames long')
        print('out is '+str(lenout)+' frames long')
