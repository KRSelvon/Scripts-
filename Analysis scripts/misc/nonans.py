def removenans(table, colindex):
    """ Function to remove nans from column in numpy table"""
    nonans = [] 
    for i in range(len(table)):
        if table[i][colindex] != 'nan':
            nonans.append(i)
    newtab = np.row_stack((table[nonans[0], :]))
    nonans.pop(0)    
    for i in nonans:
        newtab = np.row_stack((newtab[:, :],table[i, :]))
    return newtab
        
        
    
