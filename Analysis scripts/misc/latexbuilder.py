#functions to produce latex code for large sets of figures, not yet generalised 

def latexFigBuilder(filepath, molecule, analtype, txt=''):
    files = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    f = open(molecule+analtype+txt, 'w')
    f.write('\begin{figure} \n')
    f.write('\centering \n')

    for line in files:
        f.write('\hspace{0mm} \n')
        f.write('\subfloat['+line[23:37]+']{ \n')
        f.write(' '+'\includegraphics[width=75mm,height=37.5mm]{"'+line[:-4]+'".png} \n')
        f.write('  '+'\label{propgauin-} \n')
        f.write('} \n')
        print(line)
    f.write('\caption{'+molecule[0].upper()+molecule[1:]+' raw dihedral results'+'} \n')
    f.write('\end{figure}')
    f.close()
    return files

def latexFigBuilder_gauch(filepath, molecule, analtype, txt=''):
    files = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    f = open(molecule+analtype+txt, 'w')
    f.write('\begin{figure} \n')
    f.write('\centering \n')

    for line in files:
        f.write('\hspace{0mm} \n')
        f.write('\subfloat['+line[16:32]+']{ \n')
        f.write(' '+'\includegraphics[width=75mm,height=37.5mm]{"'+line[:-4]+'".png} \n')
        f.write('  '+'\label{'+molecule[:3]+line[16:32]+'} \n')
        f.write('} \n')
        print(line)
    f.write('\caption{'+molecule[0].upper()+molecule[1:]+' raw dihedral results'+'} \n')
    f.write('\end{figure}')
    f.close()
    return files


def latexFigBuilder_folding(filepath, molecule, analtype,txt=''):
    files = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    f = open(molecule+analtype+txt, 'w')
    f.write('\begin{figure} \n')
    f.write('\centering \n')

    for line in files:
        f.write('\hspace{0mm} \n')
        f.write('\subfloat['+line[10:24]+']{ \n')
        f.write(' '+'\includegraphics[width=75mm,height=37.5mm]{"'+line[:-4]+'".png} \n')
        f.write('  '+'\label{atenfold-} \n')
        f.write('} \n')
        print(line)
    f.write('\caption{'+molecule[0].upper()+molecule[1:]+' raw dihedral results'+'} \n')
    f.write('\end{figure}')
    f.close()
    return files



def latexFigBuilder_ordering(filepath, molecule, analtype, txt=''):
    files = [f for f in listdir(filepath) if isfile(join(filepath, f))]
    f = open(molecule+analtype+txt, 'w')
    f.write('\begin{figure} \n')
    f.write('\centering \n')

    for line in files:
        f.write('\hspace{0mm} \n')
        f.write('\subfloat['+line[10:24]+']{ \n')
        f.write(' '+'\includegraphics[width=75mm,height=37.5mm]{"'+line[:-4]+'".png} \n')
        f.write('  '+'\label{propord-} \n')
        f.write('} \n')
        print(line)
    f.write('\caption{'+molecule[0].upper()+molecule[1:]+' Ordering results'+'} \n')
    f.write('\end{figure}')
    f.close()
    return files

