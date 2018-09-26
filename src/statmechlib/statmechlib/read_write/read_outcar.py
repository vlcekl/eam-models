
import os
import numpy as np
import glob

def read_outcar(dataset):
    """
    Reads VASP OUTCAR file a given directory
    and returns trajectory data in a dictionary.
    
    Parameters
    ----------
    dataset : string
              directory with VASP OUTCAR file
    Returns
    -------
    traj : dictionary
           trajectory information (configuration, box, energy, forces)
    """
    
    # read configurations (box + particle coordinates)
    with open(os.path.join(dataset, 'OUTCAR'), 'r') as fc:
        xyzs = [] ; boxs = [] ; enes = [] ; temps = []
        for line in iter(fc.readline, ''):
            line = fc.readline()
            
            # box parameters
            box = np.empty((3, 3), dtype=float)
            for i in range(3):
                box[i,:] = [float(x) for x in re.findall('\S+', fc.readline())]

            # number of atoms
            line = fc.readline()
            nat = int(re.findall('\S+', fc.readline())[0])
            line = fc.readline()
            
            # atomic configuration
            xyz = np.empty((nat, 3), dtype=float)
            for i in range(nat):
                xyz[i] = [float(x) for x in re.findall('\S+', fc.readline())]
            
            boxs.append(box)
            xyzs.append(xyz)
    
    # read configurational energies
    with open(os.path.join(dataset, 'md.out'), 'r') as fe:
        enes = [] ; temps = []
        for line in iter(fe.readline, ''):
            if re.search('T=', line):
                sarr = re.findall('\S+', line)
                temps.append(float(sarr[2]))
                enes.append(float(sarr[8]))
    
    # check if the lengths of trajectory lists match
    assert len(enes) == len(xyzs), f'{dataset} XYZ and energy lenghts do not match: {len(enes)}, {len(xyzs)}'
    
    # combine trajectory data in a dictionary
    traj = {'box':boxs, 'xyz':xyzs, 'energy':enes, 'forces':[], 'temp':temps}

    return traj
