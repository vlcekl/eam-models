import os
import re
import numpy as np
import glob

def read_xdatcar(filename):
    """
    Reads VASP XDATCAR

    Parameters
    ----------
    filename: str
              full path and name of the XDATCAR

    Returns
    -------
    traj: dict
          trajectory information with keys given below
    box: list of 3x3 ndarrays
          Box dimensions
    xyz: list of natom x 3 ndarrays
          Atomic configurations
    atom_type: list of str
                atom types (names)
    atom_num: list of ints
                atom numbers for each type
    """

    with open(filename, 'r') as fc:

        xyzs = [] ; boxs = []

        for line in iter(fc.readline, ''):

            scale = float(re.findall('\S+', fc.readline())[0])
            
            # box parameters
            box = np.empty((3, 3), dtype=float)
            for i in range(3):
                box[i,:] = [float(x)*scale for x in re.findall('\S+', fc.readline())]

            # atom names
            atom_types = [name for name in re.findall('\S+', fc.readline())]

            # number of atoms
            atom_nums = [int(num) for num in re.findall('\S+', fc.readline())]
            nat = sum(atom_nums)

            line = fc.readline()
            
            # atomic configuration
            xyz = np.empty((nat, 3), dtype=float)
            for i in range(nat):
                xyz[i] = [float(x) for x in re.findall('\S+', fc.readline())]
            
            boxs.append(box)
            xyzs.append(xyz)

    traj = {'box':boxs, 'xyz':xyzs, 'atom_type':atom_types, 'atom_num':atom_nums}

    return traj

def read_poscar(filename):
    """
    Reads VASP POSCAR and CONTCAR files

    Parameters
    ----------
    filename: str
              full path and name of the POSCAR

    Returns
    -------
    traj: dict
          trajectory information with keys given below
    box: 3x3 ndarray
         Box dimensions
    xyz: natom x 3 ndarray
         Atomic configuration
    atom_type: list of str
                atom types (names)
    atom_num: list of ints
                atom numbers for each type
    """

    with open(filename, 'r') as fc:
        
        # atom names
        atom_types = [name for name in re.findall('\S+', fc.readline())]

        scale = float(re.findall('\S+', fc.readline())[0])
            
        # box parameters
        box = np.empty((3, 3), dtype=float)
        for i in range(3):
            box[i,:] = [float(x)*scale for x in re.findall('\S+', fc.readline())]

        # atom names (again, should be same as above)
        atom_types = [name for name in re.findall('\S+', fc.readline())]

        # number of atoms
        atom_nums = [int(num) for num in re.findall('\S+', fc.readline())]
        nat = sum(atom_nums)

        line = fc.readline()
            
        # atomic configuration
        xyz = np.empty((nat, 3), dtype=float)
        for i in range(nat):
            xyz[i] = [float(x) for x in re.findall('\S+', fc.readline())]
            
    traj = {'box':box, 'xyz':xyz, 'atom_type':atom_types, 'atom_num':atom_nums}

    return traj


def read_vasp_0k(dataset):
    """
    Reads configuration and energy files from a VASP energy minimization in a given directory
    and returns trajectory (1 config) data in a dictionary.
    
    Parameters
    ----------
    dataset : string
              directory with VASP MD simulation data has to contain XDATCAR and md.out files
    Returns
    -------
    traj : dictionary
           trajectory information (configuration, box, energy, forces)
    """
    
    # read configurations (box + particle coordinates)
    traj = read_poscar(os.path.join(dataset, 'rlx.CONTCAR'))
    
    # read configurational energies
    with open(os.path.join(dataset, 'rlx.md.out'), 'r') as fe:
        enes = [] ; temps = []
        for line in iter(fe.readline, ''):
            if re.search('F=', line):
                sarr = re.findall('\S+', line)
                ene_last = float(sarr[4])
        temps.append(0.0)
        enes.append(ene_last)
        
    # check if trajectory lists match
    assert len(enes) == len(traj['xyz']), f'{dataset} XYZ and energy lenghts do not match: {len(enes)}, {len(xyzs)}'
    
    # combine trajectory data in a dictionary
    traj['energy'] = enes
    traj['temp'] = temps

    return traj

def read_vasp(dataset):
    """
    Reads configuration and energy files from a VASP MD simulation in a given directory
    and returns trajectory data in a dictionary.
    
    Parameters
    ----------
    dataset : string
              directory with VASP MD simulation data, has to contain XDATCAR and md.out files
    Returns
    -------
    traj : dictionary
           trajectory information (configuration, box, energy, forces)
    """
    
    # read configurations (box + particle coordinates)
    traj = read_xdatcar(os.path.join(dataset, 'XDATCAR'))
    
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
    traj.update({'energy':enes, 'temp':temps})

    return traj


def read_weiwei(target_raw, dataset, poscar):
    """
    Parameters
    ----------
    target_raw: str
                directory with all relevant files
    dataset: str
             subdirectory with the trjectory data
    poscar: str
             POSCAR file name in target_raw directory

    Returns
    -------
    traj: dictionary
          trajectory information (configuration, box, energy, forces)
    """
    
    trajdir = target_raw + '/' + dataset

    traj = read_poscar(os.path.join(target_raw, poscar))

    # read coordinate/forces files of a trajectory
    coor_glob = os.path.join(os.path.join(target_raw, dataset), 'coor*')
    coor_files = glob.glob(coor_glob)
    fnum = lambda fname: int(re.findall('\d+$', fname)[0])
    coor_sorted = sorted(coor_files, key=fnum)

    boxs = []
    xyzs = []
    forces = []
    for fname in coor_files:
        with open(fname, 'r') as f:
            for line in iter(f.readline, ''):
                pass

    # read energy/temperature file
    fname = os.path.join(os.path.join(target_raw, dataset), dataset+'-tem-potenergy')
    with open(fname, 'r') as f:

        enes = [] ; temps = []
        for line in iter(f.readline, ''):
            pass


    # check if the lengths of trajectory lists match
    assert len(enes) == len(xyzs), f'{dataset} XYZ and energy lenghts do not match: {len(enes)}, {len(xyzs)}'
    
    # combine trajectory data in a dictionary
    traj.update({'box':boxs, 'xyz':xyzs, 'energy':enes, 'forces':forces, 'temp':temps})

    return traj
