# source /export/home/lizd20/clf/clf_jax/bin/activate

from pyscf import gto, dft, scf, ao2mo, lib, tddft,fci
from pyscf.fci import addons
import numpy as np
from pyscf.dft import numint
import scipy
from functools import reduce
import h5py
import time
import sys
#sys.path.append('/export/home/lizd20/zhaohw/tddft/')
#from SF_XTDA import SA_SF_TDA

#angles = np.arange(0,181,5)
nstates = 5
funs = ['bhandhlyp']
#fun = 'b3lyp'
def update_h5(name,sa,count,excited,vect,ground=None,ds=None):
    fname = name + '.h5'
    f = h5py.File(fname,'a')
    if ground is not None:
        f['grounds'][count] = ground
        f['Ds'][count] = ds
    f['vects'][sa,count] = np.array(vect)
    f['exciteds'][sa,count] = np.array(excited)
    f.close()
    return None 
def write_h5(R,nstates,name):
    fname = name + '.h5'
    len_r = len(R)
    f = h5py.File(fname,'w')
    f['R'] = R
    f['grounds'] = np.zeros(len_r)
    f['vects'] = np.zeros((5,len_r,nstates,5))
    f['exciteds'] = np.zeros((5,len_r,nstates))
    f['Ds'] = np.zeros((len_r,nstates))
    f.close()
    return None

def read_xyz(name):
    filename = '/export/home/lizd20/zhaohw/tddft/tadf/46_IST/xyz/'+ name
    with open(filename, 'r') as f:
        lines = f.readlines()
    atoms = []
    for line in lines[2:]:
        parts = line.split()
        if len(parts) >=4:
            atom = parts[0]
            coords = list(map(float, parts[1:4]))
            atoms.append((atom,coords))
    return atoms
mols = [str(i)+'_optimized.xyz' for i in np.arange(46)+1]
#mols = ['17_optimized.xyz','28_optimized.xyz']
#mols = ['37_optimized.xyz']
for fun in funs:
    for i,mol_i in enumerate(mols):
        atoms = read_xyz(mol_i)    
        mol = gto.M(                 
              atom = atoms,
          basis = 'aug-cc-pvdz',
          charge = 0,
          spin = 0, 
          verbose = 4,
     #     symmetry=True,
        )    
        start = time.time()
        if fun == 'hf':
            mf = scf.RHF(mol)
        else:
            mf = dft.RKS(mol)
            mf.xc = fun
    #mf.grids.level = 6
    #mf.irrep_nelec = {'A1': (4,2), 'B1':(1,1), 'B2':(1,1)}
      #  mf.level_shift = 0.1
      #  mf.grids.atom_grid = (99,590)
    #mf = mf.newton()
    #mf.max_cycle = 100
    #mf.diis_damp = 0.5
       # mf.diis_start_cycle = 3
       # mf.init_guess = 'huckel'
        #if fun != 'b3lyp':
        #data = scf.chkfile.load(f'/export/home/lizd20/zhaohw/tddft/tadf/46_IST/tda/chk_file/mol{i+1}.chk','scf')
        #mf.mo_coeff = data['mo_coeff']
        #mf.mo_occ = data['mo_occ']
        #mf.mo_energy = data['mo_energy']
        #else:
        #    mf.chkfile = f'chk_file/mol{i+1}.chk'
        mf.kernel()
    #mf.analyze()
        if not mf.converged:
            mf = mf.newton()
            mf.kernel()
        assert mf.converged
        mo, _, stable, _ = mf.stability(return_status=True)
        if stable is False:
            dm = mf.make_rdm1(mo, mf.mo_occ)
            mf.kernel(dm0=dm)
            mo, _, stable, _ = mf.stability(nroots=5,tol=1e-7,return_status=True)
            assert stable
        
        print(f'fun_{fun}')
        print(mol_i,' Singlet')
        mytd = tddft.TDDFT(mf)
        mytd.nstates = nstates
        mytd.kernel()
        mytd.analyze()

        print(mol_i,' Triplet')
        mytd2 = tddft.TDDFT(mf)
        mytd2.nstates = nstates
        mytd2.singlet = False
        mytd2.kernel()
        mytd2.analyze()

print('Done')
