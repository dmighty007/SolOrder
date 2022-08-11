from MDAnalysis.lib.nsgrid import FastNS as _FastNS
import numpy as _np
import pandas as _pd
import tqdm.auto as _tqdm
from numba import njit as _njit

@_njit
########### Calculates distance obeying PBC.......
def dist(a, b, box):
    """
    params:
        a: vector, position of 1st particle
        b: vector, position of 2nd particle
    returns:
        distance between two vectors, considering PBC
    """
    #box = box
    dx = b - a
    for i in range(3):
        if abs(dx[i]) >= box[i]/2.0:
            #print(dx[i], box[i])
            dx[i] = box[i] - abs(dx[i])
    return _np.sqrt((dx**2).sum())
@_njit 
def neighbours(i,arr):
    val = _np.unique(_np.concatenate((arr[arr[:,0] == i], arr[arr[:,1] == i])).ravel())
    mask = val != i
    return val[mask]
@_njit
def _e(i,arr,pos,global_rho,r,rsq,sigma,prefactor,cutoff,box):
    distance = []
    for item in neighbours(i,arr):
        distance.append(dist(pos[item], pos[i], box))
    rij = _np.array(distance)
    r_diff = _np.expand_dims(r, 0) - _np.expand_dims(rij, 1)
    g_m = _np.sum(_np.exp(-r_diff**2 / (2.0*sigma**2)), axis=0) / prefactor
    local_volume = 4/3 * _np.pi * cutoff**3
    rho = len(rij) / local_volume
    g_m *= global_rho / rho
    integrand = _np.where(g_m >= 1e-10, (g_m * _np.log(g_m) - g_m + 1.0) * rsq, rsq)
    int_val = -2.0 * _np.pi * rho * _np.trapz(integrand, r)
    #entropy_list.append(int_val)
    return int_val

class eFingerprint:
    def __init__(self, universe, frame = 1, include_hydrogen = False):
        
        """
        Takes MDAnalysis universe object and returns local entropy fingerprints of the oxygen atoms..
        include_hydrogen=False by default..
        
        """
        ox = universe.select_atoms("name OW").ids - 1
        _neighbours = lambda i,arr : _np.unique(_np.concatenate([arr[arr[:,0] == i], arr[arr[:,1] == i]]).ravel())
        neighbours = lambda i,arr : _neighbours(i,arr)[_neighbours(i,arr) != i]
        self.cutoff = 5.0
        self.sigma = 0.15
        self.box = universe.dimensions[:3]
        self.dimensions = universe.dimensions
        self.volume = universe.dimensions[:3][0] * universe.dimensions[:3][1] * universe.dimensions[:3][2]
        self.pos = universe.select_atoms("name OW").positions 
        self.n_atoms = universe.select_atoms("name OW").n_atoms
        self.global_rho = self.n_atoms / self.volume
        self.identifier = _pd.Series(universe.select_atoms("name OW").names)
        if include_hydrogen == True:
            self.pos = universe.atoms.positions 
            self.n_atoms = universe.atoms.n_atoms
            self.global_rho = self.n_atoms / self.volume
            self.identifier = _pd.Series(universe.atoms.names)

        #local_entropy = np.empty(n_atoms)
        nbins = int(self.cutoff / self.sigma) + 1
        self.r = _np.linspace(0.0, self.cutoff, num=nbins)
        self.rsq = self.r**2
        self.prefactor = self.rsq * (4 * _np.pi * self.global_rho * _np.sqrt(2 * _np.pi * self.sigma**2))
        self.prefactor[0] = self.prefactor[1] 
        gridsearch = _FastNS(self.cutoff, self.pos, self.dimensions, pbc=True)
        results = gridsearch.self_search()
        self.arr = results.get_pairs()
        self.ox_ids = _np.array(self.identifier[self.identifier == "OW"].index)
    def localEntropy(self):
        """
         Returns projection of entropy on local oxygen atom..
        """

        ox_num = len(self.ox_ids)
        loc_entropy = _np.zeros(ox_num)
        for i in _tqdm.trange(ox_num):
            loc_entropy[i] = (_e(self.ox_ids[i],self.arr,self.pos,self.global_rho,self.r,self.rsq,self.sigma,self.prefactor,self.cutoff, self.box))
        self.localentropy = loc_entropy
        