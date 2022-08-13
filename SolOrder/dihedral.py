#!/bin env python
import MDAnalysis as _mda
from MDAnalysis.lib.nsgrid import FastNS as _FastNS
import numpy as _np
import pandas as _pd
import tqdm.auto as _tqdm
from MDAnalysis.lib.distances import calc_angles as _calc_angles
from MDAnalysis.lib.distances import calc_bonds as _calc_bonds



class OrderParameter:
    
    def __init__(self, universe, cutoff = 4.5, frame = 0):
        """
             params :
                 universe : _mdanalysis universe object
                 cutoff : float (default = 4.5)
                 frame : int; frame number if trajectory(default 0)
        
        """
        self.cutoff = cutoff
        universe.trajectory[frame]
        self.u = universe
        self.box = universe.dimensions
        self.dimensions = universe.dimensions
        self.rids = universe.select_atoms("name OW").resids
        self.n = universe.select_atoms("name OW").n_atoms
        self.srids = _pd.Series(self.rids)
        self.maping = lambda i : self.srids[self.srids.index == i].values
        self.pos = self.u.select_atoms("name OW").positions
        self.n_atoms = self.u.select_atoms("name OW").n_atoms
        #self.self.dist = lambda a, b : _np.linalg.norm(a-b)
        self._neighbours = lambda i,arr : _np.unique(_np.concatenate([arr[arr[:,0] == i], arr[arr[:,1] == i]]).ravel())
        self.neighbours = lambda i,arr : self._neighbours(i,arr)[self._neighbours(i,arr) != i]
        gridsearch = _FastNS(self.cutoff, self.pos, self.u.dimensions, pbc=True)
        
        results = gridsearch.self_search()
        self.arr = results.get_pairs()
        
    def _filter_(self, item, arr, pos):
        nn = self.neighbours(item,arr)
        #pos = self.pos
        dist = [_calc_bonds(pos[item], pos[i], self.box) for i in nn]
        mask = _np.argsort(dist)[:4]
        return _np.array(nn)[mask], _np.array(dist)[mask]
    ########### Calculates distance obeying PBC.......
    def dist(self,a,b):
        """
        params:
            a: vector, position of 1st particle
            b: vector, position of 2nd particle
        returns:
            distance between two vectors, considering PBC
        """
        box = self.box
        dx = b - a
        for i in range(3):
            if abs(dx[i]) >= box[i]/2.0:
                #print(dx[i], box[i])
                dx[i] = box[i] - abs(dx[i])
        return _np.sqrt((dx**2).sum())
    
    ############# Positions of the hydrogens and oxygens of two adjacent molecule...
    def get_atoms(self,i,j):
        """
        params:
            i: int, resid of first oxygen molcule
            j: int, resid of second oxygen molecule
        returns:
            positions(N,3) of corresponding hydrogens and oxygens..
        """
        u = self.u
        h1a = u.select_atoms(f"resid {i} and name HW1").positions[0]
        h2a = u.select_atoms(f"resid {i} and name HW2").positions[0]
        oa = u.select_atoms(f"resid {i} and name OW").positions[0]
        h1b = u.select_atoms(f"resid {j} and name HW1").positions[0]
        h2b = u.select_atoms(f"resid {j} and name HW2").positions[0]
        ob = u.select_atoms(f"resid {j} and name OW").positions[0]
        return h1a, h2a, oa, h1b, h2b, ob
    
    ############ Choosing which hydrogens should be there in the dihedral...
    def _sort_dihedral(self, h1a, h2a, oa, h1b, h2b, ob):
        """
        params:
            takes all 6 positions of atoms of two oxygen molecules
        returns:
            position of 4 atoms that maintains furthest hydrogen criteria..
        """
        #self.dist = self.self.dist
        ##case1: dihedral atoms are : H1A-OA-OB-H1B or H1A-OA-OB-H2B
        if _np.min([self.dist(oa, h1b), _np.min([self.dist(oa, h2b), _np.min([self.dist(ob, h1a), self.dist(ob,h2a)])])]) == self.dist(ob, h2a):
            if _np.max([self.dist(oa,h1b), _np.max([self.dist(oa,h2b), _np.max([self.dist(ob,h1a), self.dist(ob,h2a)])])]) == self.dist(oa, h1b):
                dih = _np.array([h1a,oa,ob,h1b])
            else:
                dih = _np.array([h1a,oa,ob,h2b])
        ##case2: dihedral atoms are : H2A-OA-OB-H1B or H2A-OA-OB-H2B
        if _np.min([self.dist(oa, h1b), _np.min([self.dist(oa, h2b), _np.min([self.dist(ob, h1a), self.dist(ob,h2a)])])]) == self.dist(ob, h1a):
            if _np.max([self.dist(oa,h1b), _np.max([self.dist(oa,h2b), _np.max([self.dist(ob,h1a), self.dist(ob,h2a)])])]) == self.dist(oa, h1b):
                dih = _np.array([h2a,oa,ob,h1b])
            else:
                dih = _np.array([h2a,oa,ob,h2b])
        ##case3: dihedral atoms are : H2A-OA-OB-H2B OR H1A-OA-OB-H2B. 
        ##This gives the same dihedral angles as case2-b and case1-b, however, it needs         
        if _np.min([self.dist(oa, h1b), _np.min([self.dist(oa, h2b), _np.min([self.dist(ob, h1a), self.dist(ob,h2a)])])]) == self.dist(oa, h1b):
            if _np.max([self.dist(oa,h1b), _np.max([self.dist(oa,h2b), _np.max([self.dist(ob,h1a), self.dist(ob,h2a)])])]) == self.dist(ob, h1a):
                dih = _np.array([h2b,ob,oa,h1a])
            else:
                dih = _np.array([h2b,ob,oa,h2a])
        ##case4: dihedral atoms are : H2A-OA-OB-H1B OR H1A-OA-OB-H1B. 
        ##This gives the same dihedral angles as case2-a and case1-a, however, it needs        
        if _np.min([self.dist(oa, h1b), _np.min([self.dist(oa, h2b), _np.min([self.dist(ob, h1a), self.dist(ob,h2a)])])]) == self.dist(oa, h2b):
            if _np.max([self.dist(oa,h1b), _np.max([self.dist(oa,h2b), _np.max([self.dist(ob,h1a), self.dist(ob,h2a)])])]) == self.dist(ob, h1a):
                dih = _np.array([h1b,ob,oa,h1a])
            else:
                dih = _np.array([h1b,ob,oa,h2a])
        return dih
    
    ############# Calculates dihedral angles .. provided a quadralet..
    def calculate_dihedral(self, a,b,c,d):
        """
        params:
            takes 4 points in cartesian space
        returns:
            diedral angle
        """
        ab = b - a
        bc = c - b
        cd = d - c
        abc = _np.cross(ab, bc)
        bcd = _np.cross(bc, cd)
        return -_np.arccos(_np.dot(abc,bcd)/(_np.linalg.norm(abc)*_np.linalg.norm(bcd)))
    ################# Returns all dihedrals........
    def get_dihedrals(self):
        """
        returns all dihedral angles..
        """
        u = self.u
        gridsearch = _FastNS(self.cutoff, u.select_atoms("name OW").positions, u.dimensions, pbc=True)
        results = gridsearch.self_search()
        arr = results.get_pairs()
        mapped_list = []
        dihedral = []

        for i in _tqdm.trange(self.n):
            nn = self.neighbours(i,arr)
            for j in nn:
                flag1 = [i,j] not in mapped_list 
                flag2 = [j,i] not in mapped_list
                flag = flag1*flag2
                if flag:            
                    res1 = self.maping(i)[0]
                    res2 = self.maping(j)[0]
                    h1a, h2a, oa, h1b, h2b, ob = self.get_atoms(res1,res2)
                    dih = self._sort_dihedral(h1a, h2a, oa, h1b, h2b, ob)
                    dihedral.append(self.calculate_dihedral(*dih))
                    mapped_list.append([i,j])
        self.dihedrals = dihedral
    ####################### Returns individual F4 values....    
    def F4(self):
        
        """
        returns individual F4 values....
        """
        u = self.u
        gridsearch = _FastNS(self.cutoff, u.select_atoms("name OW").positions, u.dimensions, pbc=True)
        results = gridsearch.self_search()
        arr = results.get_pairs()
        counter_list = []
        f4_list = []
        for i in _tqdm.trange(self.n):
            nn = self.neighbours(i,arr)
            avg_f4 = 0
            particle_counter = 0
            for j in nn:         
                res1 = self.maping(i)[0]
                res2 = self.maping(j)[0]
                h1a, h2a, oa, h1b, h2b, ob = self.get_atoms(res1,res2)
                dih = self._sort_dihedral(h1a, h2a, oa, h1b, h2b, ob)
                #phi = _calc_dihedrals(dih[0], dih[1],dih[2],dih[3], box)
                phi = self.calculate_dihedral(dih[0], dih[1],dih[2],dih[3])
                avg_f4 += _np.cos(3*phi)
                particle_counter += 1
            counter_list.append(particle_counter)
            f4_list.append(avg_f4/particle_counter)
        self.f4 = f4_list
        
    ############# orientational tetrahedral order parameter........
    def singleOTO(self, item, arr):
        #u = self.u
        #gridsearch = _FastNS(5.0, u.select_atoms("name OW").positions, u.dimensions, pbc=True)
        pos = self.pos
        #results = gridsearch.self_search()
        #arr = self.arr
        li = self._filter_(item,arr, pos)[0]
        q = 0.0
        try:
            for i in range(3):
                for j in range(i+1, 4):
                    cos_phi = _np.cos(_calc_angles(pos[li[i]],pos[item], pos[li[j]], box = self.dimensions))
                    q += (cos_phi + 1 / 3) ** 2
            q = 1 - 3 / 8 * q
            
        except:
            print(f"Check for atom : {item}!")
            print("Specify cutoff!!")
            pass

        return q
    
    def OTO(self):
        """
         Returns Orientational tetrahedral order paramter of each oxygen atoms..
        """
        natoms = self.u.select_atoms("name OW").n_atoms
        q_array = _np.zeros(natoms)
        gridsearch = _FastNS(5.0 , self.pos, self.u.dimensions, pbc=True)
        results = gridsearch.self_search()
        arr = results.get_pairs()
        self.arr = arr
        for i in _tqdm.trange(natoms):
            q_array[i] = self.singleOTO(i,arr)
        self.tetra_orient = q_array
    ############# orientational tetrahedral order parameter........
    def singleTTO(self, item):
        #u = self.u
        #gridsearch = _FastNS(5.0, u.select_atoms("name OW").positions, u.dimensions, pbc=True)
        pos = self.pos
        #results = gridsearch.self_search()
        arr = self.arr
        dist = _np.array(self._filter_(item,arr, pos)[1])
        #q = 0.0
        try:
            r_bar = _np.mean(dist)
            #for i in range(4):
            #r_bar = _np.mean(dist)
            sqrt_dist = (dist - r_bar)**2
            a = sqrt_dist.sum()/(4* (r_bar)**2) 
            #print(sqrt_dist)
            a = 1 - (a /3)
        except:
            gridsearch = _FastNS(self.cutoff + 2.0 , self.pos, self.u.dimensions, pbc=True)
            results = gridsearch.self_search()
            arr = results.get_pairs()
            dist = _np.array(self._filter_(item,arr, pos)[1])
            r_bar = _np.mean(dist)
            #for i in range(4):
            #r_bar = _np.mean(dist)
            sqrt_dist = (dist - r_bar)**2
            a = sqrt_dist.sum()/(4* (r_bar)**2) 
            #print(sqrt_dist)
            a = 1 - (a /3)
            pass

        return a
    
    def TTO(self):
        """
         Returns Translational tetrahedral order paramter of each oxygen atoms..
        """

        natoms = self.u.select_atoms("name OW").n_atoms
        q_array = _np.zeros(natoms)
        for i in _tqdm.trange(natoms):
            q_array[i] = self.singleTTO(i)
        self.tetra_trans = q_array
    ######################################### Local Structure Index.......
    def singleLSI(self, item):
        arr = self.arr
        pos = self.pos
        try:
            newarr = _np.sort([_calc_bonds(pos[newpos], pos[item], self.dimensions) for newpos in self.neighbours(item,arr)])
            r_list = newarr[:_np.where(newarr > 3.7)[0][1]]
        except:
            gridsearch = _FastNS(self.cutoff + 2.0 , self.pos, self.u.dimensions, pbc=True)
            results = gridsearch.self_search()
            arr = results.get_pairs()
            newarr = _np.sort([_calc_bonds(pos[newpos], pos[item], self.dimensions) for newpos in self.neighbours(item,arr)])
            r_list = newarr[:_np.where(newarr > 3.7)[0][1]]
        delR = r_list[1:] - r_list[:-1]
        return ((delR - delR.mean())**2).mean()
    
    def LSI(self):
        """
         Returns Local structure index order paramter of each oxygen atoms..
        """

        lsi = _np.zeros(self.n_atoms)
        for i in _tqdm.trange(self.n_atoms):
            lsi[i] = self.singleLSI(i)
        self.lsi = lsi
        
   ####################################### Minimum Angle distribution.......
    def singleMinimumAngle(self, item):
        
        arr = self.arr
        pos = self.pos
        item2 = self._filter_(item,arr,pos)[0][0]
        res1 = self.maping(item)[0]
        res2 = self.maping(item2)[0]
        h1a, h2a, oa, h1b, h2b, ob = self.get_atoms(res1,res2)
        ang1 = _calc_angles(h1a,oa,ob, self.dimensions)
        ang2 = _calc_angles(h2a,oa,ob, self.dimensions)
        ang3 = _calc_angles(oa,ob,h1b, self.dimensions)
        ang4 = _calc_angles(oa,ob,h2b, self.dimensions)
        
        return _np.array([ang1, ang2, ang3, ang4]).min(), item2 
    
    def MinimumAngle(self):
        """
         Returns Minimum angle of each oxygen atoms with their immediate neighbours..
        """

        minAngles = _np.zeros(self.n_atoms)
        adj_list = []
        for i in _tqdm.trange(self.n_atoms):

            if i not in adj_list:
                angle, adj = self.singleMinimumAngle(i)
                adj_list.append(adj)
                minAngles[i] = angle
                minAngles[adj] = angle
            else:
                adj_list.append(i)
        self.minAngles = minAngles
