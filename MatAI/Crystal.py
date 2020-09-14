# -*- coding: utf-8 -*-
from Auxilfunction import vcross, vdot
from Gen_atom import atom
import math

def get_min_iad(cellr,iatom,jatom):
    rt = []
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                x = [i*xt for xt in cellr[0]]
                y = [i*yt for yt in cellr[1]]
                z = [i*zt for zt in cellr[2]]
                dr = [iatom[i+1]+x[i]+y[i]+z[i]-jatom[i+1] for i in range(len(x))]
                rt.append(math.sqrt(vdot(dr,dr)))
    rs = [mt for mt in rt if mt > 0.1]
    return min(rs)
                
class crystal(object):
    def __init__(self):
        self.cellr = [[1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,0.0,1.0]]
        self.atoms = []
    def set_cell(self,cr):
        self.cellr = cr
    def set_atoms(self,at):
        self.atoms = at;
        
    def get_volume(self):
        return abs(vdot(vcross(self.cellr[0],self.cellr[1]),self.cellr[2]))
    def count_atoms(self):
        catoms = {}
        atoms = [s[0] for s in self.atoms]
        for i in set(atoms):
            catoms[i] = atoms.count(i)
        return catoms
    def get_number_atoms(self):
        return len(self.atoms)
    def print_atoms(self):
        return print(self.atoms)
    def get_chemical_formula(self):
        sn = ''
        catoms = self.count_atoms()
        for k in catoms:
            sn += k
            sn += str(catoms[k])
        return sn
    def get_max_atomicNum(self):
        minn = 0
        catoms = self.count_atoms()
        for k in catoms:
            at = atom(k)
            if at.atomicNum > minn:
                minn = at.atomicNum
        return minn
    def get_iadmatrix(self):
        iadm = []
        for i in range(len(self.atoms)):
            iadm.append([])
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                iadm[i].append(get_min_iad(self.cellr,self.atoms[i],self.atoms[j]))
        return iadm
    def get_adjmatrix(self):
        iadm = self.get_iadmatrix()
        adjs = []
        for i in range(len(iadm)):
            te = [ia for ia in iadm[i] if ia !=0]
            if len(te) == 0:
                md = 0
            else:
                md = min(te)
            adj = [0 if im > 1.00001*md else 1 for im in iadm[i]] # can change
            adjs.append(adj)
        return adjs
    def get_invnmatrix(self,n=2):
        iadm = self.get_iadmatrix()
        adjs = self.get_adjmatrix()
        invns = []
        invn = [0] * len(iadm[0])
        for i in range(len(iadm)):
            invns.append(invn)
        for i in range(len(iadm)):
            for j in range(len(iadm[0])):
                if adjs[i][j] > 0 and iadm[i][j] > 0:
                    invns[i][j] = 1.0
                    for k in range(n):
                        invns[i][j] /= iadm[i][j]
        return invns
    


if __name__ == '__main__':
    cry = crystal()
    
    kt = [[2.0,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,4.0]]
            
    cry.set_cell(kt)
    for i in range(3):
        print(cry.cellr[i])
    print(cry.get_volume())
    
    sa = [['H',1,2,3],['H',2,1,3],['O',2,1,1]]
    
    cry.set_atoms(sa)
    mt = cry.count_atoms()
    print(mt)
    print('The number of atoms is: ',cry.get_number_atoms())
    for key in mt:
        print(key)
    iad = cry.get_iadmatrix()
    adj = cry.get_adjmatrix()
    invs = cry.get_invnmatrix()
    print(cry.get_chemical_formula())
    
    with open('POSCAR2','r') as xd:
        cry2 = read_poscar(xd,'VASP4')
    cry2.print_atoms()
    #print(iad)
    #print(adj)
    #print(invs)
            
    
    