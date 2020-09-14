# -*- coding: utf-8 -*-
from Crystal import crystal
from Gen_atom import atom
#from Model_input import get_input
#from Model_input import get_input
# for direct coordinate
# Read POSCAR in VASP5.x format
def read_poscar1(xd):
    row = xd.readline()
    if row == '':
        return row
    cellr = []
    row = xd.readline()
    rc = ','.join(filter(lambda x: x, row.split(' ')))
    sc = rc.split(',')
    afactor = float(sc[0])
    for i in range(3):
        temp = []
        row = xd.readline()
        rc = ','.join(filter(lambda x: x, row.split(' ')))
        sc = rc.split(',')
        for j in range(3):
            temp.append(float(sc[j]) * afactor)
        cellr.append(temp)
    nms = []
    row = xd.readline()
    rc = ','.join(filter(lambda x: x, row.split(' ')))
    sc = rc.split(',')
    for i in range(len(sc)):
        nms.append(sc[i].replace("\n",""))
    nbs = []
    row = xd.readline()
    rc = ','.join(filter(lambda x: x, row.split(' ')))
    rc = rc.replace("\n","")
    sc = rc.split(',')
    tn = 0
    for i in range(len(sc)):
        nbs.append(int(sc[i]))
        tn += int(sc[i])
    namelist = []
    for i in range(len(nms)):
        for j in range(nbs[i]):
            namelist.append(nms[i])
    atoms = []
    row = xd.readline()
    for i in range(tn):
        at = []
        at.append(namelist[i])
        row = xd.readline()
        rc = ','.join(filter(lambda x: x, row.split(' ')))
        sc = rc.split(',')
        for j in range(3):
            tm = 0
            for k in range(3):
                tm += float(sc[j]) * cellr[k][j]
            at.append(tm)
        atoms.append(at)
    cry = crystal()
    cry.set_cell(cellr)
    cry.set_atoms(atoms)
    
    return cry
# Read POSCAR in VASP4.x format
def read_poscar2(xd):
    row = xd.readline()
    if row == '':
        return row
    cellr = []
    row = xd.readline()
    rc = ','.join(filter(lambda x: x, row.split(' ')))
    sc = rc.split(',')
    afactor = float(sc[0])
    for i in range(3):
        temp = []
        row = xd.readline()
        rc = ','.join(filter(lambda x: x, row.split(' ')))
        sc = rc.split(',')
        for j in range(3):
            temp.append(float(sc[j]) * afactor)
        cellr.append(temp)
    nbs = []
    row = xd.readline()
    row = row.replace("\n","")
    rc = ','.join(filter(lambda x: x, row.split(' ')))
    sc = rc.split(',')
    tn = 0
    for i in range(len(sc)):
        nbs.append(int(sc[i])) #need change
        tn += int(sc[i])
    atoms = []
    row = xd.readline()
    for i in range(len(nbs)):
        for j in range(nbs[i]):
            at = []
            row = xd.readline()
            row = row.replace("\n","")
            rc = ','.join(filter(lambda x: x, row.split(' ')))
            sc = rc.split(',')
            for j in range(3):
                tm = 0
                for k in range(3):
                    tm += float(sc[j]) * cellr[k][j]
                at.append(tm)
            at.insert(0,sc[3])
            atoms.append(at)
    cry = crystal()
    cry.set_cell(cellr)
    cry.set_atoms(atoms)
    
    return cry
# VASP5.x and VASP4.x
def read_poscar(xd,vasp='VASP5'):
    if vasp == 'VASP4':
        return read_poscar2(xd)
    else:
        return read_poscar1(xd)
# Creat crystal object series
def deal_poscars(file_path,vasp='VASP5'):
    with open(file_path,'r') as xd:
        #row = xd.readline()
        crys = []
        cry = True
        while cry:
            cry = read_poscar(xd,vasp)
            if cry != '':
                crys.append(cry)
            #row = xd.readline()
    return crys
# Creat a crystal object
def deal_poscar(file_path,vasp='VASP5'):
    with open(file_path,'r') as xd:
        cry = read_poscar(xd,vasp)
    return cry

    
if __name__ == '__main__':
    file_path = 'POSCAR'
    
    with open(file_path,'r') as xd:
        cry1 = read_poscar(xd,'VASP4')
    #ps1 = get_input(cry1)
    #with open('POSCAR3','r') as xd:
    #    cry2 = read_poscar(xd,'VASP4')
        
    #ps2 = get_input(cry2)
    #dps = [ps1[i] - ps2[i] for i in range(len(ps1))]
    #print(dps)
    #print(ps1)
    '''
    print(cry.get_max_atomicNum())
    ca = cry.count_atoms()
    for k in ca:
        at = atom(k)
        print(at.get_atomizationEnthalpy())
    '''
    
    #crys = deal_poscars(file_path)
    #print(crys[1].get_max_atomicNum())
    #print(len(crys))
    #iad = cry.get_iadmatrix()
    #adj = cry.get_adjmatrix()
    #print(iad)
    #print(adj)
    #file_path2 = 'structrue'
    #crys = deal_poscars(file_path2)
    #print(len(crys))
    