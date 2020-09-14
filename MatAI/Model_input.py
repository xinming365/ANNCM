# -*- coding: utf-8 -*-

from ReadPOSCAR import deal_poscars
from Crystal import crystal
from Gen_atom import atom
import numpy as np

def get_property_matrix(cry):
    ppty = []
    for i in range(len(cry.atoms)):
        at = atom(cry.atoms[i][0])
        pt = at.get_property()
        ppty.append(pt)
    return ppty
def get_atomavg_property(ppty):
    atomavg = []
    for i in range(len(ppty[0])):
        ptemp = 0.0
        for j in range(len(ppty)):
            ptemp += ppty[j][i]
        ptemp /= len(ppty)
        atomavg.append(ptemp)
    return atomavg

def get_str_property(ppty,invns):
    strp = []
    nn = len(invns) * len(invns)
    for k in range(len(ppty[0])):
        ptemp = 0.0
        for i in range(len(invns)):
            for j in range(len(invns[0])):
                ptemp += abs((ppty[i][k] - ppty[j][k])) * invns[i][j]
        ptemp /= nn
        strp.append(ptemp)
    return strp

def get_str_property2(ppty,invns):
    strp = []
    nn = len(invns) * len(invns)
    for k in range(len(ppty[0])):
        ptemp = 0.0
        for i in range(len(invns)):
            for j in range(len(invns[0])):
                ptemp += abs((ppty[i][k] + ppty[j][k])) * invns[i][j]
        ptemp /= nn
        strp.append(ptemp)
    return strp

def get_volavg_property(ppty,vol):
    volavg = []
    for i in range(len(ppty[0])):
        ptemp = 0.0
        for j in range(len(ppty)):
            ptemp += ppty[j][i]
        ptemp /= vol
        volavg.append(ptemp)
    return volavg

def get_input(cry,key=0):
    #property
    ppty = get_property_matrix(cry)
    vol = cry.get_volume()
    invns = cry.get_invnmatrix()
    atomavg = get_atomavg_property(ppty)
    volavg = get_volavg_property(ppty,vol)
    strp = get_str_property(ppty,invns)
    pp = []
    for i in range(len(atomavg)):
        pp.append(atomavg[i])
    for i in range(len(volavg)):
        pp.append(volavg[i])
    for i in range(len(strp)):
        pp.append(strp[i])
    if key == 0:
        strp2 = get_str_property2(ppty,invns)
        for i in range(len(strp2)):
            pp.append(strp2[i])
    return np.array(pp)

def train_test_data(x_data,y_data,save_file = 'model.save\\test', ptrain=0.8,ptest=0.2):
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    import numpy as np
    x_train, x_test, y_train, y_test = train_test_split(x_data, y_data,train_size=ptrain,test_size=ptest, random_state=312)
    scaler = StandardScaler()
    scaler.fit(x_train)
    mean = scaler.mean_
    std = scaler.var_
    save_mean = save_file + '_mean.npy'
    save_std = save_file + '_std.npy'
    np.save(save_mean,mean)
    np.save(save_std,std)
    return x_train, x_test, y_train, y_test
    
if __name__ == '__main__':
    file_path = 'POSCAR'
    file_path2 = 'structrue'
    crys = deal_poscars(file_path)
    for i in range(len(crys)):
        pp = get_input(crys[i])
        print(pp)
        
            
    

