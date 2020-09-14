# -*- coding: utf-8 -*-
import numpy as np
import math,cmath

def computecorrelation(x,y):
    x_bar=np.mean(x)
    y_bar=np.mean(y)
    SSR=0
    Varx=0
    Vary=0
    for i in range(0,len(x)):
        SSR+=(x[i]-x_bar)*(y[i]-y_bar)
        Varx+=(x[i]-x_bar)**2
        Vary+=(y[i]-y_bar)**2
    SST=cmath.sqrt(Varx*Vary)
    return SSR/SST

def belog(x):
    logx = []
    for i in x:
        if(i > 0):
            logx.append(math.log(i,10))
        elif(i < 0):
            logx.append(0)
    return logx

def polyfot(x,y,degree):
    result={}
    coef=np.polyfit(x,y,degree)#算出各个回归系数
    result["polynomial"]=coef.tolist()
    p=np.poly1d(coef)#拟合一条线
    y_hat=p(x)
    y_bar=np.mean(y)
    SSR=np.sum((y_hat-y_bar)**2)
    SST=np.sum((y-y_bar)**2)
    result["determination"]=SSR/SST
    return result["determination"]
def rmse(x,y):
    z = []
    for i in range(len(x)):
        z.append((x[i] - y[i])**2)
    z_bar = np.mean(z)
    return np.sqrt(z_bar)
def mae(x,y):    
    z = []
    for i in range(len(x)):
        z.append(abs(x[i] - y[i]))
    return np.mean(z)
def get_stastics(y1,y2):
    mae1 = mae(y1,y2)
    rmse1 = rmse(y1,y2)
    cor = computecorrelation(y1,y2)
    print('The MAE is: ',mae1)
    print('The RMSE is: ',rmse1)
    print('The correlation parameter is: ',cor)
    return
def vdot(v1,v2):
    v3 = 0
    if(len(v1) == len(v2)):
        for i in range(len(v1)):
            v3 += v1[i] * v2[i]
    return v3

def vcross(v1,v2):
    v3 = [0,0,0]
    if(len(v1) == len(v2)):
        v3[0] = v1[1] * v2[2] - v1[2] * v2[1]
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2]
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0]
        return v3
    return v3

def zscore(mean,std,av):
    import numpy as np
    import math
    stdav = []
    for i in range(len(mean)):
        if std[i] < 1e-15:
            temp = 0
        else:
            temp = (av[i] - mean[i])/(math.sqrt(std[i]))
        stdav.append(temp)
    return np.array(stdav)

def list_split(lists,clist):
    l1 = []
    l2 = []
    nl = [ i for i in range(0,len(lists))]
    for i in clist:
        l1.append(lists[i])
    for i in nl:
        if i not in clist:
            l2.append(lists[i])
    l1 = np.array(l1)
    l2 = np.array(l2)
    return l1,l2
    
def train_test_split(x_data,y_data,train_size=0.8,test_size=0.2, random_state=0):
    import random
    random.seed(random_state)
    nx = len(x_data)
    ny = len(y_data)
    x_train = []
    x_test = []
    y_train = []
    y_test = []
    if nx == ny:
        if train_size <= 1 and test_size <= 1:
            n1 = int(nx * train_size)
            n2 = int(nx * test_size)
            nlist = [i for i in range(0,nx)]
            clist = random.sample(nlist,n1)
            x_train, x_test = list_split(x_data,clist)
            y_train, y_test = list_split(y_data,clist)
            return x_train,x_test,y_train,y_test
        else:
            nlist = [i for i in range(0,nx)]
            clist = random.sample(nlist,train_size)
            x_train, x_test = list_split(x_data,clist)
            y_train, y_test = list_split(y_data,clist)
            return x_train,x_test,y_train,y_test
    else:
        return 0


if __name__ == '__main__':
    v1 = [1,2,3]
    v2 = [4,5,6]
    v3 = vcross(v1,v2)
    print(vdot(v1,v2))
    print(v3)
    
    x_data = [[1,2],[3,4],[5,6],[7,8],[9,10]]
    y_data = [1,2,3,4,5]
    x1,x2,y1,y2 = train_test_split(x_data,y_data,train_size=0.8,test_size=0.2,random_state=23)
    print(x1)
    print(y1)
    