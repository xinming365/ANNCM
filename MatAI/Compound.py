# -*- coding: utf-8 -*-

from Gen_atom import atom

# Depart the chemical formula to atom names list and ratios list
def readComponent(comps):
    namelist = []
    numlist = []
    ccomps = comps
    while(len(ccomps) != 0):
        stemp = ccomps[1:]
        if(len(stemp) == 0):
            namelist.append(ccomps)
            numlist.append(1.0)
            break
        it = 0
        for st in stemp:
            it = it + 1
            if(st.isupper()):
                im = 0
                for mt in stemp[:it]:
                    im = im + 1
                    if(mt.isdigit()):
                        namelist.append(ccomps[0:im])
                        numlist.append(float(ccomps[im:it]))
                        #print(ccomps[0:im])
                        #print(float(ccomps[im:it]))
                        #print("st1")
                        ccomps = ccomps[it:]
                        break
                    elif(im == len(stemp[:it])):
                        namelist.append(ccomps[0:im])
                        numlist.append(1.0)
                        #print(ccomps[0:im])
                        #print(1.0)
                        #print("st2")
                        ccomps = ccomps[it:]
                        break
                break
            elif(it == len(stemp)):
                im = 0
                for mt in stemp:
                    im = im + 1
                    if(mt.isdigit()):
                        namelist.append(ccomps[0:im])
                        numlist.append(float(ccomps[im:]))
                        #print(ccomps[0:im])
                        #print(float(ccomps[im:]))
                        #print("st3")
                        ccomps = ccomps[it+1:]
                        break
                    elif(im == len(stemp)):
                        namelist.append(ccomps)
                        numlist.append(1.0)
                        #print(ccomps)
                        #print(1.0)
                        #print("st4")
                        ccomps = ccomps[it+1:]
                        break
                break
    return namelist, numlist

# Convert the chemical formula into a standard format
# for example: Fe3(Al2Cr)3 is converted to  Fe3Al6Cr3
def getFullname(comps):
    import re
    namelist = []
    numlist = []
    clist = re.split(r"\(|\)",comps)
    complist = []
    icl = len(clist)
    for i in range(icl):
        if clist[i][0].isdigit():
            cn = ""
            for j in range(len(clist[i])):
                if not clist[i][j].isalpha():
                    cn += clist[i][j]
                else:
                    break
            complist.append(clist[i][0:len(cn)])
            complist.append(clist[i][len(cn):])
        else:
            complist.append(clist[i])     
    nl = len(complist)
    ic = 0
    while ic < nl:
        nal, nul = readComponent(complist[ic])
        if ic + 1 < nl:
            if  complist[ic+1].isdigit():
                fac = float(complist[ic+1])
                nul = [fac*m for m in nul]
                ic += 2
            else:
                ic += 1
        else:
            ic += 1
        if(len(namelist) == 0):
            for k in range(len(nal)):
                namelist.append(nal[k])
                numlist.append(nul[k])
        else:
            for k in range(len(nal)):
                if nal[k] in namelist:
                    ik = namelist.index(nal[k])
                    numlist[ik] += nul[k]
                else:
                    namelist.append(nal[k])
                    numlist.append(nul[k])
    names=""
    for m in range(len(namelist)):
        names += namelist[m]
        ns = '{:g}'.format(numlist[m])
        names += ns
        #names += str(numlist[m])
    return names
# Produce the feature vector using the bag of atom method
# The dimension of the feature vector is 100
def get_atombag_vector(comps):
    namelist, numlist = readComponent(comps)
    avector = 100 * [0]
    asum = sum(numlist)
    for i in range(len(namelist)):
        at = atom(namelist[i])
        avector[at.atomicNum-1] = numlist[i]/asum
    return avector
# Produce the feature vector using the average of the atomic properties
# The dimension of the feature vector is 25
def get_atomic_vector(comps):
    namelist, numlist = readComponent(comps)
    ppty = []
    for an in namelist:
        at = atom(an)
        pt = at.get_property()
        ppty.append(pt)
    atomavg = []
    for i in range(len(ppty[0])):
        ptemp = 0.0
        for j in range(len(ppty)):
            ptemp += ppty[j][i]
        ptemp /= len(ppty)
        atomavg.append(ptemp)
    return atomavg

# Integrage above two ways to produce the feature vector
def get_component_vector(comps,ivec=0):
    cs = getFullname(comps)
    if ivec == 0:
        return get_atombag_vector(cs)
    else:
        return get_atomic_vector(cs)

if __name__ == '__main__':
    # Test the getFullname function
    sys = "Fe3(Al2Cr)3"
    av = getFullname(sys)
    print(av)
    # Test the get_component_vector function
    print(get_component_vector(sys))
    print(get_component_vector(av,1))
    

