# -*- coding: utf-8 -*-
import sys
import numpy as np

from compound_models import *
from Compound import get_component_vector

if __name__ == '__main__':

    #chems = 'MgB2'
    chems = sys.argv[1]
    x_predict = get_component_vector(chems)
    x_predict = np.array(x_predict).reshape(1,10,10,1)
    #Loading models
    print("Loading models...")
    
    com_fenergy = round(predict_model_compound_formation_energy_cnn(x_predict),3)
    com_bandgap = round(predict_model_compound_bandgap_cnn(x_predict),3)
    com_superconductivity = round(predict_model_compound_superconductivity_cnn(x_predict),2) 
    print("Results for Compound:", chems)
    print("The formation energy is :                      ",com_fenergy,"eV/atom")
    print('The band gap is:                               ',com_bandgap,'eV')
    print('The superconducting transition temperature is: ', com_superconductivity,'K')
    
