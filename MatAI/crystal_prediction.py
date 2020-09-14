# -*- coding: utf-8 -*-
import sys
from ReadPOSCAR import deal_poscars
from Model_input import get_input
from crystal_models import *


if __name__ == '__main__':
    #Read the crystal file and create crystal object 
    poscar = sys.argv[1]
    print("Reading the crystal structure...")
    cry = deal_poscars(poscar)
    #cry = deal_poscars('POSCAR','VASP4')
    chems = cry[0].get_chemical_formula()
    avectors = get_input(cry[0])
    #Loading models
    print("Loading models...")
    cry_energy = round(predict_model_crystal_energy(avectors),3)
    cry_fenergy = round( predict_model_crystal_formation_energy(avectors),3)
    cry_bandgap = round(predict_model_crystal_bandgap(avectors),3)
    
    cry_bulk_modulus = round(predict_model_crystal_bulk_modulus(avectors),2)
    cry_shear_modulus = round(predict_model_crystal_shear_modulus(avectors),2)
    cry_Cv = round(predict_model_crystal_heat_capacity_Cv(avectors),2)
    cry_Cp = round(predict_model_crystal_heat_capacity_Cp(avectors),2)
    cry_debye_T = round(predict_model_crystal_debye_temperature(avectors),2)
    cry_gruneisen = round(predict_model_crystal_gruneisen_parameter(avectors),2)
    cry_expan = round(predict_model_crystal_thermal_expansion(avectors),2)
    cry_thermcond = round(predict_model_crystal_thermal_conductivity(avectors),2)
    
    
    print("Results for Crystal: ",chems)
    print('The total energy is:                           ',cry_energy,'eV/atom')
    print("The formation energy is :                      ",cry_fenergy,"eV/atom")
    print('The band gap is:                               ',cry_bandgap,'eV')
    print("The bulk modulus is :                          ",cry_bulk_modulus,"GPa")
    print("The shear modulus is :                         ",cry_shear_modulus,"GPa")
    print('The Cv (300K) is:                              ',cry_Cv,'kB/atom')
    print('The Cp (300K)is:                               ',cry_Cp,'kB/atom')
    print('The Debye temperature is:                      ',cry_debye_T,'K')
    print('The gruneisen constant is:                     ',cry_gruneisen)
    print('The thermal expansion parameter (300K) is:     ', cry_expan,'1/K')
    print('The thermal conductivity (300K) is:            ', cry_thermcond,'10^(-5)W/m*K')