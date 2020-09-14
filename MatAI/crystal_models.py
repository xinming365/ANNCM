# -*- coding: utf-8 -*-
import numpy as np
from Auxilfunction import zscore


from keras.models import load_model

# for crystal

def predict_model_crystal_bandgap(avectors):
    mean = np.load('model.save\\crystal_bandgap_mean.npy')
    std = np.load('model.save\\crystal_bandgap_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_bandgap.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_bulk_modulus(avectors):
    mean = np.load('model.save\\crystal_bulk_modulus_mean.npy')
    std = np.load('model.save\\crystal_bulk_modulus_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_bulk_modulus.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_debye_temperature(avectors):
    mean = np.load('model.save\\crystal_debye_temperature_mean.npy')
    std = np.load('model.save\\crystal_debye_temperature_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_debye_temperature.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_energy(avectors):
    mean = np.load('model.save\\crystal_energy_mean.npy')
    std = np.load('model.save\\crystal_energy_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_energy.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_formation_energy(avectors):
    mean = np.load('model.save\\crystal_formation_energy_mean.npy')
    std = np.load('model.save\\crystal_formation_energy_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_formation_energy.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_gruneisen_parameter(avectors):
    mean = np.load('model.save\\crystal_gruneisen_parameter_mean.npy')
    std = np.load('model.save\\crystal_gruneisen_parameter_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_gruneisen_parameter.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_heat_capacity_Cp(avectors):
    mean = np.load('model.save\\crystal_heat_capacity_Cp_mean.npy')
    std = np.load('model.save\\crystal_heat_capacity_Cp_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_heat_capacity_Cp.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_heat_capacity_Cv(avectors):
    mean = np.load('model.save\\crystal_heat_capacity_Cv_mean.npy')
    std = np.load('model.save\\crystal_heat_capacity_Cv_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_heat_capacity_Cv.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_shear_modulus(avectors):
    mean = np.load('model.save\\crystal_shear_modulus_mean.npy')
    std = np.load('model.save\\crystal_shear_modulus_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_shear_modulus.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_thermal_conductivity(avectors):
    mean = np.load('model.save\\crystal_thermal_conductivity_mean.npy')
    std = np.load('model.save\\crystal_thermal_conductivity_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_thermal_conductivity.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_crystal_thermal_expansion(avectors):
    mean = np.load('model.save\\crystal_thermal_expansion_mean.npy')
    std = np.load('model.save\\crystal_thermal_expansion_std.npy')
    x_predict = zscore(mean,std,avectors)
    x_predict = np.array(x_predict).reshape(1,len(x_predict))
    model = load_model('model.save\\model_crystal_thermal_expansion.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
