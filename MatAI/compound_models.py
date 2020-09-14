# -*- coding: utf-8 -*-

from keras.models import load_model



def predict_model_compound_bandgap_cnn(x_predict):
    model = load_model('model.save/model_bandgap_cnn.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_compound_formation_energy_cnn(x_predict):
    model = load_model('model.save/model_formation_energy_cnn.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
def predict_model_compound_superconductivity_cnn(x_predict):
    model = load_model('model.save/model_superconductivity_cnn.h5')
    prediction = model.predict(x_predict)
    return prediction[0][0]
