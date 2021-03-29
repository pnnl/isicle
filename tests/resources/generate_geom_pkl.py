# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 13:23:12 2021

@author: nune558
"""
from isicle import geometry

geom = geometry.load('geom_test.smi')
geom.global_properties['load']['path'] = 'geom_test.smi'
geom.save_pickle('geom_test.pkl')