#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 15:44:21 2017

@author: litaochen
"""

import re
#import scipy.misc as misc
import matplotlib.pyplot as plt
import numpy as np
import Reader as rd
from glob import glob
#from PIL import Image

ANGLE = 0 # 0-15, 0 is front, 8 is back
OUT_RADIUS = 50 #radius of the outer circle
IN_RADIUS = 20 #radius of the inner circle
STEP = 10 #moving step both vertically and horizontally

MEAN_DIFF = 0.6 #Threshold of relative difference in mean
VAR_DIFF = 0.6 #Threshold of relative difference in varience

def ring_detect(img):
    heat_map = np.zeros(img.shape)
    backgroud = np.zeros(img.shape) + 1e-6
    
    #circle map
    xx, yy = np.mgrid[:OUT_RADIUS*2, :OUT_RADIUS*2]
    circle = (xx - OUT_RADIUS) ** 2 + (yy - OUT_RADIUS) ** 2
    roi = circle <= IN_RADIUS ** 2
    ring = np.logical_and(circle <= OUT_RADIUS ** 2, circle > IN_RADIUS ** 2)
    roi_cord = np.array(np.nonzero(roi)) - OUT_RADIUS # cord is np.adarray
    ring_cord = np.array(np.nonzero(ring)) - OUT_RADIUS
    
    #Background Threshold
    
    
    for i in list(range(OUT_RADIUS, img.shape[0]-OUT_RADIUS, STEP))+[img.shape[0]-OUT_RADIUS-1]:
        for j in list(range(OUT_RADIUS, img.shape[1]-OUT_RADIUS, STEP))+[img.shape[1]-OUT_RADIUS-1]:
            roi_cord_ = tuple(roi_cord + [[i], [j]])
            ring_cord_ = tuple(ring_cord + [[i], [j]])
            
            backgroud[roi_cord_] += 1
            
            roi_value = img[roi_cord_]
            ring_value = img[ring_cord_]
            
            mean_comp = [np.mean(roi_value), np.mean(ring_value)]
            var_comp = [np.var(roi_value), np.var(ring_value)]
            
            # camparing with given threshold 
            if np.min(mean_comp)/np.max(mean_comp) < MEAN_DIFF \
            and np.min(var_comp)/np.max(var_comp) < VAR_DIFF:
                heat_map[roi_cord_] += 1
    
    return heat_map/backgroud

#def main():
path_list = glob("./stage1_aps/*.aps")
for path in path_list[1:2]:
    data = rd.read_data(path)
    sample_id = re.search(r"(?<=s/)[^.]*", path).group(0)
    
    #we only get images from one certain angle
    img = np.rot90(data[:,:, ANGLE])
    imgplt = plt.imshow(img)
    plt.show()
    
    heat_map = ring_detect(img)
    
    plt.imshow(heat_map)
    plt.show()
    
#if __name__ == "__main__":
#    main()
    
    
    