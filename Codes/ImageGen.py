#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 18:06:04 2017

@author: litaochen
"""

import os
import re
import scipy
import numpy as np
import Reader as rd
from glob import glob
#from PIL import Image

ANGLE = 0

path_list = glob("./stage1_aps/*.aps")
img_list = []
for path in path_list:
    data = rd.read_data(path)
    sample_id = re.search(r"(?<=s/)[^.]*", path).group(0)
    
    #we only get images from one certain angle
    img = np.rot90(data[:,:, ANGLE])
    img_list.append(img)
    
#    try:
#    #    img.save( sample_id + ".png", "")
#        scipy.misc.imsave("./front/" + sample_id + ".png", img)
#    except:
#        os.makedirs("./front")
#        scipy.misc.imsave("./front/" + sample_id + ".png", img)
        
img_array = np.array(img_list)
scipy.io.savemat('Front.mat', mdict={'img_array': img_array})