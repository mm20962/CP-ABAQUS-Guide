# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:46:17 2020

@author: mi19356
"""

import numpy as np
import pandas as pd 
import os

def data_scrape(loc,cname):
    loc=os.path.join(loc,'data')
    
        
    nodes=np.array(file_read(os.path.join(loc,cname + '_nodes.inp'),',')).astype(np.float64)
    elemes=np.array(file_read(os.path.join(loc,cname + '_elems.inp'),',')).astype(np.float64)
    gbels=np.array(file_read(os.path.join(loc,cname + '_gbels.txt'),' '))
    gbels=pd.DataFrame(gbels,columns=["location","feature"],dtype='float')
    
  
    featureinfo=pd.read_csv(os.path.join(loc, cname + '.csv'),skiprows=1)
    centroid=np.asarray([featureinfo['Centroids_0'],featureinfo['Centroids_1'],featureinfo['Centroids_2']]).T
    orien=np.asarray([featureinfo['EulerAngles_0'],featureinfo['EulerAngles_1'],featureinfo['EulerAngles_2']]).T

    
    return nodes, elemes, gbels, centroid, orien

def file_read(fname,delim):
        with open(fname) as textFile:
            content_array = [line.split(delim) for line in textFile if not (line.startswith("*") or "GBManhattanDistances" in line)]
               
        return (content_array)
   