# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 12:02:24 2020

@author: Dylan Agius
"""
import numpy as np

from PyLengMorph.scrape import data_scrape


"""This function will increase the number of nodes defining the voxel
This may particularly important if you are using FFT solvers where you
dealing with voxels rather than nodes
"""

def increase_nodes(nodes,elemes):
    
   
    #increase the number of nodes by added a node between each defined node
    nodesnew=nodes[:,1:]
    
    #add nodes at every 0.5
    updatenew=[[]]*4
    
    #find maximum number of nodes
    
    findtotalmax=int(np.max(nodesnew))*2 +1
    
    column1=np.transpose(np.tile([np.arange(0,int(np.max(nodesnew[:,0]))+0.5,0.5)],int(((np.max(nodesnew[:,0])*2)+1)*((np.max(nodesnew[:,1])*2)+1))))
   
    column1=column1.reshape(np.size(column1),)
    
    column2= np.tile([i for i in np.arange(0,int(np.max(nodesnew[:,1]))+0.5,0.5) for _ in range(findtotalmax)],findtotalmax)

    column3=np.asarray([i for i in np.arange(0,int(np.max(nodesnew[:,2]))+0.5,0.5) for _ in range(int(((np.max(nodesnew[:,0])*2)+1)*((np.max(nodesnew[:,1])*2)+1))) ])

    column0=np.arange(1,np.size(column3)+1,1)
    updatenew[0]=column0
    updatenew[1]=column1
    updatenew[2]=column2
    updatenew[3]=column3
    updatenew=np.asarray(updatenew).T
    
    #update the element array to include extra nodes
    
    #start column values
    startnodes1=[1,2,3]
    maxlength=(int(np.max(nodesnew[:,0]))*2)+1
    startnodes2=[(startnodes1[0]+maxlength),(startnodes1[1]+maxlength),(startnodes1[2]+maxlength)]
    startnodes3=[(startnodes2[0]+maxlength),(startnodes2[1]+maxlength),(startnodes2[2]+maxlength)]
    
    maxlengthupdate=(maxlength*((int(np.max(nodesnew[:,1]))*2)+1))
    startnodes4=[(startnodes1[0]+maxlengthupdate),(startnodes1[1]+maxlengthupdate),(startnodes1[2]+maxlengthupdate)]
    startnodes5=[(startnodes4[0]+maxlength),(startnodes4[1]+maxlength),(startnodes4[2]+maxlength)]
    startnodes6=[(startnodes5[0]+maxlength),(startnodes5[1]+maxlength),(startnodes5[2]+maxlength)]
   
    maxlengthupdate2=(maxlengthupdate)*2
    startnodes7=[(startnodes1[0]+maxlengthupdate2),(startnodes1[1]+maxlengthupdate2),(startnodes1[2]+maxlengthupdate2)]
    startnodes8=[(startnodes7[0]+maxlength),(startnodes7[1]+maxlength),(startnodes7[2]+maxlength)]
    startnodes9=[(startnodes8[0]+maxlength),(startnodes8[1]+maxlength),(startnodes8[2]+maxlength)]
    
    totalnodearray=np.asarray(startnodes1+startnodes2+startnodes3+startnodes4+startnodes5+startnodes6+startnodes7+startnodes8+startnodes9)
    
    nodestotal=[[]]*len(elemes)
    extraval=0
    inc=0
    for i in range(0,len(elemes)):
        if i != 0 :
            if i % int(np.max(nodesnew[:,0]))==0 and i !=0 :
                if inc == int(np.max(nodesnew[:,1]))-1:
                   
  
                    extra=(maxlength*((int(np.max(nodesnew[:,1]))*2)+1)) + 2*maxlength +3
                    inc=0
                    totalnodearray=nodestotal[i-1] 
                    nodestotal[i]=np.append(nodestotal[i],nodestotal[i-1]+extra)
                else:
                    extraval=maxlength + 3
                    nodestotal[i]=np.append(nodestotal[i],nodestotal[i-1]+extraval)
                    inc += 1
            else:
               
                nodestotal[i]=np.append(nodestotal[i],nodestotal[i-1]+2)
        else:
            nodestotal[i]=totalnodearray
                
                
    nodestotal=np.asarray(nodestotal)
    elid=np.arange(1,len(elemes)+1,1).reshape(len(elemes),1)
    elemes=np.hstack((elid,nodestotal))
    
    return elemes, updatenew
    
   
    