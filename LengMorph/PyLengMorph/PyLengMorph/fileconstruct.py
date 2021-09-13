# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:59:04 2020

@author: Dylan Agius
"""

 
import numpy as np
import pandas as pd 

from PyLengMorph.scrape import data_scrape
from PyLengMorph.node_increase import increase_nodes
    

"""find grain boundary nodes corresponding to each grain"""    
    
def grainboundary(**kwargs):
    
    cname=kwargs['file']
    loc=kwargs['loc']
    
    nodes, elemes, gbels, centroid, orien = data_scrape(loc,cname)
    
   
    
    #check to see if the node increase function was used to determine
    #how many nodes per voxel
    if kwargs['nodeinc']==False :
       nodenum=8
    else:
       nodenum=27
       elemes,nodes=increase_nodes(nodes,elemes)

    #find all elements at a boundary and order according to feature id
    
    bfeature=gbels.loc[gbels['location']==0,'feature']
    indexfeature=gbels.loc[gbels['location']==0,'feature'].index.tolist()
    
    #creating an array from this info which has the element number and grain feature this belongs to 
    #respectively
    
    el_feature=np.array([list(np.asarray(indexfeature)+1),bfeature]).T
    
    #Now create an array of grains to search through the 'el_feature' array to find
    #the elements that belongs to it.  Add these elements to the row which  feature
    #the element belongs to.
    
    searchfeat=np.arange(1,np.amax(el_feature[:,1],axis=0)+1)
    
    #Search through the el_feature array for each individual value in searchfeat
    #and find when this feature appears in the array and assign the element from the
    #first column of the array.
    
    #Now elloc has the features which are the grains but are 1 value out since
    #the starting index is zero.  The rows correspond to the elements which are on the boundary
    #for that feature.
    
    elloc=[[None]*len(searchfeat)]*len(searchfeat)
    featelnode=[[None]*len(searchfeat)]*len(searchfeat)
    intersect=[[]*len(searchfeat)]*len(searchfeat)
    featurenodes=[[]*len(searchfeat)]*len(searchfeat)
    xnodes=[[]*len(searchfeat)]*len(searchfeat)
    ynodes=[[]*len(searchfeat)]*len(searchfeat)
    znodes=[[]*len(searchfeat)]*len(searchfeat)
    boundfeat=[[]*len(searchfeat)]*len(searchfeat)

 
    
    #converting dataframe to array to see if it is easier for me to use
    
    newelemes=elemes
    newnodes=nodes
    
    # finding the centroid nodal coordinates for each element
    nodecoord=[[]*len(newelemes[:,0])]*len(newelemes[:,0])
    avcoord=[[]*3]*len(newelemes[:,0])
    
    for i in range(0,len(newelemes[:,0])):
        nodecoord[i]=newnodes[(newelemes[i,1:]-1).astype(np.int64)]
        
        avcoord[i]=sum(nodecoord[i][:,1:])/nodenum
        
       # *******need to add a tool in here to decipher whether increasing nodes***
        
        
    for i in range(0,len(searchfeat)):
        
        elloc[i]=el_feature[np.where(el_feature[:,1]==searchfeat[i]),0]
        
       # Now we need to find the nodes which correspond to those elements for each feature/grain"""
        #featelnode[i]=newelemes[list((np.asarray(elloc[i])-1).astype(np.int64)),:][0,:,1:9]
        featelnode[i]=newelemes[list((np.asarray(elloc[i])-1).astype(np.int64)),:][0,:,1:]
    
    
    #Now we need to see which features have overlapping nodes with other features as these 
    #nodes will then belong to the very outer surface."""

    for i in range(0,len(searchfeat)):
        for j in range(0,len(searchfeat)):
            intersectarray=np.intersect1d(featelnode[i],featelnode[j])
            if np.size(intersectarray) > 0 :
                if i != j :
                    intersect[i]=np.append(intersect[i],intersectarray,axis=0)
                    boundfeat[i]=np.append(boundfeat[i],np.size(intersectarray)*[j+1])
                    
                    #Intersect array now has the nodes corresponding to each grain on the surface 
                    #of the grain.  The row location corresponds to the feature-1
                   
                    #now we need to find the x,y,z coordinates which correspond to these nodes
                    
                    featurenodes[i]=newnodes[np.asarray(intersect[i]-1).astype(int)]
                    xnodes[i]=featurenodes[i][:,1]
                    ynodes[i]=featurenodes[i][:,2]
                    znodes[i]=featurenodes[i][:,3]
                    
  
    
                   
    #Writing the node coordinates to separate txt files
    #firstly enumerate list to provide a dictionary depedent on the index of the list
    
    xnodelist = dict(enumerate(xnodes))
    ynodelist = dict(enumerate(ynodes))
    znodelist = dict(enumerate(znodes))
    boundfeat = dict(enumerate(boundfeat))
    
   
    #Then create a dataframe for export
    #'9999999's are added to rows missing values to ensure each row is the same dimension
    #for reading into crystal plasicity UMAT
    
    xnodedataframe= pd.DataFrame.from_dict(xnodelist, orient='index').fillna(9999999)
    ynodedataframe= pd.DataFrame.from_dict(ynodelist, orient='index').fillna(9999999)
    znodedataframe= pd.DataFrame.from_dict(znodelist, orient='index').fillna(9999999)
    boundfeatframe= pd.DataFrame.from_dict(boundfeat, orient='index').fillna(9999999)
    
    
    #saving as binary file
    #xvalues
    numRowarrx=np.array([np.size(xnodedataframe.values,0)],dtype=np.float64)
    numColarrx=np.array([np.size(xnodedataframe.values,1)],dtype=np.float64)
    xvalues=open('xvalues.bin','wb')
    numRowarrx.tofile(xvalues)
    numColarrx.tofile(xvalues)
    xnodedataframe.values.astype(np.float64).tofile(xvalues)
    xvalues.close()
    #np.savetxt('xvalues.txt', xnodedataframe.values, fmt='%f')
    
    #yvalues
    numRowarry=np.array([np.size(ynodedataframe.values,0)],dtype=np.float64)
    numColarry=np.array([np.size(ynodedataframe.values,1)],dtype=np.float64)
    yvalues=open('yvalues.bin','wb')
    numRowarry.tofile(yvalues)
    numColarry.tofile(yvalues)
    ynodedataframe.values.astype(np.float64).tofile(yvalues)
    yvalues.close()
    #np.savetxt('yvalues.txt', ynodedataframe.values, fmt='%f')
    
    #zvalues
    numRowarrz=np.array([np.size(znodedataframe.values,0)],dtype=np.float64)
    numColarrz=np.array([np.size(znodedataframe.values,1)],dtype=np.float64)
    zvalues=open('zvalues.bin','wb')
    numRowarrz.tofile(zvalues)
    numColarrz.tofile(zvalues)
    znodedataframe.values.astype(np.float64).tofile(zvalues)
    zvalues.close()
   #np.savetxt('zvalues.txt', znodedataframe.values, fmt='%f')
    
    #boundfeat
    numRowarrbf=np.array([np.size(boundfeatframe.values,0)],dtype=np.float64)
    numColarrbf=np.array([np.size(boundfeatframe.values,1)],dtype=np.float64)
    bfvalues=open('boundfeat.bin','wb')
    numRowarrbf.tofile(bfvalues)
    numColarrbf.tofile(bfvalues)
    boundfeatframe.values.astype(np.float64).tofile(bfvalues)
    bfvalues.close()
    #np.savetxt('boundfeat.txt',boundfeatframe.values,fmt='%f')

  
    elcentroid=np.array(avcoord)
    
   
    
    #saving the centroid coordinates for each voxel to text file
    #saving a binary file
    numRowarr=np.array([np.size(elcentroid,0)],dtype=np.float64)
    numColarr=np.array([np.size(elcentroid,1)],dtype=np.float64)
    centfile=open('el_centroid.bin','wb')
    numRowarr.tofile(centfile)
    numColarr.tofile(centfile)
    elcentroid.astype(np.float64).tofile(centfile)

    centfile.close()
    
    #creating a parameter file to read in array data if this is to be 
    #implemented in Abaqus
    if kwargs['abq']==True:
        paramfile= open(r"param_array.inc","w")
        
        paramfile.writelines("      PARAMETER (totalels=%d)\n" % numRowarr)
        paramfile.writelines("      PARAMETER (totalfeat=%d)\n" % numRowarrz)
        paramfile.writelines("      PARAMETER (nodeout=%d)" % numColarrz)
        
        paramfile.close()
    
    #creating an include file for fortan which contains the euler angles for 
    #each grain
    
    #creates the orientation include file in f90 format
    if kwargs['abq']==False :
    
        incfile = open(r"orien.inc","w")
        
        incfile.writelines("      orient=(/" )
        count=1
        for line in orien:
    
            for i in range(len(line)):
                if i != len(line)-1:
                    incfile.writelines("{0:.2f}".format(line[i]) + "," )
                    
                else: 
                    incfile.writelines("{0:.2f}".format(line[i]))
                    
            if count != len(orien):
        
                incfile.writelines(", &\n" )
                incfile.writelines("     " )
                count += 1
            else: 
                incfile.writelines("/)" )
    
        incfile.close()
        
    else:
        incfile = open(r"orien.inc","w")
        
        incfile.writelines("      orient=(/" )
        count=1
        for line in orien:
    
            for i in range(len(line)):
                if i != len(line)-1:
                    incfile.writelines("{0:.2f}".format(line[i]) + "," )
                    
                else: 
                    incfile.writelines("{0:.2f}".format(line[i]))
                    
            if count != len(orien):
        
                incfile.writelines(", \n" )
                incfile.writelines("     & " )
                count += 1
            else: 
                incfile.writelines("/)" )
    
        incfile.close()

    