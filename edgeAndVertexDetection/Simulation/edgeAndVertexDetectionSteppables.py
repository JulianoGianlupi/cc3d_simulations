
from PySteppables import *
import CompuCell
import sys
import numpy as np
import numpy.random as rng
class edgeAndVertexDetectionSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        
        
        
        
        for cell in self.cellList:
#             cell.targetVolume = cell.volume+rng.normal(0,10)
#             cell.lambdaVolume = 2
            pass
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        pass
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        
from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
import CompuCellSetup
# from math import *


class extraFieldsManager(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarFieldPixelNeig = self.createScalarFieldPy("Bordering Pixels")
        
        
    def start(self):
        print "extraFieldsManager: This function is called once before simulation"
        #order of pixels to look into
        self.boundaryStrategy=CompuCell.BoundaryStrategy.getInstance()
        ##2nd order
        self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromNeighborOrder(2)
    def step(self,mcs):
        
        #clears the field
        self.scalarFieldPixelNeig[:, :, :] =  0
        fourFold = {}
        threeFold = {}
        for cell in self.cellList:
            #creates/erases the cell dict for faces
            #this is a dict of dicts, as each face will have the id of the neigh
            #cell as id.
            cell.dict['faces'] = {} 
            cell.dict['edges'] = {} 

            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:

                    #creates a dict for every face
                    #a list of its pixels
                    cell.dict['faces'][str(neighbor.id)] = []
            

            #gets the boundary px for the cell
            pixelList = self.getCellBoundaryPixelList(cell)
            edgePixels =[] ## obsolite maybe
            for boundaryPixelTrackerData in pixelList:
                
                #counter to see number of neigh cells
                numberOfDifferentCellNeighbors = []
                
                #gets the actual tuple of x,y,z
                actualPixel = boundaryPixelTrackerData.pixel
                
                #goes trough neighboring pixels
                for i in xrange(self.maxNeighborIndex+1):
                    
                    #gets the data for the neigh px
                    neighborPixelData = self.boundaryStrategy.getNeighborDirect(actualPixel,i)
                    
                    #gets the cell that owns that px
                    cell2=self.cellField.get(neighborPixelData.pt)
                    

                    #checks if cell2 is actually a neigbor (don't know how it could not be, but sometimes it's not)
                    #checks if it hasn't been counted yet
                    if  (cell2 
                            and str(cell2.id) in cell.dict['faces'] 
                            and actualPixel not in cell.dict['faces'][str(cell2.id)]):
                        cell.dict['faces'][str(cell2.id)].append(actualPixel)
                    
                    #if cell2 not medium and not counted already, count it
                    if (cell2                         
                        and cell2.id not in numberOfDifferentCellNeighbors):
                        numberOfDifferentCellNeighbors.append(cell2.id)
                
                
                           
                
                #separating x, y, z
                ptx = actualPixel.x
                pty = actualPixel.y
                ptz = actualPixel.z
                
                #paints the field with the number of neighbors
                self.scalarFieldPixelNeig[ptx, pty, ptz] = len(
                                    numberOfDifferentCellNeighbors)
                
                #now lets look at vertecies
                #1st i'll crete a list of pixels to create the vertex latter
                #one of the lists will be for the 4fold pixel and another for the 3fold
                if ((len(numberOfDifferentCellNeighbors)>3) and 
                    ((actualPixel.x,actualPixel.y,actualPixel.z) not in fourFold.keys())):
                    
                    
                    orgID = sorted(numberOfDifferentCellNeighbors)
                    
                    fourFold[(actualPixel.x,actualPixel.y,actualPixel.z)] = tuple(orgID)
                if ((len(numberOfDifferentCellNeighbors)==3) and 
                    ((actualPixel.x,actualPixel.y,actualPixel.z) not in threeFold.keys())):
                    
                    
                    orgID = sorted(numberOfDifferentCellNeighbors)
                    
                    threeFold[(actualPixel.x,actualPixel.y,actualPixel.z)] = tuple(orgID)
                #####
                ###OLD METHOD
#                 if len(numberOfDifferentCellNeighbors)>3:
 
#                     vertexID =''
#                     orgID = sorted(numberOfDifferentCellNeighbors)
#                     for neigID in orgID:
#                         vertexID+=str(neigID)+','

#                     existingVIDs = list(cell.dict['edges'].keys())
#                     if len(cell.dict['edges']) == 0:
#                         print 'new edge', vertexID
#                         cell.dict['edges'][vertexID] =[]
#                         cell.dict['edges'][vertexID].append(actualPixel)

#                     else:
#                         for evID in cell.dict['edges']:
#                             if vertexID in evID:
#                                 print 'found', vertexID, 'in', evID
#                                 cell.dict['edges'][evID].append(actualPixel)

        
        #now I iterate through the marked pixels to create the vertecies
        #1st through the 4fold, as those will be the centers.
        verticies = {}
        numberOfCreatedVetex = 1
        while numberOfCreatedVetex >0:
            numberOfCreatedVetex = 0
            for vertexPx_1 in fourFold.keys():
                for vertexPx_2 in fourFold.keys(): #getting the pairs
                    if ((vertexPx_1 != vertexPx_2) and # a px with itself makes no sense
                        (fourFold[vertexPx_1] == fourFold[vertexPx_2])): 
                            #if they have the same neighboring cells they are in the same vertex
                        if fourFold[vertexPx_1] not in verticies.keys():# if this vertex hasn't been created create it
                            verticies[fourFold[vertexPx_1]] = [vertexPx_1,vertexPx_2]#save the pixels
                            numberOfCreatedVetex+=1
                        else:
                            if vertexPx_1 not in verticies[fourFold[vertexPx_1]]:
                                verticies[fourFold[vertexPx_1]].append(vertexPx_1)
                                numberOfCreatedVetex+=1
                            if vertexPx_2 not in verticies[fourFold[vertexPx_1]]:                                
                                verticies[fourFold[vertexPx_1]].append(vertexPx_2)
                                numberOfCreatedVetex+=1
            print   "joined", numberOfCreatedVetex, "4 fold pixels to vertex"
        numberOfCreatedVetex = 1
        while numberOfCreatedVetex >0:
            numberOfCreatedVetex = 0
            for vertexPx_1 in threeFold.keys():
                inExistingVertex = False
                for existingVertex in verticies.keys():
                    if threeFold[vertexPx_1] in existingVertex:
                        inExistingVertex = True
                        if vertexPx_1 not in verticies[existingVertex]:
                            verticies[existingVertex].append(vertexPx_1)
                            numberOfCreatedVetex+=1
                if not inExistingVertex:
                    for vertexPx_2 in threeFold.keys(): 
                        if ((vertexPx_1 != vertexPx_2) and # a px with itself makes no sense
                            (threeFold[vertexPx_1] == threeFold[vertexPx_2])): 
                                #if they have the same neighboring cells they are in the same vertex
                                if threeFold[vertexPx_1] not in verticies.keys():# if this vertex hasn't been created create it
                                    verticies[threeFold[vertexPx_1]] = [vertexPx_1,vertexPx_2]#save the pixels
                                    numberOfCreatedVetex+=1
                                else:
                                    if vertexPx_1 not in verticies[threeFold[vertexPx_1]]:
                                        verticies[threeFold[vertexPx_1]].append(vertexPx_1)
                                        numberOfCreatedVetex+=1
                                    if vertexPx_2 not in verticies[threeFold[vertexPx_1]]:                                
                                        verticies[threeFold[vertexPx_1]].append(vertexPx_2)
                                        numberOfCreatedVetex+=1
            print   "joined", numberOfCreatedVetex, "3 fold pixels to vertex"
#             for edgePixel in edgePixels: ## obsolite maybe
#                 for i in xrange(self.maxNeighborIndex+1):
#                     neighborPixelData = self.boundaryStrategy.getNeighborDirect(edgePixel,i)

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
 
 

from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
import CompuCellSetup



class scatterFaceVertex(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        self.pWVertex = self.addNewPlotWindow(_title='Vertex CM Position', _xAxisTitle='X coord',
                                        _yAxisTitle='Y coord', _xScaleType='linear', _yScaleType='linear')
        self.pWVertex.addPlot('Vertex_CM', _style='Dots', _color='red', _size=5)
        
        
        
        
        
    def step(self,mcs):
        self.pWVertex.eraseAllData()
        for cell in self.cellList:
            vertecies = cell.dict['edges']
            for key in vertecies:
                vCMx = 0
                vCMy = 0
                vertex = vertecies.get(key)
                for px in vertex:
                    vCMx += px.x/len(vertex)
                    vCMy += px.y/len(vertex)
                
            
                    #vCMx += int(px[0])/len(key)
                    #vCMy += int(px[1])/len(key)
                    
                
                self.pWVertex.addDataPoint("Vertex_CM", vCMx, vCMy)
            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
    
