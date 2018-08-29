
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
        for cell in self.cellList:
            #creates/erases the cell dict for faces
            #this is a dict of dicts, as each face will have the id of the neigh
            #cell as id.
            cell.dict['faces'] = {} 
            cell.dict['edges'] = {} 
#             neighborArea=0
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:
#                     neighborArea += commonSurfaceArea
                    #creates a dict for every face
                    cell.dict['faces'][str(neighbor.id)] = []
            
#             print neighborArea/cell.surface
            #gets the boundary px for the cell
            pixelList = self.getCellBoundaryPixelList(cell)
            edgePixels =[]
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
#                     print str(cell2.id) in cell.dict['faces']
                    #checks that the px is a face px not a edge px
                    #also checks if it hasn't been counted yet
                    if  (str(cell2.id) in cell.dict['faces'] 
                            and actualPixel not in cell.dict['faces'][str(cell2.id)]):
                        cell.dict['faces'][str(cell2.id)].append(actualPixel)
                    
                    #if cell2 not medium and not counted already, count it
                    if (cell2                         
                        and cell2.id not in numberOfDifferentCellNeighbors):
                        numberOfDifferentCellNeighbors.append(cell2.id)
                
                
                           
                ptx = actualPixel.x
                pty = actualPixel.y
                ptz = actualPixel.z
                self.scalarFieldPixelNeig[ptx, pty, ptz] = len(
                                    numberOfDifferentCellNeighbors)
                if len(numberOfDifferentCellNeighbors)>2:
                    
                    vertexID =''
                    orgID = sorted(numberOfDifferentCellNeighbors)
                    for neigID in orgID:
                        vertexID+=str(neigID)+','
                    
                    alreadyExists = False
                    
                    #checks if this vertex already exists or not. 
                    #Need to check if vertexID is a subset of existing vertex aswell
                    for existingVID in cell.dict['edges']:
                        if vertexID in existingVID:
                            alreadyExists = True
                            listpx = [actualPixel.x,actualPixel.y,actualPixel.z]
                            if actualPixel not in cell.dict['edges'][existingVID]:
                                cell.dict['edges'][existingVID].append(actualPixel)
                                #print actualPixel
                    if not alreadyExists:
#                         if vertexID not in cell.dict['edges']:
                        listpx = [actualPixel.x,actualPixel.y,actualPixel.z]
                        cell.dict['edges'][vertexID] = [] 
                        cell.dict['edges'][vertexID].append(actualPixel)
                        
                        alreadyExists = True
#                         elif actualPixel not in cell.dict['edges'][vertexID]:
#                             cell.dict['edges'][vertexID].append(actualPixel)
            
            for edgePixel in edgePixels:
                for i in xrange(self.maxNeighborIndex+1):
                    neighborPixelData = self.boundaryStrategy.getNeighborDirect(edgePixel,i)

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
    
