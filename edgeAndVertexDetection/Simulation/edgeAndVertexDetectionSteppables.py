
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
        self.scalarFieldVertices = self.createScalarFieldPy("Vertices")
        
        
    def start(self):
        print "extraFieldsManager: This function is called once before simulation"
        #order of pixels to look into
        self.boundaryStrategy=CompuCell.BoundaryStrategy.getInstance()
        ##2nd order
        self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromNeighborOrder(2)
        
        ##
        self.pWVertex = self.addNewPlotWindow(_title='Vertex CM Position', _xAxisTitle='X coord',
                                        _yAxisTitle='Y coord', _xScaleType='linear', _yScaleType='linear')
        self.pWVertex.addPlot('Vertex_CM', _style='Dots', _color='red', _size=5)
    def step(self,mcs):
        
        #clears the field
        self.scalarFieldPixelNeig[:, :, :] =  0
        self.scalarFieldVertices[:, :, :] =  0
        
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
                
                tuplePx = (ptx,pty,ptz)
                
                #paints the field with the number of neighbors
                self.scalarFieldPixelNeig[ptx, pty, ptz] = len(
                                    numberOfDifferentCellNeighbors)
                
                #now lets look at vertecies
                #1st i'll crete a list of pixels to create the vertex latter
                #one of the lists will be for the 4fold pixel and another for the 3fold
                if ((len(numberOfDifferentCellNeighbors)>3) and 
                    (tuplePx not in fourFold.keys())):
                    
                    
                    orgID = sorted(numberOfDifferentCellNeighbors)
                    
                    fourFold[tuplePx] = tuple(orgID)
                if ((len(numberOfDifferentCellNeighbors)==3) and 
                    (tuplePx not in threeFold.keys())):
                    
                    
                    orgID = sorted(numberOfDifferentCellNeighbors)
                    
                    threeFold[tuplePx] = tuple(orgID)
                #####
                ###OLD METHOD


        
        #now I iterate through the marked pixels to create the vertecies
        #1st through the 4fold, as those will be the centers.
        self.verticies = {}
        numberOfCreatedVetex = 1
        while numberOfCreatedVetex >0:
            numberOfCreatedVetex = 0
            for vertexPx_1 in fourFold.keys():
                for vertexPx_2 in fourFold.keys(): #getting the pairs
                    if ((vertexPx_1 != vertexPx_2) and # a px with itself makes no sense
                        (fourFold[vertexPx_1] == fourFold[vertexPx_2])): 
                        tuplePX1 = (vertexPx_1[0],vertexPx_1[1],vertexPx_1[2])
                        tuplePX2 = (vertexPx_2[0],vertexPx_2[1],vertexPx_2[2])
#                         print tuplePX1, tuplePX2
                            #if they have the same neighboring cells they are in the same vertex
                        if fourFold[vertexPx_1] not in self.verticies.keys():# if this vertex hasn't been created create it
#                             print fourFold[vertexPx_1]                           
                            self.verticies[fourFold[vertexPx_1]] = [tuplePX1,tuplePX2]#save the pixels
                            numberOfCreatedVetex+=1
                        else:
                            if vertexPx_1 not in self.verticies[fourFold[vertexPx_1]]:
                                self.verticies[fourFold[vertexPx_1]].append(tuplePX1)
                                numberOfCreatedVetex+=1
                            if vertexPx_2 not in self.verticies[fourFold[vertexPx_1]]:                                
                                self.verticies[fourFold[vertexPx_1]].append(tuplePX2)
                                numberOfCreatedVetex+=1
#             print   "joined", numberOfCreatedVetex, "4 fold pixels to vertex"
        numberOfCreatedVetex = 1
        while numberOfCreatedVetex >0:
            numberOfCreatedVetex = 0
            for vertexPx_1 in threeFold.keys():
                inExistingVertex = False
                tuplePX1 = (vertexPx_1[0],vertexPx_1[1],vertexPx_1[2])
#                 print tuplePX1
                for existingVertex in self.verticies.keys():
                    if threeFold[vertexPx_1] in existingVertex:
                        inExistingVertex = True
                        if vertexPx_1 not in self.verticies[existingVertex]:
                            self.verticies[existingVertex].append(tuplePX1)
                            numberOfCreatedVetex+=1
                if not inExistingVertex:
                    for vertexPx_2 in threeFold.keys(): 
                        tuplePX2 = (vertexPx_2[0],vertexPx_2[1],vertexPx_2[2])
#                         print tuplePX1
                        if ((vertexPx_1 != vertexPx_2) and # a px with itself makes no sense
                            (threeFold[vertexPx_1] == threeFold[vertexPx_2])): 
                                #if they have the same neighboring cells they are in the same vertex
                                if threeFold[vertexPx_1] not in self.verticies.keys():# if this vertex hasn't been created create it
                                    self.verticies[threeFold[vertexPx_1]] = [tuplePX1,tuplePX2]#save the pixels
                                    numberOfCreatedVetex+=1
                                else:
                                    if vertexPx_1 not in self.verticies[threeFold[vertexPx_1]]:
                                        self.verticies[threeFold[vertexPx_1]].append(tuplePX1)
                                        numberOfCreatedVetex+=1
                                    if vertexPx_2 not in self.verticies[threeFold[vertexPx_1]]:                                
                                        self.verticies[threeFold[vertexPx_1]].append(tuplePX2)
                                        numberOfCreatedVetex+=1
            print   "joined", numberOfCreatedVetex, "3 fold pixels to vertex"


        
        if mcs%10==0:
            self.pWVertex.eraseAllData()
            for vertexID, vertexPXs in self.verticies.iteritems():
    #             print len(vertexPXs)
                vertexCMx = 0
                vertexCMy = 0
                vertexCMz = 0
                pxColor = 0
                for j in  vertexID:
                    pxColor +=j 
                print pxColor
                for pixel in vertexPXs:
                    self.scalarFieldVertices[pixel[0], pixel[1], pixel[2]] = pxColor
                    vertexCMx += pixel[0]/len(vertexPXs)
                    vertexCMy += pixel[1]/len(vertexPXs)
                    vertexCMz += pixel[2]/len(vertexPXs)
                self.pWVertex.addDataPoint("Vertex_CM", vertexCMx, vertexCMy)
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
        pass

        
        
    def step(self,mcs):
        pass

            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
    
