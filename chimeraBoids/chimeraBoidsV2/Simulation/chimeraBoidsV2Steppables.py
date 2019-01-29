
from PySteppables import *
import CompuCell
import sys

import numpy as np
import os


global G_targetVolume_G
global G_lambdaVolume_G
global G_forceModulus_G
global G_deltaTime_G

global G_alphaBoids_G
global G_betaBoids_G
global G_gammaBoids_G
global G_deltaBoids_G


#parameter scan variables
G_targetVolume_G = 64.
G_lambdaVolume_G = 8.
G_forceModulus_G = -20
G_deltaTime_G = 10

G_alphaBoids_G = 5.5
G_betaBoids_G = 1.5
G_gammaBoids_G = .5
G_deltaBoids_G = .1




class chimeraBoidsV2Steppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        #initiating extra fields
        
        self.scalarAngleField = self.createScalarFieldCellLevelPy('angle')
        
        self.vectorIVelocityField = self.createVectorFieldCellLevelPy('Instant_Velocity')
        self.vectorTVelocityField = self.createVectorFieldCellLevelPy('dt_Velocity')
        
        self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
        
        
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = G_targetVolume_G
        self.lambdaVolume = G_lambdaVolume_G
        self.forceModulus = G_forceModulus_G
        self.deltaTime = G_deltaTime_G
        
        #boids parameters
        self.alphaBoids = G_alphaBoids_G
        self.betaBoids = G_betaBoids_G
        self.gammaBoids = G_gammaBoids_G
        self.deltaBoids = G_deltaBoids_G
        
        #set of used cluster IDs
        self.usedClusterIDs = set([0])#I think that to avoid possible errors it's best if the set is not
#                                   #empty
        #assigning parameters
        for cell in self.cellList:
            if cell:
                cell.targetVolume = self.targetVolume
                cell.lambdaVolume = self.lambdaVolume
                
                
                cell.dict['forceAngle'] = np.random.uniform(-np.pi,np.pi)
                cell.dict['angle']=cell.dict['forceAngle']
                
                cell.dict['centerMassX'] = cell.xCOM
                cell.dict['centerMassY'] = cell.yCOM
                cell.dict['centerMass'] = np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)
                
                cell.dict['positionX'] = [cell.xCOM]
                cell.dict['positionY'] = [cell.yCOM]
                
                cell.dict['velocityX_instant'] = 0.
                cell.dict['velocityY_instant'] = 0.
                
                cell.dict['velocityX_deltaT'] = 0.
                cell.dict['velocityY_deltaT'] = 0.
                
                cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                #cluster ID
                cell.dict['clusterID'] = None
    

    def positionTracking(self,mcs,cur_Cell):
        self.positionX = cur_Cell.dict['positionX']
        self.positionY = cur_Cell.dict['positionY']
        
        dx = cur_Cell.xCOM - cur_Cell.dict['centerMassX']
        dy = cur_Cell.yCOM - cur_Cell.dict['centerMassY'] 
        
        #if the cell "moved" more than half the space it whent trough the border
        if abs(dx)>0.5*self.dim.x:
            if dx < 0:
                dx +=self.dim.x
            else:
                dx -= self.dim.x
        if abs(dy)>0.5*self.dim.y:
            if dy<0:
                dy +=self.dim.y
            else:
                dy -=self.dim.y  
        
        #new absolute position is previous position + fixed dx
        cur_Cell.dict['positionX'].append(self.positionX[-1]+dx)
        cur_Cell.dict['positionY'].append(self.positionY[-1]+dy)
        
        #with a dt = 1mcs v = dx
        cur_Cell.dict['velocityX_instant'] = dx
        cur_Cell.dict['velocityY_instant'] = dy
        
        #angle of movement:
        cur_Cell.dict['angle'] = np.arctan2(dy,dx)
        
        cur_Cell.dict['centerMassX'] = cur_Cell.xCOM
        cur_Cell.dict['centerMassY'] = cur_Cell.yCOM
    
    
    def getNeighboursMeanVelocity(self,cur_Cell):
        #list of neighbours' velocities
        neigVxList = np.array([])
        neigVyList = np.array([])
        
        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cur_Cell):
            if neighbor:
                neigVxList = np.append(neigVxList,neighbor.dict['velocityX_deltaT'])
                neigVyList = np.append(neigVyList,neighbor.dict['velocityY_deltaT']) 
        
        #taking care of empty lists
        if len(neigVxList):
            neigVx = np.mean(neigVxList)
        else:       
            neigVx = 0
        if len(neigVyList):
            neigVy = np.mean(neigVyList)
        else:
            neigVy = 0
        
        return neigVx, neigVy
        
    def calculateForce(self,cV,nV,pF):
        force = self.alphaBoids*cV + self.betaBoids*nV - self.deltaBoids*pF
        return force
    
    def updateFields(self,mcs):
        
        for cCell in self.cellList:
            #updating scalar field for intantaneous angle of movement:
            self.scalarAngleField[cCell] = 180*(cCell.dict['angle']+np.pi)/np.pi
            #updating vector field for intantaneous velocity:
            self.vectorIVelocityField[cCell] = [cCell.dict['velocityX_instant'], cCell.dict['velocityX_instant'], 0]
            
            if mcs>self.deltaTime:
                self.vectorTVelocityField[cCell] = [cCell.dict['velocityX_deltaT'], cCell.dict['velocityX_deltaT'], 0]
                self.vectorForceField[cCell] = [cCell.dict['previousForceX'], cCell.dict['previousForceY'], 0]
        
#         self.vectorTVelocityField = self.createVectorFieldCellLevelPy('dt_Velocity')
        
#         self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
    
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        #############################################################
        #first cell loop to update tracking of position and velocities
        for cell in self.cellList:
            if cell:
                self.positionTracking(mcs,cell)
                
                
                
                if mcs > self.deltaTime:
                    vx = ( cell.dict['positionX'][-1] - cell.dict['positionX'][-self.deltaTime-1] )
                    vx = vx/self.deltaTime
                    
                    vy = ( cell.dict['positionY'][-1] - cell.dict['positionY'][-self.deltaTime-1] )
                    vy = vy/self.deltaTime
                    
                    cell.dict['velocityX_deltaT'] = vx
                    cell.dict['velocityY_deltaT'] = vy
        ##############################################################
        
        #################################################
        #second cell loop to calculate the boid's force 
        if mcs > self.deltaTime:
            for cell in self.cellList:
                
                cellVx = cell.dict['velocityX_deltaT']
                cellVy = cell.dict['velocityY_deltaT']
                
                #getting mean velocities of neighbours
                neigVx, neigVy = self.getNeighboursMeanVelocity(cell)
                
                #calculating new force in x and y
                forceX = self.calculateForce(cellVx,neigVx,cell.dict['previousForceX'])                
                forceY = self.calculateForce(cellVy,neigVy,cell.dict['previousForceY'])                
                
                #getting the angle of new force
                # # not sure about the implementation of noise
                cell.dict['forceAngle'] = np.arctan2(forceY,forceX) + self.gammaBoids * np.random.uniform(-np.pi,np.pi)
                
                #aplying standard force with new angle
                cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle'])  
               
                cell.dict['previousForceX'] = forceX
                cell.dict['previousForceY'] = forceY 
        #####################################################
        
        ##I need to identify the clusters. The most straitforward way I can think is
        ## simply looping the cells, checking the neighbors and assingning a colective ID. If the
        ## neighbor already has an ID use that one instead. Of course, this may lead to a single
        ## cluster having multiple IDs, so a second pass will be necessary. 
        ## there's also the issue of tracking the used IDs.
        #self.usedClusterIDs
        
        #first loop for id assignment
        for cell in self.cellList:
            if len(self.getCellNeighborDataList(cell)) == 0:
                cell.dict['clusterID'] = None
                break
            
            neigsIDs = []
            if cell.dict['clusterID'] != None:
                neigsIDs.append(cell.dict['clusterID'])
            
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor: 
                    if neighbor.dict['clusterID'] != None:
                        neigsIDs.append(neighbor.dict['clusterID'])
            
            if len(neigsIDs) == 0:
                newID = max(self.usedClusterIDs)+1
                self.usedClusterIDs.add(newID)
                cell.dict['clusterID'] = newID
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor: 
                        neighbor.dict['clusterID'] = newID
            else:
                id = min(neigsIDs)
                cell.dict['clusterID'] = id
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor: 
                        neighbor.dict['clusterID'] = id
                
        
        #second loop for id fixing
#         for 
        
        #updating extra fields
        self.updateFields(mcs)
            
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        