
from PySteppables import *
import CompuCell
import sys

import numpy as np
import os


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
        self.targetVolume = 64.
        self.lambdaVolume = 8.
        self.forceModulus = -20
        self.forceCounter = 0
        self.deltaTime = 10
        
        #boids parameters
        self.alphaBoids = 5.5
        self.betaBoids = 1.5
        self.gammaBoids = .5
        self.deltaBoids = .1
        
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
        self.updateFields(mcs)
            
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        