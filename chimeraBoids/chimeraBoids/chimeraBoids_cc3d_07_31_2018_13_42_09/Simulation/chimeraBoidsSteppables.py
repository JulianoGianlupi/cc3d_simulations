from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import numpy as np
import os

class chimeraBoidsSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField = self.createScalarFieldCellLevelPy("Angle")
        
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = 100.
        self.lambdaVolume = 8.
        self.forceModulus = -10
        self.forceCounter = 0
        self.deltaTime = 5
        #boids parameters
        self.alphaBoids = 0
        self.betaBoids = 1.
        self.gammaBoids = .5
        self.deltaBoids = 0
        
        
#         self.pWAngle = self.addNewPlotWindow(_title='Angle (rad) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
#                                         _yAxisTitle='Angle (radians)', _xScaleType='linear', _yScaleType='linear')
        
#         self.pWAngle.addPlot('Angle (rad)', _style='Dots', _color='red', _size=5)
        
        #cell initiation and dictionaries creation
        for cell in self.cellList:
            cell.targetVolume = self.targetVolume
            cell.lambdaVolume = self.lambdaVolume
            
            cell.dict['angle']=0.
            cell.dict['forceAngle'] = np.random.uniform(0.,2.*np.pi)
            
            cell.dict['centerMassX'] = [cell.xCOM]
            cell.dict['centerMassY'] = [cell.yCOM]
            cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
            
            cell.dict['velocityX'] = 0.
            cell.dict['velocityY'] = 0.
            
            cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
            cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle']) 
            
            cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
            cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
            #cell.lambdaVecZ = 0.0  
            #
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        for cell in self.cellList:
            if cell:
                self.scalarCLField[cell] = cell.dict['angle']
        #
        for cell in self.cellList:
            self.centerMassX = cell.dict['centerMassX']
            self.centerMassY = cell.dict['centerMassY']
            
            
            
            
            
            
            if len(self.centerMassX)>self.deltaTime:
                dx = cell.xCOM - self.centerMassX[-self.deltaTime]
                dy = cell.yCOM - self.centerMassY[-self.deltaTime]
                cell.dict['angle'] = np.angle(complex(dx,dy))
                
                cell.dict['velocityX'] = (cell.xCOM - self.centerMassX[-1])
                cell.dict['velocityY'] = (cell.yCOM - self.centerMassY[-1])
                
            cm = np.sqrt( (cell.xCOM)*(cell.xCOM) + (cell.yCOM)*(cell.yCOM))
            cell.dict['centerMassX'].append(cell.xCOM)            
            cell.dict['centerMassY'].append(cell.yCOM)
            cell.dict['centerMass'].append(cm)
            
        #boids force definition
        for cell in self.cellList:
            cellVx = cell.dict['velocityX']
            cellVy = cell.dict['velocityY']
            neigVxList = np.array([])
            neigVyList = np.array([])
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:
                    neigVxList = np.append(neigVxList,neighbor.dict['velocityX'])
                    neigVyList = np.append(neigVyList,neighbor.dict['velocityY'])
            neigVx = np.mean(neigVxList)
            neigVy = np.mean(neigVyList)
            
            forceX = self.alphaBoids*cellVx + self.betaBoids*neigVx - self.deltaBoids*cell.dict['previousForceX']
            forceY = self.alphaBoids*cellVy + self.betaBoids*neigVy - self.deltaBoids*cell.dict['previousForceY']
            
            cell.dict['previousForceX'] = forceX
            cell.dict['previousForceY'] = forceY
            cell.dict['forceAngle'] = np.angle(complex(forceX,forceY)) + self.gammaBoids * np.random.uniform(-np.pi,np.pi)
            
            cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
            cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        