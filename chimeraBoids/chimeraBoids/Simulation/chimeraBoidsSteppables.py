from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import numpy as np
import os
import math

class chimeraBoidsSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        #initiating the extra fields
        self.scalarCLField = self.createScalarFieldCellLevelPy("Angle")
        self.scalarVelocityField = self.createScalarFieldCellLevelPy("Speed")
        self.scalarForceField = self.createScalarFieldCellLevelPy("Force_Angle")
        self.vectorCLField = self.createVectorFieldCellLevelPy("Velocity")
        self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
        
        
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = 64.
        self.lambdaVolume = 8.
        self.forceModulus = -20
        self.forceCounter = 0
        self.deltaTime = 10
        self.maxReasonableSpeed = np.sqrt(self.dim.x*self.dim.x/4 +self.dim.y*self.dim.y/4)/self.deltaTime
        #boids parameters
        self.alphaBoids = 1.5
        self.betaBoids = 5.
        self.gammaBoids = .5
        self.deltaBoids = .1
        
        
        #plots
        
        self.pW = self.addNewPlotWindow(_title='Histogram', 
                                        _xAxisTitle='Speed',
                                        _yAxisTitle='Count',
                                        _xScaleType='linear', 
                                        _yScaleType='linear')
        self.pW.addHistogramPlot(_plotName='Hist 1', _color='green',
                                 _alpha=100)  
        
        
        
#         self.pWAngle = self.addNewPlotWindow(_title='Angle (rad) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
#                                         _yAxisTitle='Angle (radians)', _xScaleType='linear', _yScaleType='linear')
        
#         self.pWAngle.addPlot('Angle (rad)', _style='Dots', _color='red', _size=5)
        
        #cell initiation and dictionaries creation
        for cell in self.cellList:
            if cell:
                cell.targetVolume = self.targetVolume
                cell.lambdaVolume = self.lambdaVolume
                
                
                cell.dict['forceAngle'] = np.random.uniform(-np.pi,np.pi)
                cell.dict['angle']=cell.dict['forceAngle']
                
                cell.dict['centerMassX'] = [cell.xCOM]
                cell.dict['centerMassY'] = [cell.yCOM]
                cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
                
                cell.dict['positionX'] = [cell.xCOM]
                cell.dict['positionY'] = [cell.yCOM]
                
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
        
        #fields plotting
        self.cellsSpeeds = []
        for cell in self.cellList:
            if cell:
                self.scalarCLField[cell] = cell.dict['angle']/np.pi
#                 self.scalarVelocityField[cell] = np.sqrt(cell.dict['velocityX']*cell.dict['velocityX'] 
#                                             + cell.dict['velocityY']*cell.dict['velocityX'])
                self.scalarForceField[cell] = cell.dict['forceAngle']/np.pi
                self.vectorCLField[cell] = [cell.dict['velocityX'], cell.dict['velocityY'], 0]
                self.vectorForceField[cell] = [cell.lambdaVecX,cell.lambdaVecY,0]
        #
                
                speed = np.sqrt(cell.dict['velocityX']*cell.dict['velocityX'] 
                                            + cell.dict['velocityY']*cell.dict['velocityY'])
                if (float('-inf') < float(speed) < float('inf')):
                    if abs(float(speed)) < self.maxReasonableSpeed:
                        self.cellsSpeeds.append(speed)
                        self.scalarVelocityField[cell] = np.sqrt(
                                                  cell.dict['velocityX']*cell.dict['velocityX'] 
                                                + cell.dict['velocityY']*cell.dict['velocityY'])
        if mcs%50:
            self.pW.addHistogram(plot_name='Hist 1', 
                                 value_array=self.cellsSpeeds, 
                                 number_of_bins=20)
#                                  range = (np.nanmin(self.cellsSpeeds), np.nanmax(self.cellsSpeeds)))
              
        
        for cell in self.cellList:
            self.centerMassX = cell.dict['centerMassX']
            self.centerMassY = cell.dict['centerMassY']
            
            cm = np.sqrt( (cell.xCOM)*(cell.xCOM) + (cell.yCOM)*(cell.yCOM))
            cell.dict['centerMassX'].append(cell.xCOM)            
            cell.dict['centerMassY'].append(cell.yCOM)
            cell.dict['centerMass'].append(cm)
        
            
            
            
            if mcs>self.deltaTime:
                dx = cell.xCOM - self.centerMassX[-self.deltaTime]
                dy = cell.yCOM - self.centerMassY[-self.deltaTime]
                cell.dict['angle'] = np.angle(complex(dx,dy))
#                 if (self.centerMassX[-self.deltaTime] > .9*self.dim.x and
#                                             cell.xCOM < .1*self.dim.x):
#                     vlx = (cell.xCOM+self.dim.x - self.centerMassX[-self.deltaTime])/self.deltaTime                    
#                 elif (self.centerMassX[-self.deltaTime] < .1*self.dim.x and
#                                               cell.xCOM > .9*self.dim.x):
#                     vlx = (cell.xCOM-self.dim.x - self.centerMassX[-self.deltaTime])/self.deltaTime                   
#                 else:
#                     vlx = (cell.xCOM - self.centerMassX[-self.deltaTime])/self.deltaTime                                
#                 if (self.centerMassY[-self.deltaTime] > .9*self.dim.y and
#                                             cell.yCOM < .1*self.dim.y):
#                     vly = (cell.yCOM+self.dim.y - self.centerMassY[-self.deltaTime])/self.deltaTime                    
#                 elif (self.centerMassY[-self.deltaTime] < .1*self.dim.y and
#                                               cell.yCOM > .9*self.dim.y):
#                     vly = (cell.yCOM+self.dim.y - self.centerMassY[-self.deltaTime])/self.deltaTime                    
#                 else:
#                     vly = (cell.yCOM - self.centerMassY[-self.deltaTime])/self.deltaTime                    
                
                vlx = (cell.xCOM - self.centerMassX[-self.deltaTime])/self.deltaTime
                vly = (cell.yCOM - self.centerMassY[-self.deltaTime])/self.deltaTime
#                 speed = np.sqrt(vlx*vlx + vly*vly)
                
                if abs(vlx) > .5*self.dim.x/self.deltaTime:          
                    if vlx<0:
                        vlx+=self.dim.x/self.deltaTime
                    else:
                        vlx-=self.dim.x/self.deltaTime
                if abs(vly) > .5*self.dim.y/self.deltaTime:          
                    if vly<0:
                        vly+=self.dim.y/self.deltaTime
                    else:
                        vly-=self.dim.y/self.deltaTime
                
                cell.dict['velocityX'] = vlx
                cell.dict['velocityY'] = vly
#                 if abs(float(speed)) < self.maxReasonableSpeed: 
#                     cell.dict['velocityX'] = vlx
#                     cell.dict['velocityY'] = vly
                

                
                
            
        #boids force definition
        if mcs>self.deltaTime:
            for cell in self.cellList:
                cellVx = cell.dict['velocityX']
                cellVy = cell.dict['velocityY']
                neigVxList = np.array([])
                neigVyList = np.array([])
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        neigVxList = np.append(neigVxList,neighbor.dict['velocityX'])
                        neigVyList = np.append(neigVyList,neighbor.dict['velocityY'])
                if len(neigVxList):
                    neigVx = np.mean(neigVxList)
                else:       
                    neigVx = 0
                if len(neigVyList):
                    neigVy = np.mean(neigVyList)
                else:
                    neigVy = 0
                
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
        