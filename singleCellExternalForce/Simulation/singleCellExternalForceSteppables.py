from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import numpy as np
import os

class singleCellExternalForceSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField = self.createScalarFieldCellLevelPy("Angle")
        
        
        
        
        
        
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = 100.
        self.lambdaVolume = 8.
        self.forceTheta = 1.*np.pi
        self.forceModulus = [-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,-100,]#[-1.,-5.,-10.,-20.,-30.,-40.,-50.,-75.,-100.,-500.,-1000.,-5000.,-10000.,-50000.]
        self.forceCounter = 0
        self.deltaTime = 5
        
        
        
        
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
            
            cell.dict['velocityX'] = [0.]
            cell.dict['velocityY'] = [0.]
            cell.dict['velocity'] = [0.]
            
            cell.lambdaVecX = self.forceModulus[0]*np.cos(cell.dict['forceAngle']) 
            cell.lambdaVecY = self.forceModulus[0]*np.sin(cell.dict['forceAngle']) 
            #cell.lambdaVecZ = 0.0  
            #
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        for cell in self.cellList:
            self.scalarCLField[cell] = cell.dict['angle']
        
        
        
        if (mcs%100 == 0) and (mcs>0):
            self.forceCounter+=1
            t=100
            if self.forceCounter > len(self.forceModulus):
                self.stopSimulation()
            for cell in self.cellList:
                
                #calculates the avg velocity over 1000 mcs
                if self.centerMassX[-t] > cell.xCOM+.2*self.dim.x:
                    vlx = (cell.xCOM+self.dim.x - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                else:
                    vlx = (cell.xCOM - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                if self.centerMassY[-t] > cell.yCOM+.2*self.dim.y:
                    vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                else:
                    vly = (cell.yCOM - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                    
                if self.centerMassX[-t] < cell.xCOM-.2*self.dim.x:
                    vlx = (cell.xCOM+self.dim.x - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                else:
                    vlx = (cell.xCOM - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                if self.centerMassY[-t] < cell.yCOM-.2*self.dim.y:
                    vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                else:
                    vly = (cell.yCOM - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                vl = np.sqrt(vlx*vlx + vly*vly)
                cell.dict['velocity'].append(vl)
                
                #changes the force modulus
                
                #
                try:
                    cell.lambdaVecX = self.forceModulus[self.forceCounter]*np.cos(cell.dict['forceAngle']) 
                    cell.lambdaVecY = self.forceModulus[self.forceCounter]*np.sin(cell.dict['forceAngle'])    
                except:
                    self.stopSimulation()
                    break
                    
        #
        for cell in self.cellList:
            self.centerMassX = cell.dict['centerMassX']
            self.centerMassY = cell.dict['centerMassY']
            
            if len(self.centerMassX)>self.deltaTime:
                dx = cell.xCOM - self.centerMassX[-self.deltaTime]
                dy = cell.yCOM - self.centerMassY[-self.deltaTime]
                cell.dict['angle'] = np.angle(complex(dx,dy))
                
                
            cm = np.sqrt( (cell.xCOM)*(cell.xCOM) + (cell.yCOM)*(cell.yCOM))
            cell.dict['centerMassX'].append(cell.xCOM)            
            cell.dict['centerMassY'].append(cell.yCOM)
            cell.dict['centerMass'].append(cm)
            
            
            if mcs%50 == 0:
                if mcs>5:
                    vlx = (cell.xCOM - self.centerMassX[-5])/5
                    vly = (cell.yCOM - self.centerMassY[-5])/5
                    vl = np.sqrt( vlx*vlx + vly * vlx)
                
    def finish(self):
       pass
