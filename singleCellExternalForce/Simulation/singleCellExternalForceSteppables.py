from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import numpy as np
class singleCellExternalForceSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        self.targetVolume = 100.
        self.lambdaVolume = 8.
        self.forceTheta = 0.
        self.forceModulus = 5
        self.pWCM = self.addNewPlotWindow(_title='Center of Mass x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Variables', _xScaleType='linear', _yScaleType='linear')
        self.pWCM.addPlot('Center of Mass X', _style='Dots', _color='red', _size=5)
        self.pWCM.addPlot('Center of Mass Y', _style='Dots', _color='blue', _size=5)
        self.pWCM.addPlot('Center of Mass', _style='Dots', _color='white', _size=5)
        
        
        self.pWVelocity = self.addNewPlotWindow(_title='Velocity x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Velocity', _xScaleType='linear', _yScaleType='linear')
        self.pWVelocity.addPlot('Velocity X', _style='Dots', _color='red', _size=5)
        self.pWVelocity.addPlot('Velocity Y', _style='Dots', _color='blue', _size=5)
        self.pWVelocity.addPlot('Velocity ', _style='Dots', _color='white', _size=5)
        
        
        self.pWAngle = self.addNewPlotWindow(_title='Angle (rad) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Variables', _xScaleType='linear', _yScaleType='linear')
        self.pWAngle.addPlot('Angle (rad)', _style='Dots', _color='red', _size=5)
        
        
        
        
        
        self.centerMassX = []
        self.centerMassY = []
        self.centerMass = []
        
        self.velocityX = []
        self.velocityY = []
        self.velocity = []
        
        for cell in self.cellList:
            cell.targetVolume = self.targetVolume
            cell.lambdaVolume = self.lambdaVolume
            
            cell.dict['angle']=0.
            
            cell.lambdaVecX = self.forceModulus*np.cos(self.forceTheta) 
            cell.lambdaVecY = self.forceModulus*np.sin(self.forceTheta) 
            #cell.lambdaVecZ = 0.0  
#
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
        for cell in self.cellList:
            if len(self.centerMassX)>5:
                dx = cell.xCOM - self.centerMassX[-5]
                dy = cell.yCOM - self.centerMassY[-5]
                cell.dict['angle'] = np.angle(complex(dx,dy))
                self.pWAngle.addDataPoint('Angle (rad)', mcs, cell.dict['angle'])
            cm = np.sqrt( cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)
            self.centerMassX.append(cell.xCOM)
            self.centerMassY.append(cell.yCOM)
            self.centerMass.append(cm)
            self.pWCM.addDataPoint("Center of Mass X", mcs, cell.xCOM)
            self.pWCM.addDataPoint("Center of Mass Y", mcs, cell.yCOM)
            self.pWCM.addDataPoint("Center of Mass", mcs, cm)
            
        
        if len(self.centerMass) > 5:
            vlx = (self.centerMassX[-1] - self.centerMassX[-5])/4
            self.velocityX.append(vlx)
            
            vly = (self.centerMassY[-1] - self.centerMassY[-5])/4
            self.velocityY.append(vly)
            
            vl = (self.centerMass[-1] - self.centerMass[-5])/4
            self.velocity.append(vl)
            
            self.pWVelocity.addDataPoint("Velocity X", mcs, vlx)
            self.pWVelocity.addDataPoint("Velocity Y", mcs, vly)
            self.pWVelocity.addDataPoint("Velocity", mcs, vl)
        if len(self.velocity)>5:
            avgV = (vl - self.velocity[-5]) / 5
            self.pWVelocity.addDataPoint("Average Velocity", mcs, avgV)
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        