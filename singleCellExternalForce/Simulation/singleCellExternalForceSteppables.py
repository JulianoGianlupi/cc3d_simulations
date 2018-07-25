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
    def start(self):
        # any code in the start function runs before MCS=0
        self.targetVolume = 100.
        self.lambdaVolume = 8.
        self.forceTheta = 0.0*np.pi
        self.forceModulus = 100
        self.pWCM = self.addNewPlotWindow(_title='Center of Mass x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Center of Mass', _xScaleType='linear', _yScaleType='linear')
        #self.pWCM.addAutoLegend("top")
        self.pWCM.addPlot('Center of Mass X', _style='Dots', _color='red', _size=5)
        self.pWCM.addPlot('Center of Mass Y', _style='Dots', _color='blue', _size=5)
        self.pWCM.addPlot('Center of Mass', _style='Dots', _color='white', _size=5)
        
        
        self.pWVelocity = self.addNewPlotWindow(_title='Velocity x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Velocity', _xScaleType='linear', _yScaleType='linear')
        #self.pWVelocity.addAutoLegend("top")
        self.pWVelocity.addPlot('Velocity X', _style='Dots', _color='red', _size=5)
        self.pWVelocity.addPlot('Velocity Y', _style='Dots', _color='blue', _size=5)
        self.pWVelocity.addPlot('Velocity', _style='Dots', _color='white', _size=5)
        
        
        self.pWAngle = self.addNewPlotWindow(_title='Angle (rad) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Angle (radians)', _xScaleType='linear', _yScaleType='linear')
        #self.pWAngle.addAutoLegend("top")
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
            
            cell.lambdaVecX = -self.forceModulus*np.cos(self.forceTheta) 
            cell.lambdaVecY = self.forceModulus*np.sin(self.forceTheta) 
            #cell.lambdaVecZ = 0.0  
        #opening files
        
        fileDir = os.path.dirname(os.path.abspath(__file__))
        
        centerOfMassFileName = fileDir+'/comData.txt'
        self.centerOfMassFile = open(centerOfMassFileName,'w')
        self.centerOfMassFile.write('Time, x COM, y COM, COM')
        
        velocityFileName = fileDir+'/velocityData.txt'
        self.velocityFile = open(velocityFileName,'w')
        self.velocityFile.write('Time, x Vel, y Vel, Vel')
        
        angleFileName = fileDir+'/angleData.txt'
        self.angleFile = open(angleFileName, 'w')
        self.angleFile.write('Time, Angle (rad)')
#
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
        for cell in self.cellList:
            if len(self.centerMassX)>10:
                dx = cell.xCOM - self.centerMassX[-10]
                dy = cell.yCOM - self.centerMassY[-10]
                cell.dict['angle'] = np.angle(complex(dx,dy))
                self.angleFile.write('%i,%f'%(mcs,cell.dict['angle']))
            
            cm = np.sqrt( cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)
            self.centerMassX.append(cell.xCOM)
            self.centerMassY.append(cell.yCOM)
            self.centerMass.append(cm)
            self.centerOfMassFile.write('%i,%f,%f,%f'%(mcs,cell.xCOM,cell.yCOM,cm))
            
            if mcs%5 == 0:
                self.pWCM.addDataPoint("Center of Mass X", mcs, cell.xCOM)
                self.pWCM.addDataPoint("Center of Mass Y", mcs, cell.yCOM)
                self.pWCM.addDataPoint("Center of Mass", mcs, cm)
                self.pWAngle.addDataPoint('Angle (rad)', mcs, cell.dict['angle'])
            
        
        if len(self.centerMass) > 5:
            vlx = (self.centerMassX[-1] - self.centerMassX[-5])/4
            self.velocityX.append(vlx)
            
            vly = (self.centerMassY[-1] - self.centerMassY[-5])/4
            self.velocityY.append(vly)
            
            vl = (self.centerMass[-1] - self.centerMass[-5])/4
            self.velocity.append(vl)
            
            self.velocityFile.write('%i,%f,%f,%f'%(mcs,vlx,vly,vly))
            if mcs%5 ==0:
                self.pWVelocity.addDataPoint("Velocity X", mcs, vlx)
                self.pWVelocity.addDataPoint("Velocity Y", mcs, vly)
                self.pWVelocity.addDataPoint("Velocity", mcs, vl)
        
    def finish(self):
        # Finish Function gets called after the last MCS
        self.centerOfMassFile.close()
        self.velocityFile.close()
        self.angleFile.close()