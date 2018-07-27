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
        
        #constants
        self.targetVolume = 100.
        self.lambdaVolume = 8.
        self.forceTheta = 0.5*np.pi
        self.forceModulus = [-1.,-5.,-10.,-20.,-30.,-40.,-50.,-75.,-100.,-500.,-1000.,-5000.,-10000.,-50000.]
        self.forceCounter = 0
        self.deltaTime = 5
        
        
        #graphs:
        self.pWCM = self.addNewPlotWindow(_title='Center of Mass x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Center of Mass', _xScaleType='linear', _yScaleType='linear')
        #self.pWCM.addAutoLegend("top")
        self.pWCM.addPlot('Center of Mass X', _style='Dots', _color='red', _size=5)
        self.pWCM.addPlot('Center of Mass Y', _style='Dots', _color='blue', _size=5)
        self.pWCM.addPlot('Center of Mass', _style='Dots', _color='white', _size=5)
        
        
        self.pWVelocity = self.addNewPlotWindow(_title='Velocity x Force', _xAxisTitle='Force Modulus',
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
        
        self.xTurns = 0
        self.yTurns = 0
        
        #cell initiation and dictionaries creation
        for cell in self.cellList:
            cell.targetVolume = self.targetVolume
            cell.lambdaVolume = self.lambdaVolume
            
            cell.dict['angle']=0.
            
            '''
            cell.dict['centerMassX'] = [cell.xCOM]
            cell.dict['centerMassY'] = [cell.yCOM]
            cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
            '''
            cell.lambdaVecX = self.forceModulus[0]*np.cos(self.forceTheta) 
            cell.lambdaVecY = self.forceModulus[0]*np.sin(self.forceTheta) 
            #cell.lambdaVecZ = 0.0  
        
        
        
        
        #opening files
        
        fileDir = os.path.dirname(os.path.abspath(__file__))
        
        centerOfMassFileName = fileDir+'/'+'cellVolume'+str(self.targetVolume)+'forceAngle'+str(self.forceTheta)+'ForceMag'+str(self.forceModulus[0])+'comData.txt'
        self.centerOfMassFile = open(centerOfMassFileName,'w')
        self.centerOfMassFile.write('Time, x COM, y COM, COM\n')
        
        velocityFileName = fileDir+'/'+'cellVolume'+str(self.targetVolume)+'forceAngle'+str(self.forceTheta)+'ForceMag'+str(self.forceModulus[0])+'velocityData.txt'
        self.velocityFile = open(velocityFileName,'w')
        self.velocityFile.write('Time, x Vel, y Vel, Vel\n')
        
        angleFileName = fileDir+'/'+'cellVolume'+str(self.targetVolume)+'forceAngle'+str(self.forceTheta)+'ForceMag'+str(self.forceModulus[0])+'AngleData.txt'
        self.angleFile = open(angleFileName, 'w')
        self.angleFile.write('Time, Angle (rad)\n')
#
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
        
        
        
        
        
        #actual program:
        #velocity and angle calculation
        for cell in self.cellList:
            if len(self.centerMassX)>self.deltaTime:
                dx = cell.xCOM - self.centerMassX[-self.deltaTime]
                dy = cell.yCOM - self.centerMassY[-self.deltaTime]
                cell.dict['angle'] = np.angle(complex(dx,dy))
                self.angleFile.write('%i,%f\n'%(mcs,cell.dict['angle']))

            #changing parameters
            if (mcs%1000 == 0) and (mcs>0):
                self.forceCounter+=1
                
                #
                
                vlx = (cell.xCOM+self.xTurns - self.centerMassX[-1000])/1000
                self.velocityX.append(vlx)
                
                vly = (cell.yCOM+self.yTurns - self.centerMassY[-1000])/1000
                self.velocityY.append(vly)
                
                vl = np.sqrt(vlx*vlx + vly*vly)
                self.velocity.append(vl)
                
                self.velocityFile.write('%i,%f,%f,%f\n'%(mcs,vlx,vly,vly))
            
                self.pWVelocity.addDataPoint("Velocity X", -self.forceModulus[self.forceCounter-1], vlx)
                self.pWVelocity.addDataPoint("Velocity Y", -self.forceModulus[self.forceCounter-1], vly)
                self.pWVelocity.addDataPoint("Velocity", -self.forceModulus[self.forceCounter-1], vl)
            
                if self.forceCounter > len(self.forceModulus):
                    self.stopSimulation()
                #
                try:
                    cell.lambdaVecX = self.forceModulus[self.forceCounter]*np.cos(self.forceTheta) 
                    cell.lambdaVecY = self.forceModulus[self.forceCounter]*np.sin(self.forceTheta)    
                except:
                    self.stopSimulation()
                    break
            cm = np.sqrt( (cell.xCOM+self.yTurns)*(cell.xCOM+self.yTurns) + (cell.yCOM+self.yTurns)*(cell.yCOM+self.yTurns))
            self.centerMassX.append(cell.xCOM+self.xTurns)
            self.centerMassY.append(cell.yCOM+self.yTurns)
            self.centerMass.append(cm)
            self.centerOfMassFile.write('%i,%f,%f,%f\n'%(mcs,cell.xCOM,cell.yCOM,cm))
            
            if cell.xCOM < (self.centerMassX[-1] - .5*self.dim.x):
                self.yTurns += self.dim.y
                
            elif cell.xCOM > (self.centerMassX[-1] + .5*self.dim.x):
                self.yTurns -= self.dim.y
                
            if cell.yCOM < (self.centerMassY[-1] - .5*self.dim.y):
                self.yTurns += self.dim.y
                
            elif cell.yCOM > (self.centerMassY[-1] + .5*self.dim.y):
                self.yTurns -= self.dim.y
            if mcs%100 == 0:
                self.pWCM.addDataPoint("Center of Mass X", mcs, cell.xCOM)
                self.pWCM.addDataPoint("Center of Mass Y", mcs, cell.yCOM)
                self.pWCM.addDataPoint("Center of Mass", mcs, cm)
                self.pWAngle.addDataPoint('Angle (rad)', mcs, cell.dict['angle'])
        #print self.xTurns
        #print self.yTurns
        # 
        '''
        if len(self.centerMass) > self.deltaTime:
            vlx = (self.centerMassX[-1] - self.centerMassX[-self.deltaTime])/self.deltaTime
            self.velocityX.append(vlx)
            
            vly = (self.centerMassY[-1] - self.centerMassY[-self.deltaTime])/self.deltaTime
            self.velocityY.append(vly)
            
            vl = (self.centerMass[-1] - self.centerMass[-self.deltaTime])/self.deltaTime
            self.velocity.append(vl)
            
            self.velocityFile.write('%i,%f,%f,%f\n'%(mcs,vlx,vly,vly))
            if mcs%100 ==0:
                self.pWVelocity.addDataPoint("Velocity X", mcs, vlx)
                self.pWVelocity.addDataPoint("Velocity Y", mcs, vly)
                self.pWVelocity.addDataPoint("Velocity", mcs, vl)
        '''
    def finish(self):
        # Finish Function gets called after the last MCS
        self.centerOfMassFile.close()
        self.velocityFile.close()
        self.angleFile.close()