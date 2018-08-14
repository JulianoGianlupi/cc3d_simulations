from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import numpy as np
import os

#adding the python vars to be scaned

global G_targetVolume_G, G_lambda_Volume_G,G_forceAngle_G, G_forceModulus_G

G_forceAngle_G = 1
#the number given in the parameter scan is the angle in rads divided by pi multiplied by 12, 
#this is done in order to avoid a repeating decimal (might have to redo, but works for pi/6, pi/4, pi/2)
G_forceAngle_G *= np.pi/12. 

G_forceModulus_G = 1
#it's actually the energy, and f = -grad energy
G_forceModulus_G =-1.*float(G_forceModulus_G)

G_lambda_Volume_G = 1

G_targetVolume_G = 1



class singleCellExternalForceSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = G_targetVolume_G  #volume for cells
        self.lambdaVolume = G_lambda_Volume_G #how rigid cells are
        self.forceTheta = G_forceAngle_G #angle the force will be applied
        
        self.forceModulus = G_forceModulus_G 

        
        self.deltaTime = 5
        
        #turns around
        self.turnsAroundX = 0
        self.turnsAroundY = 0
        
        
        
        
        
        #cell initiation and dictionaries creation
        for cell in self.cellList:
            #asigning the volume parameters
            cell.targetVolume = self.targetVolume
            cell.lambdaVolume = self.lambdaVolume
            
            #the angle at t=0
            cell.dict['angle']=0.
            
            #writing the center of mass as a list to a cell
            #dictionary
            cell.dict['centerMassX'] = [cell.xCOM]
            cell.dict['centerMassY'] = [cell.yCOM]
            cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
            
            #writing the velocity as a list to a cell
            #dictionary
            cell.dict['velocityX'] = [0.]
            cell.dict['velocityY'] = [0.]
            cell.dict['velocity'] = [0.]
            
            #initiating the force to the cell
            cell.lambdaVecX = self.forceModulus[0]*np.cos(self.forceTheta) 
            cell.lambdaVecY = self.forceModulus[0]*np.sin(self.forceTheta) 
            #cell.lambdaVecZ = 0.0  
        
        
        
        

#
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
        for cell in self.cellList:
            #gets the previous center of mass to a variable
            self.centerMassX = cell.dict['centerMassX']
            self.centerMassY = cell.dict['centerMassY']
            
            #calculates the angle of movement
            if len(self.centerMassX)>self.deltaTime:
                dx = cell.xCOM - self.centerMassX[-self.deltaTime]
                dy = cell.yCOM - self.centerMassY[-self.deltaTime]
                cell.dict['angle'] = np.angle(complex(dx,dy))

            
            #writes the center of mass to the cell's list    
            cm = np.sqrt( (cell.xCOM)*(cell.xCOM) + (cell.yCOM)*(cell.yCOM))

            cell.dict['centerMassX'].append(cell.xCOM)            
            cell.dict['centerMassY'].append(cell.yCOM)
            cell.dict['centerMass'].append(cm)
            
    def finish(self):
       pass
from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
import CompuCellSetup


#class for writing the data from the cell
class dataWriting(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        print "dataWriting: This function is called once before simulation"
        self.fileDir = os.path.dirname(os.path.abspath(__file__))
        self.centerOfMassFile = fileDir+'/centerOfMass.dat'
        self.velocityFile = fileDir+'/velocity.dat'
        self.angleOfMovementFile = fileDir+'/angleMovement.dat'
        #headers:
        with open(self.centerOfMassFile,'a') as comDat:
            
            comDat.write('Time (mcs), xCOM, yCOM, |vec(xCOM)+vec(yCOM)|')
        with open(self.velocityFile,'a') as velDat:
            velDat.write('Time (mcs), Vx, Vy, |vec(Vx)+vec(Vy)|')
        with open(self.angleOfMovementFile,'a') as angleDat:
            angleDat.write('Time (mcs), angle (rad)')
    def step(self,mcs):
        
        for cell in self.cellList:
            with open(self.centerOfMassFile,'a') as comDat:
                #writes time, x position, y position, vec(x)+vec(y)
                comDat.write('\n%i,%f,%f,%f'%(mcs,cell.dict['centerMassX'],cell.dict['centerMassY'],cell.dict['centerMass']))
            with open(self.velocityFile,'a') as velDat:
                velDat.write('\n%i,%f,%f,%f'%(mcs, cell.dict['velocityX'],cell.dict['velocityY'],cell.dict['velocity']))
            with open(self.angleOfMovementFile,'a') as angleDat:
                angleDat.write('\n%i,%f'%(mcs,cell.dict['angle']))
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
    

from PySteppables import *
import CompuCell
import sys

from PlayerPython import *
import CompuCellSetup



class dataPloter(SteppableBasePy):
    def __init__(self,_simulator,_frequency=50):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        print "dataPloter: This function is called once before simulation"
        #graphs:
        ## center of mass graph: x in red, y in blue, x+y in white
        self.pWCM = self.addNewPlotWindow(_title='Center of Mass x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Center of Mass', _xScaleType='linear', _yScaleType='linear')
        
        self.pWCM.addPlot('Center of Mass X', _style='Dots', _color='red', _size=5)
        self.pWCM.addPlot('Center of Mass Y', _style='Dots', _color='blue', _size=5)
        self.pWCM.addPlot('Center of Mass', _style='Dots', _color='white', _size=5)
        
        ## velocity vs force graph: 
        ## x in red, y in blue, speed in white
        self.pWVelocity = self.addNewPlotWindow(_title='Velocity x Force', _xAxisTitle='Force Modulus',
                                        _yAxisTitle='Velocity', _xScaleType='linear', _yScaleType='linear')
        
        self.pWVelocity.addPlot('Velocity X', _style='Dots', _color='red', _size=5)
        self.pWVelocity.addPlot('Velocity Y', _style='Dots', _color='blue', _size=5)
        self.pWVelocity.addPlot('Velocity', _style='Dots', _color='white', _size=5)
        
        ## velocity vs time graph: 
        ## x in red, y in blue, speed in white
        self.pWVelocityTime = self.addNewPlotWindow(_title='Velocity x Time', _xAxisTitle='Time (mcs)',
                                        _yAxisTitle='Velocity', _xScaleType='linear', _yScaleType='linear')
        
        self.pWVelocityTime.addPlot('Velocity X', _style='Dots', _color='red', _size=5)
        self.pWVelocityTime.addPlot('Velocity Y', _style='Dots', _color='blue', _size=5)
        self.pWVelocityTime.addPlot('Velocity', _style='Dots', _color='white', _size=5)
        
        #angle in which cells moved graph
        self.pWAngle = self.addNewPlotWindow(_title='Angle (rad) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Angle (radians)', _xScaleType='linear', _yScaleType='linear')
        
        self.pWAngle.addPlot('Angle (rad)', _style='Dots', _color='red', _size=5)
    def step(self,mcs):
        print "dataPloter: This function is called every 50 MCS"
        for cell in self.cellList:
            self.pWCM.addDataPoint("Center of Mass X", mcs, cell.xCOM)
            self.pWCM.addDataPoint("Center of Mass Y", mcs, cell.yCOM)
            self.pWCM.addDataPoint("Center of Mass", mcs, cm)
            self.pWAngle.addDataPoint('Angle (rad)', mcs, cell.dict['angle'])
            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
    
