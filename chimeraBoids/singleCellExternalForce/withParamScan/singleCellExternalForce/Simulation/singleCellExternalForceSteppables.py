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
            
            #writing the center of mass to a cell
            #dictionary
            cell.dict['centerMassX'] = cell.xCOM
            cell.dict['centerMassY'] = cell.yCOM
#             cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
            #writing the position of the center of mass adjusted
            #for going through walls
            cell.dict['positionX'] = [cell.xCOM]
            cell.dict['positionY'] = [cell.yCOM]
            
            
            #writing the velocity  to a cell
            #dictionary
            cell.dict['velocityX'] = 0.
            cell.dict['velocityY'] = 0.
            cell.dict['velocity'] = 0.
            
            #initiating the force to the cell
            cell.lambdaVecX = self.forceModulus*np.cos(self.forceTheta) 
            cell.lambdaVecY = self.forceModulus*np.sin(self.forceTheta) 
            #cell.lambdaVecZ = 0.0  
        
        
        
        

#
#     def teleportFix(self,deltaX,deltaY):
        
#         if abs(deltaX)>self.dim.x:
#             if deltaX < 0:
#                 deltaX -=self.dim.x
#             else:
#                 deltaX += self.dim.x
#         if abs(deltaY)>self.dim.y:
#             if deltaY<0:
#                 deltaY -=self.dim.y
#             else:
#                 deltaY +=self.dim.y
#         return deltaX, deltaY
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
        for cell in self.cellList:
            #gets the previous center of mass to a variable
            self.positionX = cell.dict['positionX']
            self.positionY = cell.dict['positionY']
            
            dx = cell.xCOM - cell.dict['centerMassX']
            dy = cell.yCOM - cell.dict['centerMassY']
#             print dx, dy
#             dx, dy = teleportFix(self,dx,dy)
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
#             print dx, dy
            cell.dict['positionX'].append(self.positionX[-1]+dx)
            cell.dict['positionY'].append(self.positionY[-1]+dy)
            cell.dict['velocityX'] = dx
            cell.dict['velocityY'] = dy
            #calculates the angle of movement
            cell.dict['angle'] = np.angle(complex(dx,dy))/np.pi
#             if len(self.positionX)>self.deltaTime:
#                 angle_dx = self.positionX[-1] - self.positionX[-self.deltaTime]
#                 angle_dy = self.positionY[-1] - self.positionX[-self.deltaTime]
#                 cell.dict['angle'] = np.angle(complex(angle_dx,angle_dy))/np.pi
                

            
            

            cell.dict['centerMassX'] = cell.xCOM
            cell.dict['centerMassY'] = cell.yCOM
#             cell.dict['centerMass'].append(cm)
            
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
        fileDir = os.path.dirname(os.path.abspath(__file__))
        
        self.targetVolume = G_targetVolume_G  #volume for cells
        self.lambdaVolume = G_lambda_Volume_G #how rigid cells are
        self.forceTheta = G_forceAngle_G #angle the force will be applied
        self.forceModulus = G_forceModulus_G 
        #constant file
        self.constantsFile = fileDir+'/constants.txt'
        #data files
        self.centerOfMassFile = fileDir+'/centerOfMass.dat'
        self.velocityFile = fileDir+'/velocity.dat'
        self.angleOfMovementFile = fileDir+'/angleMovement.dat'
        self.volumeFile = fileDir+'/volume.dat'
        self.surfaceFile = fileDir +'/surface.dat'
        
        
        with open(self.constantsFile,'a') as cstFile:
            cstFile.write('target volume = ' + str(self.targetVolume) +
                          '\n lambda volume = ' + str(self.lambdaVolume) +
                          '\n force angle = ' + str(self.forceTheta) +
                          '\n force modulus = ' + str(self.forceModulus) )
        
        #headers:
        with open(self.centerOfMassFile,'a') as comDat:
            comDat.write('Time (mcs), xCOM, yCOM, |vec(xCOM)+vec(yCOM)|')
        with open(self.velocityFile,'a') as velDat:
            velDat.write('Time (mcs), Vx, Vy, |vec(Vx)+vec(Vy)|')
        with open(self.angleOfMovementFile,'a') as angleDat:
            angleDat.write('Time (mcs), angle (rad)')
        with open(self.volumeFile,'a') as volDat:
            volDat.write('Time(mcs), Cell Volume (voxels)')
        with open(self.surfaceFile,'a') as surfDat:
            surfDat.write('Time(mcs), Cell Surface (voxels)')
        #
    def step(self,mcs):
        
        for cell in self.cellList:
            self.positionX = cell.dict['positionX']
            self.positionY = cell.dict['positionY']
            with open(self.centerOfMassFile,'a') as comDat:
#                 writes time, x position, y position
                comDat.write('\n %i,%f,%f'%(mcs,self.positionX[mcs-1],self.positionY[mcs-1]))
                
            with open(self.velocityFile,'a') as velDat:
                velDat.write('\n %i,%f,%f'%(mcs, cell.dict['velocityX'],cell.dict['velocityY']))
                
            with open(self.angleOfMovementFile,'a') as angleDat:
                angleDat.write('\n %i,%f'%(mcs,cell.dict['angle']))
                
            with open(self.volumeFile,'a') as volDat:
                volDat.write('\n %i,%f'%(mcs,cell.volume))
                
            with open(self.surfaceFile,'a') as surfDat:
                surfDat.write('\n %i, %f'%(mcs,cell.surface))
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
        
        
        ## velocity vs time graph: 
        ## x in red, y in blue, speed in white
        self.pWVelocityTime = self.addNewPlotWindow(_title='Velocity x Time', _xAxisTitle='Time (mcs)',
                                        _yAxisTitle='Velocity', _xScaleType='linear', _yScaleType='linear')
        
        self.pWVelocityTime.addPlot('Velocity X', _style='Dots', _color='red', _size=5)
        self.pWVelocityTime.addPlot('Velocity Y', _style='Dots', _color='blue', _size=5)
        self.pWVelocityTime.addPlot('Velocity', _style='Dots', _color='white', _size=5)
        
        #angle in which cells moved graph
        self.pWAngle = self.addNewPlotWindow(_title='Angle (rad/pi) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
                                        _yAxisTitle='Angle (radians/pi)', _xScaleType='linear', _yScaleType='linear')
        
        self.pWAngle.addPlot('Angle (rad)', _style='Dots', _color='red', _size=5)
    def step(self,mcs):
        print "dataPloter: This function is called every 50 MCS"
        for cell in self.cellList:
            self.positionX = cell.dict['positionX']
            self.positionY = cell.dict['positionY']
            self.pWCM.addDataPoint("Center of Mass X", mcs, self.positionX[mcs])
            self.pWCM.addDataPoint("Center of Mass Y", mcs, self.positionY[mcs])
            
            self.pWAngle.addDataPoint('Angle (rad)', mcs, cell.dict['angle'])
            
            self.pWVelocityTime.addDataPoint('Velocity X',mcs,cell.dict['velocityX'])
            self.pWVelocityTime.addDataPoint('Velocity Y',mcs,cell.dict['velocityY'])
            
    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
    
