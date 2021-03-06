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
        self.targetVolume = 64. #volume for cells
        self.lambdaVolume = 8.  #how rigid cells are
        self.forceTheta = 0.*np.pi #angle the force will be applied
        #list of forces to sweep
        self.forceModulus = [-1.,-5.,-10.,-20.,-30.,-40.,-50.,-75.,-100.,-500.,-1000.,-5000.,-10000.,-50000.]

        self.forceCounter = 0
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

            #changing parameters every 100 mcs
            if (mcs%100 == 0) and (mcs>0):
                
                t=100
                #calculates the avg velocity over 1000 mcs

#                 if self.centerMassX[-t] > cell.xCOM+.2*self.dim.x:
#                     vlx = (cell.xCOM+self.dim.x - self.centerMassX[-t])/t
#                     cell.dict['velocityX'].append(vlx)
#                 else:
#                     vlx = (cell.xCOM - self.centerMassX[-t])/t
#                     cell.dict['velocityX'].append(vlx)
#                 if self.centerMassY[-t] > cell.yCOM+.2*self.dim.y:
# #                     vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
# #                     cell.dict['velocityY'].append(vly)
#                 else:
#                     vly = (cell.yCOM - self.centerMassY[-t])/t
#                     cell.dict['velocityY'].append(vly)
                    
#                 if self.centerMassX[-t] < cell.xCOM-.2*self.dim.x:
#                     vlx = (cell.xCOM+self.dim.x - self.centerMassX[-t])/t
#                     cell.dict['velocityX'].append(vlx)
#                 else:
#                     vlx = (cell.xCOM - self.centerMassX[-t])/t
#                     cell.dict['velocityX'].append(vlx)
#                 if self.centerMassY[-t] < cell.yCOM-.2*self.dim.y:
#                     vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
#                     cell.dict['velocityY'].append(vly)
#                 else:
#                     vly = (cell.yCOM - self.centerMassY[-t])/t
#                     cell.dict['velocityY'].append(vly)
                vlx = (cell.xCOM - self.centerMassX[-t])/t
                cell.dict['velocityX'].append(vlx)
                
                vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
                cell.dict['velocityY'].append(vly)
                

                #and writes it to a dictionary
                #the bunch of if statments are here to try
                #to account for cells going trough the border
                if self.centerMassX[-t] > cell.xCOM+.2*self.dim.x:
                    vlx = (cell.xCOM+self.dim.x - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                elif self.centerMassX[-t] < cell.xCOM-.2*self.dim.x:
                    vlx = (cell.xCOM+self.dim.x - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                else:
                    vlx = (cell.xCOM - self.centerMassX[-t])/t
                    cell.dict['velocityX'].append(vlx)
                if self.centerMassY[-t] > cell.yCOM+.2*self.dim.y:
                    vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                elif self.centerMassY[-t] < cell.yCOM-.2*self.dim.y:
                    vly = (cell.yCOM+self.dim.y - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                else:
                    vly = (cell.yCOM - self.centerMassY[-t])/t
                    cell.dict['velocityY'].append(vly)
                #speed

                vl = np.sqrt(vlx*vlx + vly*vly)
                
                cell.dict['velocity'].append(vl)
                
                #plots the velocities vs force
                self.pWVelocity.addDataPoint("Velocity X", -self.forceModulus[self.forceCounter], vlx)
                self.pWVelocity.addDataPoint("Velocity Y", -self.forceModulus[self.forceCounter], vly)
                self.pWVelocity.addDataPoint("Velocity", -self.forceModulus[self.forceCounter], vl)
                
                #changes the force modulus
                self.forceCounter+=1
                if self.forceCounter < len(self.forceModulus):
                    cell.lambdaVecX = self.forceModulus[self.forceCounter]*np.cos(self.forceTheta) 
                    cell.lambdaVecY = self.forceModulus[self.forceCounter]*np.sin(self.forceTheta) 
                else:#when all forces have been used stops sim
                    self.stopSimulation()
                '''if self.forceCounter > len(self.forceModulus):
                    self.stopSimulation()
                #
                try:
                    cell.lambdaVecX = self.forceModulus[self.forceCounter]*np.cos(self.forceTheta) 
                    cell.lambdaVecY = self.forceModulus[self.forceCounter]*np.sin(self.forceTheta)    
                except:
                    self.stopSimulation()
                    break

            moveX = 0
            moveY =0
            if mcs>5:
                deltaX = cell.xCOM-self.centerMassX[-5]
                deltaY = cell.yCOM-self.centerMassY[-5]
                
                if abs(deltaX) > .5*self.dim.x:
                    if deltaX < 0:
                        moveX =- self.dim.x
                    else:
                        moveX = self.dim.x
                if abs(deltaY) > .5*self.dim.y:
                    if deltaY < 0:
                        moveY = -self.dim.y
                    else:
                        moveY = self.dim.y   
                print moveX,moveY
            cm = np.sqrt( cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)

            '''
            #writes the center of mass to the cell's list    
            cm = np.sqrt( (cell.xCOM)*(cell.xCOM) + (cell.yCOM)*(cell.yCOM))

            cell.dict['centerMassX'].append(cell.xCOM)            
            cell.dict['centerMassY'].append(cell.yCOM)
            cell.dict['centerMass'].append(cm)
            
            #every 50 mcs (after 5) calculates a more
            #'instantaneous' velocity and plots it
            #also every 50 mcs it plots the angle
            if mcs%50 == 0:
                if mcs>5:
                    vlx = (cell.xCOM - self.centerMassX[-5])/5
                    vly = (cell.yCOM - self.centerMassY[-5])/5
                    vl = np.sqrt( vlx*vlx + vly * vlx)
                    self.pWVelocityTime.addDataPoint("Velocity X", mcs, vlx)
                    self.pWVelocityTime.addDataPoint("Velocity Y", mcs, vly)
                    self.pWVelocityTime.addDataPoint("Velocity", mcs, vl)
                
                
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
    
