
from PySteppables import *
from numpy import *
from math import *
import random
import CompuCell
import sys
class BoidsSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # Add attributes for boids-propulsive force (could replace with Hill function)
        # Angle of force
        # Fluctuation rate of force angle
        self.pW = self.addNewPlotWindow(_title='Position', _xAxisTitle='X Position',
                                        _yAxisTitle='Y Position', _xScaleType='linear', _yScaleType='linear')
        self.pW.addPlot('ParticlePositions', _style='Dots', _color='red', _size=2)
        self.pW2 = self.addNewPlotWindow(_title='Force', _xAxisTitle='MCS',
                                        _yAxisTitle='Force', _xScaleType='linear', _yScaleType='linear')
        self.pW2.addPlot('XForce', _style='Dots', _color='blue', _size=5)
        self.pW2.addPlot('Xvelocity', _style='Dots', _color='green', _size=5)
        
        
        
        for cell in self.cellList:
            cell.dict['propulsiveforce'] = 10
            cell.dict['angle'] = 2.0*3.14159*random.random()
            cell.dict['alpha']=0.99 #decay rate of force
            cell.dict['beta']=1.0 # conversion of CM velocity to force
            cell.dict['anglefluctuation'] = 0.1
            cell.dict['PreviousXposition'] = cell.xCOM
            cell.dict['PreviousYposition'] = cell.yCOM
            cell.lambdaVecX = -cell.dict['propulsiveforce']*cos(cell.dict['angle']) 
            cell.lambdaVecY = cell.dict['propulsiveforce']*sin(cell.dict['angle'])   
            cell.lambdaVecZ = 0.0  
            
            
    def step(self,mcs):        
        #Plot CM positions of cells
        if not(mcs%100):
            for cell in self.cellList:
                xvelocity=cell.xCOM-cell.dict['PreviousXposition']
                yvelocity=cell.yCOM-cell.dict['PreviousYposition']
                cell.lambdaVecX=cell.dict['alpha']*cell.lambdaVecX-cell.dict['beta']*(xvelocity/(3.0+abs(xvelocity)))
                cell.lambdaVecY=cell.dict['alpha']*cell.lambdaVecY-cell.dict['beta']*(yvelocity/(3.0+abs(yvelocity)))
                cell.dict['PreviousXposition'] = cell.xCOM
                cell.dict['PreviousYposition'] = cell.yCOM
                self.pW.addDataPoint('ParticlePositions', cell.xCOM, cell.yCOM) 
                self.pW2.addDataPoint('XForce', mcs, cell.lambdaVecX) 
                self.pW2.addDataPoint('Xvelocity', mcs, xvelocity) 

            
 
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        