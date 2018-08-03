
from PySteppables import *
import CompuCell
import sys
class polarityBoidsSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCLField = self.createScalarFieldCellLevelPy("Angle")
        self.scalarVelocityField = self.createScalarFieldCellLevelPy("Speed")
        self.scalarForceField = self.createScalarFieldCellLevelPy("polarity_Angle")
        self.vectorVelocityField = self.createVectorFieldCellLevelPy("Velocity")
        self.vectorPolarityField = self.createVectorFieldCellLevelPy("Polarity")
        self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
        
        
    def start(self):
        # any code in the start function runs before MCS=0
        
        #constants
        self.targetVolume = 64.
        self.lambdaVolume = 8.
        self.forceModulus = -100
        self.forceCounter = 0
        self.deltaTime = 5
        #boids parameters
        self.alphaBoids = .5
        self.betaBoids = 1.
        self.BBoids = 5
        self.nBoids = 2
        
        
#         self.pWAngle = self.addNewPlotWindow(_title='Angle (rad) x Time', _xAxisTitle='MonteCarlo Step (MCS)',
#                                         _yAxisTitle='Angle (radians)', _xScaleType='linear', _yScaleType='linear')
        
#         self.pWAngle.addPlot('Angle (rad)', _style='Dots', _color='red', _size=5)
        
        #cell initiation and dictionaries creation
        for cell in self.cellList:
            cell.targetVolume = self.targetVolume
            cell.lambdaVolume = self.lambdaVolume
            
            
            cell.dict['forceAngle'] = np.random.uniform(0.,2.*np.pi)
            cell.dict['angle']=0.
            cell.dict['polarityAngle']=cell.dict['forceAngle']
            
            cell.dict['centerMassX'] = [cell.xCOM]
            cell.dict['centerMassY'] = [cell.yCOM]
            cell.dict['centerMass'] = [np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)]
            
            cell.dict['velocityX'] = 0.
            cell.dict['velocityY'] = 0.
            
            cell.dict['polarityX'] = 0.
            cell.dict['polarityY'] = 0.
            
            cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
            cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle']) 
            
            cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
            cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
            #cell.lambdaVecZ = 0.0  
            #
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        #fields plotting
        for cell in self.cellList:
            if cell:
                self.scalarCLField[cell] = cell.dict['angle']/np.pi
                self.scalarVelocityField[cell] = np.sqrt(cell.dict['velocityX']*cell.dict['velocityX'] 
                                            + cell.dict['velocityY']*cell.dict['velocityX'])
                self.scalarForceField[cell] = cell.dict['polarityAngle']/np.pi
                self.vectorVelocityField[cell] = [cell.dict['velocityX'],cell.dict['velocityY'],0]
                self.vectorPolarityField[cell] = [cell.dict['polarityX'],cell.dict['polarityY'],0]
                self.vectorForceField[cell] = [cell.lambdaVecX,cell.lambdaVecY,0]
        #

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
                if (self.centerMassX[-self.deltaTime] > .9*self.dim.x and
                                            cell.xCOM < .1*self.dim.x):
                    vlx = (cell.xCOM+self.dim.x - self.centerMassX[-self.deltaTime])/self.deltaTime
                    
                elif (self.centerMassX[-self.deltaTime] < .1*self.dim.x and
                                              cell.xCOM > .9*self.dim.x):
                    vlx = (cell.xCOM-self.dim.x - self.centerMassX[-self.deltaTime])/self.deltaTime
                    
                else:
                    vlx = (cell.xCOM - self.centerMassX[-self.deltaTime])/self.deltaTime
                
                
                if (self.centerMassY[-self.deltaTime] > .9*self.dim.y and
                                            cell.yCOM < .1*self.dim.y):
                    vly = (cell.yCOM+self.dim.y - self.centerMassY[-self.deltaTime])/self.deltaTime
                    
                elif (self.centerMassY[-self.deltaTime] < .1*self.dim.y and
                                              cell.yCOM > .9*self.dim.y):
                    vly = (cell.yCOM+self.dim.y - self.centerMassY[-self.deltaTime])/self.deltaTime
                    
                else:
                    vly = (cell.yCOM - self.centerMassY[-self.deltaTime])/self.deltaTime
                    
                
                cell.lambdaVecX = (self.forceModulus *  cell.dict['polarityX'] * 
                                  ( cell.dict['polarityX']**self.nBoids)/ ( (self.BBoids**self.nBoids) + ( cell.dict['polarityX']**self.nBoids)))
                cell.lambdaVecY = (self.forceModulus *  cell.dict['polarityY'] * 
                                  ( cell.dict['polarityY']**self.nBoids)/ ( (self.BBoids**self.nBoids) + ( cell.dict['polarityY']**self.nBoids)))
                
                cell.dict['polarityX'] = -self.alphaBoids * cell.dict['polarityX'] + self.betaBoids * cell.dict['velocityX']
                cell.dict['polarityY'] = -self.alphaBoids * cell.dict['polarityY'] + self.betaBoids * cell.dict['velocityY']
                
                cell.dict['polarityAngle'] = np.angle(complex(cell.dict['polarityX'],cell.dict['polarityY']))
                
                cell.dict['velocityX'] = vlx
                cell.dict['velocityY'] = vly
                
              
    def finish(self):
        # Finish Function gets called after the last MCS
        pass