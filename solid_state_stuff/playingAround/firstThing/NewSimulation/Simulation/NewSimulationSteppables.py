
from PySteppables import *
import CompuCell
import sys
import numpy as np
class NewSimulationSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        self.mainLength = 10
        self.surfaceVolumeRatio = 2/self.mainLength 
        self.targetVolumeSolid1 = np.pi*self.mainLength*self.mainLength
        for cell in self.cellList:
            cell.targetVolume = self.targetVolumeSolid1
            cell.lambdaVolume = 0
            cell.targetSurface = self.surfaceVolumeRatio*self.targetVolumeSolid1
            cell.lambdaSurface = 0
            cell.Pressure = 5000
            
    def step(self,mcs):        
       pass
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        