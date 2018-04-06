
from PySteppables import *
import CompuCell
import sys
class cell_idSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # any code in the start function runs before MCS=0
        pass
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        for cell in self.cellList:
            cell.lambdaVolume = 2.0
            cell.targetVolume = 25
            cell.lambdaPressure = 2.0
            cell.Pressure = 25
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        