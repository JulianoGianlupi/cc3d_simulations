
import sys
from os import environ
from os import getcwd
import string

sys.path.append(environ["PYTHON_MODULE_PATH"])


import CompuCellSetup


sim,simthread = CompuCellSetup.getCoreSimulationObjects()
        
# add extra attributes here
        
CompuCellSetup.initializeSimulationObjects(sim,simthread)
# Definitions of additional Python-managed fields go here
        
#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()
        
from singleCellExternalForceSteppables import singleCellExternalForceSteppable
steppableInstance=singleCellExternalForceSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)
        

from singleCellExternalForceSteppables import dataWriting
instanceOfdataWriting=dataWriting(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfdataWriting)


from singleCellExternalForceSteppables import dataPloter
instanceOfdataPloter=dataPloter(_simulator=sim,_frequency=50)
steppableRegistry.registerSteppable(instanceOfdataPloter)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        