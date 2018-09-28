
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
        
from edgeAndVertexDetectionSteppables import edgeAndVertexDetectionSteppable
steppableInstance=edgeAndVertexDetectionSteppable(sim,_frequency=1)
steppableRegistry.registerSteppable(steppableInstance)
        

from edgeAndVertexDetectionSteppables import extraFieldsManager
instanceOfextraFieldsManager=extraFieldsManager(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(instanceOfextraFieldsManager)


from edgeAndVertexDetectionSteppables import scatterFaceVertex
instanceOfscatterFaceVertex=scatterFaceVertex(_simulator=sim,_frequency=10)
steppableRegistry.registerSteppable(instanceOfscatterFaceVertex)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
        
        