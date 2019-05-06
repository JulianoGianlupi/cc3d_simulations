from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys

import numpy as np
import os

global G_repetitionNumber_G

global G_targetVolume_G
global G_lambdaVolume_G
global G_forceModulus_G
global G_deltaTime_G

global G_alphaBoids_G
global G_betaBoids_G
global G_gammaBoids_G
global G_noise_G
global G_density_G

#repeat
G_repetitionNumber_G = -1 
#it's supposed to give the copy number, so if it's bnegative it's wrong


#parameter scan variables
G_targetVolume_G = 64.
G_lambdaVolume_G = 8.
G_forceModulus_G = -20 #speeds are tipically 1/100 of force.
G_deltaTime_G = 10
G_density_G = 0.20

G_alphaBoids_G = .1#5.5
G_betaBoids_G = 1.5
G_gammaBoids_G = .1
G_noise_G = .1





class chimeraBoidsV2Steppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        #initiating extra fields
        
#         self.scalarAngleField = self.createScalarFieldCellLevelPy('angle')
        
#         self.vectorIVelocityField = self.createVectorFieldCellLevelPy('Instant_Velocity')
#         self.vectorTVelocityField = self.createVectorFieldCellLevelPy('dt_Velocity')
        
#         self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
        
#         self.scalarClusterIDfield = self.createScalarFieldCellLevelPy("Cluster_ID")
        self.init_fields()
        pass
        
        
    def init_fields(self):
        self.scalarAngleField = self.createScalarFieldCellLevelPy('angle')
        
        self.vectorIVelocityField = self.createVectorFieldCellLevelPy('Instant_Velocity')
        self.vectorTVelocityField = self.createVectorFieldCellLevelPy('dt_Velocity')
        
        self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
        
        self.scalarClusterIDfield = self.createScalarFieldCellLevelPy("Cluster_ID")
    def start(self):
        # any code in the start function runs before MCS=0
        
        
        
        #constants
        self.targetVolume = G_targetVolume_G
        self.lambdaVolume = G_lambdaVolume_G
        self.forceModulus = G_forceModulus_G
        self.deltaTime = G_deltaTime_G
        
        self.density = G_density_G
        
        #boids parameters
        self.alphaBoids = G_alphaBoids_G
        self.betaBoids = G_betaBoids_G
        self.gammaBoids = G_gammaBoids_G
        self.noiseBoids = G_noise_G
        #create data folder and write initial things
        self.createWritingLoc()
        
        #set of used cluster IDs
        self.usedClusterIDs = set([0])#I think that to avoid possible errors it's best if the set is not
#                                   #empty
        self.currentIDs = set([0])
        
        
        #seeding the space
#         density = numberOfCells*self.targetVolume/(self.dim.x*self.dim.y*self.dim.z)
        
        numberOfCells = self.density * (self.dim.x*self.dim.y*self.dim.z)/self.targetVolume
        self.seedTheSpace(numberOfCells)
        
        
        
        #assigning parameters
        for cell in self.cellList:
            if cell:
                cell.targetVolume = self.targetVolume
                cell.lambdaVolume = self.lambdaVolume
                
                
                cell.dict['forceAngle'] = np.random.uniform(-np.pi,np.pi)
#                 print cell.dict['forceAngle']
                cell.dict['angle']=cell.dict['forceAngle']
                
                cell.dict['centerMassX'] = cell.xCOM
                cell.dict['centerMassY'] = cell.yCOM
                cell.dict['centerMass'] = np.sqrt(cell.xCOM*cell.xCOM + cell.yCOM*cell.yCOM)
                
                cell.dict['positionX'] = [cell.xCOM]
                cell.dict['positionY'] = [cell.yCOM]
                
                cell.dict['velocityX_instant'] = 0.
                cell.dict['velocityY_instant'] = 0.
                
                cell.dict['velocityX_deltaT'] = 0.
                cell.dict['velocityY_deltaT'] = 0.
                
                cell.dict['norm_velocityX_deltaT'] = 0.
                cell.dict['norm_velocityY_deltaT'] = 0.
                
                cell.dict['mean_neig_velX'] = 0
                cell.dict['mean_neig_velY'] = 0
                
                cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                cell.dict['orientation'] = np.zeros(2)
                #if cell.yCOM > .5*self.dim.y:
                #    cell.lambdaVecX = self.forceModulus*np.cos(.5*np.pi) 
                #    cell.lambdaVecY = self.forceModulus*np.sin(.5*np.pi) 
                #else:
                #    cell.lambdaVecX = self.forceModulus*np.cos(-.5*np.pi) 
                #    cell.lambdaVecY = self.forceModulus*np.sin(-.5*np.pi)
                
                
                
                
                
                
                cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                #cluster ID
                cell.dict['clusterID'] = None
                cell.dict['previous_clusterID'] = None
                cell.dict['color'] = None
                
                
    
    def seedTheSpace(self,numberOfCells):
        if len(self.cellList) > 0:
            print 'space already seeded'
            return
        dimX = self.dim.x
        dimY = self.dim.y
        dimZ = self.dim.z
        while numberOfCells > 0:
            x = np.random.randint(int(.01*dimX),int(.99*dimX))
            y = np.random.randint(int(.01*dimY),int(.99*dimY))
            #z = np.random.randint(int(.1*dimZ),int(.9*dimZ))
            z = 0
            pt = CompuCell.Point3D(x,y,z)
            posCell = self.cellField.get(pt) 
            if not posCell: #is medium
                newCell = self.potts.createCellG(pt)
                newCell.type = self.BOIDSA
                #print 'new cell @ ', pt
    
                self.cellField.set(pt,newCell) # to create an extension of that cell
                self.cellField[pt.x-3:pt.x+4,pt.y-3:pt.y+4,0]=newCell
                #print type(self.getCellNeighborDataList(newCell))
                numberOfCells -= 1
        return
    
    def selectForce(self,cell):
        #self.completeForce(cell)
        self.simpleForce(cell)
        
    
    def count_neighbors(self,cCell):
        n = 0
        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cCell):
            if neighbor:
                n+=1
        return n
        
    def createWritingLoc(self):
        #write location
        #instead of writig in ./parameterScan/<iteration>/Simulation/<whatever>
        #I'm writing in ./parameterScan/data/<dir name based on parameters>
        #the iteration # will be suffixed to the name of the data file
        fileDir = os.path.dirname(os.path.abspath(__file__))
        
        cc3dDir, temp = os.path.split(fileDir)
        pScanDir, temp = os.path.split(cc3dDir)
        
        dataDir = os.path.join(pScanDir,'data')
        
        if not os.path.exists(dataDir):
            os.makedirs(dataDir)
        ###
        #CHANGE THE NAME BASED ON TYPE OF SIMULATION
        #(SIMPLE FORCE, COMPLETE FORCE...)
        ###
        
        self.contactEnergy = float(self.getXMLElementValue(
                            ['Plugin', 'Name', 'Contact'], ['Energy', 'Type1', 'boidsA', 'Type2', 'boidsA']))    
        
        saveDirName = ( 'aB_' + str(self.alphaBoids) +
                        '_nB_' + str(self.noiseBoids)+
                        '_J_'+str(self.contactEnergy)+
                        '_dens_' + str(self.density))
        self.saveLoc = os.path.join(dataDir,saveDirName)
        if not os.path.exists(self.saveLoc):
            os.makedirs(self.saveLoc)
#         self.saveLoc = os.path.dirname(os.path.abspath(__file__))
#         self.saveLoc = os.path.join(self.saveLoc,'data')
#         if not os.path.exists(self.saveLoc):
#             os.makedirs(self.saveLoc)
        #str(G_repetitionNumber_G)
        self.instVelFile_name = os.path.join(self.saveLoc,'meanInstantVelocity_'+str(G_repetitionNumber_G)+'.dat')
        self.dtVelFile_name = os.path.join(self.saveLoc,'meanDeltaTVelocity_'+str(G_repetitionNumber_G)+'.dat')
        
        self.dtMeanNeighsVelsFile_name = os.path.join(self.saveLoc,'meanDeltaTNeigsVelocity_'+str(G_repetitionNumber_G)+'.dat')
        self.dtNeighsVelsFileX_name = os.path.join(self.saveLoc,'deltaTNeigsVelocityX_'+str(G_repetitionNumber_G)+'.dat')
        self.dtNeighsVelsFileY_name = os.path.join(self.saveLoc,'deltaTNeigsVelocityY_'+str(G_repetitionNumber_G)+'.dat')
        self.orderParamFile_name = os.path.join(self.saveLoc,'orderParameter_'+str(G_repetitionNumber_G)+'.dat')
        
        #writing the parameters
        with open(os.path.join(self.saveLoc,'parameters.dat'),'w+') as paramFile:
            paramFile.write('targetVolume = '+ str(self.targetVolume) +
                            '\n lambdaVolume = ' + str(self.lambdaVolume) +
                            '\n forceModulus = ' + str(self.forceModulus) +
                            '\n deltaTime = ' + str(self.deltaTime) +
                            '\n alphaBoids =' + str(self.alphaBoids) +
                            '\n betaBoids = ' + str(self.betaBoids) +
                            '\n gammaBoids = ' + str(self.gammaBoids) +
                            '\n noiseBoids = ' + str(self.noiseBoids)+
                            '\n density = '+str(self.density))
        #oppening files / writig the headers
        
        self.instVelFile = open(self.instVelFile_name,'w+')
        self.instVelFile.write('mcs,<Vx(instant)>,std,<Vy(instant)>,std\n')
        
        self.dtVelFile = open(self.dtVelFile_name,'w+')
        self.dtVelFile.write('mcs,<Vx(over dt)>,std,<Vy(over dt)>,std\n')
        
        self.dtMeanNeighsVelsFile = open(self.dtMeanNeighsVelsFile_name,'w+')
        self.dtMeanNeighsVelsFile.write('mcs,<<Vnx(over dt)>_n>_c,std,<<Vny(over dt)>_n>_c,std\n')
        
        self.orderParamFile = open(self.orderParamFile_name,'w+')
        self.orderParamFile.write('mcs, sum(v/|v|)/N\n')
        
        
#         with open(self.instVelFile,'w+') as instVel:
#             instVel.write('mcs,<Vx(instant)>,std,<Vy(instant)>,std\n')
#         with open(self.dtVelFile,'w+') as dtVel:
#             dtVel.write('mcs,<Vx(over dt)>,std,<Vy(over dt)>,std\n')
#         with open(self.dtMeanNeighsVelsFile,'w+') as dtNeighs:
#             dtNeighs.write('mcs,<<Vnx(over dt)>_n>_c,std,<<Vny(over dt)>_n>_c,std\n')  
        
#         with open(self.orderParamFile,'w+') as opf:
#             opf.write('mcs, sum(v/|v|)/N')
        
        
        #file with all the mean neighbor vel needs a header with lengh = # of cells
        headerForFullNVels = 'mcs,'
        for cell in self.cellList:
            headerForFullNVels+='<Vn>,'
        headerForFullNVels+='\n'
        
        self.dtNeighsVelsFileX = open(self.dtNeighsVelsFileX_name,'w+')
        self.dtNeighsVelsFileX.write(headerForFullNVels)
        
        self.dtNeighsVelsFileY = open(self.dtNeighsVelsFileY_name,'w+')
        self.dtNeighsVelsFileY.write(headerForFullNVels)
        
        
        
#         with open(self.dtNeighsVelsFileX, 'w+') as nvx:
#             nvx.write(headerForFullNVels)
#         with open(self.dtNeighsVelsFileY, 'w+') as nvy:
#             nvy.write(headerForFullNVels)
        # self.dtNeighsVelsFileY, self.formatFullNVels
        
    def writeData(self,mcs):
        
        #flushing data every 100mcs
        
        if mcs>0 and mcs%100==0:
            self.instVelFile.flush()
            os.fsync(self.instVelFile.fileno())
        
            self.dtVelFile.flush()
            os.fsync(self.dtVelFile.fileno())
            
            self.dtMeanNeighsVelsFile.flush()
            os.fsync(self.dtMeanNeighsVelsFile.fileno())
            
            self.orderParamFile.flush()
            os.fsync(self.orderParamFile.fileno())
            
            self.dtNeighsVelsFileX.flush()
            os.fsync(self.dtNeighsVelsFileX.fileno())
            
            self.dtNeighsVelsFileY.flush()
            os.fsync(self.dtNeighsVelsFileY.fileno())
            #os.fsync()
        
        
        
        #mean total velocity 
        
        instVelsX = []
        instVelsY = []
        
        dtVelsX = []
        dtVelsY = []
        
        dtNeighsVelsX = []
        dtNeighsVelsY = []
        
        
#         orderParam = np.zeros(2)
        orderParamX=[]
        orderParamY=[]
        for cell in self.cellList:
            instVelsX.append(cell.dict['velocityX_instant'])
            instVelsY.append(cell.dict['velocityY_instant'])
            
            if mcs>self.deltaTime:
                
                #this may or may not be called elsewhere, so, to be safe, I'm putting it here
                neigVx, neigVy = self.getNeighboursMeanVelocity(cell)
                
                
                dtVelsX.append(cell.dict['velocityX_deltaT'])
                dtVelsY.append(cell.dict['velocityY_deltaT'])
                
                vc = np.array([cell.dict['velocityX_deltaT'],cell.dict['velocityY_deltaT']])
                normVC = np.linalg.norm(vc)
                if normVC!=0:
                    orderParamX.append(cell.dict['velocityX_deltaT']/normVC)
                    orderParamY.append(cell.dict['velocityY_deltaT']/normVC)
#                     orderParam += vc/np.linalg.norm(vc)
                else:
                    orderParamX.append(0)
                    orderParamY.append(0)
                
                dtNeighsVelsX.append(cell.dict['mean_neig_velX'])
                dtNeighsVelsY.append(cell.dict['mean_neig_velY'])
        
        orderParamX = np.mean(orderParamX)
        orderParamX_std = np.std(orderParamX)
        
        orderParamY = np.mean(orderParamY)
        orderParamY_std = np.std(orderParamY)
        
        
        
        orderParam = np.linalg.norm((orderParamX,orderParamY))
        
        
        #writting instant velocity
        mivx = np.mean(instVelsX)
        eivx = np.std(instVelsX)
        
        mivy = np.mean(instVelsY)
        eivy = np.std(instVelsY)
        self.instVelFile.write('%i,%f,%f,%f,%f\n'%(mcs, mivx,eivx, mivy,eivy))
        
#         with open(self.instVelFile, 'a+') as ivf:
#             mivx = np.mean(instVelsX)
#             eivx = np.std(instVelsX)
            
#             mivy = np.mean(instVelsY)
#             eivy = np.std(instVelsY)
            
#             ivf.write('%i,%f,%f,%f,%f\n'%(mcs, mivx,eivx, mivy,eivy))
#         self.instVelFile = os.path.join(self.saveLoc,'meanInstantVelocity.dat')
#         self.dtVelFile = os.path.join(self.saveLoc,'meanDeltaTVelocity.dat')
#         self.dtMeanNeighsVelsFile = os.path.join(self.saveLoc,'meanDeltaTNeigsVelocity.dat')
        if mcs>self.deltaTime:
            
            
            #writting dt velocity
            mtvx = np.mean(dtVelsX)
            etvx = np.std(dtVelsX)
            
            mtvy = np.mean(dtVelsY)
            etvy = np.std(dtVelsY)
            
            self.dtVelFile.write('%i,%f,%f,%f,%f\n'%(mcs, mtvx,etvx, mtvy,etvy))
            
#             with open(self.dtVelFile, 'a+') as tvf:
#                 mtvx = np.mean(dtVelsX)
#                 etvx = np.std(dtVelsX)
                
#                 mtvy = np.mean(dtVelsY)
#                 etvy = np.std(dtVelsY)
                
#                 tvf.write('%i,%f,%f,%f,%f\n'%(mcs, mtvx,etvx, mtvy,etvy))

            #writting neighbor mean vel
            mtvnx = np.mean(dtNeighsVelsX)
            etvnx = np.std(dtNeighsVelsX)
            
            mtvny = np.mean(dtNeighsVelsX)
            etvny = np.std(dtNeighsVelsX)
            
            self.dtMeanNeighsVelsFile.write('%i,%f,%f,%f,%f\n'%(mcs,mtvnx,etvnx,mtvny,etvny))
            
            
            
#             with open(self.dtMeanNeighsVelsFile,'a+') as tnvf:
#                 mtvnx = np.mean(dtNeighsVelsX)
#                 etvnx = np.std(dtNeighsVelsX)
                
#                 mtvny = np.mean(dtNeighsVelsX)
#                 etvny = np.std(dtNeighsVelsX)
                
#                 tnvf.write('%i,%f,%f,%f,%f\n'%(mcs,mtvnx,etvnx,mtvny,etvny))
            # self.dtNeighsVelsFileY, self.formatFullNVels
            
            #writting all of neigh vels
            self.dtNeighsVelsFileX.write('%i,'%(mcs))
            for v in dtNeighsVelsX:
                self.dtNeighsVelsFileX.write('%f,'%(v))
            self.dtNeighsVelsFileX.write('\n')
            
            self.dtNeighsVelsFileY.write('%i,'%(mcs))
            for v in dtNeighsVelsY:
                self.dtNeighsVelsFileY.write('%f,'%(v))
            self.dtNeighsVelsFileY.write('\n')  
            
            
#             with open(self.dtNeighsVelsFileX,'a+') as tnvxf:
#                 tnvxf.write('%i,'%(mcs))
#                 for v in dtNeighsVelsX:
#                     tnvxf.write('%f,'%(v))
#                 tnvxf.write('\n')
#             with open(self.dtNeighsVelsFileY,'a+') as tnvyf:
#                 tnvyf.write('%i,'%(mcs))
#                 for v in dtNeighsVelsY:
#                     tnvyf.write('%f,'%(v))
#                 tnvyf.write('\n')    

            #writting the order parameter
            self.orderParamFile.write('%i,%f\n'%(mcs,orderParam))
            
#             with open(self.orderParamFile,'a+') as opf:
#                 opf.write('%i,%f'%(mcs,orderParam))
            
    def positionTracking(self,mcs,cur_Cell):
        self.positionX = cur_Cell.dict['positionX']
        self.positionY = cur_Cell.dict['positionY']
        
        dx = cur_Cell.xCOM - cur_Cell.dict['centerMassX']
        dy = cur_Cell.yCOM - cur_Cell.dict['centerMassY'] 
        
        #if the cell "moved" more than half the space it whent trough the border
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
        
        #new absolute position is previous position + fixed dx
        cur_Cell.dict['positionX'].append(self.positionX[-1]+dx)
        cur_Cell.dict['positionY'].append(self.positionY[-1]+dy)
        
        #with a dt = 1mcs v = dx
        cur_Cell.dict['velocityX_instant'] = dx
        cur_Cell.dict['velocityY_instant'] = dy
        
        #angle of movement:
        cur_Cell.dict['angle'] = np.arctan2(dy,dx)
        
        cur_Cell.dict['centerMassX'] = cur_Cell.xCOM
        cur_Cell.dict['centerMassY'] = cur_Cell.yCOM
    
    
    def getNeighboursMeanVelocity(self,cur_Cell):
        #list of neighbours' velocities
        neigVxList = np.array([])
        neigVyList = np.array([])
        #for whatever reason len(self.getCellNeighborDataList()) doesn't work on karst....
        numberNeighs = self.count_neighbors(cur_Cell)
        if numberNeighs <1:
            #no neighbors, return 0. 
            cur_Cell.dict['mean_neig_velX'] = 0
            cur_Cell.dict['mean_neig_velY'] = 0
            return 0, 0
        
        #KEEP! AS THIS IS THE PROPER WAY, GODDAMN KARST
#         if len(self.getCellNeighborDataList(cur_Cell)) <=1:
#             #no neighbors, return 0. (medium is counted as neighbor, so <=1 not ==0)
#             cur_Cell.dict['mean_neig_velX'] = 0
#             cur_Cell.dict['mean_neig_velY'] = 0
#             return 0, 0
        
        
        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cur_Cell):
            if neighbor:
                neigVxList = np.append(neigVxList,neighbor.dict['velocityX_deltaT'])
                neigVyList = np.append(neigVyList,neighbor.dict['velocityY_deltaT']) 
        
        #taking care of empty lists, just in case something goes wrong
        if len(neigVxList):
            neigVx = np.mean(neigVxList)
            cur_Cell.dict['mean_neig_velX'] = neigVx
        else:  
            cur_Cell.dict['mean_neig_velX'] = 0      
            neigVx = 0
        if len(neigVyList):
            neigVy = np.mean(neigVyList)
            cur_Cell.dict['mean_neig_velY'] = neigVy
        else:
            cur_Cell.dict['mean_neig_velY'] = 0
            neigVy = 0
        
        return neigVx, neigVy
    
    def getNormNeighboursMeanVelocity(self,cur_Cell):
        #list of neighbours' velocities
        
        #for whatever reason len(self.getCellNeighborDataList()) doesn't work on karst....
        numberNeighs = self.count_neighbors(cur_Cell)
        if numberNeighs <1:
            #no neighbors, return 0. 
            return 0, 0
        
        neigVxList = np.array([])
        neigVyList = np.array([])
        #KEEP! AS THIS IS THE PROPER WAY, GODDAMN KARST
#         if len(self.getCellNeighborDataList(cur_Cell)) <=1:
#             #no neighbors, return 0. (medium is counted as neighbor, so <=1 not ==0)
#             cur_Cell.dict['mean_neig_velX'] = 0
#             cur_Cell.dict['mean_neig_velY'] = 0
#             return 0, 0
        
        
        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cur_Cell):
            if neighbor:
                neigVxList = np.append(neigVxList,neighbor.dict['norm_velocityX_deltaT'])
                neigVyList = np.append(neigVyList,neighbor.dict['norm_velocityY_deltaT']) 
        
        #taking care of empty lists, just in case something goes wrong
        if len(neigVxList):
            neigVx = np.mean(neigVxList)
        else:  
            neigVx = 0
        if len(neigVyList):
            neigVy = np.mean(neigVyList)
        else:
            neigVy = 0
        
        return neigVx, neigVy
    
    
    
    def completeForce(self,cell): # F = a * V + b * <Vn> - d * F(t-dt)  all unitary vectors
        
        #normalized velocity
        cellVx = cell.dict['norm_velocityX_deltaT']
        cellVy = cell.dict['norm_velocityY_deltaT']
        V = np.array([cellVx,cellVy])
        
        #getting mean normalized velocities of neighbours
        neigV = np.zeros(2)
        neigV = self.getNormNeighboursMeanVelocity(cell)
        neigV = np.array(neigV)
        
        #previous force normalized
        pForce = np.array([cell.dict['previousForceX']/self.forceModulus,
                           cell.dict['previousForceY']/self.forceModulus])
        
        #noise vector (normalized)
        noisevec = np.random.normal(size=2)
        normn = np.linalg.norm(noisevec)
        noisevec = noisevec/normn
        
        #new force
        force = np.zeros(2)
        force[:] = self.alphaBoids*V[:] + self.betaBoids*neigV[:] - self.gammaBoids*pForce[:] + self.noiseBoids * noisevec[:]
        normForce = np.linalg.norm(force)
        force = force / normForce #also normalized
        
        #updating the force
        cell.lambdaVecX = self.forceModulus*force[0] 
        cell.lambdaVecY = self.forceModulus*force[1]  
        
        #updating cell dictionaries
        cell.dict['forceAngle'] = np.arctan2(force[0],force[1])
        cell.dict['previousForceX'] = self.forceModulus*force[0]
        cell.dict['previousForceY'] = self.forceModulus*force[1]
        
        
        
    def completeForce_old(self,cell): # F = a * V + b * <Vn> - d * F(t-dt)
        #need to vectorize this, make it neater
        cellVx = cell.dict['velocityX_deltaT']
        cellVy = cell.dict['velocityY_deltaT']
        
        #getting mean velocities of neighbours
        neigVx, neigVy = self.getNeighboursMeanVelocity(cell)
        
        #calculating new force in x and y
        forceX = self.alphaBoids*cellVx + self.betaBoids*neigVx - self.gammaBoids*cell.dict['previousForceX']
        forceY = self.alphaBoids*cellVy + self.betaBoids*neigVy - self.gammaBoids*cell.dict['previousForceY']
        
        #getting the angle of new force
        cell.dict['forceAngle'] = np.arctan2(forceY,forceX) + self.noiseBoids * np.random.uniform(-np.pi,np.pi)
        
        #aplying standard force with new angle
        cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
        cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle'])  
        
        cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
        cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle'])  
    
    def simpleForce(self,cCell): # d orientation /dt = -A orientation + (1-A) velocity + noise
        #orientation is a vector
        dt = 1 #this calculation will happen every mcs
        
        #velocity vector (normalized)
        v = np.array([cCell.dict['velocityX_deltaT'],cCell.dict['velocityY_deltaT']])
        normV = np.linalg.norm(v)
        v = v/normV
        
        #noise vector (normalized)
        nv = np.random.normal(size=2)
        normn = np.linalg.norm(nv)
        nv = nv/normn
        
        # -A orientation + (1-A) velocity
        oneAlpha = 1- self.alphaBoids
        variation = - self.alphaBoids * cCell.dict['orientation'] + oneAlpha*v
        
        #new orientation (normalized)
        cCell.dict['orientation'] += dt*variation + self.noiseBoids * nv
        normOrient = np.linalg.norm(cCell.dict['orientation'])
        cCell.dict['orientation'] = cCell.dict['orientation']/normOrient
        
        cCell.lambdaVecX = self.forceModulus*cCell.dict['orientation'][0] 
        cCell.lambdaVecY = self.forceModulus*cCell.dict['orientation'][1]  
        
        cCell.dict['forceAngle'] = np.arctan2(cCell.dict['orientation'][0],cCell.dict['orientation'][1])
        cCell.dict['previousForceX'] = self.forceModulus*cCell.dict['orientation'][0]
        cCell.dict['previousForceY'] = self.forceModulus*cCell.dict['orientation'][1]
    
    
    def assignClusterIDs(self,curr_Cell):
        #print type(self.getCellNeighborDataList(curr_Cell))
        #print len(self.getCellNeighborDataList(curr_Cell))
        #for whatever reason len(self.getCellNeighborDataList()) doesn't work on karst....
        
        numberNeighs = self.count_neighbors(curr_Cell)
        if numberNeighs < 1: 

            curr_Cell.dict['clusterID'] = None

            return
        #KEEP! AS THIS IS THE PROPER WAY, GODDAMN KARST
#         if len(self.getCellNeighborDataList(curr_Cell)) <= 1: #medium is counted!!! 

#             curr_Cell.dict['clusterID'] = None

#             return
        
        neigsIDs = []
        if curr_Cell.dict['clusterID'] != None:
            neigsIDs.append(curr_Cell.dict['clusterID'])
        
        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(curr_Cell):
            if neighbor: 
                if neighbor.dict['clusterID'] != None:
                    neigsIDs.append(neighbor.dict['clusterID'])
        
        if len(neigsIDs) == 0:
            newID = max(self.usedClusterIDs)+1
            #self.currentIDs
            
            
            #print 'new ID:', newID
            self.usedClusterIDs.add(newID)
            curr_Cell.dict['clusterID'] = newID
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(curr_Cell):
                if neighbor: 
                    neighbor.dict['clusterID'] = newID
        else:
            id = min(neigsIDs)
            curr_Cell.dict['clusterID'] = id
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(curr_Cell):
                if neighbor: 
                    neighbor.dict['clusterID'] = id
        
        if len(set(neigsIDs)) > 1 and sorted(neigsIDs) not in self.equivalentIDs:
            self.equivalentIDs.append(sorted(neigsIDs))
    
    def clusterBreakingDetection(self):
        
        ## there's the issue of a cluster breaking in
        ## two clusters, detecting that and spliting the ID.
        ## 
        ## I could keep track of the neighbors of every cell. Then, if the neighborhood
        ## has changed check if the cell that left still has some neighbors. If this is so
        ## then I can just assign a new ID to those that left. I'll do this check after the reasignment.
        ## I'll have to delete the list of equivalencies before entering this function.
        ## I also need a dictionary for the neighbors, probably the neighbors IDs. 
        ## Also, I only care if the neigborhood has diminished.
        ##
        ## well, that didn't work
        ##
        ## Maybe I could reset the whole cluster ID thing every MCS. To keep history I could have the old 
        ## ID in a dict. Then assing that one over the new one. But before that do the whole reassingDoubleIDs
        ## and then have some kind of method for the separation. something like, if 2 cells have the same old ID 
        ## only one of them gets to keep it.
        ## 
        ## right now it works, but it takes a couple of mcs to work. Well, latest test has it working imediatly.
        ## so probably there was some glitch before. idk
        
        self.equivalencyToOld  = set([])
        
        for cell in self.cellList:
            if (cell.dict['previous_clusterID'] != None and 
                cell.dict['clusterID'] != None):
                #
                self.equivalencyToOld.add((cell.dict['previous_clusterID'],cell.dict['clusterID']))
        toDiscard = set([])
        for pairA in self.equivalencyToOld:
            for pairB in self.equivalencyToOld:
                if (pairA[0] == pairB[0] and
                    pairA[1] != pairB[1]):
                    #
                    toDiscard.add(pairB)
        
        for d in toDiscard:
            self.equivalencyToOld.discard(d)
        
        #get the old ids that will be used again
        oldInUse = [ids[0] for ids in self.equivalencyToOld]
        
        #set of current ids that will conflict with old ids
        #ie, after reassigning cells that will conflict with that reasignment,
        #having a new id that has the same value that an in use old id.
        repeatedNew = set([])
        for d in toDiscard:
            if d[1] in oldInUse:
                
                #if there's a conflict generate a new id for the conflicting cell
                newID = max(self.usedClusterIDs)+1
                repeatedNew.add((newID,d[1]))
        
        for cell in self.cellList:
            #solve conflict
            self.reassingDoubleIDs(cell,repeatedNew)
        
        for cell in self.cellList:
            self.reassingDoubleIDs(cell,self.equivalencyToOld)
    
    def reassingDoubleIDs(self,curr_Cell,eqs):
        if curr_Cell.dict['clusterID'] == None:
            return#
        for e in eqs:
                for notSmalestID in e[1:]:
                    if curr_Cell.dict['clusterID'] == notSmalestID:
                        curr_Cell.dict['clusterID'] = e[0]
                        #print 'reasigned', notSmalestID, 'to', e[0]
                        return
    
    def coloringByID(self,inUseClusterIDs):
        
        color = 0
        if len(inUseClusterIDs)>0:
            maxID = max(inUseClusterIDs)
            inUseClusterIDs = sorted(list(set(inUseClusterIDs)))
            interval = float(maxID)/len(inUseClusterIDs)
            #print maxID, len(inUseClusterIDs), interval
            for cell in self.cellList:
#                 print color
                for i, id in enumerate(inUseClusterIDs):
                    if cell.dict['clusterID'] == id:
                        color = .5*((i+1)*interval)+.5*maxID
                        cell.dict['color'] = color
                        break
    
    def updateFields(self,mcs):
        
        ##for better colors in the cluster IDs I'll set them to be a % of the max ID
        ##(previouslly the ID itself was used as color
        for cCell in self.cellList:
            
            #updating scalar field for cluster coloring.
            if cCell.dict['color'] != None:
                self.scalarClusterIDfield[cCell] = cCell.dict['color']
            else:
                self.scalarClusterIDfield[cCell] = 0
            
            #updating scalar field for intantaneous angle of movement:
            self.scalarAngleField[cCell] = 180*(cCell.dict['angle']+np.pi)/np.pi
            #updating vector field for intantaneous velocity:
            self.vectorIVelocityField[cCell] = [cCell.dict['velocityX_instant'], cCell.dict['velocityY_instant'], 0]
            
            if mcs>self.deltaTime:
                self.vectorTVelocityField[cCell] = [cCell.dict['velocityX_deltaT'], cCell.dict['velocityY_deltaT'], 0]
                self.vectorForceField[cCell] = [-cCell.dict['previousForceX'], -cCell.dict['previousForceY'], 0]
        
#         self.vectorTVelocityField = self.createVectorFieldCellLevelPy('dt_Velocity')
        
#         self.vectorForceField = self.createVectorFieldCellLevelPy("Force")
    
    
        
    
    
    
    
    
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
        #####################################################
        
        ##I need to identify the clusters. The most straitforward way I can think is
        ## simply looping the cells, checking the neighbors and assingning a colective ID. If the
        ## neighbor already has an ID use that one instead. This is what's done in the 1st loop.
        ## Of course, this may lead to a single
        ## cluster having multiple IDs, so a second pass will be necessary. 
        ## there's also the issue of tracking the used IDs (witch will only really be an issue if 
        ## the number of used IDs gets to be close to 2**32, I guess).
        ########
        
        #first cell loop to update tracking of position, velocities, id assignment
        
        self.equivalentIDs = []
        #first loop for id assignment

        for cell in self.cellList:
            cell.dict['color'] = None
            self.assignClusterIDs(cell)
            self.positionTracking(mcs,cell)
            
            if mcs == 5:
                cell.lambdaVolume = self.lambdaVolume
            
            if mcs > self.deltaTime:
                vx = ( cell.dict['positionX'][-1] - cell.dict['positionX'][-self.deltaTime-1] )
                vx = vx/self.deltaTime
                
                vy = ( cell.dict['positionY'][-1] - cell.dict['positionY'][-self.deltaTime-1] )
                vy = vy/self.deltaTime
                
                cell.dict['velocityX_deltaT'] = vx
                cell.dict['velocityY_deltaT'] = vy
                
                v = np.array([vx,vy])
                nv = np.linalg.norm(v)
                v = v/nv
                
                cell.dict['norm_velocityX_deltaT'] = v[0]
                cell.dict['norm_velocityY_deltaT'] = v[1]
                
        
        for cell in self.cellList:
            self.reassingDoubleIDs(cell,self.equivalentIDs)
        
        
        #################################################
        
        ## there's the issue of a cluster breaking in
        ## two clusters, detecting that and spliting the ID.
        ## 
        ## I could keep track of the neighbors of every cell. Then, if the neighborhood
        ## has changed check if the cell that left still has some neighbors. If this is so
        ## then I can just assign a new ID to those that left. I'll do this check after the reasignment.
        ## I'll have to delete the list of equivalencies before entering this function.
        ## I also need a dictionary for the neighbors, probably the neighbors IDs. 
        ## Also, I only care if the neigborhood has diminished.
        ##
        ## well, that didn't work
        ##
        ## Maybe I could reset the whole cluster ID thing every MCS. To keep history I could have the old 
        ## ID in a dict. Then assing that one over the new one. But before that do the whole reassingDoubleIDs
        ## and then have some kind of method for the separation. something like, if 2 cells have the same old ID 
        ## only one of them gets to keep it.
        ## 
        ## right now it works, but it takes a couple of mcs to work. Well, latest test has it working imediatly.
        ## so probably there was some glitch before. idk
        
        
        self.clusterBreakingDetection()
        
        #################################################
        # cell loop to calculate the boid's force 
        if mcs > self.deltaTime:
            for cell in self.cellList:
                self.selectForce(cell)
        
        #################################################
        
        self.writeData(mcs)
        
        #################################################
        inUseClusterIDs = []##will use this to assign colors
        for cell in self.cellList:
            if cell.dict['clusterID'] != None:
                inUseClusterIDs.append(cell.dict['clusterID'])
        
        self.coloringByID(inUseClusterIDs)
        self.currentIDs = set(inUseClusterIDs)
        #updating extra fields
        self.updateFields(mcs)
        
        

        #################################################
        #cleanup
        #print '!!!!!!!!!!!!'
        #print self.currentIDs
        #print '!!!!!!!!!!!!'
        self.oldIDs = self.currentIDs.copy()
        self.usedClusterIDs = None
        self.usedClusterIDs = set([0])
        self.currentIDs = None
        self.currentIDs = set([0])
        for cell in self.cellList:
            cell.dict['previous_clusterID'] = cell.dict['clusterID']
            cell.dict['clusterID'] = None
        
        
        
        
        
        
    def finish(self):
        # Finish Function gets called after the last MCS
        self.instVelFile.close()
        
        self.dtVelFile.close()
        
        self.dtMeanNeighsVelsFile.close()
        
        self.orderParamFile.close()
        
        self.dtNeighsVelsFileX.close()
        
        self.dtNeighsVelsFileY.close()
        