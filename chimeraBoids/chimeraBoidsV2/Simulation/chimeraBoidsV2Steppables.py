from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys

import numpy as np
import os


global G_targetVolume_G
global G_lambdaVolume_G
global G_forceModulus_G
global G_deltaTime_G

global G_alphaBoids_G
global G_betaBoids_G
global G_gammaBoids_G
global G_noise_G


#parameter scan variables
G_targetVolume_G = 64.
G_lambdaVolume_G = 8.
G_forceModulus_G = -20 #speeds are tipically 1/100 of force.
G_deltaTime_G = 10

G_alphaBoids_G = .1#5.5
G_betaBoids_G = 1.5
G_gammaBoids_G = .1
G_noise_G = .1




class chimeraBoidsV2Steppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        #initiating extra fields
        
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
                
                
    def selectForce(self,cell):
        #self.completeForce(cell)
        self.simpleForce(cell)
        
    def createWritingLoc(self):
        #write location
        
        self.saveLoc = os.path.dirname(os.path.abspath(__file__))
        self.saveLoc = os.path.join(self.saveLoc,'data')
        if not os.path.exists(self.saveLoc):
            os.makedirs(self.saveLoc)
        
        self.instVelFile = os.path.join(self.saveLoc,'meanInstantVelocity.dat')
        self.dtVelFile = os.path.join(self.saveLoc,'meanDeltaTVelocity.dat')
        self.dtNeighsVelsFile = os.path.join(self.saveLoc,'meanDeltaTNeigsVelocity.dat')
        
        #writing the parameters
        with open(os.path.join(self.saveLoc,'parameters.dat'),'w+') as paramFile:
            paramFile.write('targetVolume = '+ str(self.targetVolume) +
                            '\n lambdaVolume = ' + str(self.lambdaVolume) +
                            '\n forceModulus = ' + str(self.forceModulus) +
                            '\n deltaTime = ' + str(self.deltaTime) +
                            '\n alphaBoids =' + str(self.alphaBoids) +
                            '\n betaBoids = ' + str(self.betaBoids) +
                            '\n gammaBoids = ' + str(self.gammaBoids) +
                            '\n noiseBoids = ' + str(self.noiseBoids))
        #writig the headers
        with open(self.instVelFile,'w+') as instVel:
            instVel.write('mcs,<Vx(instant)>,std,<Vy(instant)>,std\n')
        with open(self.dtVelFile,'w+') as dtVel:
            dtVel.write('mcs,<Vx(over dt)>,std,<Vy(over dt)>,std\n')
        with open(self.dtNeighsVelsFile,'w+') as dtNeighs:
            dtNeighs.write('mcs,<<Vnx(over dt)>_n>_c,std,<<Vny(over dt)>_n>_c,std\n')    
        
        
    def writeData(self,mcs):
        #mean total velocity 
        
        instVelsX = []
        instVelsY = []
        
        dtVelsX = []
        dtVelsY = []
        
        dtNeighsVelsX = []
        dtNeighsVelsY = []
        
        for cell in self.cellList:
            instVelsX.append(cell.dict['velocityX_instant'])
            instVelsY.append(cell.dict['velocityY_instant'])
            
            if mcs>self.deltaTime:
                dtVelsX.append(cell.dict['velocityX_deltaT'])
                dtVelsY.append(cell.dict['velocityY_deltaT'])
                
                dtNeighsVelsX.append(cell.dict['mean_neig_velX'])
                dtNeighsVelsY.append(cell.dict['mean_neig_velY'])
        
        with open(self.instVelFile, 'a+') as ivf:
            mivx = np.mean(instVelsX)
            eivx = np.std(instVelsX)
            
            mivy = np.mean(instVelsY)
            eivy = np.std(instVelsY)
            
            ivf.write('%i,%f,%f,%f,%f\n'%(mcs, mivx,eivx, mivy,eivy))
#         self.instVelFile = os.path.join(self.saveLoc,'meanInstantVelocity.dat')
#         self.dtVelFile = os.path.join(self.saveLoc,'meanDeltaTVelocity.dat')
#         self.dtNeighsVelsFile = os.path.join(self.saveLoc,'meanDeltaTNeigsVelocity.dat')
        if mcs>self.deltaTime:
            with open(self.dtVelFile, 'a+') as tvf:
                mtvx = np.mean(dtVelsX)
                etvx = np.std(dtVelsX)
                
                mtvy = np.mean(dtVelsY)
                etvy = np.std(dtVelsY)
                
                tvf.write('%i,%f,%f,%f,%f\n'%(mcs, mtvx,etvx, mtvy,etvy))
            
            with open(self.dtNeighsVelsFile,'a+') as tnvf:
                mtvnx = np.mean(dtNeighsVelsX)
                etvnx = np.std(dtNeighsVelsX)
                
                mtvny = np.mean(dtNeighsVelsX)
                etvny = np.std(dtNeighsVelsX)
                
                tnvf.write('%i,%f,%f,%f,%f\n'%(mcs,mtvnx,etvnx,mtvny,etvny))
        
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
        
        if len(self.getCellNeighborDataList(cur_Cell)) <=1:
            #no neighbors, return 0. (medium is counted as neighbor, so <=1 not ==0)
            cur_Cell.dict['mean_neig_velX'] = 0
            cur_Cell.dict['mean_neig_velY'] = 0
            return 0, 0
        
        
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
        
    def completeForce(self,cell): # F = a * V + b * <Vn> - d * F(t-dt)
        cellVx = cell.dict['velocityX_deltaT']
        cellVy = cell.dict['velocityY_deltaT']
        
        #getting mean velocities of neighbours
        neigVx, neigVy = self.getNeighboursMeanVelocity(cell)
        
        #calculating new force in x and y
        forceX = self.alphaBoids*cellVx + self.betaBoids*neigVx - self.noiseBoids*cell.dict['previousForceX']
        forceY = self.alphaBoids*cellVy + self.betaBoids*neigVy - self.noiseBoids*cell.dict['previousForceY']
        
        #getting the angle of new force
        cell.dict['forceAngle'] = np.arctan2(forceY,forceX) + self.gammaBoids * np.random.uniform(-np.pi,np.pi)
        
        #aplying standard force with new angle
        cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
        cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle'])  
        
        cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
        cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle'])  
    
    def simpleForce(self,cCell): # d orientation /dt = -A orientation + (1-A) velocity
        #self.deltaTime
        dt = 1 #this calculation will happen every mcs
        v = np.array([cCell.dict['velocityX_deltaT'],cCell.dict['velocityY_deltaT']])
        normV = np.linalg.norm(v)
        v = v/normV
        oneAlpha = 1- self.alphaBoids
        variation = - self.alphaBoids * cCell.dict['orientation'] + oneAlpha*v
        
        #orientation is a vector
        cCell.dict['orientation'] += dt*variation
        normOrient = np.linalg.norm(cCell.dict['orientation'])
        cCell.dict['orientation'] = cCell.dict['orientation']/normOrient
        
        cCell.lambdaVecX = self.forceModulus*cCell.dict['orientation'][0] 
        cCell.lambdaVecY = self.forceModulus*cCell.dict['orientation'][1]  
        
        cCell.dict['forceAngle'] = np.arctan2(cCell.dict['orientation'][0],cCell.dict['orientation'][1])
        cCell.dict['previousForceX'] = self.forceModulus*cCell.dict['orientation'][0]
        cCell.dict['previousForceY'] = self.forceModulus*cCell.dict['orientation'][1]
    
    
    def assignClusterIDs(self,curr_Cell):
        #print len(self.getCellNeighborDataList(curr_Cell))
        if len(self.getCellNeighborDataList(curr_Cell)) <= 1: #medium is counted!!! 

            curr_Cell.dict['clusterID'] = None

            return
        
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
            
            
            print 'new ID:', newID
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
            print maxID, len(inUseClusterIDs), interval
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
                self.vectorForceField[cCell] = [cCell.dict['previousForceX'], cCell.dict['previousForceY'], 0]
        
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
            
            
            
            if mcs > self.deltaTime:
                vx = ( cell.dict['positionX'][-1] - cell.dict['positionX'][-self.deltaTime-1] )
                vx = vx/self.deltaTime
                
                vy = ( cell.dict['positionY'][-1] - cell.dict['positionY'][-self.deltaTime-1] )
                vy = vy/self.deltaTime
                
                cell.dict['velocityX_deltaT'] = vx
                cell.dict['velocityY_deltaT'] = vy
        
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
        print '!!!!!!!!!!!!'
        print self.currentIDs
        print '!!!!!!!!!!!!'
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
        pass
        