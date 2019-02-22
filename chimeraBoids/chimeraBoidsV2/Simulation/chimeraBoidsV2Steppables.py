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
global G_deltaBoids_G


#parameter scan variables
G_targetVolume_G = 64.
G_lambdaVolume_G = 8.
G_forceModulus_G = -20
G_deltaTime_G = 10

G_alphaBoids_G = 5.5
G_betaBoids_G = 1.5
G_gammaBoids_G = .5
G_deltaBoids_G = .1




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
        self.deltaBoids = G_deltaBoids_G
        
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
                
                cell.dict['previousForceX'] = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.dict['previousForceY'] = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                
                
                cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle']) 
                
                #cluster ID
                cell.dict['clusterID'] = None
                cell.dict['color'] = None
                
                #neighbor tracking:
                cell.dict['true_breakAway'] = False
                cell.dict['previous_neighbors'] = []
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        cell.dict['previous_neighbors'].append(neighbor.id)
    

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
        
        for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cur_Cell):
            if neighbor:
                neigVxList = np.append(neigVxList,neighbor.dict['velocityX_deltaT'])
                neigVyList = np.append(neigVyList,neighbor.dict['velocityY_deltaT']) 
        
        #taking care of empty lists
        if len(neigVxList):
            neigVx = np.mean(neigVxList)
        else:       
            neigVx = 0
        if len(neigVyList):
            neigVy = np.mean(neigVyList)
        else:
            neigVy = 0
        
        return neigVx, neigVy
        
    def calculateForce(self,cV,nV,pF):
        force = self.alphaBoids*cV + self.betaBoids*nV - self.deltaBoids*pF
        return force
    
    
    
    
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
    
    
    def reassingDoubleIDs(self,curr_Cell,eqs):
        if curr_Cell.dict['clusterID'] == None:
            return#
        for e in eqs:
                for notSmalestID in e[1:]:
                    if curr_Cell.dict['clusterID'] == notSmalestID:
                        curr_Cell.dict['clusterID'] = e[0]
                        print 'reasigned', notSmalestID, 'to', e[0]
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
        
        #############################################################
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
        ##############################################################
        
        #####################################################
        
        ##I need to identify the clusters. The most straitforward way I can think is
        ## simply looping the cells, checking the neighbors and assingning a colective ID. If the
        ## neighbor already has an ID use that one instead. This is what's done in the 1st loop.
        ## Of course, this may lead to a single
        ## cluster having multiple IDs, so a second pass will be necessary. 
        ## there's also the issue of tracking the used IDs (witch will only really be an issue if 
        ## the number of used IDs gets to be close to 2**32, I guess).
        ########
        ## All of the above I believe is solved. However there's the issue of a cluster breaking in
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
        
        
        
        for cell in self.cellList:
            self.reassingDoubleIDs(cell,self.equivalentIDs)
        
        '''
        completeBreakAways = set([])
        #trueBreakAways = set([])
        for cell in self.cellList:
            if (len(cell.dict['previous_neighbors']) == 0 or #if the cell didn't have neighs dont bother
                len(self.getCellNeighborDataList(cell)) <= 1):# or #if the cell is alone now don't bother
                #len(cell.dict['previous_neighbors']) <= len(currentNeigs)):# or#if the n of neigs increased dont bother
                #cell.id in completeBreakAways):
                
                cell.dict['true_breakAway'] = False
                break
            
            #comparing the neig before and now
            currentNeigs = []
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:  
                    currentNeigs.append(neighbor.id)
            if len(cell.dict['previous_neighbors']) <= len(currentNeigs):
                break
            #getting the cells that broke away
            # setA - setB = set of elements in A that are not in B
            breakAways = list(set(cell.dict['previous_neighbors']) - set(currentNeigs))
            for bA in self.cellList:                
                if bA.id in breakAways: #this cell is a break away
                    stillConnected = False
                    for neighbor, commonSurfaceArea in self.getCellNeighborDataList(bA):
                        if neighbor:   
                            if neighbor.id in currentNeigs: #there's still a connection
                                stillConnected = True
                                cell.dict['true_breakAway'] = False
                else:
                    break
                if stillConnected == False:#if there's no direct conection
                    cell.dict['true_breakAway'] = True
                    #trueBreakAways.add(bA.id)
            completeBreakAways.update(breakAways)
        
        self.divisionID = set([max(self.usedClusterIDs)])
        self.divisionEquivalentIDs = []
        for cell in self.cellList:
            if cell.dict['true_breakAway'] == False:
                break
            else:
                neigsNewIDs = []
                if cell.dict['clusterID'] in self.divisionID:
                    neigsNewIDs.append(cell.dict['clusterID'])
                
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        if neighbor.dict['clusterID'] in self.divisionID:
                            neigsNewIDs.append(neighbor.dict['clusterID'])
                
                #if no neigh has a new id
                if len(neigsNewIDs) == 0:
                    #create a new id and use that
                    newID = max(self.divisionID)+1
                    self.divisionID.add(newID)
                    cell.dict['clusterID'] = newID
                    for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                        if neighbor:
                            neighbor.dict['clusterID'] = newID
                else:
                    #use the tiniest of the neig ids
                    id = min(neigsNewIDs)
                    cell.dict['clusterID'] = id
                    for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                        if neighbor:
                            neighbor.dict['clusterID'] = id
                #update the equivalencies
                if len(set(neigsNewIDs)) > 1 and sorted(neigsNewIDs) not in self.divisionEquivalentIDs:
                    self.divisionEquivalentIDs.append(neigsNewIDs)
        
        #reasign duplicates
        for cell in self.cellList:
            self.reassingDoubleIDs(cell,self.divisionEquivalentIDs)
            cell.dict['previous_neighbors']=[]
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        cell.dict['previous_neighbors'].append(neighbor.id)
        
        
        
        
        self.usedClusterIDs.update(self.divisionID)    
        '''
        
        '''
        #find breakAways
        breakAways = []
        
        for cell in self.cellList:
            foundBreakAways = []
            if len(cell.dict['previous_neighbors']) == 0:
                break
            if len(self.getCellNeighborDataList(cell)) <= 1:
                break            
            if len(cell.dict['previous_neighbors']) < len(currentNeigs):
                break
            
            currentNeigs = []
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:  
                    currentNeigs.append(neighbor.id)
            
            
            # setA - setB = set of elements in A that are not in B
            breakAway = list(set(cell.dict['previous_neighbors']) - set(currentNeigs))
            for b in breakAway:
                if b not in foundBreakAways:
                    foundBreakAways.append(b)
        
            breakAways.extend(foundBreakAways)
        breakAways = list(set(breakAways))
        
        
        #reasign breakAways
        
        self.divisionID = set([max(self.usedClusterIDs)])
        self.divisionEquivalentIDs = []
        for cell in self.cellList:
            if len(self.getCellNeighborDataList(cell)) <= 1:
                break
            if cell.id in breakAways:
                neigsNewIDs = []
                if cell.dict['clusterID'] in self.divisionID:
                    neigsNewIDs.append(cell.dict['clusterID'])
                
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        if neighbor.dict['clusterID'] in self.divisionID:
                            neigsNewIDs.append(neighbor.dict['clusterID'])
            else:
                break
            
            if len(neigsNewIDs) == 0:
                newID = max(self.divisionID)+1
                self.divisionID.add(newID)
                cell.dict['clusterID'] = newID
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        neighbor.dict['clusterID'] = newID
            else:
                id = min(neigsNewIDs)
                cell.dict['clusterID'] = id
                for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        neighbor.dict['clusterID'] = id
            if len(set(neigsNewIDs)) > 1 and sorted(neigsNewIDs) not in self.divisionEquivalentIDs:
                self.divisionEquivalentIDs.append(neigsNewIDs)
        
        
        for cell in self.cellList:
            self.reassingDoubleIDs(cell,self.divisionEquivalentIDs)
            cell.dict['previous_neighbors']=[]
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                    if neighbor:
                        cell.dict['previous_neighbors'].append(neighbor.id)
        
        
        self.usedClusterIDs.update(self.divisionID)
        '''
        #################################################
        
        
        # cell loop to calculate the boid's force 
        if mcs > self.deltaTime:
            for cell in self.cellList:
                
                cellVx = cell.dict['velocityX_deltaT']
                cellVy = cell.dict['velocityY_deltaT']
                
                #getting mean velocities of neighbours
                neigVx, neigVy = self.getNeighboursMeanVelocity(cell)
                
                #calculating new force in x and y
                forceX = self.calculateForce(cellVx,neigVx,cell.dict['previousForceX'])                
                forceY = self.calculateForce(cellVy,neigVy,cell.dict['previousForceY'])                
                
                #getting the angle of new force
                # # not sure about the implementation of noise
                cell.dict['forceAngle'] = np.arctan2(forceY,forceX) + self.gammaBoids * np.random.uniform(-np.pi,np.pi)
                
                #aplying standard force with new angle
                cell.lambdaVecX = self.forceModulus*np.cos(cell.dict['forceAngle']) 
                cell.lambdaVecY = self.forceModulus*np.sin(cell.dict['forceAngle'])  
               
                cell.dict['previousForceX'] = forceX
                cell.dict['previousForceY'] = forceY 
#        
        
       
        
        
        inUseClusterIDs = []##will use this to assign colors
        for cell in self.cellList:
            if cell.dict['clusterID'] != None:
                inUseClusterIDs.append(cell.dict['clusterID'])
        
        self.coloringByID(inUseClusterIDs)
        self.currentIDs = set(inUseClusterIDs)
        #updating extra fields
        self.updateFields(mcs)
        
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        