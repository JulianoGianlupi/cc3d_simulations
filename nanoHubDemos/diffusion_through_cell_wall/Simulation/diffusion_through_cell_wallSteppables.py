
from PySteppables import *
import CompuCell
import sys
class diffusion_through_cell_wallSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        self.wallThickness = 3
        self.borderSideLenght = 32-3
        
        self.relaxTime = 50
        
    def start(self):
        # any code in the start function runs before MCS=0
        #selecting grid type (hex/brick)
        self.selectGrid(type = 'hex')
        
        #create a source
        pt = CompuCell.Point3D(int(0.9*self.dim.x),105,0)
        newSource = self.potts.createCellG(pt)
        newSource.type = self.SOURCE
        self.cellField[int(pt.x-2):pt.x,
                       int(pt.y-2):pt.y,
                       0]=newSource
        #create a sink
        pt = CompuCell.Point3D(self.dim.x-5,self.dim.y-5,0)
        newSink = self.potts.createCellG(pt)
        newSink.type = self.SINK
        self.cellField[int(pt.x-2):pt.x,
                       int(pt.y-2):pt.y,
                       0]=newSink
        
        
    def selectGrid(self,type = 'brick'):
        if type == 'brick':
            self.brickGrid()
        elif type == 'hex':
            self.hexGrid() 
        else:
            raise Exception('Not a valid grid type.\n Select type = "hex" or type = "brick"')
            self.stopSimulation()
            
            
    def hexGrid(self):
        c = self.borderSideLenght + self.wallThickness
        initialCenter = (c,c)
        center = list(initialCenter)
        z=int(0)
        i=1
        pt =CompuCell.Point3D(0,0,0)
        newBorder = self.potts.createCellG(pt)
        newBorder.type = self.BOUNDARY
        while center[1] + self.borderSideLenght + self.wallThickness < self.dim.y:
            while center[0] + self.borderSideLenght + self.wallThickness < self.dim.x:
                x = center[0]
                y = center[1]
                pt = CompuCell.Point3D(x,y,z)
                newCell = self.potts.createCellG(pt)
                newCell.type = self.CELL
                newCell.targetVolume = 1.5*np.sqrt(3)*self.borderSideLenght**2
                newCell.lambdaVolume = 50
                self.cellField[int(pt.x-0.6*self.borderSideLenght):int(pt.x+0.6*self.borderSideLenght),
                               int(pt.y-0.6*self.borderSideLenght):int(pt.y+0.6*self.borderSideLenght),
                               0]=newCell
                #top/bottom left wall
                x0 = int(self.borderSideLenght * np.cos(np.pi*150./180.))        
                a = (1-np.sin(np.pi*150/180))/np.cos(np.pi*150./180.)
                b = self.borderSideLenght * np.sin(np.pi*150./180.)
                xThick = self.wallThickness*np.cos(np.pi*60./180.)
                yThick = self.wallThickness*np.sin(np.pi*120./180.)
                for xi in range(x0,0):
                    #print xi
                    x = int(center[0]+xi-xThick)
                    yi =( b - a * (xi - x0))
                    #print yi
                    y = int(center[1]+yi)
                    pt = CompuCell.Point3D(x,y,z)
                    #print pt
                    #newBorder = self.potts.createCellG(pt)
                    #newBorder.type = self.BOUNDARY
                    self.cellField.set(pt,newBorder) # to create an extension of that cell
                    self.cellField[pt.x:int(pt.x+xThick),
                                   pt.y:int(pt.y+yThick),
                                   0]=newBorder
                    y= int(center[1]-yi)
                    pt = CompuCell.Point3D(x,y,z)
                    #newBorder = self.potts.createCellG(pt)
                    #newBorder.type = self.BOUNDARY
                    self.cellField.set(pt,newBorder) # to create an extension of that cell
                    self.cellField[pt.x:int(pt.x+xThick),
                                   pt.y:int(pt.y+yThick),
                                   0]=newBorder
                
                
                
                #top/bottom right wall
                x0 = 0
                xf = int(self.borderSideLenght * np.cos(np.pi*30./180.))
                a = (np.sin(np.pi*30./180.) -1)/np.cos(np.pi*30./180.)
                xThick = self.wallThickness*np.cos(np.pi*60./180.)
                yThick = self.wallThickness*np.sin(np.pi*60./180.)
                for xi in range(x0,xf):
                    x = int(center[0] + xi)
                    yi = self.borderSideLenght + a *xi
                    y = int(center[1] + yi)
                    pt = CompuCell.Point3D(x,y,z)
                    #newBorder = self.potts.createCellG(pt)
                    #newBorder.type = self.BOUNDARY
                    self.cellField.set(pt,newBorder) # to create an extension of that cell
                    self.cellField[pt.x:int(pt.x + xThick),
                                   pt.y:int(pt.y + yThick),
                                   pt.z] = newBorder
                    
                    y= int(center[1]-yi)
                    pt = CompuCell.Point3D(x,y,z)
                    #newBorder = self.potts.createCellG(pt)
                    #newBorder.type = self.BOUNDARY
                    self.cellField.set(pt,newBorder) # to create an extension of that cell
                    self.cellField[pt.x:int(pt.x+xThick),
                                   pt.y:int(pt.y+yThick),
                                   0]=newBorder
                
                #righ side wall
                
                x0 = self.borderSideLenght*np.cos(np.pi*330./180.)
                y0 = self.borderSideLenght*np.sin(np.pi*330./180.)
                
                x = int(center[0]+x0)
                y = int(center[1]+y0)
                
                pt = CompuCell.Point3D(x,y,z)
                #newBorder = self.potts.createCellG(pt)
                #newBorder.type = self.BOUNDARY
                self.cellField.set(pt,newBorder)
                self.cellField.set(pt,newBorder)
                self.cellField[pt.x:int(pt.x+self.wallThickness),
                               pt.y:int(pt.y+self.borderSideLenght+self.wallThickness),
                               0]=newBorder
                
                #left side wall
                x = int(center[0]-x0)
                pt = CompuCell.Point3D(x,y,z)
                #newBorder = self.potts.createCellG(pt)
                #newBorder.type = self.BOUNDARY
                self.cellField.set(pt,newBorder)
                self.cellField.set(pt,newBorder)
                self.cellField[int(pt.x-self.wallThickness):pt.x,
                               pt.y:int(pt.y+self.borderSideLenght+self.wallThickness),
                               0]=newBorder
                
                center[0] += int(2*x0 + self.wallThickness)
            center[1] += int(1.5*self.borderSideLenght + self.wallThickness)
            center[0] = int(initialCenter[0] + (x0+0.5*self.wallThickness)*(i%2))
            i+=1
        
    def brickGrid(self):
        c = self.borderSideLenght + self.wallThickness
        initialCenter = (c,c)
        center = list(initialCenter)
        i=1
        while center[1]+0.5*self.borderSideLenght+self.wallThickness <self.dim.y:
            print i
            while center[0]+0.5*self.borderSideLenght+self.wallThickness <self.dim.x:
                x = center[0] 
                y = center[1]
                z = 0
                
                #seed cell at the center
                pt = CompuCell.Point3D(x,y,z)
                newCell = self.potts.createCellG(pt)
                newCell.type = self.CELL
                newCell.targetVolume = self.borderSideLenght**2
                newCell.lambdaVolume = 8
                self.cellField.set(pt,newCell) # to create an extension of that cell
                self.cellField[int(pt.x-0.45*self.borderSideLenght):int(pt.x+0.45*self.borderSideLenght),
                               int(pt.y-0.45*self.borderSideLenght):int(pt.y+0.45*self.borderSideLenght),
                               0]=newCell
                
                #side wall on the left:        
                
                pt = CompuCell.Point3D(int(x- 0.5*self.borderSideLenght),y,z)
                newBorder = self.potts.createCellG(pt)
                newBorder.type = self.BOUNDARY  #Boundary
                self.cellField.set(pt,newBorder) # to create an extension of that cell
                self.cellField[pt.x-self.wallThickness:pt.x,
                               int(pt.y-0.5*(self.borderSideLenght+self.wallThickness)):int(pt.y+0.5*(self.borderSideLenght+self.wallThickness)),
                               0]=newBorder
                
                #side wall on the right
                pt = CompuCell.Point3D(int(x+ 0.5*self.borderSideLenght),y,z)
                newBorder = self.potts.createCellG(pt)
                newBorder.type = self.BOUNDARY  #Boundary
                self.cellField.set(pt,newBorder) # to create an extension of that cell
                self.cellField[pt.x:pt.x+self.wallThickness,
                               int(pt.y-0.5*(self.borderSideLenght+self.wallThickness)):int(pt.y+0.5*(self.borderSideLenght+self.wallThickness)),
                               0]=newBorder
                
                #bottom wall
                pt = CompuCell.Point3D(x,int(y-0.5*self.borderSideLenght),z)
                newBorder = self.potts.createCellG(pt)
                newBorder.type = self.BOUNDARY  #Boundary
                self.cellField.set(pt,newBorder) # to create an extension of that cell
                self.cellField[int(pt.x-0.5*(self.borderSideLenght+self.wallThickness)):int(pt.x+0.5*(self.borderSideLenght+self.wallThickness)),
                               pt.y-self.wallThickness:pt.y,
                               0]=newBorder
                #top wall
                pt = CompuCell.Point3D(x,int(y+0.5*self.borderSideLenght),z)
                newBorder = self.potts.createCellG(pt)
                newBorder.type = self.BOUNDARY  #Boundary
                self.cellField.set(pt,newBorder) # to create an extension of that cell
                self.cellField[int(pt.x-0.5*(self.borderSideLenght+self.wallThickness)):int(pt.x+0.5*(self.borderSideLenght+self.wallThickness)),
                               pt.y:pt.y+self.wallThickness,
                               0]=newBorder
                #update center
                print center
                center[0] += self.borderSideLenght + self.wallThickness
                
            center[1] += self.borderSideLenght + self.wallThickness
            
            center[0] = int(initialCenter[0] + (0.5*self.borderSideLenght)*(i%2))
            i+=1
    
    def step(self,mcs):        
        #type here the code that will run every _frequency MCS
        
        
#         pixelList = self.getCellBoundaryPixelList(cell)
#         for boundaryPixelTrackerData in pixelList:
#             print "pixel of cell id=", cell.id, " type:", cell.type, " = ", boundaryPixelTrackerData.pixel, " number of pixels=", pixelList.numberOfPixels()
        
        for cell in self.cellListByType(self.CELL):
            print cell.id
            for neighbor, commonSurfaceArea in self.getCellNeighborDataList(cell):
                if neighbor:
                    print "neighbor.id", neighbor.id, " commonSurfaceArea=", commonSurfaceArea
                else:
                    print "Medium commonSurfaceArea=", commonSurfaceArea
            
            
        #after relaxTime
        #for each cell iterate the boundary, draw circles around each pixel, find the nearest voxel
        #that is not boundary and not itself. That will be a pixel-pair, each cell will have a dict of 
        #pixel-pairs. Also have cell-pairs.
        #
        #for the transport, 
        # d C1/dt = P_channel * (C1-C2)*(P1-P2)
        # each cell-pair will have a P_channel. The eq above transfers quantity for each pixel-pair
        if mcs == self.relaxTime:
            for cell in self.cellListByType(self.CELL):
                cell.lambdaVolume = 8
        pass
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        