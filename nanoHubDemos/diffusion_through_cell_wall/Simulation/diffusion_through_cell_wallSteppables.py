
from PySteppables import *
import CompuCell
import sys
class diffusion_through_cell_wallSteppable(SteppableBasePy):

    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
        self.wallThickness = 8
        self.borderSideLenght = 128-8
        
    def start(self):
        # any code in the start function runs before MCS=0
        self.brickGrid()
        
        pt = CompuCell.Point3D(55,55,0)
        newSource = self.potts.createCellG(pt)
        newSource.type = self.SOURCE
        self.cellField[pt.x-25:pt.x+2,
                       pt.y-25:pt.y+2,
                       0]=newSource
        
        
        #hex grid:
        
#         c = self.borderSideLenght + self.wallThickness
#         initialCenter = (c,c)
#         center = list(initialCenter)
#         z=0
#         for x in range(int(center[0]-self.borderSideLenght),int(center[0]+self.borderSideLenght)):
#             for y in range(int(center[1]-self.borderSideLenght),int(center[1]+self.borderSideLenght)):
#                 r = np.sqrt(x**2 + y**2)
#                 theta = np.arctan2(x,y)
#                 s = np.sin( theta - (np.pi/3)*np.floor( 1+ 3*theta/np.pi ) )
#                 h = .5*np.sqrt(3)/s
                
#                 if r > h and r <= h + self.wallThickness:
#                     pt = CompuCell.Point3D(x,y,z)
#                     newBorder = self.potts.createCellG(pt)
#                     newBorder.type = self.BOUNDARY
                
                
                
                
        
#         #side wall on the left:
        
#         pt = CompuCell.Point3D(int(x- 0.5*self.borderSideLenght),y,z)
#         newBorder = self.potts.createCellG(pt)
#         newBorder.type = self.BOUNDARY  #Boundary
#         self.cellField.set(pt,newBorder) # to create an extension of that cell
#         self.cellField[pt.x-self.wallThickness:pt.x,
#                        int(pt.y-0.5*self.borderSideLenght):int(pt.y+0.5*self.borderSideLenght),
#                        0]=newBorder
        
#         #side wall on the right
#         pt = CompuCell.Point3D(int(x+ 0.5*self.borderSideLenght),y,z)
#         newBorder = self.potts.createCellG(pt)
#         newBorder.type = self.BOUNDARY  #Boundary
#         self.cellField.set(pt,newBorder) # to create an extension of that cell
#         self.cellField[pt.x:pt.x+self.wallThickness,
#                        int(pt.y-0.5*self.borderSideLenght):int(pt.y+0.5*self.borderSideLenght),
#                        0]=newBorder
        
#         #top left wall
#         wallCenter_x = int(center[0] + self.borderSideLenght * np.cos(np.pi*240./180.))
#         wallCenter_y = int(center[1] + self.borderSideLenght * np.sin(np.pi*240./180.))
        
#         pt = CompuCell.Point3D(wallCenter_x,wallCenter_y,z)
#         newBorder = self.potts.createCellG(pt)
#         newBorder.type = self.BOUNDARY  #Boundary
        
#         z=0
#         hex_x = []
#         hex_y = []
#         for i in range(6):
#             deg = 60 * i - 30
#             theta = deg * np.pi/180
#             x = center[0] + self.borderSideLenght * np.cos(theta)
#             y = center[1] + self.borderSideLenght * np.sin(theta)
#             hex_x.append(x)
#             hex_y.append(y)
        
#         for i,p in enumerate(hex_x):
#             xp = hex_x[i-1]
#             yp = hex_y[i-1]
            
#             xc = hex_x[i]
#             yc = hex_y[i] 
            
            
        
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
        pass
    def finish(self):
        # Finish Function gets called after the last MCS
        pass
        