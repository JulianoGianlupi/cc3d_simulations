#-----------------------------------------------------------------------------------------
#CC3D RECIPES ============================================================================
#  ++  Main Python File  ++  

EXTRA: #==================================================================================
  # to list all global variables defined in python
  try:
    Old
  except NameError:
    pass
  else:
    for key in keys:
      del globals()[key]
  global Old, keys; Old={}; keys=[]
  for key in globals():
    Old[key]=globals()[key]
  #
  global a,b,c,etc...
  #
  # inside def configureSimulation(sim):
  Globals=cc3d.ElementCC3D("Globals")
    for key in sorted(globals(), key=str.lower):
      if key not in Old:
        Globals.ElementCC3D(key,{},globals()[key])
        keys.append(key)

#
#
BASIC: #==================================================================================
  #
  Metadata =cc3d.ElementCC3D("Metadata")
    Metadata.ElementCC3D("VirtualProcessingUnits",{"ThreadsPerVPU":2},2)
    Metadata.ElementCC3D("DebugOutputFrequency",{},0)
  #
  Potts =cc3d.ElementCC3D("Potts")
    Potts.ElementCC3D("Dimensions",{"x":70,"y":150,"z":70})
    Potts.ElementCC3D("Steps",{},10000)
    Potts.ElementCC3D("Temperature",{},40)
    Potts.ElementCC3D("NeighborOrder",{},7)
    Potts.ElementCC3D("Flip2DimRatio",{},1)
    Potts.ElementCC3D("LatticeType",{},"Hexagonal")
    Potts.ElementCC3D("Boundary_x",{},"Periodic")
    Potts.ElementCC3D("Boundary_y",{},"Periodic")
#
#
PLUGINS: #================================================================================
  #
  PlaySet =cc3d.ElementCC3D("Plugin",{"Name":"PlayerSettings"})
    PlaySet.ElementCC3D("Cell", {"Type":"0", "Color":"#000000"}) #black
    PlaySet.ElementCC3D("Cell", {"Type":"1", "Color":"#FF0000"}) #red
    PlaySet.ElementCC3D("Cell", {"Type":"2", "Color":"#00FF00"}) #green
    PlaySet.ElementCC3D("Cell", {"Type":"3", "Color":"#0000FF"}) #blue
    PlaySet.ElementCC3D("Cell", {"Type":"4", "Color":"#CCCCCC"}) #grey
    PlaySet.ElementCC3D("Cell", {"Type":"5", "Color":"#FFFF99"}) #light yellow
    PlaySet.ElementCC3D("TypesInvisibleIn3D",{"Types":"0,3,4"})
    PlaySet.ElementCC3D("VisualControl", {"ScreenshotFrequency":20, "ScreenUpdateFrequency":10})
    PlaySet.ElementCC3D("MainWindow",{"Projection":"2D", "XZProj":20})
    PlaySet.ElementCC3D("NewWindow", {"Projection":"3D", "WindowNumber":1, "CameraClippingRange":"0.13 129.9","CameraDistance":56.7, "CameraViewUp":"0.01 0.99 0.02"})
    PlaySet.ElementCC3D("NewWindow", {"Projection":"3D", "WindowNumber":2, "CameraFocalPoint":"-16 18 61","CameraPos":"-19 17. 66."})
    
    View3D=playSet.ElementCC3D("View3D")
    View3D.ElementCC3D("CameraPosition",{"x":1,"y":1,"z":1})
  #
  VolumeLocalFlex =cc3d.ElementCC3D("Plugin",{"Name":"VolumeLocalFlex"})
  #
  Volume =cc3d.ElementCC3D("Plugin",{"Name":"Volume"})
    Volume.ElementCC3D("VolumeEnergyParameters",{"CellType":"CELL_TYPE", "LambdaVolume":2, "TargetVolume":25})
    Volume.ElementCC3D("TargetVolume",{},25)
    Volume.ElementCC3D("LambdaVolume",{},2)
  #
  SurfaceLocalFlex =cc3d.ElementCC3D("Plugin",{"Name":"SurfaceLocalFlex"})
  #
  Surface =cc3d.ElementCC3D("Plugin",{"Name":"Surface"})
    Surface.ElementCC3D("SurfaceEnergyParameters",{"CellType":"CELL_TYPE", "LambdaSurface":2, "TargetSurface":25})
    Surface.ElementCC3D("TargetSurface",{},25)
    Surface.ElementCC3D("LambdaSurface",{},2)
  #
  ClusterSurface =cc3d.ElementCC3D("Plugin",{"Name":"ClusterSurface"})
    ClusterSurface.ElementCC3D("TargetClusterSurface":80})
    ClusterSurface.ElementCC3D("LambdaClusterSurface":2.0})
  #
  CenterOfMass =cc3d.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
  #
  NeighborTracker =cc3d.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
  #
  MomentOfInertia =cc3d.ElementCC3D("Plugin",{"Name":"MomentOfInertia"})
  #
  PixelTracker =cc3d.ElementCC3D("Plugin",{"Name":"PixelTracker"})
  #
  BoundaryPixelTracker =cc3d.ElementCC3D("Plugin",{"Name":"BoundaryPixelTracker"})
    BoundaryPixelTracker.ElementCC3D("NeighborOrder",{},1)
  #
  LengthConstraint =cc3d.ElementCC3D("Plugin",{"Name":"LengthConstraint"})
    LengthConstraint.ElementCC3D("LengthEnergyParameters",{"CellType":"CELL_TYPE","LambdaLength":2.0,"TargetLength":25, ,"MinorTargetLength":5})
  #
  ExternalPotentialLocalFlex =cc3d.ElementCC3D("Plugin",{"Name":"ExternalPotentialLocalFlex"})
    ExternalPotentialLocalFlex.ElementCC3D("Algorithm",{},"PixelBased") #CenterOfMassBased
  #
  ExternalPotential =cc3d.ElementCC3D("Plugin",{"Name":"ExternalPotential"})
    ExternalPotential.ElementCC3D("Algorithm",{},"PixelBased") #CenterOfMassBased
    ExternalPotential.ElementCC3D("ExternalPotentialParameters",{"CellType":"CELL_TYPE", "x":"-0.5", "y":"0.0", "z":"0.0"})
  #
  Connectivity =cc3d.ElementCC3D("Plugin",{"Name":"Connectivity"})
    Connectivity.ElementCC3D("Penalty",{},"10000000")
  #
  ConnectivityGlobal =cc3d.ElementCC3D("Plugin",{"Name":"ConnectivityGlobal"})
    ConnectivityGlobal.ElementCC3D("DoNotPrecheckConnectivity")
    ConnectivityGlobal.ElementCC3D("Penalty",{"Type":"CELL_TYPE"},1000000)
  #
  FocalPointPlasticity =cc3d.ElementCC3D("Plugin",{"Name":"FocalPointPlasticity"})
    FocalPointPlasticity.ElementCC3D("Local")
    Parameters=FocalPointPlasticity.ElementCC3D("Parameters",{"Type1":"CELL_TYPE","Type2":"CELL_TYPE"})
      Parameters.ElementCC3D("Lambda",{},10)
      Parameters.ElementCC3D("ActivationEnergy",{},-50)
      Parameters.ElementCC3D("TargetDistance",{},7)
      Parameters.ElementCC3D("MaxDistance",{},20)
      Parameters.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":1},1)
    ParametersIn=FocalPointPlasticity.ElementCC3D("InternalParameters",{"Type1":"CELL_TYPE","Type2":"CELL_TYPE"})
      ParametersIn.ElementCC3D("Lambda",{},10)
      ParametersIn.ElementCC3D("ActivationEnergy",{},-50)
      ParametersIn.ElementCC3D("TargetDistance",{},7)
      ParametersIn.ElementCC3D("MaxDistance",{},20)
      ParametersIn.ElementCC3D("MaxNumberOfJunctions",{"NeighborOrder":1},1)
    FocalPointPlasticity.ElementCC3D("NeighborOrder",{},1)
  #
  Secretion =cc3d.ElementCC3D("Plugin",{"Name":"Secretion"})
    Field=Secretion.ElementCC3D("Field",{"Name":"FIELD"})
      Field.ElementCC3D("Secretion",{"Type":"CELL_TYPES"},1.0)
  #
  Chemotaxis =cc3d.ElementCC3D("Plugin",{"Name":"Chemotaxis"})
    Field=Chemotaxis.ElementCC3D("ChemicalField",{"Name":"FIELD", "Source":"PDE_SOLVER"})
      Field.ElementCC3D("ChemotaxisByType",{"Type":"CELL_TYPE", "ChemotactTowards":"CELL_TYPES", "Lambda":1.0})
      Field.ElementCC3D("ChemotaxisByType",{"Type":"CELL_TYPE", "ChemotactTowards":"CELL_TYPES", "Lambda":1.0, "SaturationCoef":100.0})
      Field.ElementCC3D("ChemotaxisByType",{"Type":"CELL_TYPE", "ChemotactTowards":"CELL_TYPES", "Lambda":1.0, "SaturationLinearCoef":10.1})
#
#
STEPPABLES: #=============================================================================
  #
  PIFInitializer =cc3d.ElementCC3D("Steppable",{"Type":"PIFInitializer"})
    PIFInitializer.ElementCC3D("PIFName",{},"PLEASE_PUT_PROPER_FILE_NAME_HERE")
  #
  BlobInitializer =cc3d.ElementCC3D("Steppable",{"Type":"BlobInitializer"})
    Region=BlobInitializer.ElementCC3D("Region")
      Region.ElementCC3D("Center",{"x":50,"y":50,"z":0})
      Region.ElementCC3D("Radius",{},20)
      Region.ElementCC3D("Gap",{},0)
      Region.ElementCC3D("Width",{},5)
      Region.ElementCC3D("Types",{},"CELL_TYPES")
  #
  UniformInitializer =cc3d.ElementCC3D("Steppable",{"Type":"UniformInitializer"})
    Region=UniformInitializer.ElementCC3D("Region")
      Region.ElementCC3D("BoxMin", {"x":1,  "y":1,  "z":0})
      Region.ElementCC3D("BoxMax", {"x":100, "y":100, "z":1})
      Region.ElementCC3D("Gap", {},1)
      Region.ElementCC3D("Width",{},6)
      Region.ElementCC3D("Types",{},"CELL_TYPES")
  #
  DiffusionSolverFE =cc3d.ElementCC3D("Steppable",{"Type":"DiffusionSolverFE"})
    Field=DiffusionSolverFE.ElementCC3D("DiffusionField")
      Data=Field.ElementCC3D("DiffusionData")
        Data.ElementCC3D("FieldName",{},"FIELD")
        Data.ElementCC3D("GlobalDiffusionConstant",{},0.1)
        Data.ElementCC3D("GlobalDecayConstant",{},1e-05)
        Data.ElementCC3D("InitialConcentrationExpression",{},"x*y")
        Data.ElementCC3D("ConcentrationFileName",{},"NAME_OF_THE_FILE.txt")
        Data.ElementCC3D("DiffusionCoefficient",{"CellType":"CELL_TYPE"},0.1)
        Data.ElementCC3D("DiffusionCoefficient",{"CellType":"CELL_TYPE"},0.1)
        Data.ElementCC3D("DecayCoefficient",{"CellType":"CELL_TYPE"},0.0001)
        Data.ElementCC3D("DecayCoefficient",{"CellType":"CELL_TYPE"},0.0001)
      Secretion=Field.ElementCC3D("SecretionData")
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.1)
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.1)
        Secretion.ElementCC3D("SecretionOnContact",{"SecreteOnContactWith":"CELL_TYPES","Type":"CELL_TYPE"},0.2)
        Secretion.ElementCC3D("ConstantConcentration",{"Type":"CELL_TYPE"},0.1)
      BoundaryConditions=Field.ElementCC3D("BoundaryConditions")
        PlaneX=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"X"})
          PlaneX.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":10.0})
          PlaneX.ElementCC3D("ConstantValue",{"PlanePosition":"Max","Value":5.0})
          PlaneX.ElementCC3D("Periodic")
          PlaneX.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":10.0})
        PlaneY=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"Y"})
          PlaneY.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":10.0})
          PlaneY.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":5.0})
          PlaneY.ElementCC3D("Periodic")
          PlaneY.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":10.0})
  #
  FlexibleDiffusionSolverFE =cc3d.ElementCC3D("Steppable",{"Type":"FlexibleDiffusionSolverFE"})
    Field=FlexibleDiffusionSolverFE.ElementCC3D("DiffusionField")
      Data=Field.ElementCC3D("DiffusionData")
        Data.ElementCC3D("FieldName",{},"FIELD")
        Data.ElementCC3D("DiffusionConstant",{},0.1)
        Data.ElementCC3D("DecayConstant",{},1e-05)
        Data.ElementCC3D("DoNotDiffuseTo",{},"CELL_TYPES")
        Data.ElementCC3D("DoNotDecayIn",{},"CELL_TYPES")
        Data.ElementCC3D("InitialConcentrationExpression",{},"x*y")
        Data.ElementCC3D("ConcentrationFileName",{},"NAME_OF_THE_FILE.txt")
        Data.ElementCC3D("ExtraTimesPerMCS",{},0)
        Data.ElementCC3D("DeltaX",{},1.0)
        Data.ElementCC3D("DeltaT",{},1.0)
    BoundaryConditions=Field.ElementCC3D("BoundaryConditions")
        PlaneX=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"X"})
          PlaneX.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":10.0})
          PlaneX.ElementCC3D("ConstantValue",{"PlanePosition":"Max","Value":5.0})
          PlaneX.ElementCC3D("Periodic")
          PlaneX.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":10.0})
        PlaneY=BoundaryConditionsElmnt.ElementCC3D("Plane",{"Axis":"Y"})
          PlaneY.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":10.0})
          PlaneY.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":5.0})
          PlaneY.ElementCC3D("Periodic")
          PlaneY.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":10.0})
  #
  SteadyStateDiffusionSolver2D =cc3d.ElementCC3D("Steppable",{"Type":"SteadyStateDiffusionSolver2D"})
    Field=SteadyStateDiffusionSolver2D.ElementCC3D("DiffusionField")
      Data=Field.ElementCC3D("DiffusionData")
        Data.ElementCC3D("FieldName",{},"FIELD")
        Data.ElementCC3D("DiffusionConstant",{},"1.0")
        Data.ElementCC3D("DecayConstant",{},"1e-05")
        Data.ElementCC3D("InitialConcentrationExpression",{},"x*y")
      Secretion=Field.ElementCC3D("SecretionData")
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.1)
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.2)
      BoundaryConditions=Field.ElementCC3D("BoundaryConditions")
        PlaneX=BoundaryConditions.ElementCC3D("Plane",{"Axis":"X"})
          PlaneX.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":10.0})
          PlaneX.ElementCC3D("ConstantValue",{"PlanePosition":"Max","Value":5.0})
          PlaneX.ElementCC3D("Periodic")
          PlaneX.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":10.0})
        PlaneY=BoundaryConditions.ElementCC3D("Plane",{"Axis":"Y"})
          PlaneY.ElementCC3D("ConstantDerivative",{"PlanePosition":"Min","Value":10.0})
          PlaneY.ElementCC3D("ConstantDerivative",{"PlanePosition":"Max","Value":5.0})
          PlaneY.ElementCC3D("Periodic")
          PlaneY.ElementCC3D("ConstantValue",{"PlanePosition":"Min","Value":10.0})
  #
  ReactionDiffusionSolverFE=cc3d.ElementCC3D("Steppable",{"Type":"ReactionDiffusionSolverFE"})
    Field=ReactionDiffusionSolverFE.ElementCC3D("DiffusionField")
      Data=Field.ElementCC3D("DiffusionData")
        Data.ElementCC3D("FieldName",{},"FIELD")
        Data.ElementCC3D("DiffusionConstant",{},0.0011)
        Data.ElementCC3D("AdditionalTerm",{},"FORMULA")
        Data.ElementCC3D("DecayConstant",{},0.1)
        Data.ElementCC3D("DoNotDecayIn",{},"CELL_TYPES")
      Secretion=Field.ElementCC3D("SecretionData")
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.1)
    Field2=ReactionDiffusionSolverFE.ElementCC3D("DiffusionField")
      Data=Field2.ElementCC3D("DiffusionData")
        Data.ElementCC3D("FieldName",{},"FIELD2")
        Data.ElementCC3D("DiffusionConstant",{},0.1)
        Data.ElementCC3D("AdditionalTerm",{},"FORMULA")
        Data.ElementCC3D("DecayConstant",{},1.0)
        Data.ElementCC3D("DoNotDecayIn",{},"CELL_TYPES")
      Secretion=Field2.ElementCC3D("SecretionData")
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.1)
    Field3=ReactionDiffusionSolverFE.ElementCC3D("DiffusionField")
      Data=Field3.ElementCC3D("DiffusionData")
        Data.ElementCC3D("FieldName",{},"FIELD3")
        Data.ElementCC3D("DiffusionConstant",{},0.1)
        Data.ElementCC3D("AdditionalTerm",{},"FORMULA")
      Secretion=Field3.ElementCC3D("SecretionData")
        Secretion.ElementCC3D("Secretion",{"Type":"CELL_TYPE"},0.1)
  #
  Multiple Calls:
    PDEcaller =cc3d.ElementCC3D("Steppable",{"Type":"PDESolverCaller"})
    PDEcaller.ElementCC3D("CallPDE",{"PDESolverName":"FlexibleDiffusionSolverFE", "ExtraTimesPerMC":9})

  #
#
#
PYTHON MAIN FILE TEMPLATE: #==============================================================
  
  import sys,time
  from os import environ
  from os import getcwd
  import string
  sys.path.append(environ["PYTHON_MODULE_PATH"])
  sys.path.append(environ["SWIG_LIB_INSTALL_DIR"])
  
  try:
    Old
  except NameError:
    pass
  else:
    for key in keys:
      del globals()[key]
  global Old, keys; Old={}; keys=[]
  for key in globals():
    Old[key]=globals()[key]
  
  global A, B, C
  
  def configureSimulation(sim):
    import CompuCellSetup
    from XMLUtils import ElementCC3D
    cc3d=ElementCC3D("CompuCell3D")
    
    Globals=cc3d.ElementCC3D("Globals")
    for key in sorted(globals(), key=str.lower):
      if key not in Old:
        Globals.ElementCC3D(key,{},globals()[key])
        keys.append(key)
    
    md=cc3d.ElementCC3D("Metadata")
    md.ElementCC3D("VirtualProcessingUnits",{"ThreadsPerVPU":2},2)
    md.ElementCC3D("DebugOutputFrequency",{},0)
    
    potts=cc3d.ElementCC3D("Potts")
    potts.ElementCC3D("Dimensions",{"x":Lx,"y":Ly,"z":Lz})
    potts.ElementCC3D("Steps",{},Time)
    potts.ElementCC3D("Temperature",{},int(T))
    potts.ElementCC3D("NeighborOrder",{},2)
    
    #CELL TYPES:
    cellType=cc3d.ElementCC3D("Plugin",{"Name":"CellType"})
    cellType.ElementCC3D("CellType", {"TypeName":"Medium","TypeId":"0"})
    cellType.ElementCC3D("CellType", {"TypeName":"cyto",  "TypeId":"1"})
    
    #CELL COLORS:
    playSet=cc3d.ElementCC3D("Plugin",{"Name":"PlayerSettings"})
    playSet.ElementCC3D("Cell", {"Type":"0", "Color":"#000000"}) #black
    playSet.ElementCC3D("Cell", {"Type":"1", "Color":"#FFFFFF"}) #white
     
    #CONTACT ENERGIES:
    contact=cc3d.ElementCC3D("Plugin",{"Name":"Contact"})
    contact.ElementCC3D("Energy", {"Type1":"Medium", "Type2":"Medium"},0)
    contact.ElementCC3D("Energy", {"Type1":"Medium", "Type2":"cyto"  },10)
    #
    contact.ElementCC3D("Energy", {"Type1":"cyto", "Type2":"cyto"},10)
    #-neighbor order
    contact.ElementCC3D("NeighborOrder",{},nOrder)
    
    CompuCellSetup.setSimulationXMLDescription(cc3d)
  
  import CompuCellSetup
  sim,simthread = CompuCellSetup.getCoreSimulationObjects()
  configureSimulation(sim)
  
  import CompuCell
  # Create extra player fields here or add attributes
  pyAttributeAdder,dictAdder=CompuCellSetup.attachDictionaryToCells(sim)
  CompuCellSetup.initializeSimulationObjects(sim,simthread)
  
  #Add Python steppables here
  steppableRegistry=CompuCellSetup.getSteppableRegistry()
  
  from Project_Step import Steppable
  steppable=Steppable(_simulator=sim,_frequency=1)
  steppableRegistry.registerSteppable(steppable)
  
  CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
  ##sys.exit()
#
#
PYTHON STEPPABLE FILE TEMPLATE: #=========================================================
  import sys
  import os
  from PySteppables import *
  from PySteppablesExamples import MitosisSteppableBase
  from PySteppablesExamples import MitosisSteppableClustersBase
  import CompuCell
  import CompuCellSetup
  from PlayerPython import *
  from math import *
  from random import *
  from copy import deepcopy
  import time
  
  class Steppable(SteppableBasePy): #or MitosisSteppableBase or MitosisSteppableClustersBase
    def __init__(self,_simulator,_frequency):
      SteppableBasePy.__init__(self,_simulator,_frequency)
      
    def start(self):
      pass
    
    def step(self,mcs):
      pass
      
    def finish(self):
      pass
#
#
PYTHON STEPPABLES: #=======================================================================
  #
  File Management:
    import os
    import inspect
     
    sourceFile = os.path.abspath(__file__)
    sourceFile = inspect.getfile(inspect.currentframe())
     
    SourceDir = os.path.dirname(os.path.abspath(__file__))
    SourceDir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    
    screenDir = CompuCellSetup.getScreenshotDirectoryName()
    
    source=os.path.abspath( __file__ )
    if (source[-1]=="c"): source=source[:-1]
    os.system("copy "+source[:source.rfind("_")]+".py "+screenDir)
    os.system("copy "+source+" "+screenDir)
  #
  Loading/Saving files:
    File=open(FileName,'w')#=write  | 'a')=append  | 'r')=read  | 'r+')=read+write
    
    File.write("%s%d%f\n" % ("string", integer, real) )
    
    File.read() #read entire file
    File.readline() #read one line  | File.readline() again reads next line
    for line in File:
      #line = current line in file
      
    File.close()
  #
  Saving PIFFs:
      import CompuCellSetup
      Dir=CompuCellSetup.getScreenshotDirectoryName()
      FileName=Dir+"/PiffFile_"+str(mcs)+".piff"
      print mcs, " -==// Saving simulation piff \\==- ", FileName
      File=open(FileName,'w')
      
      for cell in self.cellList:
        name=self.typeIdTypeNameDict[cell.type]
        id=str(cell.id)
        pixelList=self.getCellPixelList(cell)
        for pixel in pixelList:
          x=pixel.pixel.x; y=pixel.pixel.y; z=pixel.pixel.z
          File.write("%s\n" % (id+" "+name+(" "+str(x))*2+(" "+str(y))*2+(" "+str(z))*2) )
      File.close()
      
      File.write("%s\n" % ("Include Clusters") )
      for cell in self.cellList:
        name=self.typeIdTypeNameDict[cell.type]
        id=str(cell.id)
        cId=srt(cell.clusterId)
        pixelList=self.getCellPixelList(cell)
        for pixel in pixelList:
          x=pixel.pixel.x; y=pixel.pixel.y; z=pixel.pixel.z
          File.write("%s\n" % (cId+" "+id+" "+name+(" "+str(x))*2+(" "+str(y))*2+(" "+str(z))*2) )
      File.close()
  #
  Temperature Manipulation:
    pottsXMLData=self.simulator.getCC3DModuleData("Potts")
    temperatureElement=pottsXMLData.getFirstElement("Temperature")
    currentT=float(temperatureElement.getText())
    temperatureElement.updateElementValue(str(newT))
    self.simulator.updateCC3DModule(pottsXMLData)
  #
  Contact Energy Manipulation:
    from XMLUtils import dictionaryToMapStrStr as d2mss
    
    contactXMLData=self.simulator.getCC3DModuleData("Plugin","Contact")
    Type1name=self.typeIdTypeNameDict[cell.type]
    cellNeighborList=self.getCellNeighbors(cell) # generates list of neighbors of cell 'cell'
    for neighbor in cellNeighborList:
      if neighbor.neighborAddress:
        Type2name=self.typeIdTypeNameDict[neighbor.neighborAddress.type]
        contactEnergy=contactXMLData.getFirstElement("Energy",d2mss({"Type1":Type1name,"Type2":Type2name}))
        if not contactEnergy:
          contactEnergy=contactXMLData.getFirstElement("Energy",d2mss({"Type1":Type2name,"Type2":Type1name}))
        EnergyValue=float(contactEnergy.getText())
        contactEnergy.updateElementValue(str(NEW_VALUE))
  #
  Mitosis:
    from PySteppablesExamples import MitosisSteppableBase
    
    cells2div=[]
    for cell in self.cellList:
      if cell.volume>50:
        cells2div.append(cell)
    for cell in cells2div:
      self.divideCellRandomOrientation(cell)
      self.divideCellOrientationVectorBased(cell,1,0,0)
      self.divideCellAlongMajorAxis(cell) #  |#---#|
      self.divideCellAlongMinorAxis(cell) #  |##|##|
      
    def updateAttributes(self):
      parentCell=self.mitosisSteppable.parentCell
      childCell=self.mitosisSteppable.childCell
      childCell.type=parentCell.type
      ...
      bionetAPI.copyBionetworkFromParent(parentCell,childCell)
      ...
      parentCellDict=CompuCell.getPyAttrib(parentCell)
      childCellDict=CompuCell.getPyAttrib(childCell)
      for key in parentCellDict:
        if (key!="Bionetwork"): 
          childCellDict[key] = deepcopy(parentCellDict[key])
  #
  Mitosis Cluster:
    from PySteppablesExamples import MitosisSteppableClustersBase
    
    cells2div=[]
    for cell in self.cellListByType(self.TYPE):
      if (cell.volume>1000):
        cells2div.append(cell)
    for cell in cells2div:
      self.divideClusterRandomOrientation(cell.clusterId)
      self.divideClusterOrientationVectorBased(cell.clusterId,1,0,0)
      self.divideClusterAlongMajorAxis(cell.clusterId)
      self.divideClusterAlongMinorAxis(cell.clusterId)
    
    def updateAttributes(self):
      childCell = self.mitosisSteppable.childCell
      parentCell = self.mitosisSteppable.parentCell
      ...
      compListChild=self.inventory.getClusterCells(childCell.clusterId)
      compListParent=self.inventory.getClusterCells(parentCell.clusterId)
      for cell in compListChild:
        ...
      for cell in compListParent:
        ...
      for cell in compListChild+compListParent:
        ...
  #
  Visiting Neighbor Pixels:
    self.boundaryStrategy=CompuCell.BoundaryStrategy.getInstance()
    self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromDepth(D)
    self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromNeighborOrder(N)
    #  D 2D/3D  |   N  | #Pixels 2D/3D
    #  1  /  1  |   1  |    4 / 6
    # 1.5 / 1.5 |   2  |    8 / 18
    #  2  / 1.8 |   3  |   12 / 26
    # 2.3 /  2  |   4  |   20 / 32
    # 2.9 / 2.3 |   5  |   24 / 56
    #  3  / 2.5 |   6  |   28 / 80
    # 3.2 / 2.9 |   7  |   36 / 92
    # 3.7 /  3  |   8  |   44 / 122
    #     / 3.2 |   9  |   48 / 146
    #     / 3.5 |  10  |      / 170
    #     / 3.7 |  11  |      / 178
    #     / 3.8 |  12  |      / 202
    #     /     |  13  |      / 250
    
    pt=CompuCell.Point3D(x,y,z)
    for i in xrange(self.maxNeighborIndex+1): 
      pN=self.boundaryStrategy.getNeighborDirect(pt,i) #pt = original pixel
      cell2=self.cellField.get(pN.pt) #pN.pt = neighbor pixel
  #
  FPP Links Manipulation:
    for fpp in self.getFocalPointPlasticityDataList(cell):
      cell2=fpp.neighborAddress
      targetDistance=fpp.targetDistance
      lambdaDistance=fpp.lambdaDistance
      self.focalPointPlasticityPlugin.createFocalPointPlasticityLink(cell,cell2,lambda,targetDistance,maxDistance)
      self.focalPointPlasticityPlugin.setFocalPointPlasticityParameters(cell,cell2,lambda,targetDistance,maxDistance)
      self.focalPointPlasticityPlugin.deleteFocalPointPlasticityLink(cell,cell2)
    for fpp in self.getInternalFocalPointPlasticityDataList(cell):
      cell2=fpp.neighborAddress
      targetDistance=fpp.targetDistance
      lambdaDistance=fpp.lambdaDistance
      self.focalPointPlasticityPlugin.createInternalFocalPointPlasticityLink(cell,cell2,lambda,targetDistance,maxDistance)
      self.focalPointPlasticityPlugin.setInternalFocalPointPlasticityParameters(cell,cell2,lambda,targetDistance,maxDistance)
      self.focalPointPlasticityPlugin.deleteInternalFocalPointPlasticityLink(cell,cell2)
    Anchor:
      anchorId=self.fppPlugin.createAnchor(cell,lambda,targetDistance,maxDistance,x,y,z)
      self.fppPlugin.deleteAnchor(cell,anchorId)
      self.fppPlugin.setAnchorParameters(cell,anchorId,lambda,targetDistance,maxDistance,x,y,z)
  #
  Ploting Histograms:
    self.pH=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
    #Plot Title - properties
    self.pH.setTitle("TITLE")
    # properties of x,y axes
    self.pH.setXAxisTitle("X_TITLE")
    self.pH.setYAxisTitle("Y_TITLE")
    # Add histogram plot
    self.pH.addHistPlot("NAME",_r=255,_g=0,_b=0,_alpha=160) #alpha: 0=transparent, 255=opaque 
    self.pH.addAutoLegend("top")
    self.pH.addGrid()
    
    #create histogram
    from numpy import *
    L=[list of values]
    (Hist,Bin)=histogram(L,bins=10)
    Hist,Bin = self.HistList(L)
    D={keys:values}
    Hist,Bin=self.HistBins(D,nBins)
    Hist=[list of values] #y axis
    Bin =[list of bins]   #x axis   len(Bin)=len(Hist)+1
    
    self.pH.eraseAllData()
    self.pH.setHistogramColor("NAME",_r=255,_g=0,_b=255,_alpha=160)
    self.pH.addHistPlotData("NAME",Hist,Bin)
    self.pH.showAllHistPlots()
    
    def HistList(self,L): #histogram from list with fixed bin size of 1 unit
      Hist=[0]*max(L); Bin=[]
      for i in range(max(L)):
        Hist[i]+=L.count(i+1)
        Bin.append(i+1)
      Bin.append(max(L)+1)
      return Hist, Bin
    
    def HistBins(self,L,nBins): #Histogram from dictionary
      Bin=[]; Hist=[0]*nBins
      maxK=max(L.keys()); minK=min(L.keys())
      k=(maxK-minK)/(nBins-1.)
      b=(nBins-1.)/(maxK-minK)
      for i in range(nBins+1):
        Bin.append(minK + i*k)
      for key,val in L.items():
        bin=int((key-minK)*b)
        Hist[bin]+=val
      return Hist,Bin
  #
  Plotting Bars:
    self.pW.setBarPlotView()
    self.pW.setTitle("TITLE")
    self.pW.setXAxisTitle("X_TITLE")
    self.pW.setYAxisTitle("Y_TITLE")
    
    Ly=[list of y values]    # height of each bar
    Lx=[list of x locations] # location of each bar  |  len(Ly)=len(Lx)
    self.pW.addBarPlotData(Ly,Lx,width) # width = width of the bar
    self.pW.showAllBarCurvePlots()
  #
  Ploting Graphs:
    #CREATING GRAPH (full)
    self.pW=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
    # Plot Title - properties
    self.pW.setTitle("Average Volume And Surface")
    self.pW.setTitleSize(12)
    self.pW.setTitleColor("Green")
    # plot background
    self.pW.setPlotBackgroundColor("orange")
    # properties of x axis
    self.pW.setXAxisTitle("MonteCarlo Step (MCS)")
    self.pW.setXAxisLogScale()
    self.pW.setXAxisTitleSize(10)
    self.pW.setXAxisTitleColor("blue")
    self.pW.setXAxisScale(low,high)
    # properties of y axis
    self.pW.setYAxisTitle("Variables")
    self.pW.setYAxisLogScale()
    self.pW.setYAxisTitleSize(10)
    self.pW.setYAxisTitleColor("red")
    self.pW.setYAxisScale(low,high)
    # add plot
    self.pW.addPlot("MVol",_style='Dots') #NoCurve,Lines,Sticks,Steps,Dots
    self.pW.changePlotProperty("MVol","LineWidth",5)
    self.pW.changePlotProperty("MVol","LineColor","red")
    # add plot
    self.pW.addPlot("MSur",_style='Steps')
    self.pW.changePlotProperty("MSur","LineWidth",1)
    self.pW.changePlotProperty("MSur","LineColor","green")
    # extra
    self.pW.addGrid()
    self.pW.addAutoLegend("top")
    
    #CREATING THE GRAPH (short)
    self.pW=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
    self.pW.setTitle("TITLE")
    self.pW.setXAxisTitle("X_AXIS")
    self.pW.setXAxisTitle("Y_AXIS")
    self.pW.addPlot("DATA",_style='Dots') #NoCurve,Lines,Sticks,Steps,Dots
    self.pW.changePlotProperty("DATA","LineWidth",5)
    self.pW.changePlotProperty("DATA","LineColor","red")
    self.pW.addGrid()
    self.pW.addAutoLegend("top")
    
    #UPDATING GRAPH
    self.pW.setXAxisScale(min,max)
    self.pW.setYAxisScale(min,max)
    
    self.pW.eraseAllData()
    self.pW.eraseData("MVol")
    self.pW.eraseData("MSur")
    
    self.pW.addDataPoint("MVol",mcs,volume)
    self.pW.addDataPoint("MSur",mcs,surface)
    
    self.pW.showAllPlots()
    self.pW.showPlot("MVol")
    self.pW.showPlot("MSur")
  #
  Saving Plots:
    qwtPlotWidget=self.pW.getQWTPLotWidget()
    Size=qwtPlotWidget.size()
    self.pW.savePlotAsPNG(fileName,Size.width(),Size.height()) #screen
    self.pW.savePlotAsPNG(fileName) #default 400x400
    self.pW.savePlotAsPNG(fileName,1000,500) #custom
  #
  Moments Of Inertia (cell):
    cell.iXX
    cell.iYY
    cell.iZZ
    cell.iXY; cell.iXZ; cell.iYZ
    
    cell.ecc
    
    axes=self.momentOfInertiaPlugin.getSemiaxes(cell)
   "minorAxis =",axes[0]
   "majorAxis =",axes[2]
   "medianAxis=",axes[1]
  #
  Calculating Tissue Moments of Inertia 2D:
    from numpy import *
    Ixx=0; Iyy=0; Ixy=0
    vol=0; x=0; y=0
    for cell in self.cellList:
      vol+=cell.volume
      x+=cell.xCM; y+=cell.yCM
    xCM=float(x)/vol; yCM=float(y)/vol
    for cell in self.cellList:
      pixelList=CellPixelList(self.pixelTrackerPlugin,cell)
      for pixel in pixelList:
        x=(pixel.pixel.x-xCM); y=(pixel.pixel.y-yCM)
        Ixx+=y*y; Iyy+=x*x; Ixy-=x*y
    A=array( [(Ixx, Ixy),(Ixy, Iyy)] )
    B=linalg.eig(A)
    minorAxis=sorted(B[0])[0]
    majorAxis=sorted(B[0])[1]
  #
  Calculating Tissue Moments of Inertia 3D:
    Ixx=0; Iyy=0; Izz=0
    Ixy=0; Ixz=0; Iyz=0
    vol=0; x=0; y=0; z=0
    for cell in self.cellList:
      vol+=cell.volume
      x+=cell.xCM; y+=cell.yCM; z+=cell.zCM
    xCM=float(x)/vol; yCM=float(y)/vol; zCM=float(z)/vol
    for cell in self.cellList:
      pixelList=CellPixelList(self.pixelTrackerPlugin,cell)
      for pixelData in pixelList:
        pt=pixelData.pixel
        x=(pt.x-xCM); y=(pt.y-yCM); z=(pt.z-zCM)
        x2=x*x; y2=y*y; z2=z*z
        Ixx+=y2+z2; Iyy+=x2+z2; Izz+=x2+y2
        Ixy-=x*y; Ixz-=x*z; Iyz-=y*z
    A=array( [(Ixx, Ixy, Ixz), (Ixy, Iyy, Iyz), (Ixz, Iyz, Izz)] )
    B=linalg.eig(A)
    minorAxis=sorted(B[0])[0]
    medianAxis=sorted(B[0])[1]
    majorAxis=sorted(B[0])[2]
  #
  Getting Cell Type Info:
    #Cell type-id and type-name:
    self.typeIdTypeNameDict # {0:"Medium", 1:"CELL_TYPE_1", ...}
    self.typeIdTypeNameDict[cell.type] # current cell's type name
    len(self.typeIdTypeNameDict) # number of cell types + medium
    
    #Cell from id:
    cell=self.inventory.attemptFetchingCellById(id)
    
    #Distance between cells:
    D = self.distanceBetweenCells(cell,cell2)       # D = sqrt(dx**2+dy**2+dz**2)
    V = self.distanceVectorBetweenCells(cell,cell2) # [dx,dy,dz]
  #
  Access/Modify Cell Lattice:
    pt=CompuCell.Point3D(x,y,z) # defines a lattice vector
    
    medium=CompuCell.getMediumCell() #get Medium cell
    
    cell=self.cellField.get(pt) # get cell that is on pt [self.cellField.get(pt.x,pt.y,pt.z)]
    cell=self.cellField[x,y,0]
    
    # to create a brand new cell
    newcell=self.potts.createCellG(pt)
    newcell.type=1 # don’t forget to assign a type to the new cell
    cell=self.potts.createCell()
    
    self.cellField.set(pt,cell) # to create an extension of that cell
    self.cellField[x:x+4,y:y+4,0]=cell
  #
  SBML solver:
    #Loading/adding SBML
    self.addSBMLToCellTypes(_modelFile='',_modelName='',_types=[],_stepSize=1.0,_initialConditions={})
    
    self.addSBMLToCell(_modelFile='',_modelName='',_cell=None,_stepSize=1.0,_initialConditions={},_coreModelName='',_modelPathNormalized='')
    
    self.addSBMLToCellIds(_modelFile='',_modelName='',_ids=[],_stepSize=1.0,_initialConditions={})
    
    self.addFreeFloatingSBML(_modelFile='',_modelName='',_stepSize=1.0,_initialConditions={})
    
    #Removing/deleting SBML
    self.deleteSBMLFromCellTypes(_modelName='',_types=[])
    
    self.deleteSBMLFromCell(_modelName='',_cell=None)
    
    self.deleteSBMLFromCellIds(_modelName='',_ids=[])
    
    self.deleteFreeFloatingSBML(_modelName='')
    
    #Check if SBML exists in a cell
    self.getSBMLSimulator(_modelName='',_cell=None) #supress _cell entry to check free floating SBML
    
    #TimeStep SBML
    self.timestepSBML()  #time step all models
    
    self.timestepCellSBML() #time step all models associated with cells
    
    self.timestepFreeFloatingSBML() #time step all free floating models
    
    #Change time steps
    self.setStepSizeForCellTypes( _modelName='',_types=[],_stepSize=1.0)
    
    self.setStepSizeForCell(_modelName='',_cell=None,_stepSize=1.0)
    
    self.setStepSizeForCellIds(_modelName='',_ids=[],_stepSize=1.0)
    
    self.setStepSizeForFreeFloatingSBML(_modelName='',_stepSize=1.0)
    
    #Get concentrations/parameters values
    state = self.getSBMLState(_modelName='',_cell=None) #returns dictionary with all values
    
    value = self.getSBMLValue(_modelName='',_valueName='',_cell=None) #return value of specific parameter/concentrarion
    
    #Set concentrations/parameters values
    self.setSBMLState(_modelName='',_cell=None,_state={}) #modify all parameters/concentrarions
    
    self.setSBMLValue(_modelName='',_valueName='',_value=0.0,_cell=None) #modify only 1parameter/concentrarion
    
    #Copy SBML from one cell to another:
    self.copySBMLs(_fromCell,_toCell,_sbmlNames=[])
  #
  Secretion:
    Secretor=self.getFieldSecretor("FIELDNAME") # you may reuse secretor for many cells. Simply define it outside the loop
    
    Secretor.secreteInsideCell(cell,SecretionConstant)
    Secretor.secreteInsideCellAtBoundary(cell,SecretionConstant)
    Secretor.secreteInsideCellAtCOM(cell,SecretionConstant)
    Secretor.secreteOutsideCellAtBoundary(cell,SecretionConstant)
  #
  Controlling time:
    self.setMaxMCS(100000)
    
    self.stopSimulation()
  #
  Controlling Lattice:
    
    self.resizeAndShiftLattice(_newSize=(X,Y,Z), _shiftVec=(VX,VY,VZ))
  #
  Cluster Boundary Pixels:
    for pixel in self.getClusterBoundary(cell.clusterId):
      ...
    
    def getClusterBoundary(self,clusterId):
      self.boundaryStrategy=CompuCell.BoundaryStrategy.getInstance()
      self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromNeighborOrder(1)
      L=[]; compList=self.inventory.getClusterCells(clusterId)
      for cell in compList: #going over all compartments
        pixelList=self.getCellBoundaryPixelList(cell) 
        for bPixel in pixelList: #going over all boundary pixels of each compartment
          for i in xrange(self.maxNeighborIndex+1): #looping over boundary pixels' neighbors
            pN=self.boundaryStrategy.getNeighborDirect(bPixel.pixel,i)
            cell=self.cellField[pN.pt.x,pN.pt.y,pN.pt.z] #getting cell at each neighbor pixel
            if (not cell or cell.clusterId!=clusterId): #Check if neighbor pixel belong to a different cluster
              L.append(bPixel.pixel)
              break
      return L
  #
  Getting Outside Boundary Pixels:
    self.boundaryStrategy=CompuCell.BoundaryStrategy.getInstance()
    self.maxNeighborIndex=self.boundaryStrategy.getMaxNeighborIndexFromNeighborOrder(1)

    pixelList=self.getCellBoundaryPixelList(cell)
    L=[]
    for bPixel in pixelList:
      for i in xrange(self.maxNeighborIndex+1): 
        pN=self.boundaryStrategy.getNeighborDirect(bPixel.pixel,i)
        x,y,z=pN.pt.x,pN.pt.y,pN.pt.z
        if ([x,y,z] not in L):
          cell2=self.cellField[x,y,z]
          if (not cell2 or cell2.id!=cell.id):
            L.append([x,y,z])
    return L
  #
  Cutting a region of the cells (Wound Infliction) 2D
    def cut(self,x0,y0,dx,dy):
        medium=CompuCell.getMediumCell() #get medium
        deleteCells = []
        for xx in range(int(dx)):
            x = int(x0-dx/2 + xx + .5)
            for yy in range(int(dy)):
                y = int(y0-dy/2 + yy + .5)
                for z in range(self.dim.z):
                    cell=self.cellField[x,y,z]
                    if cell:
                        cell.targetVolume-=1      #decrease cell target volume by 1
                        if (-1<cell.targetVolume<0): deleteCells.append(cell)
                    self.cellField[x,y,z]=medium  #replacing cell pixel by medium
        for cell in deleteCells:  #making sure cells are deleted
            self.deleteCell(cell)
#
#
##