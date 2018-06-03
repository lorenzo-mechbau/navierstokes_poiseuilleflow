#> \file
#> \author Chris Bradley
#> \brief This example program solves Poiseuille flow in a multi-block tube
#>        using OpenCMISS.
#>

#================================================================================================================================
#  Start Program
#================================================================================================================================

# Import the libraries (OpenCMISS,python,numpy,scipy)
import math,numpy,csv,time,sys,os,pdb
from opencmiss.iron import iron

# Set the number of elements in the mesh
numberOfSquareElements = 2
numberOfArmElements = 1
numberOfLengthElements = 3
# Set the Reynolds number
Re = 1000
# Set the maximum inlet velocity
maxInletFlow = 1.0
# Time stepping parameters
startTime = 0.0
stopTime  = 2.0
timeStep  = 0.1

# Override with command line arguments if need be
print sys.argv
if len(sys.argv) > 1:
    if len(sys.argv) > 9:
        sys.exit('Error: too many arguments - currently only accepting 8 options: numberOfSquareElements numberOfArmElements numberOfLengthElements Re maxInletFlow startTime stopTime timeStep')        
    numberOfSquareElements = int(sys.argv[1])
    if len(sys.argv) > 2:
        numberOfArmElements = int(sys.argv[2])
    if len(sys.argv) > 3:
        numberOfLengthElements = int(sys.argv[3])
    if len(sys.argv) > 4:
        Re = float(sys.argv[4])
    if len(sys.argv) > 5:
        maxInletFlow = float(sys.argv[5])
    if len(sys.argv) > 6:
        startTime = float(sys.argv[6])
    if len(sys.argv) > 7:
        stopTime = float(sys.argv[7])
    if len(sys.argv) > 8:
        timeStep = float(sys.argv[8])

LINEAR = 1
QUADRATIC = 2

pipeRadius = 1.0
lengthSize = 3.0
squareSizeRatio = 0.500

fluidVelocityInterpolation = QUADRATIC
fluidPressureInterpolation = LINEAR

## Inlet velocity parameters
A = maxInletFlow/2.0
B = maxInletFlow/2.0
C = stopTime/5.0

fluidDensity = 1.0
fluidDynamicViscosity = fluidDensity*maxInletFlow*(2.0*pipeRadius)/Re

fluidPInit = 0.0

#RBS = False
RBS = True

outputFrequency = 1 # Result output frequency

# Output flags
fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.NONE
#fluidEquationsSetOutputType = iron.EquationsSetOutputTypes.PROGRESS
fluidEquationsOutputType = iron.EquationsOutputTypes.NONE
#fluidEquationsOutputType = iron.EquationsOutputTypes.TIMING
#fluidEquationsOutputType = iron.EquationsOutputTypes.MATRIX
#fluidEquationsOutputType = iron.EquationsOutputTypes.ELEMENT_MATRIX
fluidDynamicSolverOutputType = iron.SolverOutputTypes.NONE
#fluidDynamicSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidDynamicSolverOutputType = iron.SolverOutputTypes.MATRIX
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.NONE
fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MONITOR
#fluidNonlinearSolverOutputType = iron.SolverOutputTypes.MATRIX
fluidLinearSolverOutputType = iron.SolverOutputTypes.NONE
#fluidLinearSolverOutputType = iron.SolverOutputTypes.PROGRESS
#fluidLinearSolverOutputType = iron.SolverOutputTypes.MATRIX

# Set solver parameters
fluidDynamicSolverTheta    = [0.5]
nonlinearMaximumIterations      = 100000000 #default: 100000
nonlinearRelativeTolerance      = 1.0E-9   #default: 1.0E-05
nonlinearSolutionTolerance      = 1.0E-9   #default: 1.0E-05
nonlinearAbsoluteTolerance      = 1.0E-8    #default: 1.0E-10
nonlinearMaxFunctionEvaluations = 10000
nonlinearLinesearchAlpha        = 1.0
linearMaximumIterations      = 100000000 #default: 100000
linearRelativeTolerance      = 1.0E-6    #default: 1.0E-05
linearAbsoluteTolerance      = 1.0E-6    #default: 1.0E-10
linearDivergenceTolerance    = 1.0E5     #default: 1.0E5
linearRestartValue           = 30        #default: 30

progressDiagnostics = True
debug = False

#================================================================================================================================
#  Should not need to change anything below here.
#================================================================================================================================

if numberOfLengthElements == 0:
    numberOfDimensions = 2
else:
    numberOfDimensions = 3

fluidCoordinateSystemUserNumber = 1
  
fluidRegionUserNumber = 1

fluidVelocityBasisUserNumber = 1
fluidPressureBasisUserNumber = 2

fluidMeshUserNumber = 1
  
fluidDecompositionUserNumber = 1
  
fluidGeometricFieldUserNumber = 1
fluidEquationsSetFieldUserNumber = 2
fluidDependentFieldUserNumber = 3
fluidMaterialsFieldUserNumber = 4
fluidIndependentFieldUserNumber = 5
bcCellMLModelsFieldUserNumber = 6
bcCellMLStateFieldUserNumber = 7
bcCellMLParametersFieldUserNumber = 8
bcCellMLIntermediateFieldUserNumber = 9

fluidEquationsSetUserNumber  = 1

bcCellMLUserNumber = 1
  
fluidProblemUserNumber = 1
 
#================================================================================================================================
#  Initialise OpenCMISS
#================================================================================================================================

worldRegion = iron.Region()
iron.Context.WorldRegionGet(worldRegion)

# Set the OpenCMISS random seed so that we can test this example by using the
# same parallel decomposition
numberOfRandomSeeds = iron.Context.RandomSeedsSizeGet()
randomSeeds = [0]*numberOfRandomSeeds
randomSeeds[0] = 100
iron.Context.RandomSeedsSet(randomSeeds)

# Get the computational nodes info
computationEnvironment = iron.ComputationEnvironment()
iron.Context.ComputationEnvironmentGet(computationEnvironment)
numberOfComputationalNodes = computationEnvironment.NumberOfWorldNodesGet()
computationalNodeNumber = computationEnvironment.WorldNodeNumberGet()

#================================================================================================================================
#  Coordinate Systems
#================================================================================================================================

if (progressDiagnostics):
    print(' ')
    print('Coordinate systems ...')

# Create a RC coordinate system for the fluid region
fluidCoordinateSystem = iron.CoordinateSystem()
fluidCoordinateSystem.CreateStart(fluidCoordinateSystemUserNumber,iron.Context)
fluidCoordinateSystem.DimensionSet(numberOfDimensions)
fluidCoordinateSystem.CreateFinish()

if (progressDiagnostics):
    print('Coordinate systems ... Done')
  
#================================================================================================================================
#  Regions
#================================================================================================================================

if (progressDiagnostics):
    print('Regions ...')

# Create a fluid region
fluidRegion = iron.Region()
fluidRegion.CreateStart(fluidRegionUserNumber,worldRegion)
fluidRegion.label = 'FluidRegion'
fluidRegion.coordinateSystem = fluidCoordinateSystem
fluidRegion.CreateFinish()

if (progressDiagnostics):
    print('Regions ... Done')

#================================================================================================================================
#  Bases
#================================================================================================================================

if (progressDiagnostics):
    print('Basis functions ...')
    
numberOfNodesXi = fluidVelocityInterpolation+1
numberOfGaussXi = fluidVelocityInterpolation+1

fluidVelocityBasis = iron.Basis()
fluidVelocityBasis.CreateStart(fluidVelocityBasisUserNumber,iron.Context)
fluidVelocityBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
fluidVelocityBasis.numberOfXi = numberOfDimensions
if (fluidVelocityInterpolation == LINEAR):
    fluidVelocityBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions
elif (fluidVelocityInterpolation == QUADRATIC):
    fluidVelocityBasis.interpolationXi = [iron.BasisInterpolationSpecifications.QUADRATIC_LAGRANGE]*numberOfDimensions
else:
    print('Invalid velocity interpolation')
    exit()
fluidVelocityBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*numberOfDimensions
fluidVelocityBasis.CreateFinish()

fluidPressureBasis = iron.Basis()
fluidPressureBasis.CreateStart(fluidPressureBasisUserNumber,iron.Context)
fluidPressureBasis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
fluidPressureBasis.numberOfXi = numberOfDimensions
fluidPressureBasis.interpolationXi = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*numberOfDimensions
fluidPressureBasis.quadratureNumberOfGaussXi = [numberOfGaussXi]*numberOfDimensions
fluidPressureBasis.CreateFinish()

numberOfLocalNodes = numberOfNodesXi*numberOfNodesXi
if (numberOfDimensions==3):
    numberOfLocalNodes = numberOfLocalNodes*numberOfNodesXi
localNodeIdx000 = 0
localNodeIdx100 = numberOfNodesXi-1
localNodeIdx010 = numberOfNodesXi*(numberOfNodesXi-1)
localNodeIdx110 = numberOfNodesXi*numberOfNodesXi-1
if (numberOfDimensions==3):
    localNodeIdx001 = numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
    localNodeIdx101 = numberOfNodesXi-1+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
    localNodeIdx011 = numberOfNodesXi*(numberOfNodesXi-1)+numberOfNodesXi*numberOfNodesXi*(numberOfNodesXi-1)
    localNodeIdx111 = numberOfLocalNodes-1

if (progressDiagnostics):
    print('Basis functions ... Done')
  
#================================================================================================================================
#  Mesh
#================================================================================================================================

numberOfNodesPerBlock = numberOfSquareElements*(numberOfNodesXi-1)*(numberOfArmElements*(numberOfNodesXi-1)+1)
numberOfElementsPerBlock = numberOfSquareElements*numberOfArmElements
numberOfNodesPerLength = 4*numberOfNodesPerBlock+ \
                     (numberOfSquareElements*(numberOfNodesXi-1)-1)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
numberOfElementsPerLength = 4*numberOfElementsPerBlock+numberOfSquareElements*numberOfSquareElements
numberOfFluidNodes = numberOfNodesPerLength
numberOfFluidElements = numberOfElementsPerLength
if (numberOfDimensions==3):
    numberOfFluidNodes = numberOfFluidNodes*(numberOfLengthElements*(numberOfNodesXi-1)+1)
    numberOfFluidElements = numberOfFluidElements*numberOfLengthElements

if (debug):
    print('  Mesh Parameters:')
    print('    numberOfSquareElements: %d' % (numberOfSquareElements))
    print('    numberOfArmElements: %d' % (numberOfArmElements))
    if (numberOfDimensions==3):
        print('    numberOfLengthElements: %d' % (numberOfLengthElements))
    print('    numberOfNodesXi: %d' % (numberOfNodesXi))
    print('    numberOfNodesPerBlock: %d' % (numberOfNodesPerBlock))
    print('    numberOfElementPerBlock: %d' % (numberOfElementsPerBlock))
    print('    numberOfNodesPerLength: %d' % (numberOfNodesPerLength))
    print('    numberOfElementsPerLength: %d' % (numberOfElementsPerLength))
    print('    numberOfFluidNodes: %d' % (numberOfFluidNodes))
    print('    numberOfFluidElements: %d' % (numberOfFluidElements))
        
fluidNodes = iron.Nodes()
fluidNodes.CreateStart(fluidRegion,numberOfFluidNodes)
fluidNodes.CreateFinish()

fluidMesh = iron.Mesh()
fluidMesh.CreateStart(fluidMeshUserNumber,fluidRegion,numberOfDimensions)
fluidMesh.NumberOfElementsSet(numberOfFluidElements)
fluidMesh.NumberOfComponentsSet(2)

fluidVelocityElements = iron.MeshElements()
fluidVelocityElements.CreateStart(fluidMesh,1,fluidVelocityBasis)
fluidPressureElements = iron.MeshElements()
fluidPressureElements.CreateStart(fluidMesh,2,fluidPressureBasis)

if (debug):
    print('  Fluid Elements:')
for zElementIdx in range(1,max(numberOfLengthElements+1,2)):
    #Handle the arm blocks first
    previousBlock = 4
    for blockIdx in range(1,5):
        for yElementIdx in range(1,numberOfArmElements+1):
            for xElementIdx in range(1,numberOfSquareElements+1):
                localNodes = [0]*numberOfLocalNodes
                elementNumber = xElementIdx+(yElementIdx-1)*numberOfSquareElements+(blockIdx-1)*numberOfSquareElements*numberOfArmElements+\
                                (zElementIdx-1)*numberOfElementsPerLength
                if (xElementIdx == 1):
                    localNodes[localNodeIdx000] = (previousBlock-1)*numberOfNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)+ \
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                 (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx100] = (blockIdx-1)*numberOfNodesPerBlock+numberOfNodesXi-1+\
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                                                 (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                else:
                    localNodes[localNodeIdx000] = (blockIdx-1)*numberOfNodesPerBlock+(xElementIdx-2)*(numberOfNodesXi-1)+(numberOfNodesXi-2)+1+\
                                                 (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSquareElements*(numberOfNodesXi-1))+\
                                                 (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+numberOfNodesXi-1
                localNodes[localNodeIdx010] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1)*(numberOfNodesXi-1)
                localNodes[localNodeIdx110] = localNodes[localNodeIdx100] + numberOfSquareElements*(numberOfNodesXi-1)*(numberOfNodesXi-1)
                if (numberOfDimensions == 3):
                    localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfNodesPerLength*(numberOfNodesXi-1)
                if(fluidVelocityInterpolation == QUADRATIC):
                    localNodes[1] = localNodes[localNodeIdx100] - 1
                    localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1)
                    localNodes[4] = localNodes[1] + numberOfSquareElements*(numberOfNodesXi-1)
                    localNodes[5] = localNodes[4] + 1
                    localNodes[7] = localNodes[localNodeIdx110] - 1
                    if (numberOfDimensions == 3):
                        localNodes[9] = localNodes[0]+numberOfNodesPerLength
                        localNodes[10] = localNodes[1]+numberOfNodesPerLength
                        localNodes[11] = localNodes[2]+numberOfNodesPerLength
                        localNodes[12] = localNodes[3]+numberOfNodesPerLength
                        localNodes[13] = localNodes[4]+numberOfNodesPerLength
                        localNodes[14] = localNodes[5]+numberOfNodesPerLength
                        localNodes[15] = localNodes[6]+numberOfNodesPerLength
                        localNodes[16] = localNodes[7]+numberOfNodesPerLength
                        localNodes[17] = localNodes[8]+numberOfNodesPerLength
                        localNodes[19] = localNodes[10]+numberOfNodesPerLength
                        localNodes[21] = localNodes[12]+numberOfNodesPerLength
                        localNodes[22] = localNodes[13]+numberOfNodesPerLength
                        localNodes[23] = localNodes[14]+numberOfNodesPerLength
                        localNodes[25] = localNodes[16]+numberOfNodesPerLength
                if (numberOfDimensions == 2):
                    linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110]]
                    if (debug):
                        print('    Element %8d; Nodes: %8d, %8d, %8d, %8d' % (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3]))
                        if (fluidVelocityInterpolation==QUADRATIC):
                            print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                else:
                    linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                                   localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
                    if (debug):
                        print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                        if (fluidVelocityInterpolation==QUADRATIC):
                            print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                            print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                            print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                fluidPressureElements.NodesSet(elementNumber,linearNodes)
                fluidVelocityElements.NodesSet(elementNumber,localNodes)
        previousBlock = blockIdx
    #Handle the square block
    if (numberOfSquareElements==1):
        elementNumber = elementNumber + 1
        localNodes[localNodeIdx000] = 3*numberOfNodesPerBlock+\
                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx100] = 4*numberOfNodesPerBlock+\
                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx010] = 2*numberOfNodesPerBlock+\
                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
        localNodes[localNodeIdx110] = numberOfNodesPerBlock+\
                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
        if(fluidVelocityInterpolation == QUADRATIC):
            localNodes[1] = localNodes[localNodeIdx100] - 1
            localNodes[3] = localNodes[localNodeIdx000] - 1
            localNodes[4] = localNodes[localNodeIdx100] + 1
            localNodes[5] = localNodes[localNodeIdx110] - 1
            localNodes[7] = localNodes[localNodeIdx010] - 1
        if (numberOfDimensions == 2):
            linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110]]
            if (debug):
                print('    Element %8d; Nodes: %8d, %8d, %8d, %8d' % (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3]))
            if (fluidVelocityInterpolation==QUADRATIC):
                print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
        else:
            localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfNodesPerLength*(numberOfNodesXi-1)
            localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfNodesPerLength*(numberOfNodesXi-1)
            localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfNodesPerLength*(numberOfNodesXi-1)
            localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfNodesPerLength*(numberOfNodesXi-1)
            linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                           localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
            if (fluidVelocityInterpolation == QUADRATIC):
                localNodes[9] = localNodes[0]+numberOfNodesPerLength
                localNodes[10] = localNodes[1]+numberOfNodesPerLength
                localNodes[11] = localNodes[2]+numberOfNodesPerLength
                localNodes[12] = localNodes[3]+numberOfNodesPerLength
                localNodes[13] = localNodes[4]+numberOfNodesPerLength
                localNodes[14] = localNodes[5]+numberOfNodesPerLength
                localNodes[15] = localNodes[6]+numberOfNodesPerLength
                localNodes[16] = localNodes[7]+numberOfNodesPerLength
                localNodes[17] = localNodes[8]+numberOfNodesPerLength
                localNodes[19] = localNodes[10]+numberOfNodesPerLength
                localNodes[21] = localNodes[12]+numberOfNodesPerLength
                localNodes[22] = localNodes[13]+numberOfNodesPerLength
                localNodes[23] = localNodes[14]+numberOfNodesPerLength
                localNodes[25] = localNodes[16]+numberOfNodesPerLength
            if (debug):
                print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                      (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                if (fluidVelocityInterpolation==QUADRATIC):
                    print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                    print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                    print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                          (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                                
        fluidPressureElements.NodesSet(elementNumber,linearNodes)
        fluidVelocityElements.NodesSet(elementNumber,localNodes)
    else:
        for yElementIdx in range(1,numberOfSquareElements+1):
            for xElementIdx in range(1,numberOfSquareElements+1):
                localNodes = [0]*numberOfLocalNodes
                elementNumber = 4*numberOfElementsPerBlock+xElementIdx+(yElementIdx-1)*numberOfSquareElements+\
                                (zElementIdx-1)*numberOfElementsPerLength
                if (yElementIdx == 1):
                    if (xElementIdx == 1):
                        #Bottom-left
                        localNodes[localNodeIdx000] = 3*numberOfNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 3*numberOfNodesPerBlock+numberOfArmElements*(numberOfNodesXi-1)*\
                                                      numberOfSquareElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 3*numberOfNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-2)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements*(numberOfNodesXi-1)
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Bottom-right
                        localNodes[localNodeIdx000] = 4*numberOfNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)-(numberOfNodesXi-2)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx110] - 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                        elif(fluidVelocityInterpolation == cubic):
                            print("Not implemented.")
                            exit()
                    else:
                        #Bottom
                        localNodes[localNodeIdx000] = 3*numberOfNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-2)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]+(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                        elif(fluidVelocityInterpolation == cubic):
                            print("Not implemented.")
                            exit()
                elif (yElementIdx == numberOfSquareElements):
                    if (xElementIdx == 1):
                        #Top-left
                        localNodes[localNodeIdx000] = 2*numberOfNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 2*numberOfNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = 2*numberOfNodesPerBlock-(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[1] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] + 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Top-right
                        localNodes[localNodeIdx000] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+\
                                                      (numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+\
                                                      (numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = numberOfNodesPerBlock+numberOfSquareElements*(numberOfNodesXi-1)*\
                                                      numberOfArmElements*(numberOfNodesXi-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = numberOfNodesPerBlock+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx110] - 1
                            localNodes[7] = localNodes[localNodeIdx010] - 1
                    else:
                        #Top
                        localNodes[localNodeIdx000] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((numberOfSquareElements-1)*(numberOfNodesXi-1)-1)+\
                                                      (xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = 2*numberOfNodesPerBlock-(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]-(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx010] - 1
                else:
                    if (xElementIdx == 1):
                        #Left
                        localNodes[localNodeIdx000] = 3*numberOfNodesPerBlock-(yElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]-(numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx100]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx100] - 1
                            localNodes[3] = localNodes[localNodeIdx000] - 1
                            localNodes[4] = localNodes[localNodeIdx110] - numberOfSquareElements*(numberOfNodesXi-1) 
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx110] - 1
                    elif (xElementIdx == numberOfSquareElements):
                        #Right
                        localNodes[localNodeIdx000] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(numberOfSquareElements-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = numberOfSquareElements*(numberOfNodesXi-1)*numberOfArmElements*(numberOfNodesXi-1)+\
                                                      (yElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx100]+(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx010] - numberOfSquareElements*(numberOfNodesXi-1) + 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[localNodeIdx100] + 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                    else:
                        #Middle
                        localNodes[localNodeIdx000] = 4*numberOfNodesPerBlock+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      ((yElementIdx-1)*(numberOfNodesXi-1)-1)+(xElementIdx-1)*(numberOfNodesXi-1)+\
                                                      (zElementIdx-1)*numberOfNodesPerLength*(numberOfNodesXi-1)
                        localNodes[localNodeIdx100] = localNodes[localNodeIdx000]+(numberOfNodesXi-1)
                        localNodes[localNodeIdx010] = localNodes[localNodeIdx000]+(numberOfSquareElements*(numberOfNodesXi-1)-1)*\
                                                      (numberOfNodesXi-1)
                        localNodes[localNodeIdx110] = localNodes[localNodeIdx010]+(numberOfNodesXi-1)
                        if(fluidVelocityInterpolation == QUADRATIC):
                            localNodes[1] = localNodes[localNodeIdx000] + 1
                            localNodes[3] = localNodes[localNodeIdx000] + numberOfSquareElements*(numberOfNodesXi-1) - 1
                            localNodes[4] = localNodes[3] + 1
                            localNodes[5] = localNodes[4] + 1
                            localNodes[7] = localNodes[localNodeIdx010] + 1
                if (numberOfDimensions == 2):
                    linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110]]
                    if (debug):
                        print('    Element %8d; Nodes: %8d, %8d, %8d, %8d' % (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3]))
                    if (fluidVelocityInterpolation==QUADRATIC):
                        print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                              (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                else:
                    localNodes[localNodeIdx001] = localNodes[localNodeIdx000]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx101] = localNodes[localNodeIdx100]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx011] = localNodes[localNodeIdx010]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    localNodes[localNodeIdx111] = localNodes[localNodeIdx110]+numberOfNodesPerLength*(numberOfNodesXi-1)
                    linearNodes = [localNodes[localNodeIdx000],localNodes[localNodeIdx100],localNodes[localNodeIdx010],localNodes[localNodeIdx110], \
                                   localNodes[localNodeIdx001],localNodes[localNodeIdx101],localNodes[localNodeIdx011],localNodes[localNodeIdx111]]
                    if (fluidVelocityInterpolation == QUADRATIC):
                        localNodes[9] = localNodes[0]+numberOfNodesPerLength
                        localNodes[10] = localNodes[1]+numberOfNodesPerLength
                        localNodes[11] = localNodes[2]+numberOfNodesPerLength
                        localNodes[12] = localNodes[3]+numberOfNodesPerLength
                        localNodes[13] = localNodes[4]+numberOfNodesPerLength
                        localNodes[14] = localNodes[5]+numberOfNodesPerLength
                        localNodes[15] = localNodes[6]+numberOfNodesPerLength
                        localNodes[16] = localNodes[7]+numberOfNodesPerLength
                        localNodes[17] = localNodes[8]+numberOfNodesPerLength
                        localNodes[19] = localNodes[10]+numberOfNodesPerLength
                        localNodes[21] = localNodes[12]+numberOfNodesPerLength
                        localNodes[22] = localNodes[13]+numberOfNodesPerLength
                        localNodes[23] = localNodes[14]+numberOfNodesPerLength
                        localNodes[25] = localNodes[16]+numberOfNodesPerLength
                        if (debug):
                            print('    Element %8d; Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                  (elementNumber,linearNodes[0],linearNodes[1],linearNodes[2],linearNodes[3],linearNodes[4],linearNodes[5],linearNodes[6],linearNodes[7]))
                            if (fluidVelocityInterpolation==QUADRATIC):
                                print('                      Nodes: %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                      (localNodes[0],localNodes[1],localNodes[2],localNodes[3],localNodes[4],localNodes[5],localNodes[6],localNodes[7],localNodes[8]))
                                print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                      (localNodes[9],localNodes[10],localNodes[11],localNodes[12],localNodes[13],localNodes[14],localNodes[15],localNodes[16],localNodes[17]))
                                print('                             %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d, %8d' % \
                                      (localNodes[18],localNodes[19],localNodes[20],localNodes[21],localNodes[22],localNodes[23],localNodes[24],localNodes[25],localNodes[26]))
                                
                fluidPressureElements.NodesSet(elementNumber,linearNodes)
                fluidVelocityElements.NodesSet(elementNumber,localNodes)

fluidVelocityElements.CreateFinish()
fluidPressureElements.CreateFinish()

fluidMesh.CreateFinish()

if (progressDiagnostics):
    print('Meshes ... Done')    

#================================================================================================================================
#  Decomposition
#================================================================================================================================

if (progressDiagnostics):
    print('Decomposition ...')
    
# Create a decomposition for the fluid mesh
fluidDecomposition = iron.Decomposition()
fluidDecomposition.CreateStart(fluidDecompositionUserNumber,fluidMesh)
fluidDecomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
fluidDecomposition.NumberOfDomainsSet(numberOfComputationalNodes)
fluidDecomposition.CalculateFacesSet(True)
fluidDecomposition.CreateFinish()

if (progressDiagnostics):
    print('Decomposition ... Done')
    
#================================================================================================================================
#  Geometric Field
#================================================================================================================================

if (progressDiagnostics):
    print('Geometric Field ...')
    
# Start to create a default (geometric) field on the fluid region
fluidGeometricField = iron.Field()
fluidGeometricField.CreateStart(fluidGeometricFieldUserNumber,fluidRegion)
# Set the decomposition to use
fluidGeometricField.DecompositionSet(fluidDecomposition)
# Set the scaling to use
fluidGeometricField.ScalingTypeSet(iron.FieldScalingTypes.NONE)
fluidGeometricField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidGeometry')
# Set the domain to be used by the field components.
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
if (numberOfDimensions == 3):
    fluidGeometricField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
# Finish creating the second field
fluidGeometricField.CreateFinish()

if (progressDiagnostics):
    print('Geometric Field ... Done')
    
if (progressDiagnostics):
    print('Geometric Parameters ...')

armSize = (1.0-squareSizeRatio)*pipeRadius
squareSize = 2.0*(pipeRadius-armSize)/math.sqrt(2.0)

for zNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
    #Handle the arm blocks first
    previousBlock = 4
    for blockIdx in range(1,5):
        for yNodeIdx in range(1,numberOfArmElements*(numberOfNodesXi-1)+2):
            for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
                nodeNumber = (blockIdx-1)*numberOfNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)+\
                             (zNodeIdx-1)*numberOfNodesPerLength
                nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
                if (nodeDomain == computationalNodeNumber):
                    if (yNodeIdx == numberOfArmElements*(numberOfNodesXi-1)+1):
                        #On the square
                        if (blockIdx == 1):
                            xPosition = squareSize/2.0
                            yPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize/2.0
                        elif (blockIdx == 2):
                            xPosition = squareSize/2.0-xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                            yPosition = squareSize/2.0
                        elif (blockIdx == 3):
                            xPosition = -squareSize/2.0
                            yPosition = squareSize/2.0-xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))
                        elif (blockIdx == 4):
                            xPosition = xNodeIdx*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize/2.0
                            yPosition = -squareSize/2.0
                    else:
                        #In the arm
                        #Work out the r, theta position
                        theta = xNodeIdx*math.pi/(2*numberOfSquareElements*(numberOfNodesXi-1))+(blockIdx-1)*math.pi/2.0-math.pi/4.0
                        armRadius = armSize+(math.sqrt(2.0)-1.0)*squareSize*math.cos(2.0*theta)*math.cos(2.0*theta)/2.0
                        squareRadius = squareSize/2.0+(math.sqrt(2.0)-1.0)*squareSize*math.sin(2.0*theta)*math.sin(2.0*theta)/2.0
                        radius = (numberOfArmElements*(numberOfNodesXi-1)-yNodeIdx+1)*armRadius/ \
                                 (numberOfArmElements*(numberOfNodesXi-1))+squareRadius
                        xPosition = radius*math.cos(theta)
                        yPosition = radius*math.sin(theta)
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
                    if (numberOfDimensions == 2):
                        if (debug):
                            print('      Node        %d:' % (nodeNumber))
                            print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))
                    else:
                        zPosition = (zNodeIdx-1)*lengthSize/(numberOfLengthElements*(numberOfNodesXi-1))
                        fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                     1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
                        if (debug):
                            print('      Node        %d:' % (nodeNumber))
                            print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition,yPosition,zPosition))

    #Now handle square
    for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = 4*numberOfNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)+\
                         (zNodeIdx-1)*numberOfNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                xPosition = (xNodeIdx-1)*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize/2.0
                yPosition = (yNodeIdx-1)*squareSize/(numberOfSquareElements*(numberOfNodesXi-1))-squareSize/2.0
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,xPosition)
                fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                             1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,2,yPosition)
                if (numberOfDimensions == 2):
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Position         = [ %.2f, %.2f ]' % (xPosition,yPosition))                                     
                else:
                    zPosition = (zNodeIdx-1)*lengthSize/(numberOfLengthElements*(numberOfNodesXi-1))
                    fluidGeometricField.ParameterSetUpdateNodeDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                                 1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,3,zPosition)
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Position         = [ %.2f, %.2f, %.2f ]' % (xPosition,yPosition,zPosition))
                        
# Update fields            
fluidGeometricField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidGeometricField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Geometric Parameters ... Done')

#================================================================================================================================
#  Equations Set
#================================================================================================================================

if (progressDiagnostics):
    print('Equations Sets ...')

# Create the equations set for the fluid region - Navier-Stokes
fluidEquationsSetField = iron.Field()
fluidEquationsSet = iron.EquationsSet()
if RBS:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
else:
    fluidEquationsSetSpecification = [iron.EquationsSetClasses.FLUID_MECHANICS,
                                      iron.EquationsSetTypes.NAVIER_STOKES_EQUATION,
                                      iron.EquationsSetSubtypes.TRANSIENT_NAVIER_STOKES]
fluidEquationsSet.CreateStart(fluidEquationsSetUserNumber,fluidRegion,fluidGeometricField,
                              fluidEquationsSetSpecification,fluidEquationsSetFieldUserNumber,
                              fluidEquationsSetField)
fluidEquationsSet.OutputTypeSet(fluidEquationsSetOutputType)
fluidEquationsSet.CreateFinish()

if RBS:
    # Set max CFL number (default 1.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,2,1.0E20)
    # Set time increment (default 0.0)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,3,timeStep)
    # Set stabilisation type (default 1.0 = RBS)
    fluidEquationsSetField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U1,
                                                       iron.FieldParameterSetTypes.VALUES,4,1.0)

if (progressDiagnostics):
    print('Equations Sets ... Done')


#================================================================================================================================
#  Dependent Field
#================================================================================================================================

if (progressDiagnostics):
    print('Dependent Fields ...')

# Create the equations set dependent field variables for dynamic Navier-Stokes
fluidDependentField = iron.Field()
fluidEquationsSet.DependentCreateStart(fluidDependentFieldUserNumber,fluidDependentField)
fluidDependentField.VariableLabelSet(iron.FieldVariableTypes.U,'FluidDependent')
# Set the mesh component to be used by the field components.
for componentIdx in range(1,numberOfDimensions+1):
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,componentIdx,1)
    fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,componentIdx,1)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.U,numberOfDimensions+1,2)
fluidDependentField.ComponentMeshComponentSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,2)
# fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.U,numberOfDimensions+1,iron.FieldInterpolationTypes.NODE_BASED)
# fluidDependentField.ComponentInterpolationSet(iron.FieldVariableTypes.DELUDELN,numberOfDimensions+1,iron.FieldInterpolationTypes.NODE_BASED)
# Finish the equations set dependent field variables
fluidEquationsSet.DependentCreateFinish()

# Initialise the fluid dependent field
for componentIdx in range(1,numberOfDimensions+1):
    fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,componentIdx,0.0)
# Initialise pressure component
fluidDependentField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                numberOfDimensions+1,fluidPInit)
     
fluidDependentField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
fluidDependentField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)

if (progressDiagnostics):
    print('Dependent Fields ... Done')
     
#================================================================================================================================
#  Materials Field
#================================================================================================================================

if (progressDiagnostics):
    print('Materials Fields ...')

# Create the equations set materials field variables for dynamic Navier-Stokes
fluidMaterialsField = iron.Field()
fluidEquationsSet.MaterialsCreateStart(fluidMaterialsFieldUserNumber,fluidMaterialsField)
# Finish the equations set materials field variables
fluidEquationsSet.MaterialsCreateFinish()
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,fluidDynamicViscosity)
fluidMaterialsField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,fluidDensity)

if (progressDiagnostics):
    print('Materials Fields ... Done')
    
#================================================================================================================================
#  Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Equations ...')

# Fluid equations 
fluidEquations = iron.Equations()
fluidEquationsSet.EquationsCreateStart(fluidEquations)
fluidEquations.sparsityType = iron.EquationsSparsityTypes.SPARSE
fluidEquations.outputType = fluidEquationsOutputType
fluidEquationsSet.EquationsCreateFinish()

if (progressDiagnostics):
    print('Equations ... Done')

#================================================================================================================================
#  CellML
#================================================================================================================================

if (progressDiagnostics):
    print('CellML ...')

# Create CellML equations for the temporal boundary conditions
bcCellML = iron.CellML()
bcCellML.CreateStart(bcCellMLUserNumber,fluidRegion)
bcCellMLIdx = bcCellML.ModelImport("input/poiseuilleinlet.cellml")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/pipeRadius")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/length")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/dynamicViscosity")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/A")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/B")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/C")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/x")
bcCellML.VariableSetAsKnown(bcCellMLIdx,"main/y")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inletx")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inlety")
bcCellML.VariableSetAsWanted(bcCellMLIdx,"main/inletz")
bcCellML.CreateFinish()

# Create CellML <--> OpenCMISS field maps
bcCellML.FieldMapsCreateStart()
# Map geometric field to x and y
bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/x",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidGeometricField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/y",iron.FieldParameterSetTypes.VALUES)
# Map fluid velocity to ensure dependent field isn't cleared when the velocities are copied back
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/inletx",iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES,
	                        bcCellMLIdx,"main/inlety",iron.FieldParameterSetTypes.VALUES)
if (numberOfDimensions==3):
    bcCellML.CreateFieldToCellMLMap(fluidDependentField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES,
	                            bcCellMLIdx,"main/inletz",iron.FieldParameterSetTypes.VALUES)
# Map inletx, inlety and inletz to dependent field
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inletx",iron.FieldParameterSetTypes.VALUES,
	                        fluidDependentField,iron.FieldVariableTypes.U,1,iron.FieldParameterSetTypes.VALUES)
bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inlety",iron.FieldParameterSetTypes.VALUES,
	                        fluidDependentField,iron.FieldVariableTypes.U,2,iron.FieldParameterSetTypes.VALUES)
if (numberOfDimensions==3):
    bcCellML.CreateCellMLToFieldMap(bcCellMLIdx,"main/inletz",iron.FieldParameterSetTypes.VALUES,
                                    fluidDependentField,iron.FieldVariableTypes.U,3,iron.FieldParameterSetTypes.VALUES)
bcCellML.FieldMapsCreateFinish()

# Create the CellML models field
bcCellMLModelsField = iron.Field()
bcCellML.ModelsFieldCreateStart(bcCellMLModelsFieldUserNumber,bcCellMLModelsField)
bcCellMLModelsField.VariableLabelSet(iron.FieldVariableTypes.U,"BCModelMap")
bcCellML.ModelsFieldCreateFinish()

# Only evaluate BC on inlet nodes
bcCellMLModelsField.ComponentValuesInitialiseIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,0)
if (debug):
    print('  CellML Boundary Conditions:')
    print('    Inlet Model Set:')
for blockIdx in range(1,5):
    for yNodeIdx in range(2,numberOfArmElements*(numberOfNodesXi-1)+2):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                               1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
    for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        nodeNumber = 4*numberOfNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            bcCellMLModelsField.ParameterSetUpdateNodeIntg(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,
                                                           1,iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV,nodeNumber,1,1)
            if (debug):
                print('      Node        %d:' % (nodeNumber))

bcCellMLModelsField.ParameterSetUpdateStart(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
bcCellMLModelsField.ParameterSetUpdateFinish(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES)
                
# Create the CellML state field
bcCellMLStateField = iron.Field()
bcCellML.StateFieldCreateStart(bcCellMLStateFieldUserNumber,bcCellMLStateField)
bcCellMLStateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCState")
bcCellML.StateFieldCreateFinish()

# Create the CellML parameters field
bcCellMLParametersField = iron.Field()
bcCellML.ParametersFieldCreateStart(bcCellMLParametersFieldUserNumber,bcCellMLParametersField)
bcCellMLParametersField.VariableLabelSet(iron.FieldVariableTypes.U,"BCParameters")
bcCellML.ParametersFieldCreateFinish()

# Get the component numbers
pipeRadiusComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/pipeRadius")
lengthComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/length")
dynamicViscosityComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/dynamicViscosity")
AComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/A")
BComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/B")
CComponentNumber = bcCellML.FieldComponentGet(bcCellMLIdx,iron.CellMLFieldTypes.PARAMETERS,"main/C")
# Set up the parameters field
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    pipeRadiusComponentNumber,pipeRadius)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    lengthComponentNumber,lengthSize)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    dynamicViscosityComponentNumber,fluidDynamicViscosity)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    AComponentNumber,A)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    BComponentNumber,B)
bcCellMLParametersField.ComponentValuesInitialiseDP(iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES, \
                                                    CComponentNumber,C)

# Create the CELL intermediate field
bcCellMLIntermediateField = iron.Field()
bcCellML.IntermediateFieldCreateStart(bcCellMLIntermediateFieldUserNumber,bcCellMLIntermediateField)
bcCellMLIntermediateField.VariableLabelSet(iron.FieldVariableTypes.U,"BCIntermediate")
bcCellML.IntermediateFieldCreateFinish()

if (progressDiagnostics):
    print('CellML ... Done')

#================================================================================================================================
#  Problem
#================================================================================================================================

if (progressDiagnostics):
    print('Problems ...')

# Create a fluid problem
fluidProblem = iron.Problem()
if RBS:
    fluidProblemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                 iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                 iron.ProblemSubtypes.TRANSIENT_RBS_NAVIER_STOKES]
else:
    fluidProblemSpecification = [iron.ProblemClasses.FLUID_MECHANICS,
                                 iron.ProblemTypes.NAVIER_STOKES_EQUATION,
                                 iron.ProblemSubtypes.TRANSIENT_NAVIER_STOKES]
fluidProblem.CreateStart(fluidProblemUserNumber,iron.Context,fluidProblemSpecification)
fluidProblem.CreateFinish()

if (progressDiagnostics):
    print('Problems ... Done')

#================================================================================================================================
#  Control Loop
#================================================================================================================================

if (progressDiagnostics):
    print('Control Loops ...')

# Create the fluid problem control loop
fluidControlLoop = iron.ControlLoop()
fluidProblem.ControlLoopCreateStart()
fluidProblem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE],fluidControlLoop)
fluidControlLoop.LabelSet('TimeLoop')
fluidControlLoop.TimesSet(startTime,stopTime,timeStep)
fluidControlLoop.TimeOutputSet(outputFrequency)
fluidProblem.ControlLoopCreateFinish()

if (progressDiagnostics):
    print('Control Loops ... Done')

#================================================================================================================================
#  Solvers
#================================================================================================================================

if (progressDiagnostics):
    print('Solvers ...')

# Create the problem solver
bcCellMLEvaluationSolver = iron.Solver()
fluidDynamicSolver = iron.Solver()
fluidNonlinearSolver = iron.Solver()
fluidLinearSolver = iron.Solver()
movingMeshLinearSolver = iron.Solver()

fluidProblem.SolversCreateStart()
# Solvers for a Navier Stokes problem
# Get the BC CellML solver
fluidProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],1,bcCellMLEvaluationSolver)
bcCellMLEvaluationSolver.outputType = iron.SolverOutputTypes.PROGRESS
# Get the dynamic solver
fluidProblem.SolverGet([iron.ControlLoopIdentifiers.NODE],2,fluidDynamicSolver)
fluidDynamicSolver.OutputTypeSet(fluidDynamicSolverOutputType)
fluidDynamicSolver.DynamicThetaSet(fluidDynamicSolverTheta)
# Get the dynamic nonlinear solver
fluidDynamicSolver.DynamicNonlinearSolverGet(fluidNonlinearSolver)
fluidNonlinearSolver.NewtonLineSearchTypeSet(iron.NewtonLineSearchTypes.LINEAR)
fluidNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.EQUATIONS) #(.FD/EQUATIONS)
#fluidNonlinearSolver.NewtonJacobianCalculationTypeSet(iron.JacobianCalculationTypes.FD) #(.FD/EQUATIONS)
fluidNonlinearSolver.NewtonMaximumFunctionEvaluationsSet(nonlinearMaxFunctionEvaluations)
fluidNonlinearSolver.OutputTypeSet(fluidNonlinearSolverOutputType)
fluidNonlinearSolver.NewtonAbsoluteToleranceSet(nonlinearAbsoluteTolerance)
fluidNonlinearSolver.NewtonMaximumIterationsSet(nonlinearMaximumIterations)
fluidNonlinearSolver.NewtonRelativeToleranceSet(nonlinearRelativeTolerance)
fluidNonlinearSolver.NewtonLineSearchAlphaSet(nonlinearLinesearchAlpha)
# Get the dynamic nonlinear linear solver
fluidNonlinearSolver.NewtonLinearSolverGet(fluidLinearSolver)
#fluidLinearSolver.LinearTypeSet(iron.LinearSolverTypes.ITERATIVE)
#fluidLinearSolver.LinearIterativeMaximumIterationsSet(linearMaximumIterations)
#fluidLinearSolver.LinearIterativeDivergenceToleranceSet(linearDivergenceTolerance)
#fluidLinearSolver.LinearIterativeRelativeToleranceSet(linearRelativeTolerance)
#fluidLinearSolver.LinearIterativeAbsoluteToleranceSet(linearAbsoluteTolerance)
fluidLinearSolver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
fluidLinearSolver.OutputTypeSet(fluidLinearSolverOutputType)
# Finish the creation of the problem solver
fluidProblem.SolversCreateFinish()

if (progressDiagnostics):
    print('Solvers ... Done')

#================================================================================================================================
#  CellML Equations
#================================================================================================================================

if (progressDiagnostics):
    print('CellML Equations ...')

# Create CellML equations and add BC equations to the solver
bcEquations = iron.CellMLEquations()
fluidProblem.CellMLEquationsCreateStart()
bcCellMLEvaluationSolver.CellMLEquationsGet(bcEquations)
bcEquationsIndex = bcEquations.CellMLAdd(bcCellML)
fluidProblem.CellMLEquationsCreateFinish()

if (progressDiagnostics):
    print('CellML Equations ... Done')

#================================================================================================================================
#  Solver Equations
#================================================================================================================================

if (progressDiagnostics):
    print('Solver Equations ...')

# Start the creation of the fluid problem solver equations
fluidProblem.SolverEquationsCreateStart()
# Get the fluid dynamic solver equations
fluidSolverEquations = iron.SolverEquations()
fluidDynamicSolver.SolverEquationsGet(fluidSolverEquations)
fluidSolverEquations.sparsityType = iron.SolverEquationsSparsityTypes.SPARSE
fluidEquationsSetIndex = fluidSolverEquations.EquationsSetAdd(fluidEquationsSet)
# Finish the creation of the fluid problem solver equations
fluidProblem.SolverEquationsCreateFinish()

if (progressDiagnostics):
    print('Solver Equations ...')

#================================================================================================================================
#  Boundary Conditions
#================================================================================================================================

if (progressDiagnostics):
    print('Boundary Conditions ...')

# Start the creation of the fluid boundary conditions
fluidBoundaryConditions = iron.BoundaryConditions()
fluidSolverEquations.BoundaryConditionsCreateStart(fluidBoundaryConditions)

if (debug):
    print('  Fluid Boundary Conditions:')
    print('    Inlet Boundary conditions:')
# Set inlet boundary conditions on the left hand edge
for blockIdx in range(1,5):
    for yNodeIdx in range(2,numberOfArmElements*(numberOfNodesXi-1)+2):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfNodesPerBlock+xNodeIdx+(yNodeIdx-1)*numberOfSquareElements*(numberOfNodesXi-1)
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                if (numberOfDimensions==2):
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Velocity         = [ %.2f, %.2f ]' % (inletVelocity,0.0))                 
                else:
                    fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                    nodeNumber,3,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Velocity         = [ %.2f, %.2f, %.2f ]' % (0.0,0.0,0.0))                 
for yNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
    for xNodeIdx in range(2,numberOfSquareElements*(numberOfNodesXi-1)+1):
        nodeNumber = 4*numberOfNodesPerBlock+(xNodeIdx-1)+(yNodeIdx-2)*(numberOfSquareElements*(numberOfNodesXi-1)-1)
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
        if (nodeDomain == computationalNodeNumber):
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumber,1,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumber,2,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
            if (numberOfDimensions==2):
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))                 
            else:
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,3,iron.BoundaryConditionsTypes.FIXED_INLET,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Velocity         = [ %.2f, %.2f, %.2f ]' % (0.0,0.0,0.0))                 
if (debug):
    print('    Wall Boundary conditions:')
# Set no slip boundary conditions on the wall
for zNodeIdx in range(1,numberOfLengthElements*(numberOfNodesXi-1)+2):
    for blockIdx in range(1,5):
        for xNodeIdx in range(1,numberOfSquareElements*(numberOfNodesXi-1)+1):
            nodeNumber = (blockIdx-1)*numberOfNodesPerBlock+xNodeIdx+(zNodeIdx-1)*numberOfNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,1)
            if (nodeDomain == computationalNodeNumber):
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,1,iron.BoundaryConditionsTypes.FIXED,0.0)
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,2,iron.BoundaryConditionsTypes.FIXED,0.0)
                if (numberOfDimensions==2):
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Velocity         = [ %.2f, %.2f ]' % (0.0,0.0))                 
                else:
                    fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                    iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                    nodeNumber,3,iron.BoundaryConditionsTypes.FIXED,0.0)
                    if (debug):
                        print('      Node        %d:' % (nodeNumber))
                        print('         Velocity         = [ %.2f, %.2f, %.2f ]' % (0.0,0.0,0.0))                 
if (debug):
    print('    No Pressure Boundary conditions:')
# Set no pressure boundary conditions on the outlet
for blockIdx in range(1,5):
    for yElementIdx in range(2,numberOfArmElements+2):
        for xElementIdx in range(1,numberOfSquareElements+1):
            nodeNumber = (blockIdx-1)*numberOfNodesPerBlock+xElementIdx*(numberOfNodesXi-1)+\
                         (yElementIdx-1)*(numberOfNodesXi-1)*numberOfSquareElements*(numberOfNodesXi-1)+ \
                         numberOfLengthElements*(numberOfNodesXi-1)*numberOfNodesPerLength
            nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
            if (nodeDomain == computationalNodeNumber):
                fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                                iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                                nodeNumber,numberOfDimensions+1,iron.BoundaryConditionsTypes.FIXED,0.0)
                if (debug):
                    print('      Node        %d:' % (nodeNumber))
                    print('         Pressure         = %.2f' % (0.0))                 
for yElementIdx in range(1,numberOfSquareElements):
    for xElementIdx in range(1,numberOfSquareElements):
        nodeNumber = 4*numberOfNodesPerBlock+xElementIdx*(numberOfNodesXi-1)+numberOfSquareElements*(numberOfNodesXi-1)-1+\
                     (yElementIdx-1)*(numberOfNodesXi-1)*(numberOfSquareElements*(numberOfNodesXi-1)-1)+ \
                     numberOfLengthElements*(numberOfNodesXi-1)*numberOfNodesPerLength
        nodeDomain = fluidDecomposition.NodeDomainGet(nodeNumber,2)
        if (nodeDomain == computationalNodeNumber):
            fluidBoundaryConditions.SetNode(fluidDependentField,iron.FieldVariableTypes.U,1, \
                                            iron.GlobalDerivativeConstants.NO_GLOBAL_DERIV, \
                                            nodeNumber,numberOfDimensions+1,iron.BoundaryConditionsTypes.FIXED,0.0)
            if (debug):
                print('      Node        %d:' % (nodeNumber))
                print('         Pressure         = %.2f' % (0.0))                 

# Finish fluid boundary conditions
fluidSolverEquations.BoundaryConditionsCreateFinish()

if (progressDiagnostics):
    print('Boundary Conditions ... Done')

#================================================================================================================================
#  Run Solvers
#================================================================================================================================

#quit()

# Export results
fields = iron.Fields()
fields.CreateRegion(fluidRegion)
fields.NodesExport("Poiseuille","FORTRAN")
fields.ElementsExport("Poiseuille","FORTRAN")
fields.Finalise()

# Solve the problem
print('Solving problem...')
start = time.time()
fluidProblem.Solve()
end = time.time()
elapsed = end - start
print('Calculation Time = %3.4f' %elapsed)
print('Problem solved!')
print('#')

