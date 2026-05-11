from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import part
import material
import assembly
import step
import mesh
import job
import odbAccess
import sys

# =============================================================================
# PARAMETERS
# =============================================================================
L = 10.0
w = 1.0
E = 200E9
nu = 0.3
meshSize = 0.5

# =============================================================================
# MATERIAL CREATION
# =============================================================================
mdb.models['Model-1'].Material(name='Steel')
mdb.models['Model-1'].materials['Steel'].Elastic(table=((E, nu), ))

# =============================================================================
# PART CREATION
# =============================================================================
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=20.0)
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0, 0), point2=(L, w))
p = mdb.models['Model-1'].Part(name='Bar', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Bar']
p.BaseSolidExtrude(sketch=s, depth=w)
s.unsetPrimaryObject()

# Section assignment
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
    material='Steel', thickness=None)
p = mdb.models['Model-1'].parts['Bar']
c = p.cells
cells = c.getByBoundingBox(-0.1, -0.1, -0.1, L+0.1, w+0.1, w+0.1)
region = p.Set(cells=cells, name='All')
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

# Mesh
p.seedPart(size=meshSize, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

# =============================================================================
# ASSEMBLY
# =============================================================================
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Bar']
a.Instance(name='Bar-1', part=p, dependent=ON)

# =============================================================================
# STEP 1 SETUP (STATIC STEP)
# =============================================================================
mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial', 
    nlgeom=ON, initialInc=0.1, minInc=1e-09, maxNumInc=1000)

# Boundary conditions
a = mdb.models['Model-1'].rootAssembly
v = a.instances['Bar-1'].vertices
verts1 = v.getByBoundingBox(-0.1, -0.1, -0.1, 0.1, w+0.1, w+0.1)
region = a.Set(vertices=verts1, name='Fixed')
mdb.models['Model-1'].DisplacementBC(name='Fix', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)

# Load
f = a.instances['Bar-1'].faces
faces = f.getByBoundingBox(L-0.1, -0.1, -0.1, L+0.1, w+0.1, w+0.1)
region = a.Surface(side1Faces=faces, name='LoadFace')
mdb.models['Model-1'].Pressure(name='Load-1', createStepName='Step-1', 
    region=region, magnitude=-1e7)

# Output requests
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=('U',))

# =============================================================================
# RUN JOB 1
# =============================================================================
job_name = 'Job1'
mdb.Job(name=job_name, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
    contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
    resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, numDomains=1)

my_job = mdb.jobs[job_name]
my_job.submit()
my_job.waitForCompletion()

# =============================================================================
# CHECK DEFORMED COORDINATES
# =============================================================================
odb = odbAccess.openOdb(path=job_name + '.odb')
instance = odb.rootAssembly.instances['BAR-1']
last_frame = odb.steps['Step-1'].frames[-1]
U = last_frame.fieldOutputs['U']

# Build displacement dictionary
disp_dict = {v.nodeLabel: v.data for v in U.values}

# Check: Did right end move past X = 10.5?
max_x_deformed = 0
for node in instance.nodes:
    x_original = node.coordinates[0]
    
    if node.label in disp_dict:
        x_displacement = disp_dict[node.label][0]
        x_deformed = x_original + x_displacement
        
        if x_deformed > max_x_deformed:
            max_x_deformed = x_deformed

print("Maximum deformed X coordinate: {0}".format(max_x_deformed))

odb.close()

# =============================================================================
# CONDITIONAL STEP 2
# =============================================================================
THRESHOLD_X = 10.01

if max_x_deformed <= THRESHOLD_X:
    print("Did not deform past {0}! Running Step-2...".format(THRESHOLD_X))
    
    # Add Step 2
    mdb.models['Model-1'].StaticStep(name='Step-2', 
        previous='Step-1', nlgeom=ON, initialInc=0.1, minInc=1e-09, 
        maxNumInc=1000)
    
    # Increase load for Step-2
    mdb.models['Model-1'].loads['Load-1'].setValuesInStep(
        stepName='Step-2', magnitude=-8e7)
    
    # Run with both steps
    job_name_2 = 'Job2'
    mdb.Job(name=job_name_2, model='Model-1', description='', 
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, 
        queue=None, memory=90, memoryUnits=PERCENTAGE, 
        getMemoryFromAnalysis=True, nodalOutputPrecision=SINGLE, 
        echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, 
        historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numDomains=1)
    
    my_job2 = mdb.jobs[job_name_2]
    my_job2.submit()
    my_job2.waitForCompletion()
    
    print("Step-2 complete!")
else:
    print("Exceeded threshold. No Step-2.")


sys.exit(0)