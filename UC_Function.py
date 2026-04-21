# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *
from abaqusConstants import *
import __main__
import numpy as np
import time

def UCModel(L, w_f, E1_FRP, E2_FRP, nu12_FRP, G12_FRP, G13_FRP, G23_FRP, rho_FRP, rho_m, C10_m, D1_m, rho_t, E_t, nu_t, t_t, t_FRP, layup, meshSize, prestress, uz_pull, cpus, job_id):
    """
    UC Model function for Abaqus simulation.
    
    Returns:
        radius (float): Fitted cylinder radius
        has_inflection (bool): Whether the geometry has inflection points (not tristable if True)
    """
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import odbAccess
    
    # -------------------------------------------------------------------------------------
    # MATERIAL CREATION
    # -------------------------------------------------------------------------------------
   
   # FRP - composite lamina
    mdb.models['Model-1'].Material(name='FRP')
    mdb.models['Model-1'].materials['FRP'].Density(table=((rho_FRP, ), ))
    mdb.models['Model-1'].materials['FRP'].Elastic(type=LAMINA, table=((E1_FRP, 
        E2_FRP, nu12_FRP, G12_FRP, G13_FRP, G23_FRP), ))
    # membrane - hyperelastic, neo hookean
    mdb.models['Model-1'].Material(name='Membrane')
    mdb.models['Model-1'].materials['Membrane'].Density(table=((rho_m, ), ))
    mdb.models['Model-1'].materials['Membrane'].Hyperelastic(
        materialType=ISOTROPIC, testData=OFF, type=NEO_HOOKE, 
        volumetricResponse=VOLUMETRIC_DATA, table=((C10_m, D1_m), ))
    # tape - linear elastic
    mdb.models['Model-1'].Material(name='Tape')
    mdb.models['Model-1'].materials['Tape'].Density(table=((rho_t, ), ))
    mdb.models['Model-1'].materials['Tape'].Elastic(table=((E_t, nu_t), ))
    
    # -------------------------------------------------------------------------------------
    # FRP FRAME PART
    # -------------------------------------------------------------------------------------
    
    # create frame geometry
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(-0.5*L, -0.5*L), point2=(0.5*L, 0.5*L))
    s.rectangle(point1=(-0.5*L + w_f, -0.5*L + w_f), point2=(0.5*L - w_f, 0.5*L - w_f))
    p = mdb.models['Model-1'].Part(name='frame', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['frame']
    p.BaseShell(sketch=s)
    
    # partition frame
    f, e = p.faces, p.edges
    edges = e.getByBoundingBox(0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchUpEdge=edges[0], 
        sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=86.63, gridSpacing=2.16, transform=t)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    s1.Line(point1=(-0.5*L + w_f, -0.5*L), point2=(-0.5*L + w_f, -0.5*L + w_f))
    s1.Line(point1=(-0.5*L, -0.5*L + w_f), point2=(-0.5*L + w_f, -0.5*L + w_f))
    s1.Line(point1=(-0.5*L, 0.5*L - w_f), point2=(-0.5*L + w_f, 0.5*L - w_f))
    s1.Line(point1=(-0.5*L + w_f, 0.5*L), point2=(-0.5*L + w_f, 0.5*L - w_f))
    
    s1.Line(point1=(0.5*L - w_f, -0.5*L), point2=(0.5*L - w_f, -0.5*L + w_f))
    s1.Line(point1=(0.5*L, -0.5*L + w_f), point2=(0.5*L - w_f, -0.5*L + w_f))
    s1.Line(point1=(0.5*L, 0.5*L - w_f), point2=(0.5*L - w_f, 0.5*L - w_f))
    s1.Line(point1=(0.5*L - w_f, 0.5*L), point2=(0.5*L - w_f, 0.5*L - w_f))
    
    p = mdb.models['Model-1'].parts['frame']
    f = p.faces
    e = p.edges
    edges = e.getByBoundingBox(0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    p.PartitionFaceBySketch(sketchUpEdge=edges[0], faces=f[0], sketch=s1)
    
    # assign composite layups
    # note: this is done in the global coordinate system and so the midplane strains and curvatures are also output in the global coordinate system

    # Bottom Left Corner
    faces = f.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, -0.5*L + w_f + 0.1, -0.5*L + w_f + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Bottom_Left', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    layup_BL = layup + [angle + 90 for angle in layup]
    for i, angle in enumerate(layup_BL):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
            
    # Bottom Right Corner
    faces = f.getByBoundingBox(0.5*L - w_f - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, -0.5*L + w_f + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Bottom_Right', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    layup_BR = [angle + 90 for angle in layup] + layup 
    for i, angle in enumerate(layup_BR):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # Top Left Corner
    faces = f.getByBoundingBox(-0.5*L - 0.1, 0.5*L - w_f - 0.1, -0.1, -0.5*L + w_f + 0.1, 0.5*L + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Top_Left', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    layup_TL = [angle + 90 for angle in layup] + layup 
    for i, angle in enumerate(layup_TL):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
            
    # Top Right Corner
    faces = f.getByBoundingBox(0.5*L - w_f - 0.1, 0.5*L - w_f - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Top_Right', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    layup_TR = layup + [angle + 90 for angle in layup]
    for i, angle in enumerate(layup_TR):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # Bottom
    faces = f.getByBoundingBox(-0.5*L + w_f - 0.1, -0.5*L - 0.1, -0.1, 0.5*L -w_f + 0.1, -0.5*L + w_f + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Bottom', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    for i, angle in enumerate(layup):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
            
    # Top
    faces = f.getByBoundingBox(-0.5*L + w_f - 0.1, 0.5*L - w_f - 0.1, -0.1, 0.5*L -w_f + 0.1, 0.5*L + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Top', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    for i, angle in enumerate(layup):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # Left
    faces = f.getByBoundingBox(-0.5*L - 0.1, -0.5*L + w_f - 0.1, -0.1, -0.5*L + w_f + 0.1, 0.5*L - w_f + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Left', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    for i, angle in enumerate(layup):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle+90, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
            
    # Right
    faces = f.getByBoundingBox(0.5*L - w_f - 0.1, -0.5*L + w_f - 0.1, -0.1, 0.5*L + 0.1, 0.5*L - w_f + 0.1, 0.1)
    region = regionToolset.Region(faces=faces)
    
    compositeLayup = mdb.models['Model-1'].parts['frame'].CompositeLayup(
        name='Right', description='', elementType=SHELL, 
        offsetType=MIDDLE_SURFACE, symmetric=False, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    compositeLayup.CompositePly(suppressed=False, plyName='Tape', region=region, 
        material='Tape', thicknessType=SPECIFY_THICKNESS, thickness=t_t, 
        orientationType=SPECIFY_ORIENT, orientationValue=0.0, 
        additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        axis=AXIS_3, angle=0.0, numIntPoints=3)
    
    # create corner overlap layup
    for i, angle in enumerate(layup):
        ply_name = f'Ply-{i+1}'
    
        compositeLayup.CompositePly(suppressed=False, plyName=ply_name, region=region, 
            material='FRP', thicknessType=SPECIFY_THICKNESS, thickness=t_FRP, 
            orientationType=SPECIFY_ORIENT, orientationValue=angle+90, 
            additionalRotationType=ROTATION_NONE, additionalRotationField='', 
            axis=AXIS_3, angle=0.0, numIntPoints=3)
            
    # Mesh Frame Part
    p.seedPart(size=meshSize, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=S4, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['frame']
    f = p.faces
    faces = f.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p.generateMesh()
    
    # -------------------------------------------------------------------------------------
    # MEMBRANE PART
    # -------------------------------------------------------------------------------------
    
    # create membrane part
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(-0.5*L, -0.5*L), point2=(0.5*L, 0.5*L))
    p = mdb.models['Model-1'].Part(name='membrane', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['membrane']
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['membrane']

    # create partition
    f, e = p.faces, p.edges
    edges = e.getByBoundingBox(0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    t = p.MakeSketchTransform(sketchPlane=f[0], sketchUpEdge=edges[0], 
        sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=84.85, gridSpacing=2.12, transform=t)
    s1.setPrimaryObject(option=SUPERIMPOSE)
    
    p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
    
    s1.Line(point1=(-0.5*L + w_f, -0.5*L), point2=(-0.5*L + w_f, 0.5*L))
    s1.Line(point1=(0.5*L - w_f, -0.5*L), point2=(0.5*L - w_f, 0.5*L))
    s1.Line(point1=(-0.5*L, -0.5*L + w_f), point2=(0.5*L, -0.5*L + w_f))
    s1.Line(point1=(-0.5*L, 0.5*L - w_f), point2=(0.5*L, 0.5*L - w_f))

    
    p.PartitionFaceBySketch(sketchUpEdge=edges[0], faces=f[0], sketch=s1)
    s1.unsetPrimaryObject()
    
    # create and assign membrane section
    mdb.models['Model-1'].MembraneSection(name='Membrane', material='Membrane', 
        thicknessType=UNIFORM, thickness=0.025, thicknessField='', 
        poissonDefinition=DEFAULT)
    
    p = mdb.models['Model-1'].parts['membrane']
    f = p.faces
    faces = f.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    region = p.Set(faces=faces, name='Set-1')
    p.SectionAssignment(region=region, sectionName='Membrane', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    
    # mesh membrane using membrane elements, reduced integration
    p.seedPart(size=meshSize, deviationFactor=0.1, minSizeFactor=0.1)
    elemType1 = mesh.ElemType(elemCode=M3D4R, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=M3D3, elemLibrary=STANDARD)

    faces = f.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    pickedRegions =(faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))

    p.generateMesh()
    
    # -------------------------------------------------------------------------------------
    # ASSEMBLY CREATION
    # -------------------------------------------------------------------------------------
    # create instances
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['frame']
    a.Instance(name='frame-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['membrane']
    a.Instance(name='membrane-1', part=p, dependent=ON)
    
    # Tie constraint between frame and membrane
    a = mdb.models['Model-1'].rootAssembly

    s1 = a.instances['frame-1'].faces
    side2Faces1 = s1.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    region1=a.Surface(side2Faces=side2Faces1, name='frame_bottom_surf')
    s1 = a.instances['membrane-1'].faces
    side1Faces1 = s1.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    region2=a.Surface(side1Faces=side1Faces1, name='membrane_top_surf')
    mdb.models['Model-1'].Tie(name='Tie', main=region1, secondary=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, 
        thickness=ON)
        
    # membrane prestretch - equal biaxial prestress
    f = a.instances['membrane-1'].faces
    faces = f.getByBoundingBox(-0.5*L - 0.1, -0.5*L - 0.1, -0.1, 0.5*L + 0.1, 0.5*L + 0.1, 0.1)
    region = a.Set(faces=faces, name='Set-1')
    mdb.models['Model-1'].Stress(name='Biaxial_prestress', region=region, 
        distributionType=UNIFORM, sigma11=prestress, sigma22=prestress, sigma12=0.0, 
        sigma33=None, sigma13=None, sigma23=None)
        
    # -------------------------------------------------------------------------------------
    # ANALYSIS SETUP
    # -------------------------------------------------------------------------------------
    # NOTE: step setup (incrementation, damping, etc) is currently hard coded
    
    # initial BCs
    a = mdb.models['Model-1'].rootAssembly
    v = a.instances['frame-1'].vertices
    verts1 = v.getByBoundingBox(0.5*L - w_f - 0.1, 0.5*L - 0.1, -0.1, 0.5*L - w_f + 0.1, 0.5*L + 0.1, 0.1)
    verts2 = v.getByBoundingBox(0.5*L - 0.1, 0.5*L - w_f - 0.1, -0.1, 0.5*L + 0.1, 0.5*L - w_f + 0.1, 0.1)
    region = a.Set(vertices=verts1+verts2, name='fixedZNodes')
    mdb.models['Model-1'].DisplacementBC(name='FixZ', createStepName='Initial', 
        region=region, u1=UNSET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)

    verts1 = v.getByBoundingBox(-0.5*L - 0.1, -0.5*L + w_f - 0.1, -0.1, -0.5*L + 0.1, -0.5*L + w_f + 0.1, 0.1)
    verts2 = v.getByBoundingBox(-0.5*L + w_f - 0.1, -0.5*L - 0.1, -0.1, -0.5*L + w_f + 0.1, -0.5*L + 0.1, 0.1)
    region = a.Set(vertices=verts1 + verts2, name='pinnedZNodes')
    mdb.models['Model-1'].DisplacementBC(name='PinnedBC', createStepName='Initial', 
        region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=SET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
   
    # Shape forming step, implicit dynamic quasi-static NLGEOM
    mdb.models['Model-1'].ImplicitDynamicsStep(name='ShapeForming', 
        previous='Initial', maxNumInc=10000, application=QUASI_STATIC, 
        initialInc=0.1, minInc=1e-09, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, 
        initialConditions=OFF, nlgeom=ON)

    # pull corners to transition to cylindrical step, static general, automatic stabilization
    mdb.models['Model-1'].StaticStep(name='pullCorners', previous='ShapeForming', 
        maxNumInc=10000, stabilizationMagnitude=2e-06, 
        stabilizationMethod=DISSIPATED_ENERGY_FRACTION, 
        continueDampingFactors=False, adaptiveDampingRatio=0.01, 
        initialInc=0.01, minInc=1e-09)
        
    verts1 = v.getByBoundingBox(-0.5*L + w_f - 0.1, -0.5*L + w_f - 0.1, -0.1, -0.5*L + w_f + 0.1, -0.5*L + w_f + 0.1, 0.1)
    verts2 = v.getByBoundingBox(0.5*L - w_f - 0.1, 0.5*L - w_f - 0.1, -0.1, 0.5*L - w_f + 0.1, 0.5*L - w_f + 0.1, 0.1)
    region = a.Set(vertices=verts1+verts2, name='pullingNodes')
    
    mdb.models['Model-1'].DisplacementBC(name='pullCorners', 
        createStepName='pullCorners', region=region, u1=UNSET, u2=UNSET, 
        u3=uz_pull, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
        distributionType=UNIFORM, fieldName='', localCsys=None)
        
    # release step to find stable shape, implicit dynamic quasi-static
    mdb.models['Model-1'].ImplicitDynamicsStep(name='Release', 
        previous='pullCorners', maxNumInc=10000, application=QUASI_STATIC, 
        initialInc=0.1, minInc=1e-09, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, 
        initialConditions=OFF, nlgeom=ON)
        
    mdb.models['Model-1'].boundaryConditions['pullCorners'].deactivate('Release')
    
    # output requests
    mdb.models['Model-1'].FieldOutputRequest(name='F-Output-2', 
        createStepName='ShapeForming', variables=('S', 'SE', 'U'), frequency=1)
    mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-2', 
        createStepName='ShapeForming', variables=('ALLAE', 'ALLIE', 'ALLKE', 
        'ALLSE', 'ALLSD', 'ALLWK', 'ETOTAL'), frequency=1)
    
    # create analysis job
    job_name = f'Job_{job_id}' #make unique job name based on prestress level
    mdb.Job(name=job_name, model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=cpus, numDomains=cpus, numGPUs=0)
        
    # run job and wait for completion
    my_job = mdb.jobs[job_name]
    start_time = time.time()
    my_job.submit()
    my_job.waitForCompletion()
    elapsed_time = time.time() - start_time
    
    # -------------------------------------------------------------------------------------
    # EXTRACT RESULTS
    # -------------------------------------------------------------------------------------
    
    # -------------------------------------------------------------------------------------
    # JOB STATUS CHECK
    # -------------------------------------------------------------------------------------
    
    # open odb
    odb_path = job_name + '.odb'
    odb = odbAccess.openOdb(path=odb_path)
    
    # Check if simulation reached completion in the 'Release' step
    if 'Release' in odb.steps and len(odb.steps['Release'].frames) > 0:
        release_step = odb.steps['Release']
        last_frame = release_step.frames[-1]
        simulation_complete = last_frame.frameValue >= 1.0
    else: 
        simulation_complete = False 
    
    # -------------------------------------------------------------------------------------
    # DATA EXTRACTION: ONLY IF SIMULATION CONVERGED
    # -------------------------------------------------------------------------------------
    if simulation_complete:
        # -------------------------------------------------------------------------------------
        # ENERGY CHECK FOR VALIDITY OF RESULTS
        # -------------------------------------------------------------------------------------
        
        # Get history region (whole model)
        step_key = odb.steps.keys()[-1]
        hist_region = odb.steps[step_key].historyRegions['Assembly ASSEMBLY']
        
        # Extract energy data
        kinetic_energy = hist_region.historyOutputs['ALLKE'].data
        strain_energy = hist_region.historyOutputs['ALLSE'].data
        artificial_energy = hist_region.historyOutputs['ALLAE'].data
        static_dissipation = hist_region.historyOutputs['ALLSD'].data
        
        # Convert to arrays
        ke_values = np.array([point[1] for point in kinetic_energy])
        se_values = np.array([point[1] for point in strain_energy])
        ae_values = np.array([point[1] for point in artificial_energy])
        sd_values = np.array([point[1] for point in static_dissipation])
        
        # Calculate ratios at each time step
        # Avoid division by zero
        se_nonzero = np.where(se_values != 0, se_values, 1e-10)
        
        ke_se_ratio = ke_values / se_nonzero
        ae_se_ratio = ae_values / se_nonzero
        sd_se_ratio = sd_values / se_nonzero
        
        # Find maximum ratios
        max_ke_ratio = np.max(ke_se_ratio)
        max_ae_ratio = np.max(ae_se_ratio)
        max_sd_ratio = np.max(sd_se_ratio)
        
        # -------------------------------------------------------------------------------------
        # TRISTABILITY CHECK
        # -------------------------------------------------------------------------------------
        # Currently done by looking for inflection points on frame strips
            
        # Get last frame from Release step
        last_frame = odb.steps['Release'].frames[-1]

        # Get displacement field for specific instance
        instance = odb.rootAssembly.instances['FRAME-1']
        U = last_frame.fieldOutputs['U'].getSubset(region=instance)

        # Check for inflection points at 4 boundaries
        conditions = [
            ('x', -0.5*L),
            ('x', 0.5*L),
            ('y', -0.5*L),
            ('y', 0.5*L)
        ]

        has_inflection = False

        for coord_name, coord_value in conditions:
            coord_idx = 0 if coord_name == 'x' else 1
            
            # Find nodes at this boundary and extract deformed coordinates
            deformed_coords = []
            for value in U.values:
                node = instance.nodes[value.nodeLabel - 1]  # nodeLabel is 1-indexed
                if abs(node.coordinates[coord_idx] - coord_value) < 1e-2:
                    deformed = [node.coordinates[i] + value.data[i] for i in range(3)]
                    deformed_coords.append(deformed)
            
            # Check for inflection if we have enough points
            if len(deformed_coords) > 3:
                deformed_coords = np.array(deformed_coords)
                
                # Sort by z-coordinate
                sort_idx = np.argsort(deformed_coords[:, 2])
                sorted_coords = deformed_coords[sort_idx]
                
                # Get the other coordinate
                z = sorted_coords[:, 2]
                other = sorted_coords[:, 1] if coord_name == 'x' else sorted_coords[:, 0]
                
                # Compute second derivative
                d_dz = np.gradient(other, z)
                d2_dz2 = np.gradient(d_dz, z)
                
                # Find inflection points
                sign_changes = np.diff(np.sign(d2_dz2))
                if np.any(sign_changes != 0):
                    has_inflection = True
                    break
        
        # -------------------------------------------------------------------------------------
        # CYLINDER FIT
        # -------------------------------------------------------------------------------------
        # cylinder of best fit - this is done even if geometry not tristable
        
        # Extract deformed coordinates at last frame of simulation
        last_step_key = odb.steps.keys()[-1]
        last_step = odb.steps[last_step_key]
        last_frame = last_step.frames[-1]
        
        # Get displacement field
        displacement_field = last_frame.fieldOutputs['U']
        
        # Get all instances in the assembly
        assembly = odb.rootAssembly
        
        # Store coordinates for cylinder fitting
        deformed_coords = []
        
        # Loop through all displacement values
        for value in displacement_field.values:
            node_label = value.nodeLabel
            
            # Get the instance name and node
            # The value.instance gives us the instance this node belongs to
            instance = value.instance
            
            # Access the node from the instance
            node = instance.nodes[node_label - 1]  # Node labels start at 1, index at 0
            
            # Get original coordinates
            x_orig = node.coordinates[0]
            y_orig = node.coordinates[1]
            z_orig = node.coordinates[2]
            
            # Get displacements
            u1 = value.data[0]
            u2 = value.data[1]
            u3 = value.data[2]
            
            # Calculate deformed coordinates
            x_def = x_orig + u1
            y_def = y_orig + u2
            z_def = z_orig + u3
            
            # Store for cylinder fitting
            deformed_coords.append([x_def, y_def, z_def])
        
        # deformed_coords is assumed to exist and have shape (N, 3)
        points = np.asarray(deformed_coords, dtype=float)
        N = points.shape[0]

        # Axis direction via covariance PCA 
        centroid = points.mean(axis=0)
        X = points - centroid
        cov = np.dot(X.T, X) / N

        eigvals, eigvecs = np.linalg.eigh(cov)
        axis_dir = eigvecs[:, np.argmax(eigvals)]
        axis_dir /= np.linalg.norm(axis_dir)

        # Build orthonormal basis for cross-section plane
        # Choose any vector not parallel to axis
        ref = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(ref, axis_dir)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])

        e1 = np.cross(axis_dir, ref)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(axis_dir, e1)
        
        # Project points onto plane normal to axis
        proj = points - np.outer(np.dot(points - centroid, axis_dir), axis_dir)

        x2d = np.dot(proj - centroid, e1)
        y2d = np.dot(proj - centroid, e2)

        # Algebraic circle fit (Taubin-style, NumPy-only)
        x = x2d
        y = y2d

        x_m = x.mean()
        y_m = y.mean()

        u = x - x_m
        v = y - y_m

        Su2 = np.sum(u**2)
        Sv2 = np.sum(v**2)
        Suv = np.sum(u*v)
        Su3 = np.sum(u**3)
        Sv3 = np.sum(v**3)
        Su2v = np.sum(u**2 * v)
        Suv2 = np.sum(u * v**2)

        A = np.array([[Su2, Suv],
                      [Suv, Sv2]])

        B = 0.5 * np.array([Su3 + Suv2,
                            Sv3 + Su2v])

        uc, vc = np.linalg.solve(A, B)

        xc = x_m + uc
        yc = y_m + vc

        radius = np.sqrt((u - uc)**2 + (v - vc)**2).mean()
        
        axis_point = centroid + xc * e1 + yc * e2

        # Goodness-of-fit
        dist = np.sqrt((x - xc)**2 + (y - yc)**2)
        residuals = dist - radius

        rmse = np.sqrt(np.mean(residuals**2))
        mae = np.mean(np.abs(residuals))
        max_error = np.max(np.abs(residuals))
    
    else:
        # If simulation did not complete, return None values
        radius = None
        has_inflection = None
    
    # Close ODB
    odb.close()
    
    # Return the two requested values
    return radius, has_inflection