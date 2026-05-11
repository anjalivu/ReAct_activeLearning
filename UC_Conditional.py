# -*- coding: mbcs -*-
# Do not delete the following import lines
from tkinter import OFF, ON
from unittest.mock import DEFAULT

from abaqus import *
from abaqusConstants import *
from helper_functions import *
import __main__
import numpy as np
import time

def UC_conditional(L, w_f, E1_FRP, E2_FRP, nu12_FRP, G12_FRP, G13_FRP, G23_FRP, rho_FRP, rho_m, C10_m, D1_m, rho_t, E_t, nu_t, t_t, t_FRP, layup, meshSize, prestress, uz_pull, cpus, job_name='Job-1'):
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
    import os
    import glob
    
    # =========================================================================
    # CLEANUP SECTION
    # =========================================================================
    print("="*60)
    print("CLEANING UP FROM PREVIOUS RUNS")
    print("="*60)
    
    # 1. Clean up files
    print("Deleting old job files...")
    patterns = [f'{job_name}.*', f'{job_name}_pullCorners.*']
    files_deleted = 0
    
    for pattern in patterns:
        for filepath in glob.glob(pattern):
            try:
                # Skip .cae files (ABAQUS model database)
                if not filepath.endswith('.cae'):
                    os.remove(filepath)
                    files_deleted += 1
            except Exception as e:
                print(f"Warning: Could not delete {filepath}: {str(e)}")
    
    print(f"Deleted {files_deleted} old files")
    
    # 2. Clean up jobs
    print("Deleting old jobs...")
    if job_name in mdb.jobs.keys():
        del mdb.jobs[job_name]
        print(f"Deleted job: {job_name}")
    
    if job_name + '_pullCorners' in mdb.jobs.keys():
        del mdb.jobs[job_name + '_pullCorners']
        print(f"Deleted job: {job_name}_pullCorners")
    
    # 3. Clean up model
    try:
        del mdb.models['Model-1']
        print("Deleted existing Model-1")
    except KeyError:
        print("Model-1 did not exist (first run)")
        
    mdb.Model(name='Model-1', modelType=STANDARD_EXPLICIT)
    
    print("Cleanup complete. Fresh Model-1 created")
    
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
        
    print('Model creation complete. Starting analysis...')
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
    
    print('Initial boundary conditions applied. Starting shape forming step...')

    # Shape forming step, implicit dynamic quasi-static NLGEOM
    mdb.models['Model-1'].ImplicitDynamicsStep(name='ShapeForming', 
        previous='Initial', maxNumInc=10000, application=QUASI_STATIC, 
        initialInc=0.1, minInc=1e-09, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, 
        initialConditions=OFF, nlgeom=ON)
    






    # # output requests
    # mdb.models['Model-1'].FieldOutputRequest(name='F-Output-2', 
    #     createStepName='ShapeForming', variables=('S', 'SE', 'U'), frequency=1)
    # mdb.models['Model-1'].HistoryOutputRequest(name='H-Output-2', 
    #     createStepName='ShapeForming', variables=('ALLAE', 'ALLIE', 'ALLKE', 
    #     'ALLSE', 'ALLSD', 'ALLWK', 'ETOTAL'), frequency=1)
    
    # obtain displacement
    mdb.models['Model-1'].FieldOutputRequest(name='Step-1-Displacement', 
            createStepName='ShapeForming', 
            variables=('U',),
              frequency=1)
    
    # create analysis job
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
    odb1 = odbAccess.openOdb(path=odb_path)
    odb_curr = odb1 

    odb2_exists = False
    has_inflection = None
    radius = None
    pull_corners_needed = False


    # Check if simulation reached completion in the 'ShapeForming' step
    if 'ShapeForming' in odb1.steps and len(odb1.steps['ShapeForming'].frames) > 0:
        shapeforming_step = odb1.steps['ShapeForming']
        last_frame = shapeforming_step.frames[-1]
        simulation_complete = last_frame.frameValue >= 1.0
        
    else: 
        simulation_complete = False 



    print(f"Step 1: Simulation complete: {simulation_complete}")
    # -------------------------------------------------------------------------------------
    # DATA EXTRACTION: ONLY IF SIMULATION CONVERGED
    # -------------------------------------------------------------------------------------
    if simulation_complete:
        # -------------------------------------------------------------------------------------
        # TRISTABILITY CHECK
        # -------------------------------------------------------------------------------------
            
        # Get last frame from Shapeforming step
        last_frame = odb1.steps['ShapeForming'].frames[-1]

        # Get displacement field for specific instance
        instance = odb1.rootAssembly.instances['FRAME-1']
        U = last_frame.fieldOutputs['U'].getSubset(region=instance)

        # Check for inflection points at 4 boundaries
        conditions = [
            ('x', -0.5*L),
            ('x', 0.5*L),
            ('y', -0.5*L),
            ('y', 0.5*L)
        ]

        has_inflection = check_inflection(U, conditions, instance)

        print(f"Step 1: Inflection point detected: {has_inflection}")
        


        if has_inflection:
            print("Performing pull_corners step")
            pull_corners_needed = True


            # add step 2
            # pull corners to transition to cylindrical step, implicit dynamic quasi-static
            mdb.models['Model-1'].ImplicitDynamicsStep(name='pullCorners', 
                previous='ShapeForming', maxNumInc=10000, application=QUASI_STATIC, 
                initialInc=0.01, minInc=1e-09, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, #CHANGED THIS FROM 0.1 TO 0.05 TO GO FAST
                initialConditions=OFF, nlgeom=ON)
                
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

            job_name2 = job_name + '_pullCorners'

            
            # obtain displacement
            mdb.models['Model-1'].FieldOutputRequest(name='Step-1-Displacement', 
                    createStepName='ShapeForming', 
                    variables=('U',),
                    frequency=1)
            
            # create analysis job
            mdb.Job(name=job_name2, model='Model-1', description='', type=ANALYSIS, 
                atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
                memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
                modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
                multiprocessingMode=DEFAULT, numCpus=cpus, numDomains=cpus, numGPUs=0)

            # run job and wait for completion
            my_job2 = mdb.jobs[job_name2]
            start_time = time.time()
            my_job2.submit()
            my_job2.waitForCompletion()
            elapsed_time = time.time() - start_time


            odb_path2 = job_name2 + '.odb'
            odb2 = odbAccess.openOdb(path=odb_path2)
            odb_curr = odb2
            odb2_exists = True


            # Get last frame from Release step
            last_frame = odb2.steps['Release'].frames[-1]

            # Get displacement field for specific instance
            instance = odb2.rootAssembly.instances['FRAME-1']
            U = last_frame.fieldOutputs['U'].getSubset(region=instance)

            # Check for inflection points at 4 boundaries
            conditions = [
                ('x', -0.5*L),
                ('x', 0.5*L),
                ('y', -0.5*L),
                ('y', 0.5*L)
            ]

            has_inflection = check_inflection(U, conditions, instance)

            print(f"Step 2: Inflection point detected: {has_inflection}")
            
        # -------------------------------------------------------------------------------------
        # CYLINDER FIT
        # -------------------------------------------------------------------------------------
        # Extract deformed coordinates at last frame of simulation



    if not has_inflection and simulation_complete:
        print("Performing cylinder fit on final shape")

        last_step_key = odb_curr.steps.keys()[-1]
        last_step = odb_curr.steps[last_step_key]
        last_frame = last_step.frames[-1]
        
        # Get displacement field
        displacement_field = last_frame.fieldOutputs['U']
        
        # Get all instances in the assembly
        assembly = odb_curr.rootAssembly
        
        
        # Loop through all displacement values
        radius, dist, residuals = fit_cylinder(displacement_field)
    
    else:
        # If simulation did not complete, return None values
        radius = None
    
    # Close ODBs
    odb1.close()

    if odb2_exists:
        odb2.close()
    
    # Return the two requested values
    return radius, has_inflection, pull_corners_needed