import numpy as np
from scipy.special import gamma as Gamma
import moorpy as mp
import copy
import pandas as pd
import yaml
from calcLineProps import getMoorpyLineType
import os



def spread_mooring(moor_dict, tol = 0.01, maxIter = 500, no_fail = True, finite_difference = False, includeBody = True):
    """
    Credit : Serag-Eldin Abdelmoteleb
    Creates a moorpy mooring system object with axis-symmetric, multi-segmented line configuration for multi-column platform.
    (N.B. The platform body type in this system is set to fixed in position so that no information about platform's mass or 
     buoyancy is required. The body type can be changed in subsequent analysis once the mass and buoyancy information are added.)

    Parameters
    ----------
    moor_dict : dictionary
        A dictionary with the following keys:
            nCols: int
                number of outer columns of the platform.
            legsPerCol: int
                number of mooring legs attached to a single column
            colAngles: list of floats
                oreintation of outer columns in degrees
            spreadAngle: float
                angle on which the legs in a single column are spread in degrees (0 if 1 leg per column)
            depth: float
                seabed depth
            rCol: float
                distance between an outer column center and platfomr center
            Dcol: float
                column diameter
            lF2A: float
                fairlead to anchor distance
            zFair: float
                fairlead vertical position w.r.t SWL

            NOTE : The line lengths for multi segmented lines has to start from anchor to the fairlead - This is the order which is 
            followed the the function.

            lineLengths: list of floats
                line lengths in meters
            segLengths: list of floats
                maximum length of the segment lengths in meters used for discretization
            lineDiameters: list of floats
                line nominal diameters in meters
            lineTypes: list of strings
                line material type (must match the name key of one of the types in the lineTypes input)
            materialDicts: list of dictionaries
                line type material dictionaries each with the following keys:
                    name: (str) type name. e.g. chain, wire etc.
                    massden: (float) linear density in kg/m
                    EA: (float) extensional stiffness in N/m
                    d: (float) volume equivalent diameter in m
                    MBL: (float) minimum breaking load in N
    includeBody : the body is not included in moordyn file if it is coupled

    Returns
    -------
    conv : bool
        moorpy convergence flag.
    ms : moorpy.System
        A moorpy mooring System object.
    """
    # TODO: make line discretization customizable
    # Read parameters from input dictionary
    n_legs = moor_dict['legsPerCol']
    col_angles = np.radians(moor_dict['colAngles'])
    if n_legs > 1:
        spread_angle =  np.radians(moor_dict['spreadAngle'])
    else:
        spread_angle = 0
    depth = moor_dict['depth']
    r_col = moor_dict['rCol']
    D_col = moor_dict['Dcol']
    l_f2a = moor_dict['lF2A']
    z_fair = moor_dict['zFair']
    line_lengths = moor_dict['lineLengths'] 
    seg_lengths = moor_dict['segLengths']
    line_types = moor_dict['lineTypes']
    line_types_names = moor_dict['lineTypesNames']
    L_tot = np.sum(line_lengths)
    point_masses = moor_dict['point_masses'] # a list of intermediate clump masses (must be 1 less than the number of lines)
    point_volumes = moor_dict['point_volumes'] # a list of intermediate buoy volumes (must be 1 less than the number of lines)
    
    # Check input data
    n_lines = len(line_lengths)
    np.testing.assert_equal(n_lines, len(line_types_names))
    np.testing.assert_equal(n_lines, len(seg_lengths))
    np.testing.assert_equal(n_lines-1, len(point_masses))
    np.testing.assert_equal(n_lines-1, len(point_volumes))

    # Create mooring system
    ms = mp.System(depth = depth)
    
    # Add line types
    ms.lineTypes = moor_dict['lineTypes']
    
    if includeBody:
        # Add a fixed body (to be changed to floating in subsequent analyses)
        ms.addBody(-1, r6 = [0]*6) #coupled body

    # Add connection points and lines
    for col, col_angle in enumerate(col_angles): # loop over platform's columns
        for n in range(n_legs): # loop over lines per column
            ## set line angle for evaluation of fairlead position    
            if n_legs > 1:
                line_angle = col_angle - spread_angle/2 + n*spread_angle/(n_legs-1) 
            else:
                line_angle = col_angle
            
            r_fairlead = np.array([r_col*np.cos(col_angle) + D_col/2*np.cos(line_angle), # fairlead x-coord
                                   r_col*np.sin(col_angle) + D_col/2*np.sin(line_angle), # fairlead y-coord
                                   z_fair])            # fairlead z-coord
            
            r_anchor = np.array([r_col*np.cos(col_angle) + (D_col/2+l_f2a)*np.cos(line_angle), # anchor x-coord
                                 r_col*np.sin(col_angle) + (D_col/2+l_f2a)*np.sin(line_angle), # anchor y-coord
                                 -depth]) # anchor z-coord
            
            ## add anchor point
            ms.addPoint(1, np.around(r_anchor,4))
            
            ## Add intermediate connection points (if any)
            r_0 = r_fairlead

            if n_lines >1:
                for line_len, line_type_name, seg_len, m, v in zip(line_lengths[:-1],line_types_names[:-1], seg_lengths[:-1], point_masses, point_volumes):
                    lt = line_types[line_type_name]
                    r_0 = r_0 + line_len/L_tot*(r_anchor-r_fairlead)
                    n_segs = max(2, int(line_len/seg_len)) #ensure that there is at least 2 segments
                    ms.addPoint(0, np.around(r_0,4), m = m, v = v)
                    ms.addLine(np.around(line_len,4), lt, nSegs=n_segs, pointA = len(ms.pointList)-1, pointB = len(ms.pointList))

            ## add fairlead point
            if includeBody:
                ms.addPoint(1, np.around(r_fairlead,4)) #fixed point 
                ms.bodyList[0].attachPoint(pointID = len(ms.pointList), rAttach = ms.pointList[len(ms.pointList)-1].r)
            else:
                ms.addPoint(-1, np.around(r_fairlead,4)) #fairlead point as coupled with no body
            
            n_segs = max(2, int(line_lengths[-1]/seg_lengths[-1]))
            ms.addLine(np.around(line_lengths[-1],4), line_types[line_types_names[-1]],nSegs = n_segs, pointA = len(ms.pointList)-1, pointB = len(ms.pointList))
    
    ms.initialize()
    conv = ms.solveEquilibrium(tol = tol, no_fail = no_fail, maxIter = maxIter, finite_difference = finite_difference)
    
    return conv,ms



def genMoorpySys2(morringYamlFile, bdyYamlFile):

    moorData = yaml.safe_load(open(morringYamlFile))   
    linTypList = moorData['linTypList'] #hold the line types
    moor_dict = moorData['moor_dict'] #holds the geometry info

    moorpyLtypeDict = {}
    for linTyp in linTypList:
        res = getMoorpyLineType(**linTyp)
        moorpyLtypeDict[res['name']] = res
    moor_dict['lineTypes'] = moorpyLtypeDict


    conv, ms = spread_mooring(moor_dict, tol = 0.01, maxIter = 500, no_fail = True, finite_difference = False, includeBody = True)

    bdyData = yaml.safe_load(open(bdyYamlFile))
    ms.bodyList[0].rCG = bdyData['r_cg']
    ms.bodyList[0].m = bdyData['m']
    ms.bodyList[0].v = bdyData['v']
    ms.bodyList[0].AWP = bdyData['AWP']
    ms.bodyList[0].rM = bdyData['rM'] #TODO: add horizontal location of cb from hydrodyn - this is uncorrected for the location
    ms.bodyList[0].type = -1 #coupled
    ms.initialize()
    
    return ms



