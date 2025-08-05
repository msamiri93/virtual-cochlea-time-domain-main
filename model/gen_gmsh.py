#!C:\Users\Mohammad Shokrian\AppData\Local\Programs\Python\Python39\python.exe
import argparse
import numpy as np
import gmsh

import scipy.io

def addNodes(outputList, lines, boundary):
    tags = []
    for line in lines: 
        node_tags, _, _ = gmsh.model.mesh.getNodes(1, line, True)
        tags.extend(list(set(list(node_tags))))
    outputList[boundary] = list(set(list(tags)))

def generate_mesh(args):

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.model.add("Model")

    # Parameters
    Outer_bound = args.BM_length-args.Chamber_height+args.Heli_gap
    FLD_height = args.Chamber_height-args.OOC_radi
    
    # Points
    gmsh.model.geo.addPoint(0, args.OOC_radi, 0, 1.0, 1)
    gmsh.model.geo.addPoint(0, -args.OOC_radi, 0, 1.0, 2)
    gmsh.model.geo.addPoint(0, args.Chamber_height, 0, 1.0, 3)
    gmsh.model.geo.addPoint(0, -args.Chamber_height, 0, 1.0, 4)
    
    gmsh.model.geo.addPoint(args.BM_length, args.OOC_radi, 0, 1.0, 5)
    gmsh.model.geo.addPoint(args.BM_length, -args.OOC_radi, 0, 1.0, 6)
    gmsh.model.geo.addPoint(Outer_bound, args.Chamber_height, 0, 1.0, 7)
    gmsh.model.geo.addPoint(Outer_bound, -args.Chamber_height, 0, 1.0, 8)

    gmsh.model.geo.addPoint(Outer_bound, args.OOC_radi, 0, 1.0, 9)
    gmsh.model.geo.addPoint(Outer_bound, 0, 0, 1.0, 10)

    gmsh.model.geo.addPoint(Outer_bound, -args.OOC_radi, 0, 1.0, 11)
    
    gmsh.model.geo.addPoint(args.BM_length, 0, 0, 1.0, 12)
    
    gmsh.model.geo.addPoint(Outer_bound + args.Chamber_height*np.cos(args.Split_angle), -args.Chamber_height*np.sin(args.Split_angle), 0, 1.0, 13)
    gmsh.model.geo.addPoint(Outer_bound + args.Chamber_height*np.cos(args.Split_angle), args.Chamber_height*np.sin(args.Split_angle), 0, 1.0, 14)
    gmsh.model.geo.addPoint(args.BM_length + args.Heli_gap, 0, 0, 1.0, 15)
    
    # Lines
    gmsh.model.geo.addLine(9, 7, 1)
    gmsh.model.geo.addLine(9, 5, 2)
    gmsh.model.geo.addLine(11, 6, 3)
    gmsh.model.geo.addLine(11, 8, 4)
    gmsh.model.geo.addLine(9, 1, 5)
    gmsh.model.geo.addLine(1, 3, 6)
    gmsh.model.geo.addLine(2, 4, 7)
    gmsh.model.geo.addLine(2, 11, 8)
    gmsh.model.geo.addLine(8, 4, 9)
    gmsh.model.geo.addLine(3, 7, 10)

    gmsh.model.geo.addCircleArc(7, 10, 14, 11)
    gmsh.model.geo.addCircleArc(14, 10, 13, 12)
    gmsh.model.geo.addCircleArc(13, 10, 8, 13)
    gmsh.model.geo.addCircleArc(6, 12, 5, 14)
    gmsh.model.geo.addCircleArc(6, 15, 13, 15)
    gmsh.model.geo.addCircleArc(5, 15, 14, 16)
    
    # Curve Loops and Plane Surfaces
    gmsh.model.geo.addCurveLoop([1, 11, -16, -2], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop([16, 12, -15, 14], 2)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addCurveLoop([13, -4, 3, 15], 3)
    gmsh.model.geo.addPlaneSurface([3], 3) 
    gmsh.model.geo.addCurveLoop([9, -7, 8, 4], 4)
    gmsh.model.geo.addPlaneSurface([4], 4)  
    gmsh.model.geo.addCurveLoop([5, 6, 10, -1], 5)
    gmsh.model.geo.addPlaneSurface([5], 5)

    # Transfinite Curves and Surfaces
    for curve in [5,8,9,10]: # main length
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, int(Outer_bound/args.dz)+1, "Progression", 1)

    for curve in [2,3,13,11]: # apical length
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, int((args.Chamber_height-args.Heli_gap)/args.dz)+1, "Progression", 1)
    
    for curve in [4,1,7,6]: # main depth
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, int(FLD_height/args.dh)+1, "Progression", args.depth_scaling)

    for curve in [15,16]: # apical depth
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, int(FLD_height/args.dh)+1, "Progression", args.heli_scaling)

    for curve in [14,12]: # helicotrema
        gmsh.model.geo.mesh.setTransfiniteCurve(curve, int(np.pi*args.OOC_radi/(args.dz*0.666))+1, "Progression", 1)

    
    gmsh.model.geo.mesh.setTransfiniteSurface(1)
    gmsh.model.geo.mesh.setTransfiniteSurface(2)
    gmsh.model.geo.mesh.setTransfiniteSurface(3)
    gmsh.model.geo.mesh.setTransfiniteSurface(4)
    gmsh.model.geo.mesh.setTransfiniteSurface(5)

    gmsh.model.geo.mesh.setRecombine(2, 1)
    gmsh.model.geo.mesh.setRecombine(2, 2)
    gmsh.model.geo.mesh.setRecombine(2, 3)
    gmsh.model.geo.mesh.setRecombine(2, 4)
    gmsh.model.geo.mesh.setRecombine(2, 5)

    gmsh.model.geo.synchronize()  
    gmsh.model.removeEntities([(0, 10), (0, 12), (0, 15)])

    # Synchronize and Generate Mesh
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(2)
 
    outputList = {}
    addNodes(outputList, [2, 5], 'TOP_SURFACE')
    addNodes(outputList, [3, 8], 'BOTTOM_SURFACE')
    addNodes(outputList, [6], 'OW_SURFACE')
    addNodes(outputList, [7], 'RW_SURFACE')

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()

    # Store nodes in an array (node_id, [x, y, z])
    nodes = [ list(node_coords[i*3:(i+1)*3]) for i in range(len(node_tags))]

    
    elementIDs, nodeIds = gmsh.model.mesh.getElementsByType(3)

    nodeIds = list(np.array(nodeIds).reshape((len(elementIDs), 4)))

    elements = [list(nodeids) for (eid, nodeids) in zip(elementIDs, nodeIds)]

    mexport = {
        'boundaries': outputList,
        'nodes': nodes,
        'elements': elements,
    }
    scipy.io.savemat(f'{args.file_name}.mat', mexport)

    #gmsh.fltk.run()
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser();
    parser.add_argument('--file_name',type=str,required=True);
    parser.add_argument('--depth_scaling',type=float,default=1);
    parser.add_argument('--heli_scaling',type=float,default=1);
    parser.add_argument('--dz',type=float,default=10);
    parser.add_argument('--dh',type=float,default=10);
    parser.add_argument('--BM_length',type=float,default=12000);
    parser.add_argument('--Chamber_height',type=float,default=300);
    parser.add_argument('--OOC_radi',type=float,default=25);
    parser.add_argument('--Heli_gap',type=float,default=100);
    parser.add_argument('--Split_angle',type=float,default=np.pi/9);
    args = parser.parse_args();

    generate_mesh(args)