# This file contains the functions reading in the data of a VTU file.
# Author: Dave May

#======================= import essential libraries===========================
import numpy as np
import vtk as vtk
#=========================== Animation of field data ==============================
def vtu_extract_element_connectivity(fname):
    
    """
    Loads a vtu file, extracts the "connectivity" field and
    converts the data into a NumPy array (int64).
    The element-vertex map is returned as a 2D ndarray.
    It is assumed that all cells are of the same VTK type.
    The "offset" field in the vtu file is processed and checked to ensure that 
    all cells are the same VTK type (e.g. triangle).
    
    Parameters:
    -----------
    fname : string
      Filename of a VTU file      
    Returns:
    --------
    elmap : ndarray, shape = (ncells, npoints_per_cell)
        Element-vertex map associated with VTK cells, with dimensions (nCells, npoints_per_cell)
    """  

    # Read the source vtu file.
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update() # Needed because of GetScalarRange
    output = reader.GetOutput()

    cells = output.GetCells()

    offset_per_cell = cells.IsHomogeneous()
    if offset_per_cell <= 0:
        raise ValueError('Require all cells have the same size (homogeneous).')

    nCells = cells.GetNumberOfCells()
    array = cells.GetConnectivityArray()
    elmap = np.asarray(array ,dtype=np.int64)
    elmap = elmap.reshape(nCells, offset_per_cell)
    return elmap

def vtu_extract_fields(fname, extract_coor=True):
    """
    Extract point fields, cell fields and point coordiantes from a vtu file.
    
    Parameters:
    -----------
    extract_coor : boolean
        Flag to indicate whether you want the coordinates to be extracted.
    
    Returns:
    --------
    point_field : dict
      All point fields found. Key is the name provided in the VTU file.
    cell_field : dict
      All cell fields found. Key is the name provided in the VTU file.
    coor_field : dict
      Coordinates of the mesh. `coor_field` be empty if extract_coor = False.
      Key used is "coor".
    """
    
    point_field = {}
    cell_field = {}
    coor_field = {}
    
    # Read the source vtu file.
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update() # Needed because of GetScalarRange
    output = reader.GetOutput()

    points = output.GetPoints()
    array = points.GetData()
    xyz = np.asarray(array, dtype=np.float64)
    if extract_coor:
        coor_field["coor"] = xyz
    
    pointFields = output.GetPointData()
    
    nfields = pointFields.GetNumberOfArrays()
    print('#point fields', nfields)
    for f in range(nfields):
        print('  pointFields.name', pointFields.GetArrayName(f))
        #print('pointFields.ndof', pointFields.GetNumberOfComponents(f))
        array = pointFields.GetArray(f)
        pf = np.asarray(array, dtype=np.float64)
        point_field[ pointFields.GetArrayName(f) ] = pf
        
    #cells = output.GetCells()
    cellFields = output.GetCellData()
    nfields = cellFields.GetNumberOfArrays()
    print('#cell fields', nfields)
    for f in range(nfields):
        print('  cellFields.name', cellFields.GetArrayName(f))
        array = cellFields.GetArray(f)
        cf = np.asarray(array, dtype=np.float64)
        cell_field[ cellFields.GetArrayName(f) ] = cf
        
    return point_field, cell_field, coor_field