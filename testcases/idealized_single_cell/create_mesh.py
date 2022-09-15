from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def create_mesh():

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    fileGrid = Dataset("grid_in.nc","w",format="NETCDF3_CLASSIC")

    fileGrid.on_a_sphere = "NO"

    nCells = 1
    nVertices = 4
    vertexDegree = 4

    fileGrid.createDimension("nCells", nCells)
    fileGrid.createDimension("nVertices", nVertices)
    fileGrid.createDimension("vertexDegree", vertexDegree)

    xCell = np.zeros(nCells)
    yCell = np.zeros(nCells)
    zCell = np.zeros(nCells)

    xVertex = np.zeros(nVertices)
    yVertex = np.zeros(nVertices)
    zVertex = np.zeros(nVertices)
    cellsOnVertex = np.zeros((nVertices,vertexDegree),dtype="i")
    cellsOnVertex[:] = -1

    xVertex[0] = -0.5
    yVertex[0] = -0.5
    cellsOnVertex[0,0] = 0

    xVertex[1] =  0.5
    yVertex[1] = -0.5
    cellsOnVertex[1,0] = 0

    xVertex[2] =  0.5
    yVertex[2] =  0.5
    cellsOnVertex[2,0] = 0

    xVertex[3] = -0.5
    yVertex[3] =  0.5
    cellsOnVertex[3,0] = 0

    var = fileGrid.createVariable("xCell","d",dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileGrid.createVariable("yCell","d",dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileGrid.createVariable("zCell","d",dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileGrid.createVariable("xVertex","d",dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileGrid.createVariable("yVertex","d",dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileGrid.createVariable("zVertex","d",dimensions=["nVertices"])
    var[:] = zVertex[:]

    cellsOnVertex[:] += 1
    var = fileGrid.createVariable("cellsOnVertex","i",dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileGrid.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.nc grid.nc" %(mpas_tools_dir))

    fileGrid = Dataset("grid.nc","a")

    fileGrid.variables["latCell"][:]   = math.radians(80.0)
    fileGrid.variables["latVertex"][:] = math.radians(80.0)
    fileGrid.variables["latEdge"][:]   = math.radians(80.0)

    fileGrid.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_mesh()
