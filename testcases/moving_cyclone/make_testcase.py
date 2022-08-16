from netCDF4 import Dataset
import math
import os
import netCDF4
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import matplotlib.cm as cmm
from mpl_toolkits.axes_grid1 import make_axes_locatable

lx = 512000.0
ly = 512000.0

#-------------------------------------------------------------------------------

def create_grid_quad(gridfileOut, dc):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    dx = dc
    dy = dc

    nx = int(lx / dx)
    ny = int(ly / dy)

    vertexDegree = 4

    fileGrid = Dataset("grid_in.nc","w",format="NETCDF3_CLASSIC")

    fileGrid.on_a_sphere = "NO"

    nCells = nx * ny
    nVertices = (nx+1) * (ny+1)

    fileGrid.createDimension("nCells", nCells)
    fileGrid.createDimension("nVertices", nVertices)
    fileGrid.createDimension("vertexDegree", vertexDegree)

    xCell = np.zeros(nCells)
    yCell = np.zeros(nCells)
    zCell = np.zeros(nCells)

    for j in range(0,ny):
        for i in range(0,nx):
            iCell = i + j * nx
            xCell[iCell] = (float(i) + 0.5) * dc
            yCell[iCell] = (float(j) + 0.5) * dc

    xVertex = np.zeros(nVertices)
    yVertex = np.zeros(nVertices)
    zVertex = np.zeros(nVertices)
    cellsOnVertex = np.zeros((nVertices,vertexDegree),dtype="i")

    for j in range(0,ny+1):
        for i in range(0,nx+1):
            iVertex = i + j * (nx+1)
            xVertex[iVertex] = float(i) * dc
            yVertex[iVertex] = float(j) * dc

            ic = i - 1 ; jc = j - 1
            iCell = ic + jc * nx
            if (ic >= 0 and jc >= 0):
                cellsOnVertex[iVertex,0] = iCell
            else:
                cellsOnVertex[iVertex,0] = -1

            ic = i ; jc = j - 1
            iCell = ic + jc * nx
            if (ic < nx and jc >= 0):
                cellsOnVertex[iVertex,1] = iCell
            else:
                cellsOnVertex[iVertex,1] = -1

            ic = i ; jc = j
            iCell = ic + jc * nx
            if (ic < nx and jc < ny):
                cellsOnVertex[iVertex,2] = iCell
            else:
                cellsOnVertex[iVertex,2] = -1

            ic = i - 1 ; jc = j
            iCell = ic + jc * nx
            if (ic >= 0 and jc < ny):
                cellsOnVertex[iVertex,3] = iCell
            else:
                cellsOnVertex[iVertex,3] = -1

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

    cellsOnVertex[:] = cellsOnVertex[:] + 1
    var = fileGrid.createVariable("cellsOnVertex","i",dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileGrid.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.nc %s" %(mpas_tools_dir,gridfileOut))

#-------------------------------------------------------------------------------

def create_grid_hex(gridfileOut, dc):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    dx = dc
    dy = dc * math.sqrt(0.75)

    nx = int(lx / dx)
    if (nx % 2 != 0):
        nx += 1
    ny = int(ly / dy)
    if (ny % 2 != 0):
        ny += 1

    fileout = open("namelist.input","w")

    line = "&periodic_grid\n"
    fileout.write(line)

    line = "   nx = %i,\n" %(nx+2)
    fileout.write(line)

    line = "   ny = %i,\n" %(ny+2)
    fileout.write(line)

    line = "   dc = %f,\n" %(dc)
    fileout.write(line)

    line = "   nVertLevels = 1,\n"
    fileout.write(line)

    line = "   nTracers = 1,\n"
    fileout.write(line)

    line = "   nproc = 1,\n"
    fileout.write(line)

    line = "/\n"
    fileout.write(line)

    fileout.close()

    os.system("rm grid.nc")
    os.system("%s/mesh_tools/periodic_hex/periodic_grid" %(mpas_tools_dir))
    os.system("python %s/mesh_tools/periodic_hex/mark_periodic_boundaries_for_culling.py -f grid.nc" %(mpas_tools_dir))
    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid.nc %s" %(mpas_tools_dir,gridfileOut))

#-------------------------------------------------------------------------------

def create_forcing(gridFilename, forcingFilename, forcingPlotFilename):

    fileGrid = Dataset(gridFilename, "r")

    nCells = len(fileGrid.dimensions["nCells"])

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]

    fileGrid.close()


    simulationDurationInDays = 4
    hoursInDay = 24

    nTimes = hoursInDay * simulationDurationInDays

    uAirVelocity = np.zeros((nTimes,nCells))
    vAirVelocity = np.zeros((nTimes,nCells))
    airSpeed = np.zeros((nTimes,nCells))

    uOceanVelocity = np.zeros((nTimes,nCells))
    vOceanVelocity = np.zeros((nTimes,nCells))
    oceanSpeed = np.zeros((nTimes,nCells))

    xtimes = []

    vOceanMax = 0.01
    vAirMax = 30.0 / math.e

    alpha = math.radians(-90.0-18.0)
    salpha = math.sin(alpha)
    calpha = math.cos(alpha)

    year = 1
    month = 1
    day = 1
    hour = 0
    minute = 0
    second = 0

    for iTime in range(0,nTimes):
        for iCell in range(0,nCells):

            time = float(iTime) # hours

            mx = 0.5*lx + 0.1*lx * (time / 24.0)
            my = 0.5*ly + 0.1*ly * (time / 24.0)
            r = math.sqrt(math.pow((mx - xCell[iCell]), 2) +
                          math.pow((my - yCell[iCell]), 2))
            s = (math.e / 100.0) * math.exp(-(r / (100.0 * 1000.0)))

            uAirVelocity[iTime,iCell] = s * vAirMax * ( calpha * (xCell[iCell]-mx) + salpha * (yCell[iCell]-my)) * 0.001
            vAirVelocity[iTime,iCell] = s * vAirMax * (-salpha * (xCell[iCell]-mx) + calpha * (yCell[iCell]-my)) * 0.001
            airSpeed[iTime,iCell] = math.sqrt(math.pow(uAirVelocity[iTime,iCell],2) +
                                              math.pow(vAirVelocity[iTime,iCell],2))

            uOceanVelocity[iTime,iCell] = vOceanMax * ( (2.0*yCell[iCell] - lx) / lx)
            vOceanVelocity[iTime,iCell] = vOceanMax * (-(2.0*xCell[iCell] - ly) / ly)
            oceanSpeed[iTime,iCell] = math.sqrt(math.pow(uOceanVelocity[iTime,iCell],2) +
                                                math.pow(vOceanVelocity[iTime,iCell],2))

        timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
        xtimes.append(timeStr)

        hour += 1
        if hour > 23:
            hour = 0
            day += 1


    fileForcing = Dataset(forcingFilename, "w", format="NETCDF3_CLASSIC")

    fileForcing.createDimension("nCells", nCells)
    fileForcing.createDimension("StrLen", 64)
    fileForcing.createDimension("Time", None)

    var = fileForcing.createVariable("xtime", "c", dimensions=["Time","StrLen"])
    for iTime in range(0,nTimes):
        var[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        var[iTime,19:] = " "*45

    var = fileForcing.createVariable("uAirVelocity","d",dimensions=["Time","nCells"])
    var[:] = uAirVelocity[:]

    var = fileForcing.createVariable("vAirVelocity","d",dimensions=["Time","nCells"])
    var[:] = vAirVelocity[:]

    var = fileForcing.createVariable("uOceanVelocity","d",dimensions=["Time","nCells"])
    var[:] = uOceanVelocity[:]

    var = fileForcing.createVariable("vOceanVelocity","d",dimensions=["Time","nCells"])
    var[:] = vOceanVelocity[:]

    fileForcing.close()


    cm = 1/2.54  # centimeters in inches
    plt.rc('font', family="Times New Roman")
    plt.rc('mathtext',fontset="stix")
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    patches = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex]*0.001,yVertex[iVertex]*0.001))
        patches.append(Polygon(vertices,True))

    stride = 10

    cmap = matplotlib.cm.viridis

    fig, axes = plt.subplots(2, 2, figsize=(15*cm,13*cm))

    iTime = int(nTimes / 2)

    # u air velocity
    pc1 = PatchCollection(patches, cmap=cmap)
    pc1.set_array(uAirVelocity[iTime,:])
    axes[0,0].add_collection(pc1)

    divider = make_axes_locatable(axes[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc1, cax=cax, orientation='vertical')

    axes[0,0].set_aspect("equal")
    axes[0,0].set_xlabel("x (km)")
    axes[0,0].set_ylabel("y (km)")
    axes[0,0].set_title("u air velocity (m/s)")
    axes[0,0].set_xlim((np.amin(xVertex)*0.001,np.amax(xVertex)*0.001))
    axes[0,0].set_ylim((np.amin(yVertex)*0.001,np.amax(yVertex)*0.001))

    # v air velocity
    pc2 = PatchCollection(patches, cmap=cmap)
    pc2.set_array(vAirVelocity[iTime,:])
    axes[0,1].add_collection(pc2)

    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc2, cax=cax, orientation='vertical')

    axes[0,1].set_aspect("equal")
    axes[0,1].set_xlabel("x (km)")
    axes[0,1].set_ylabel("y (km)")
    axes[0,1].set_title("v air velocity (m/s)")
    axes[0,1].set_xlim((np.amin(xVertex)*0.001,np.amax(xVertex)*0.001))
    axes[0,1].set_ylim((np.amin(yVertex)*0.001,np.amax(yVertex)*0.001))

    # air velocity
    pc3 = PatchCollection(patches, cmap=cmap)
    pc3.set_array(airSpeed[iTime,:])
    axes[1,0].add_collection(pc3)
    axes[1,0].quiver(xCell[::stride]*0.001, yCell[::stride]*0.001, uAirVelocity[iTime,::stride], vAirVelocity[iTime,::stride])

    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc3, cax=cax, orientation='vertical')

    axes[1,0].set_aspect("equal")
    axes[1,0].set_xlabel("x (km)")
    axes[1,0].set_ylabel("y (km)")
    axes[1,0].set_title("Air velocity (m/s)")
    axes[1,0].set_xlim((np.amin(xVertex)*0.001,np.amax(xVertex)*0.001))
    axes[1,0].set_ylim((np.amin(yVertex)*0.001,np.amax(yVertex)*0.001))

    # ocean velocity
    pc4 = PatchCollection(patches, cmap=cmap)
    pc4.set_array(oceanSpeed[iTime,:])
    axes[1,1].add_collection(pc4)
    axes[1,1].quiver(xCell[::stride]*0.001, yCell[::stride]*0.001, uOceanVelocity[iTime,::stride], vOceanVelocity[iTime,::stride])

    divider = make_axes_locatable(axes[1,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc4, cax=cax, orientation='vertical')

    axes[1,1].set_aspect("equal")
    axes[1,1].set_xlabel("x (km)")
    axes[1,1].set_ylabel("y (km)")
    axes[1,1].set_title("Ocean velocity (m/s)")
    axes[1,1].set_xlim((np.amin(xVertex)*0.001,np.amax(xVertex)*0.001))
    axes[1,1].set_ylim((np.amin(yVertex)*0.001,np.amax(yVertex)*0.001))

    plt.tight_layout()
    plt.savefig(forcingPlotFilename,dpi=300)

#-------------------------------------------------------------------------------

def create_ic(gridFilename, icFilename):

    nCategories = 5

    fileGrid = Dataset(gridFilename, "r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]

    fileGrid.close()

    iceAreaCategory = np.zeros((nCells,nCategories,1))
    iceVolumeCategory = np.zeros((nCells,nCategories,1))

    for iCell in range(0,nCells):
        iceAreaCategory[iCell,0,0] = 1.0
        iceVolumeCategory[iCell,0,0] = 0.3 + 0.005 * (math.sin(6e-5*xCell[iCell]) + math.sin(3e-5*yCell[iCell]))

    fVertex = 1.46e-4

    fileIC = Dataset(icFilename,"w",format="NETCDF3_CLASSIC")

    fileIC.createDimension("nCells",nCells)
    fileIC.createDimension("nVertices",nVertices)
    fileIC.createDimension("nCategories",nCategories)
    fileIC.createDimension("ONE",1)

    var = fileIC.createVariable("iceAreaCategory","d",dimensions=["nCells","nCategories","ONE"])
    var[:] = iceAreaCategory[:]

    var = fileIC.createVariable("iceVolumeCategory","d",dimensions=["nCells","nCategories","ONE"])
    var[:] = iceVolumeCategory[:]

    var = fileIC.createVariable("fVertex","d",dimensions=["nVertices"])
    var[:] = fVertex

    fileIC.close()

#-------------------------------------------------------------------------------

def make_testcase():

    create_grid_quad("grid_moving_cyclone_quad_8km.nc", 8000.0)
    create_grid_quad("grid_moving_cyclone_quad_4km.nc", 4000.0)
    create_grid_quad("grid_moving_cyclone_quad_2km.nc", 2000.0)

    create_forcing("grid_moving_cyclone_quad_8km.nc", "forcing_moving_cyclone_quad_8km.nc", "forcing_moving_cyclone_quad_8km.png")
    create_forcing("grid_moving_cyclone_quad_4km.nc", "forcing_moving_cyclone_quad_4km.nc", "forcing_moving_cyclone_quad_4km.png")
    create_forcing("grid_moving_cyclone_quad_2km.nc", "forcing_moving_cyclone_quad_2km.nc", "forcing_moving_cyclone_quad_2km.png")

    create_ic("grid_moving_cyclone_quad_8km.nc", "ic_moving_cyclone_quad_8km.nc")
    create_ic("grid_moving_cyclone_quad_4km.nc", "ic_moving_cyclone_quad_4km.nc")
    create_ic("grid_moving_cyclone_quad_2km.nc", "ic_moving_cyclone_quad_2km.nc")

    create_grid_hex("grid_moving_cyclone_hex_8km.nc", 8000.0)
    create_grid_hex("grid_moving_cyclone_hex_4km.nc", 4000.0)
    create_grid_hex("grid_moving_cyclone_hex_2km.nc", 2000.0)

    create_forcing("grid_moving_cyclone_hex_8km.nc", "forcing_moving_cyclone_hex_8km.nc", "forcing_moving_cyclone_hex_8km.png")
    create_forcing("grid_moving_cyclone_hex_4km.nc", "forcing_moving_cyclone_hex_4km.nc", "forcing_moving_cyclone_hex_4km.png")
    create_forcing("grid_moving_cyclone_hex_2km.nc", "forcing_moving_cyclone_hex_2km.nc", "forcing_moving_cyclone_hex_2km.png")

    create_ic("grid_moving_cyclone_hex_8km.nc", "ic_moving_cyclone_hex_8km.nc")
    create_ic("grid_moving_cyclone_hex_4km.nc", "ic_moving_cyclone_hex_4km.nc")
    create_ic("grid_moving_cyclone_hex_2km.nc", "ic_moving_cyclone_hex_2km.nc")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    make_testcase()
