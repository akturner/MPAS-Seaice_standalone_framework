from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
import matplotlib.cm as mcm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean

#-------------------------------------------------------------------------------

def plot_testcase(gridType,resolution):

    filenameIn = "output_%s_%s/output.0001.nc" %(gridType,resolution)
    fileIn = Dataset(filenameIn,"r")

    nCells = len(fileIn.dimensions["nCells"])
    nVertices = len(fileIn.dimensions["nVertices"])
    vertexDegree = len(fileIn.dimensions["vertexDegree"])

    nEdgesOnCell = fileIn.variables["nEdgesOnCell"][:]

    verticesOnCell = fileIn.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    cellsOnVertex = fileIn.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    interiorVertex = fileIn.variables["interiorVertex"][:]

    xCell = fileIn.variables["xCell"][:]
    yCell = fileIn.variables["yCell"][:]

    xVertex = fileIn.variables["xVertex"][:]
    yVertex = fileIn.variables["yVertex"][:]

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    iceAreaCell = fileIn.variables["iceAreaCell"][:]
    iceVolumeCell = fileIn.variables["iceVolumeCell"][:]
    shear = fileIn.variables["shear"][:]
    uVelocity = fileIn.variables["uVelocity"][:]
    vVelocity = fileIn.variables["vVelocity"][:]

    speed = np.sqrt(np.add(np.power(uVelocity,2),np.power(vVelocity,2)))

    fileIn.close()


    cmap = mcm.viridis



    patchesCells = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex]*0.001,yVertex[iVertex]*0.001))
        patchesCells.append(Polygon(vertices,True))

    patchesVertices = []
    verticesIndex = []
    for iVertex in range(0,nVertices):
        if (interiorVertex[-1,iVertex]):
            verticesIndex.append(iVertex)
            vertices = []
            for iCellOnVertex in range(0,vertexDegree):
                iCell = cellsOnVertex[iVertex,iCellOnVertex]
                vertices.append((xCell[iCell]*0.001,yCell[iCell]*0.001))
            patchesVertices.append(Polygon(vertices,True))



    cm = 1/2.54  # centimeters in inches
    plt.rc('font',family="Times New Roman")
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

    stride = 10

    fig, axes = plt.subplots(2,2,figsize=(15*cm,13*cm))

    # velocity
    pc1 = PatchCollection(patchesVertices, cmap=cmap)
    pc1.set_array(speed[-1,verticesIndex])
    pc1.set_clim([0.0,None])
    axes[0,0].add_collection(pc1)
    axes[0,0].quiver(xVertex[::stride]*0.001, yVertex[::stride]*0.001, uVelocity[-1,::stride], vVelocity[-1,::stride])

    divider = make_axes_locatable(axes[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc1, cax=cax, orientation='vertical')

    # shear
    norm = matplotlib.colors.LogNorm(vmin=1e-8,vmax=1e-4)
    #norm = matplotlib.colors.LogNorm()
    pc2 = PatchCollection(patchesCells, cmap="Greys", norm=norm)
    shear[:] /= (100.0 * 86400.0) # % per day to s^-1
    #shear = np.where(shear < 1e-8, 1e-8, shear)
    pc2.set_array(shear[-1,:])
    #pc2.set_clim([0.0,None])
    axes[0,1].add_collection(pc2)

    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc2, cax=cax, orientation='vertical')

    # ice concentration
    norm = matplotlib.colors.Normalize(vmin=0.75,vmax=1.0)
    pc3 = PatchCollection(patchesCells, cmap=cmocean.cm.ice, norm=norm)
    pc3.set_array(iceAreaCell[-1,:])
    axes[1,0].add_collection(pc3)

    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc3, cax=cax, orientation='vertical')

    # ice volume
    pc4 = PatchCollection(patchesCells, cmap=cmap)
    pc4.set_array(iceVolumeCell[-1,:])
    axes[1,1].add_collection(pc4)

    divider = make_axes_locatable(axes[1,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pc4, cax=cax, orientation='vertical')




    axes[0,0].set_xlim((xMin*0.001,xMax*0.001))
    axes[0,0].set_ylim((yMin*0.001,yMax*0.001))
    axes[0,0].set_title("Velocity (m/s)")
    axes[0,0].set_aspect("equal")
    axes[0,0].set_xlabel("x (km)")
    axes[0,0].set_ylabel("y (km)")

    axes[0,1].set_xlim((xMin*0.001,xMax*0.001))
    axes[0,1].set_ylim((yMin*0.001,yMax*0.001))
    axes[0,1].set_title("Shear (s^-1)")
    axes[0,1].set_aspect("equal")
    axes[0,1].set_xlabel("x (km)")
    axes[0,1].set_ylabel("y (km)")

    axes[1,0].set_xlim((xMin*0.001,xMax*0.001))
    axes[1,0].set_ylim((yMin*0.001,yMax*0.001))
    axes[1,0].set_title("Concentration (-)")
    axes[1,0].set_aspect("equal")
    axes[1,0].set_xlabel("x (km)")
    axes[1,0].set_ylabel("y (km)")

    axes[1,1].set_xlim((xMin*0.001,xMax*0.001))
    axes[1,1].set_ylim((yMin*0.001,yMax*0.001))
    axes[1,1].set_title("Thickness (m)")
    axes[1,1].set_aspect("equal")
    axes[1,1].set_xlabel("x (km)")
    axes[1,1].set_ylabel("y (km)")

    plt.tight_layout()
    plt.savefig("moving_cyclone_%s_%s.png" %(gridType,resolution),dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-m', dest='gridType', choices=["quad","hex"], help='')
    parser.add_argument('-r', dest='resolution', choices=["2km","4km","8km"], help='')

    args = parser.parse_args()

    plot_testcase(args.gridType,args.resolution)
