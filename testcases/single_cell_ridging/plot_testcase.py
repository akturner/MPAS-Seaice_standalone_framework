import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def plot_testcase():

    # read in data
    filenameIn = "output/output.2000.nc"

    filein = Dataset(filenameIn,"r")

    nTimes = len(filein.dimensions["Time"])
    nIceLayers = len(filein.dimensions["nIceLayers"])

    ridgeConvergence = filein.variables["ridgeConvergence"][:]
    ridgeShear = filein.variables["ridgeShear"][:]

    iceAreaCell = filein.variables["iceAreaCell"][:]
    iceVolumeCell = filein.variables["iceVolumeCell"][:]

    iceAreaCategory = filein.variables["iceAreaCategory"][:]
    iceVolumeCategory = filein.variables["iceVolumeCategory"][:]

    surfaceTemperature = filein.variables["surfaceTemperature"][:]


    filein.close()

    iceThicknessCell = np.zeros(nTimes)
    for iTime in range(0,nTimes):
        if (iceAreaCell[iTime,0] > 0.0):
            iceThicknessCell[iTime] = iceVolumeCell[iTime,0] / iceAreaCell[iTime,0]


    # plot
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

    fig, axes = plt.subplots(3,2,figsize=(15*cm,15*cm))

    # ridgeConvergence, ridgeShear
    axes[0,0].plot(ridgeConvergence[:,0], label="conv.")
    axes[0,0].plot(ridgeShear[:,0], label="shear")

    axes[0,0].legend()
    axes[0,0].set_ylabel(r'ridgeConvergence/Shear ($\mathrm{s}^{-1}$)')

    # ice area, volume, thickness
    axes[0,1].plot(iceAreaCell[:,0], label="ice conc.", color="k")
    axis2 = axes[0,1].twinx()
    axis2.plot(iceVolumeCell[:,0], label="ice volume", color="b", lw=2)
    axis2.plot(iceThicknessCell, label="ice thickness", color="g")

    axes[0,1].legend()
    axes[0,1].set_ylabel("iceAreaCell (-)")

    axis2.legend()
    axis2.set_ylabel("iceVolumeCell/iceThicknessCell (m)")

    # iceAreaCategory
    axes[1,0].plot(iceAreaCategory[:,0,0,0]-iceAreaCategory[0,0,0,0], label="Cat. 1")
    axes[1,0].plot(iceAreaCategory[:,0,1,0]-iceAreaCategory[0,0,1,0], label="Cat. 2")
    axes[1,0].plot(iceAreaCategory[:,0,2,0]-iceAreaCategory[0,0,2,0], label="Cat. 3")
    axes[1,0].plot(iceAreaCategory[:,0,3,0]-iceAreaCategory[0,0,3,0], label="Cat. 4")
    axes[1,0].plot(iceAreaCategory[:,0,4,0]-iceAreaCategory[0,0,4,0], label="Cat. 5")

    axes[1,0].legend()
    axes[1,0].set_ylabel(r'$\Delta$iceAreaCategory (-)')

    # iceVolumeCategory
    axes[1,1].plot(iceVolumeCategory[:,0,0,0]-iceVolumeCategory[0,0,0,0], label="Cat. 1")
    axes[1,1].plot(iceVolumeCategory[:,0,1,0]-iceVolumeCategory[0,0,1,0], label="Cat. 2")
    axes[1,1].plot(iceVolumeCategory[:,0,2,0]-iceVolumeCategory[0,0,2,0], label="Cat. 3")
    axes[1,1].plot(iceVolumeCategory[:,0,3,0]-iceVolumeCategory[0,0,3,0], label="Cat. 4")
    axes[1,1].plot(iceVolumeCategory[:,0,4,0]-iceVolumeCategory[0,0,4,0], label="Cat. 5")

    axes[1,1].legend()
    axes[1,1].set_ylabel(r'$\Delta$iceVolumeCategory (m)')

    # surface temperature
    axes[2,0].plot(surfaceTemperature[:,0,0,0]-surfaceTemperature[0,0,0,0], label="Cat. 1")
    axes[2,0].plot(surfaceTemperature[:,0,1,0]-surfaceTemperature[0,0,1,0], label="Cat. 2")
    axes[2,0].plot(surfaceTemperature[:,0,2,0]-surfaceTemperature[0,0,2,0], label="Cat. 3")
    axes[2,0].plot(surfaceTemperature[:,0,3,0]-surfaceTemperature[0,0,3,0], label="Cat. 4")
    axes[2,0].plot(surfaceTemperature[:,0,4,0]-surfaceTemperature[0,0,4,0], label="Cat. 5")

    axes[2,0].legend()
    axes[2,0].set_ylabel(r'$\Delta$surfaceTemperature ($\mathrm{C}^\circ$)')

    axes[2,1].axis('off')

    plt.tight_layout()
    plt.savefig("ridging.png",dpi=400)


#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase()
