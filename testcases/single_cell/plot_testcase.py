from netCDF4 import Dataset
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def plot_testcase():

    filein = Dataset("./output/output.2000.nc","r")

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

    fig, axis = plt.subplots(figsize=(8*cm,7*cm))

    surfaceTemperatureCell = filein.variables["surfaceTemperatureCell"][:]

    axis.plot(surfaceTemperatureCell,color="green")
    axis.set_ylabel("Temperature (C)")
    axis.set_xlabel("Time step")
    axis.set_xlim(0,4700)
    axis.set_ylim(None,0)
    axis.set_title("MPAS_Seaice single cell")

    axis2 = axis.twinx()

    iceVolumeCell = filein.variables["iceVolumeCell"][:]
    snowVolumeCell = filein.variables["snowVolumeCell"][:]

    axis2.plot(iceVolumeCell,color="red")
    axis2.plot(snowVolumeCell,color="blue")
    axis2.set_ylabel("Thickness (m)")
    axis2.set_ylim(0,None)

    plt.tight_layout()
    plt.savefig("single_cell_thickness.eps")
    plt.savefig("single_cell_thickness.png",dpi=300)
    plt.cla()
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8*cm,7*cm))

    shortwaveDown = filein.variables["shortwaveDown"][:]
    longwaveDown = filein.variables["longwaveDown"][:]
    shortwaveVisibleDirectDown = filein.variables["shortwaveVisibleDirectDown"][:]
    shortwaveVisibleDiffuseDown = filein.variables["shortwaveVisibleDiffuseDown"][:]
    shortwaveIRDirectDown = filein.variables["shortwaveIRDirectDown"][:]
    shortwaveIRDiffuseDown = filein.variables["shortwaveIRDiffuseDown"][:]

    axis.plot(shortwaveDown, label="shortwaveDown", lw=0.2)
    axis.plot(longwaveDown, label="longwaveDown", lw=0.2)
    axis.plot(shortwaveVisibleDirectDown, label="shortwaveVisibleDirectDown", lw=0.2)
    axis.plot(shortwaveVisibleDiffuseDown, label="shortwaveVisibleDiffuseDown", lw=0.2)
    axis.plot(shortwaveIRDirectDown, label="shortwaveIRDirectDown", lw=0.2)
    axis.plot(shortwaveIRDiffuseDown, label="shortwaveIRDiffuseDown", lw=0.2)

    axis.set_ylabel("Flux (W/m2)")
    axis.legend()

    plt.tight_layout()
    plt.savefig("single_cell_shortwave_flux.eps")
    plt.savefig("single_cell_shortwave_flux.png",dpi=300)
    plt.cla()
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(8*cm,7*cm))

    albedoVisibleDirectCell = filein.variables["albedoVisibleDirectCell"][:]
    albedoVisibleDiffuseCell = filein.variables["albedoVisibleDiffuseCell"][:]
    albedoIRDirectCell = filein.variables["albedoIRDirectCell"][:]
    albedoIRDiffuseCell = filein.variables["albedoIRDiffuseCell"][:]
    bareIceAlbedoCell = filein.variables["bareIceAlbedoCell"][:]
    snowAlbedoCell = filein.variables["snowAlbedoCell"][:]
    pondAlbedoCell = filein.variables["pondAlbedoCell"][:]

    axis.plot(albedoVisibleDirectCell, label="albedoVisibleDirectCell", lw=0.2)
    axis.plot(albedoVisibleDiffuseCell, label="albedoVisibleDiffuseCell", lw=0.2)
    axis.plot(albedoIRDirectCell, label="albedoIRDirectCell", lw=0.2)
    axis.plot(albedoIRDiffuseCell, label="albedoIRDiffuseCell", lw=0.2)
    axis.plot(bareIceAlbedoCell, label="bareIceAlbedoCell", lw=0.2)
    axis.plot(snowAlbedoCell, label="snowAlbedoCell", lw=0.2)
    axis.plot(pondAlbedoCell, label="pondAlbedoCell", lw=0.2)

    axis.set_ylabel("Albedo")
    axis.legend()

    plt.tight_layout()
    plt.savefig("single_cell_albedo.eps")
    plt.savefig("single_cell_albedo.png",dpi=300)
    plt.cla()
    plt.close(fig)

    filein.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase()
