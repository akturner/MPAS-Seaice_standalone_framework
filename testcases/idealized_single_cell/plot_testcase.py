from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import math
import argparse

#-------------------------------------------------------------------------------

def plot_subtestcase(process):

    filein = Dataset("./output_%s/output.2000.nc" %(process),"r")

    nCategories = len(filein.dimensions["nCategories"])
    nTimes = len(filein.dimensions["Time"])
    nIceLayers = len(filein.dimensions["nIceLayers"])
    nSnowLayers = len(filein.dimensions["nSnowLayers"])

    iceAreaCategory = filein.variables["iceAreaCategory"][:,0,:,0]

    iceAreaCell = filein.variables["iceAreaCell"][:,0]
    iceVolumeCell = filein.variables["iceVolumeCell"][:,0]
    snowVolumeCell = filein.variables["snowVolumeCell"][:,0]

    surfaceHeatFlux = filein.variables["surfaceHeatFlux"][:,0]
    surfaceShortwaveFluxCategory = filein.variables["surfaceShortwaveFlux"][:,0,:]
    longwaveDown = filein.variables["longwaveDown"][:,0]
    longwaveUp = filein.variables["longwaveUp"][:,0]
    sensibleHeatFlux = filein.variables["sensibleHeatFlux"][:,0]
    latentHeatFlux = filein.variables["latentHeatFlux"][:,0]

    surfaceIceMelt = filein.variables["surfaceIceMelt"][:,0]
    basalIceMelt = filein.variables["basalIceMelt"][:,0]
    lateralIceMelt = filein.variables["lateralIceMelt"][:,0]
    snowMelt = filein.variables["snowMelt"][:,0]
    congelation = filein.variables["congelation"][:,0]
    snowiceFormation = filein.variables["snowiceFormation"][:,0]
    frazilFormation = filein.variables["frazilFormation"][:,0]

    surfaceTemperatureCell = filein.variables["surfaceTemperatureCell"][:]
    iceTemperature = filein.variables["iceTemperature"][:,0,0,:]
    snowTemperature = filein.variables["snowTemperature"][:,0,0,:]

    seaSurfaceTemperature = filein.variables["seaSurfaceTemperature"][:,0]
    seaFreezingTemperature = filein.variables["seaFreezingTemperature"][:,0]
    seaSurfaceSalinity = filein.variables["seaSurfaceSalinity"][:,0]
    oceanMixedLayerDepth = filein.variables["oceanMixedLayerDepth"][:,0]

    oceanFreshWaterFlux = filein.variables["oceanFreshWaterFlux"][:,0]
    oceanSaltFlux = filein.variables["oceanSaltFlux"][:,0]
    oceanHeatFlux = filein.variables["oceanHeatFlux"][:,0]
    oceanShortwaveFlux = filein.variables["oceanShortwaveFlux"][:,0]

    filein.close()

    surfaceShortwaveFlux = np.zeros(nTimes)
    oceanFluxTemperature = np.zeros(nTimes)
    for iCategory in range(0,nCategories):
        surfaceShortwaveFlux[:] += surfaceShortwaveFluxCategory[:,iCategory] * iceAreaCategory[:,iCategory]

    seaiceDensitySeaWater = 1026.0
    seaiceSeaWaterSpecificHeat = 4218.0
    oceanFluxTemperature = np.zeros(nTimes)
    for iTime in range(0,nTimes):
        if (surfaceIceMelt[iTime] > 0.0):
            if (oceanFreshWaterFlux[iTime] > 0.0):
                oceanFluxTemperature[iTime] = oceanHeatFlux[iTime] / (oceanFreshWaterFlux[iTime] * seaiceDensitySeaWater * seaiceSeaWaterSpecificHeat)
            else:
                oceanFluxTemperature[iTime] = 0.0
        else:
            oceanFluxTemperature[iTime] = math.nan

    emissivity = 1.0

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

    # ice area and volume
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    axis.plot(iceAreaCell,color="green",label="iceAreaCell",lw=0.5)
    axis.set_ylabel("Area")
    axis.set_xlabel("Time step")
    axis.set_ylim(0,1)
    axis.set_xlim(0,None)
    axis.set_title("Ice area and volumes")
    axis.legend(loc='center left', bbox_to_anchor=(1.2, 0.3))
    axis.yaxis.label.set_color("green")

    axis2 = axis.twinx()

    axis2.plot(iceVolumeCell,color="red",label="iceVolumeCell",lw=0.5)
    axis2.plot(snowVolumeCell,color="blue",label="snowVolumeCell",lw=0.5)
    axis2.set_ylabel("Volume (m)")
    axis2.set_ylim(0,None)
    axis2.legend(loc='center left', bbox_to_anchor=(1.2, 0.7))
    axis2.yaxis.label.set_color("red")

    plt.tight_layout()
    plt.savefig("ice_area_and_volume_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

    # surface heat balance
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    axis.plot(surfaceHeatFlux,label="surfaceHeatFlux",color="black",lw=0.5) # fsurfn
    axis.plot(surfaceShortwaveFlux,label="surfaceShortwaveFlux",lw=0.5) # fswsfc
    axis.plot(longwaveDown*emissivity,label="longwaveDown Abs.",lw=0.5) # flwdabs = emissivity * flw
    axis.plot(longwaveUp,label="longwaveUp",lw=0.5) # flwoutn
    axis.plot(sensibleHeatFlux,label="sensibleHeatFlux",lw=0.5) # fsensn
    axis.plot(latentHeatFlux,label="latentHeatFlux",lw=0.5) # flatn

    axis.set_ylabel("Surface heat flux")
    axis.set_xlabel("Time step")
    axis.set_xlim(0,None)
    axis.set_title("Surface heat balance")
    axis.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("surface_heat_balance_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

    # melt growth rates
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    axis.plot(surfaceIceMelt,label="surfaceIceMelt",lw=0.5)
    axis.plot(basalIceMelt,label="basalIceMelt",lw=0.5)
    axis.plot(lateralIceMelt,label="lateralIceMelt",lw=0.5)
    axis.plot(snowMelt,label="snowMelt",lw=0.5)
    axis.plot(congelation,label="congelation",lw=0.5)
    axis.plot(snowiceFormation,label="snowiceFormation",lw=0.5)
    axis.plot(frazilFormation,label="frazilFormation",lw=0.5)

    axis.set_ylabel("Melt growth rates")
    axis.set_xlabel("Time step")
    axis.set_xlim(0,None)
    axis.set_title("Melt growth rates")
    axis.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("melt_growth_rates_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

    # temperatures
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    axis.plot(surfaceTemperatureCell,label="Surface Temperature",color="black",lw=0.5)
    for iSnowLayer in range(0,nSnowLayers):
        axis.plot(snowTemperature[:,iSnowLayer],label="Snow Temperature %i" %(iSnowLayer),lw=0.5)
    for iIceLayer in range(0,nIceLayers):
        axis.plot(iceTemperature[:,iIceLayer],label="Ice Temperature %i" %(iIceLayer),lw=0.5)
    axis.plot(seaSurfaceTemperature,label="seaSurfaceTemperature",lw=0.5)
    axis.plot(seaFreezingTemperature,label="seaFreezingTemperature",lw=0.5)

    axis.set_ylabel(r'Temperature ($^\circ$ C)')
    axis.set_xlabel("Time step")
    axis.set_xlim(0,None)
    axis.set_title("Temperatures")
    axis.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("temperatures_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

    # ocean temperatures
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    axis.plot(seaSurfaceTemperature,label="seaSurfaceTemperature",lw=0.5)
    axis.plot(seaFreezingTemperature,label="seaFreezingTemperature",lw=0.5)

    axis.set_ylabel(r'Temperature ($^\circ$ C)')
    axis.set_xlabel("Time step")
    axis.set_xlim(0,None)
    axis.set_title("Temperatures")
    axis.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("sea_temperatures_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

    # ocean
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    color1 = plt.cm.jet(0)
    color2 = plt.cm.jet(0.33)
    color3 = plt.cm.jet(0.66)
    color4 = plt.cm.jet(0.99)

    par1 = axis.twinx()
    par2 = axis.twinx()

    p1, = axis.plot(seaSurfaceTemperature,label="seaSurfaceTemperature",color=color1)
    p2, = axis.plot(oceanFluxTemperature,label="oceanFluxTemperature",color=color2)
    p3, = par1.plot(seaSurfaceSalinity,label="seaSurfaceSalinity",color=color3)
    p4, = par2.plot(oceanMixedLayerDepth,label="oceanMixedLayerDepth",color=color4)

    axis.set_xlabel("Time step")
    axis.set_ylabel(r'Temperature ($^\circ$ C)')
    par1.set_ylabel("Sea Surface Salinity (ppt)")
    par2.set_ylabel("Mixed layer depth (m)")

    par2.spines['right'].set_position(('outward', 60))
    par2.xaxis.set_ticks([])

    axis.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p3.get_color())
    par2.yaxis.label.set_color(p4.get_color())

    lns = [p1, p2, p3, p4]
    axis.legend(handles=lns, loc='best')

    plt.tight_layout()
    plt.savefig("ocean_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

    # ocean fluxes
    fig, axis = plt.subplots(figsize=(12*cm,7*cm))

    color1 = plt.cm.jet(0)
    color2 = plt.cm.jet(0.33)
    color3 = plt.cm.jet(0.66)
    color4 = plt.cm.jet(0.99)

    par1 = axis.twinx()
    par2 = axis.twinx()

    p1, = axis.plot(oceanHeatFlux,label="oceanHeatFlux",color=color1)
    p2, = axis.plot(oceanShortwaveFlux,label="oceanShortwaveFlux",color=color2)
    p3, = par1.plot(oceanFreshWaterFlux,label="oceanFreshWaterFlux",color=color3, lw=2.5)
    p4, = par2.plot(oceanSaltFlux,label="oceanSaltFlux",color=color4)

    axis.set_xlabel("Time step")
    axis.set_ylabel(r'Heat flux ($\mathrm{W} \mathrm{m}^{-2}$)')
    par1.set_ylabel(r'Fresh water flux ($\mathrm{kg} \mathrm{m}^{-2}$)')
    par2.set_ylabel(r'Salt flux ($\mathrm{kg} \mathrm{m}^{-2}$)')

    par2.spines['right'].set_position(('outward', 60))
    par2.xaxis.set_ticks([])

    axis.yaxis.label.set_color(p1.get_color())
    par1.yaxis.label.set_color(p3.get_color())
    par2.yaxis.label.set_color(p4.get_color())

    lns = [p1, p2, p3, p4]
    axis.legend(handles=lns, loc='best')

    plt.tight_layout()
    plt.savefig("ocean_fluxes_%s.png" %(process),dpi=300)
    plt.cla()
    plt.close(fig)

#-------------------------------------------------------------------------------

def plot_testcase(process=None):

    if (process is None):
        processes = ["congelation",
                     "surface_ice_melt",
                     "frazil"]
    else:
        processes = [process]

    for process in processes:
        plot_subtestcase(process)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-p', dest='process', default=None, help='')

    args = parser.parse_args()

    plot_testcase(args.process)
