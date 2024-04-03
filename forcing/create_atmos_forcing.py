from __future__ import print_function
from netCDF4 import Dataset
import netCDF4
import numpy as np
import os
import subprocess
import sys
import glob
import configparser
from create_forcing import create_scrip_grid_file, create_T62_remap_file, get_mpas_grid_info, create_scrip_file_MPAS, write_scrip_in_file, create_output_times, get_remapping_data

#-------------------------------------------------------------------------------

def create_forcing(\
        filenameOutTemplate, \
        varnameInput, \
        yearStart, \
        yearStop, \
        filenameInputTemplate, \
        varnameOutput, \
        inputTimesPerYear, \
        yearStartData, \
        remapMatrix, \
        dstGridSize):

    # loop over years
    for year in range(yearStart,yearStop+1):

        print("  Year: %i of %i to %i" %(year, yearStart, yearStop))

        # create output file
        filenameOut = filenameOutTemplate.replace("$Y",str(year))
        fileForcing = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

        # dimensions
        nCells = fileForcing.createDimension("nCells",dstGridSize)
        StrLen = fileForcing.createDimension("StrLen",64)
        Time   = fileForcing.createDimension("Time",)

        # time
        xtimes = create_output_times(inputTimesPerYear, year)
        varXtime = fileForcing.createVariable("xtime","c",dimensions=["Time","StrLen"])
        for iTime in range(0,inputTimesPerYear):
            varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
            varXtime[iTime,19:] = " "*45

        # loop over variables
        for iVariable in range(0,len(varnameInput)):

            print("    Variable: %s to %s" %(varnameInput[iVariable], varnameOutput[iVariable]))

            # open input file
            filenameInput = sorted(glob.glob(filenameInputTemplate[iVariable].replace("$Y",str(year))))[0]
            fileInput = Dataset(filenameInput,"r")
            arrayIn = fileInput.variables[varnameInput[iVariable]][:]
            fileInput.close()

            arrayOut = np.zeros((inputTimesPerYear,dstGridSize))

            # loop over times
            for iTime in range(0,inputTimesPerYear):

                arrayInTime = arrayIn[iTime,:,:].flatten()
                arrayOut[iTime,:] = remapMatrix.dot(arrayInTime)

            # output variable to netcdf file
            var = fileForcing.createVariable(varnameOutput[iVariable],"d",dimensions=["Time","nCells"])
            var[:] = arrayOut[:]

        # close forcing file
        fileForcing.close()

#-------------------------------------------------------------------------------

def write_scrip_in_file(srcTitle):

    scripFile = open("scrip_in","w")

    scripFile.write("&remapInputs\n")
    scripFile.write("    num_maps = 1\n")
    scripFile.write("    gridFile1 = 'remap_grid_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    gridFile2 = 'remap_grid_MPAS_tmp.nc'\n")
    scripFile.write("    interpFile1 = 'remap_%s_to_MPAS_tmp.nc'\n" %(srcTitle))
    scripFile.write("    interpFile2 = 'remap_MPAS_to_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    mapName1 = '%s to MPAS bilinear mapping'\n" %(srcTitle))
    scripFile.write("    mapName2 = 'MPAS to %s bilinear mapping'\n" %(srcTitle))
    scripFile.write("    mapMethod = 'bilinear'\n")
    scripFile.write("    normalizeOpt = 'fracArea'\n")
    scripFile.write("    outputFormat = 'scrip'\n")
    scripFile.write("    restrict_type = 'latitude'\n")
    scripFile.write("    num_srch_bins = 90 \n")
    scripFile.write("    luse_grid1_area = .false.\n")
    scripFile.write("    luse_grid2_area = .false.\n")
    scripFile.write("/\n")

    scripFile.close()

#-------------------------------------------------------------------------------

def perform_remapping(\
        filenameMPASGrid, \
        outputDir, \
        startYear, \
        endYear, \
        dataDirSixHourly, \
        dataDirMonthly):

    # create MPAS scrip grid file
    print("create_scrip_file_MPAS")
    scripGridFilename  = "remap_grid_MPAS_tmp.nc"
    dstGridSize = create_scrip_file_MPAS(filenameMPASGrid, scripGridFilename)

    # create T62 remapping file
    print("create_T62_remap_file")
    scripT62Filename = "remap_grid_T62_tmp.nc"
    filenames = sorted(glob.glob(dataDirSixHourly+"/t_10/t_10.*.nc"))
    filenameT62 = filenames[0]
    srcGridSize = create_T62_remap_file(scripT62Filename, "T62", filenameT62)

    # create input scrip file
    print("write_scrip_in_file")
    write_scrip_in_file("T62")

    # run scrip to generate weights
    print("ESMF_RegridWeightGen")
    process = subprocess.Popen(["ESMF_RegridWeightGen","--source","remap_grid_T62_tmp.nc","--destination","remap_grid_MPAS_tmp.nc","--weight","remap_T62_to_MPAS_tmp.nc","--method","bilinear","--weight_only"])
    process.wait()
    if (process.returncode != 0):
        print("ESMF_RegridWeightGen error: ", process.returncode)
        exit(1);

    # get remapping weights
    print("get_remapping_data")
    filenameRemapping = "remap_T62_to_MPAS_tmp.nc"
    remapMatrix = get_remapping_data(filenameRemapping, srcGridSize, dstGridSize)

    print("create_forcing six hourly")
    # combined six hourly file
    create_forcing(\
       outputDir+"/LYq_six_hourly.$Y.nc", \
       ["T_10_MOD","Q_10_MOD","U_10_MOD","V_10_MOD"], \
       startYear, \
       endYear, \
       [dataDirSixHourly+"/t_10/t_10.$Y.*.nc", \
        dataDirSixHourly+"/q_10/q_10.$Y.*.nc", \
        dataDirSixHourly+"/u_10/u_10.$Y.*.nc", \
        dataDirSixHourly+"/v_10/v_10.$Y.*.nc"], \
       ["airTemperature", \
        "airSpecificHumidity", \
        "uAirVelocity", \
        "vAirVelocity"], \
       1460, \
       1948, \
       remapMatrix, \
       dstGridSize)

    print("create_forcing monthly")
    # combined monthly file
    create_forcing(\
       outputDir+"/LYq_monthly.nc", \
       ["cldf","prec"], \
       0, \
       0, \
       [dataDirMonthly+"/cldf.omip.nc", \
        dataDirMonthly+"/prec.nmyr.nc"], \
       ["cloudFraction", \
        "rainfallRate"], \
       12, \
       0, \
       remapMatrix, \
       dstGridSize)

#-------------------------------------------------------------------------------

'''
create_atmos_forcing.py
=======================

Usage
-----

This script creates atmospheric forcing using six hourly CORE-II data and
monthly AOMIP climatologies.

Usage: python create_atmos_forcing.py configFilename

where configFilename is a python config file with the following example format:

[forcing_generation]
filenameMPASGrid = /location/of/MPAS/grid
outputDir = /location/to/put/output/forcing
startYear = 1948
endYear = 2007
dataDirSixHourly = /location/of/CORE-II/data
dataDirMonthly = /location/of/AOMIP/climatologies

SCRIP
-----

This script requires the SCRIP package to be installed.
SCRIP is a software package which computes addresses and weights for remapping
and interpolating fields between grids in spherical coordinates. It can be
obtained from https://github.com/SCRIP-Project/SCRIP

CORE-II data
------------

Six-hourly air temperature, velocity and specific humidity comes from CORE-II.
Data files can be obtained from
https://data1.gfdl.noaa.gov/nomads/forms/core/COREv2/CIAF_v2.html
To generate forcing for a given year YYYY, the following files are required:
${dataDirSixHourly}/t_10/t_10.YYYY.*.nc
${dataDirSixHourly}/q_10/q_10.YYYY.*.nc
${dataDirSixHourly}/u_10/u_10.YYYY.*.nc
${dataDirSixHourly}/v_10/v_10.YYYY.*.nc
where ${dataDirSixHourly} is the local location of the six hourly data.

AOMIP climatologies
-------------------

Monthly climatologies of cloudiness and precipitation comes from AOMIP.
The following data files are required:
${dataDirMonthly}/cldf.omip.nc
${dataDirMonthly}/prec.nmyr.nc
where ${dataDirMonthly} is the local location of the monthly data.
These files can be obtained from:
https://web.lcrc.anl.gov/public/e3sm/mpas_standalonedata/mpas-seaice/forcing/
MPAS-Seaice_clim_data.tar.gz
'''

if (len(sys.argv) != 2):
    print("Usage: python create_atmos_forcing.py configFilename")
    sys.exit()

config = configparser.ConfigParser()
config.read(sys.argv[1])

filenameMPASGrid = config.get   ('forcing_generation','filenameMPASGrid')
outputDir        = config.get   ('forcing_generation','outputDir')
startYear        = config.getint('forcing_generation','startYear')
endYear          = config.getint('forcing_generation','endYear')
dataDirSixHourly = config.get   ('forcing_generation','dataDirSixHourly')
dataDirMonthly   = config.get   ('forcing_generation','dataDirMonthly')

perform_remapping(\
        filenameMPASGrid, \
        outputDir, \
        startYear, \
        endYear, \
        dataDirSixHourly, \
        dataDirMonthly)
