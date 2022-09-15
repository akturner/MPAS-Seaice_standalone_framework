from netCDF4 import Dataset
import numpy as np
import netCDF4

#-------------------------------------------------------------------------------

def create_output_times(inputTimesPerYear, year):

    daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]

    xtimes = []

    if (inputTimesPerYear == 1460):

        minute = 0
        second = 0

        for iMonth in range(0,12):
            for iDay in range(0,daysInMonth[iMonth]):
                for iSixHours in range(0,4):

                    month = iMonth + 1
                    day = iDay + 1
                    hour = (iSixHours + 1) * 6

                    timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
                    xtimes.append(timeStr)

    elif (inputTimesPerYear == 12):

        day = 15
        hour = 0
        minute = 0
        second = 0

        for iMonth in range(0,12):

            month = iMonth + 1

            timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
            xtimes.append(timeStr)

    return xtimes

#-------------------------------------------------------------------------------

def create_forcing_subtest(
        forcingName,
        airTemperature,
        airSpecificHumidity,
        uAirVelocity,
        vAirVelocity,
        cloudFraction,
        rainfallRate,
        seaSurfaceSalinity,
        seaSurfaceTemperature,
        uOceanVelocity,
        vOceanVelocity,
        seaSurfaceTiltU,
        seaSurfaceTiltV,
        oceanMixedLayerDepth,
        oceanHeatFluxConvergence):

    # six hourly atmos forcing
    fileOut = Dataset("atmosphere_forcing_six_hourly.%s.2000.nc" %(forcingName),"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",1)
    fileOut.createDimension("StrLen",64)
    fileOut.createDimension("Time",None)

    # time
    xtimes = create_output_times(1460, 2000)
    varXtime = fileOut.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,1460):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45

    airTemperatureVar = fileOut.createVariable("airTemperature","d",dimensions=["Time", "nCells"])
    airSpecificHumidityVar = fileOut.createVariable("airSpecificHumidity","d",dimensions=["Time", "nCells"])
    uAirVelocityVar = fileOut.createVariable("uAirVelocity","d",dimensions=["Time", "nCells"])
    vAirVelocityVar = fileOut.createVariable("vAirVelocity","d",dimensions=["Time", "nCells"])

    airTemperatureVar[:] = airTemperature
    airSpecificHumidityVar[:] = airSpecificHumidity
    uAirVelocityVar[:] = uAirVelocity
    vAirVelocityVar[:] = vAirVelocity

    fileOut.close()

    # monthly atmos forcing
    fileOut = Dataset("atmosphere_forcing_monthly.%s.nc" %(forcingName),"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",1)
    fileOut.createDimension("StrLen",64)
    fileOut.createDimension("Time",None)

    # time
    xtimes = create_output_times(12, 0)
    varXtime = fileOut.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,12):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45

    cloudFractionVar = fileOut.createVariable("cloudFraction","d",dimensions=["Time", "nCells"])
    rainfallRateVar = fileOut.createVariable("rainfallRate","d",dimensions=["Time", "nCells"])

    cloudFractionVar[:] = cloudFraction
    rainfallRateVar[:] = rainfallRate

    fileOut.close()

    # monthly ocean forcing
    fileOut = Dataset("ocean_forcing_monthly.%s.nc" %(forcingName),"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",1)
    fileOut.createDimension("StrLen",64)
    fileOut.createDimension("Time",None)

    # time
    xtimes = create_output_times(12, 0)
    varXtime = fileOut.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,12):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(xtimes[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45

    cloudFractionVar = fileOut.createVariable("cloudFraction","d",dimensions=["Time", "nCells"])
    rainfallRateVar = fileOut.createVariable("rainfallRate","d",dimensions=["Time", "nCells"])

    seaSurfaceSalinityVar = fileOut.createVariable("seaSurfaceSalinity","d",dimensions=["Time", "nCells"])
    seaSurfaceTemperatureVar = fileOut.createVariable("seaSurfaceTemperature","d",dimensions=["Time", "nCells"])
    uOceanVelocityVar = fileOut.createVariable("uOceanVelocity","d",dimensions=["Time", "nCells"])
    vOceanVelocityVar = fileOut.createVariable("vOceanVelocity","d",dimensions=["Time", "nCells"])
    seaSurfaceTiltUVar = fileOut.createVariable("seaSurfaceTiltU","d",dimensions=["Time", "nCells"])
    seaSurfaceTiltVVar = fileOut.createVariable("seaSurfaceTiltV","d",dimensions=["Time", "nCells"])
    oceanMixedLayerDepthVar = fileOut.createVariable("oceanMixedLayerDepth","d",dimensions=["Time", "nCells"])
    oceanHeatFluxConvergenceVar = fileOut.createVariable("oceanHeatFluxConvergence","d",dimensions=["Time", "nCells"])

    seaSurfaceSalinityVar[:] = seaSurfaceSalinity
    seaSurfaceTemperatureVar[:] = seaSurfaceTemperature
    uOceanVelocityVar[:] = uOceanVelocity
    vOceanVelocityVar[:] = vOceanVelocity
    seaSurfaceTiltUVar[:] = seaSurfaceTiltU
    seaSurfaceTiltVVar[:] = seaSurfaceTiltV
    oceanMixedLayerDepthVar[:] = oceanMixedLayerDepth
    oceanHeatFluxConvergenceVar[:] = oceanHeatFluxConvergence

    fileOut.close()

#-------------------------------------------------------------------------------

def create_forcing():

    processes = ["congelation",
                 "surface_ice_melt",
                 "frazil"]

    args = {"congelation":
            {"airTemperature":273.15 - 30.0, # Kelvin
             "airSpecificHumidity":0.0001,
             "uAirVelocity":0.0,
             "vAirVelocity":0.0,
             "cloudFraction":1.0,
             "rainfallRate":0.0,
             "seaSurfaceSalinity":34.0,
             "seaSurfaceTemperature":-1.8,
             "uOceanVelocity":0.0,
             "vOceanVelocity":0.0,
             "seaSurfaceTiltU":0.0,
             "seaSurfaceTiltV":0.0,
             "oceanMixedLayerDepth":50.0,
             "oceanHeatFluxConvergence":0.0},
            "surface_ice_melt":
            {"airTemperature":273.15 + 30.0, # Kelvin
             "airSpecificHumidity":0.0001,
             "uAirVelocity":0.0,
             "vAirVelocity":0.0,
             "cloudFraction":1.0,
             "rainfallRate":0.0,
             "seaSurfaceSalinity":34.0,
             "seaSurfaceTemperature":-1.8,
             "uOceanVelocity":0.0,
             "vOceanVelocity":0.0,
             "seaSurfaceTiltU":0.0,
             "seaSurfaceTiltV":0.0,
             "oceanMixedLayerDepth":50.0,
             "oceanHeatFluxConvergence":0.0},
            "frazil":
            {"airTemperature":273.15 - 30.0, # Kelvin
             "airSpecificHumidity":0.0001,
             "uAirVelocity":0.0,
             "vAirVelocity":0.0,
             "cloudFraction":1.0,
             "rainfallRate":0.0,
             "seaSurfaceSalinity":34.0,
             "seaSurfaceTemperature":-1.8,
             "uOceanVelocity":0.0,
             "vOceanVelocity":0.0,
             "seaSurfaceTiltU":0.0,
             "seaSurfaceTiltV":0.0,
             "oceanMixedLayerDepth":50.0,
             "oceanHeatFluxConvergence":0.0}}

    for process in processes:
        create_forcing_subtest(
            process,
            args[process]["airTemperature"],
            args[process]["airSpecificHumidity"],
            args[process]["uAirVelocity"],
            args[process]["vAirVelocity"],
            args[process]["cloudFraction"],
            args[process]["rainfallRate"],
            args[process]["seaSurfaceSalinity"],
            args[process]["seaSurfaceTemperature"],
            args[process]["uOceanVelocity"],
            args[process]["vOceanVelocity"],
            args[process]["seaSurfaceTiltU"],
            args[process]["seaSurfaceTiltV"],
            args[process]["oceanMixedLayerDepth"],
            args[process]["oceanHeatFluxConvergence"])

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_forcing()
