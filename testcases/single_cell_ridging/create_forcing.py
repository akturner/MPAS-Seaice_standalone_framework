from netCDF4 import Dataset
import netCDF4
import numpy as np
from datetime import datetime, timedelta
import argparse

#-------------------------------------------------------------------------------

def create_forcing_file(filenameout, nTimes, timesStr, ridgeConvergence, ridgeShear):

    fileout = Dataset(filenameout,"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",1)
    fileout.createDimension("Time",None)
    fileout.createDimension("StrLen",64)

    varXtime = fileout.createVariable("xtime","c",dimensions=["Time","StrLen"])
    for iTime in range(0,nTimes):
        varXtime[iTime,0:19] = netCDF4.stringtochar(np.array(timesStr[iTime], 'S19'))
        varXtime[iTime,19:] = " "*45

    varRidgeConvergence = fileout.createVariable("ridgeConvergence","d",dimensions=["Time","nCells"])
    varRidgeConvergence[:] = ridgeConvergence[:]

    varRidgeShear = fileout.createVariable("ridgeShear","d",dimensions=["Time","nCells"])
    varRidgeShear[:] = ridgeShear[:]

    fileout.close()

#-------------------------------------------------------------------------------

def create_forcing(filenameout):

    nTimes = 8760

    ridgeConvergence = np.zeros(nTimes)
    ridgeShear = np.zeros(nTimes)

    timesStr = []
    time0 = datetime.fromisoformat('0001-01-01T00:00:00')
    timeDelta = timedelta(hours=1)
    time = time0
    for iTime in range(0,nTimes):
        timeSinceStart = float((time - time0).total_seconds())

        if (iTime > 1000 and
            iTime < 2000):

            ridgeConvergence[iTime] = 1e-8

        if (iTime > 3000 and
            iTime < 4000):

            ridgeConvergence[iTime] = -1e-8

        if (iTime > 5000 and
            iTime < 6000):

            ridgeShear[iTime] = 1e-8

        if (iTime > 7000 and
            iTime < 8000):

            ridgeConvergence[iTime] = 1e-8
            ridgeShear[iTime] = 1e-8

        timesStr.append(time.strftime("%04Y-%m-%d_%H:%M:%S"))
        time = time + timeDelta

    create_forcing_file(filenameout, nTimes, timesStr, ridgeConvergence, ridgeShear)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-o', dest='filenameOut', required=True, help='')

    args = parser.parse_args()

    create_forcing(args.filenameOut)
