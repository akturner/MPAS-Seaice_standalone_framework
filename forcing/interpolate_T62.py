from netCDF4 import Dataset
import subprocess
import argparse
from create_forcing import create_scrip_file_MPAS, create_T62_remap_file, get_remapping_data
import numpy as np

#-------------------------------------------------------------------------------

def interpolate_T62(filenameMPASGrid,
                    filenameInput,
                    filenameOut,
                    varname,
                    latVarName,
                    lonVarName,
                    index):

    # create MPAS scrip grid file
    print("create_scrip_file_MPAS")
    scripGridFilename  = "remap_grid_MPAS_tmp.nc"
    dstGridSize = create_scrip_file_MPAS(filenameMPASGrid, scripGridFilename)

    # create T62 remapping file
    print("create_T62_remap_file")
    scripT62Filename = "remap_grid_T62_tmp.nc"
    srcGridSize = create_T62_remap_file(scripT62Filename, "T62", filenameInput, latVarName, lonVarName)

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

    # input file
    fileInput = Dataset(filenameInput,"r")
    varIn = fileInput.variables[varname]
    if (len(varIn.dimensions) == 2):
        arrayIn = varIn[:]
    elif (len(varIn.dimensions) == 3):
        arrayIn = varIn[index,:]
    else:
        raise Exception("Unexpected number of dimensions in T62 file")
    fileInput.close()

    # interpolate
    arrayOut = np.zeros(dstGridSize)
    arrayOut = remapMatrix.dot(arrayIn.flatten())

    # output
    fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

    fileOut.createDimension("nCells",dstGridSize)

    var = fileOut.createVariable(varname,"d",dimensions=["nCells"])
    var[:] = arrayOut[:]

    fileOut.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Interpolate from T62 to MPAS')
    parser.add_argument('-m', dest='filenameMPASGrid', required=True, help="MPAS grid filename")
    parser.add_argument('-i', dest='filenameInput',    required=True, help="T62 input filename")
    parser.add_argument('-o', dest='filenameOut',      required=True, help="Output filename")
    parser.add_argument('-v', dest='varname',          required=True, help="Variable name to interpolate")
    parser.add_argument('--latName', dest='latVarName', default="lat", help="Name of latitude varname in T62 file")
    parser.add_argument('--lonName', dest='lonVarName', default="lon", help="Name of longitude varname in T62 file")
    parser.add_argument('--index',   dest='index', type=int, default=0, help="Slice index of input variable")

    args = parser.parse_args()

    interpolate_T62(args.filenameMPASGrid,
                    args.filenameInput,
                    args.filenameOut,
                    args.varname,
                    args.latVarName,
                    args.lonVarName,
                    args.index)
