from netCDF4 import Dataset
import argparse
import numpy as np

#-------------------------------------------------------------------------------

def create_input_itd(filenameout):

    nCells = 1
    nCategories = 5

    categoryThicknessLimits = np.zeros(nCategories+1)
    categoryThicknessLimits[0] = 0.0
    categoryThicknessLimits[1] = 0.6
    categoryThicknessLimits[2] = 1.4
    categoryThicknessLimits[3] = 2.4
    categoryThicknessLimits[4] = 3.6
    categoryThicknessLimits[5] = 1e6

    initialCategoryIceArea = np.zeros(nCategories)
    #initialCategoryIceArea[0] = 0.05
    #initialCategoryIceArea[1] = 0.1
    #initialCategoryIceArea[2] = 0.3
    #initialCategoryIceArea[3] = 0.35
    #initialCategoryIceArea[4] = 0.2
    initialCategoryIceArea[:] = 0.2

    surfaceTemperature = np.zeros((nCells,nCategories,1))
    surfaceTemperature[:,0,0] = -1.0
    surfaceTemperature[:,1,0] = -2.0
    surfaceTemperature[:,2,0] = -3.0
    surfaceTemperature[:,3,0] = -4.0
    surfaceTemperature[:,4,0] = -5.0

    fileout = Dataset(filenameout,"w",format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nCategories",nCategories)
    fileout.createDimension("nCategoriesP1",nCategories+1)
    fileout.createDimension("ONE",1)

    categoryThicknessLimitsVar = fileout.createVariable("categoryThicknessLimits","d",dimensions=["nCategoriesP1"])
    categoryThicknessLimitsVar[:] = categoryThicknessLimits[:]

    initialCategoryIceAreaVar = fileout.createVariable("initialCategoryIceArea","d",dimensions=["nCategories"])
    initialCategoryIceAreaVar[:] = initialCategoryIceArea[:]

    surfaceTemperatureVar = fileout.createVariable("surfaceTemperature","d",dimensions=["nCells","nCategories","ONE"])
    surfaceTemperatureVar[:] = surfaceTemperature[:]

    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-o', dest='filenameOut', required=True, help='')

    args = parser.parse_args()

    create_input_itd(args.filenameOut)
