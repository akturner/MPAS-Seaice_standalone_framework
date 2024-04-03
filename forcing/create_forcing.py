from netCDF4 import Dataset
import numpy as np
import scipy.sparse as sparse

#-------------------------------------------------------------------------------

def create_scrip_grid_file(filenameScrip, nGridSize, nGridCorners, gridRank, gridDims, gridCenterLat, gridCenterLon, gridImask, gridCornerLat, gridCornerLon, title):

    fileScrip = Dataset(filenameScrip,"w",format="NETCDF3_CLASSIC")

    # dimensions
    fileScrip.createDimension("grid_size", nGridSize)
    fileScrip.createDimension("grid_corners", nGridCorners)
    fileScrip.createDimension("grid_rank", gridRank)

    # variables
    var = fileScrip.createVariable("grid_dims","i",dimensions=["grid_rank"])
    var[:] = gridDims[:]

    var = fileScrip.createVariable("grid_center_lat","d",dimensions=["grid_size"])
    var.units = "radians"
    var[:] = gridCenterLat[:]

    var = fileScrip.createVariable("grid_center_lon","d",dimensions=["grid_size"])
    var.units = "radians"
    var[:] = gridCenterLon[:]

    var = fileScrip.createVariable("grid_imask","i",dimensions=["grid_size"])
    var[:] = gridImask[:]

    var = fileScrip.createVariable("grid_corner_lat","d",dimensions=["grid_size","grid_corners"])
    var.units = "radians"
    var[:] = gridCornerLat[:]

    var = fileScrip.createVariable("grid_corner_lon","d",dimensions=["grid_size","grid_corners"])
    var.units = "radians"
    var[:] = gridCornerLon[:]

    # attributes
    fileScrip.title = title

    fileScrip.close()

#-------------------------------------------------------------------------------

def create_T62_remap_file(filenameScrip, title, filenameT62, latVarName="LAT", lonVarName="LON"):

    nLat = 94
    nLon = 192

    fileT62 = Dataset(filenameT62,"r")

    LATin = fileT62.variables[latVarName][:]
    LONin = fileT62.variables[lonVarName][:]

    fileT62.close()

    LAT = np.zeros(nLat+2)
    LON = np.zeros(nLon+2)

    LAT[1:-1] = LATin[:]
    LON[1:-1] = LONin[:]

    # center
    latCenter = np.zeros((nLat,nLon))
    lonCenter = np.zeros((nLat,nLon))

    for iLon in range(0,nLon):
       latCenter[:,iLon] = LATin[:]

    for iLat in range(0,nLat):
       lonCenter[iLat,:] = LONin[:]

    latCenter = np.radians(latCenter)
    lonCenter = np.radians(lonCenter)

    # corners
    latCorner = np.zeros((nLat,nLon,4))
    lonCorner = np.zeros((nLat,nLon,4))

    for iLon in range(0,nLon):

        iLon2 = iLon + 1

        lonCorner[:,iLon,0] = 0.5 * (LON[iLon2] + LON[iLon2-1])
        lonCorner[:,iLon,1] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,2] = 0.5 * (LON[iLon2] + LON[iLon2+1])
        lonCorner[:,iLon,3] = 0.5 * (LON[iLon2] + LON[iLon2-1])

    lonCorner = (np.where(lonCorner < 0.0, lonCorner + 360.0, lonCorner))
    lonCorner = np.radians(lonCorner)

    for iLat in range(0,nLat):

        iLat2 = iLat + 1

        latCorner[iLat,:,0] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,1] = 0.5 * (LAT[iLat2] + LAT[iLat2-1])
        latCorner[iLat,:,2] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])
        latCorner[iLat,:,3] = 0.5 * (LAT[iLat2] + LAT[iLat2+1])

    latCorner = np.radians(latCorner)

    # create file
    nGridSize = nLat * nLon
    nGridCorners = 4
    gridRank = 2
    gridDims = np.array([nLon,nLat])
    gridImask = np.ones(nGridSize,dtype="i")

    latCornerScrip = np.zeros((nLat*nLon,4))
    lonCornerScrip = np.zeros((nLat*nLon,4))

    for iLat in range(0,nLat):
        for iLon in range(0,nLon):
            for iCorner in range(0,4):
                ij = iLat * nLon + iLon
                latCornerScrip[ij,iCorner] = latCorner[iLat,iLon,iCorner]
                lonCornerScrip[ij,iCorner] = lonCorner[iLat,iLon,iCorner]

    create_scrip_grid_file(filenameScrip, nGridSize, nGridCorners, gridRank, gridDims, latCenter.flatten(), lonCenter.flatten(), gridImask, latCornerScrip, lonCornerScrip, title)

    return nGridSize

#-------------------------------------------------------------------------------

def get_mpas_grid_info(filenameMPASGrid):

    fileGrid = Dataset(filenameMPASGrid,"r")

    nCells    = len(fileGrid.dimensions["nCells"])
    maxEdges  = len(fileGrid.dimensions["maxEdges"])
    nVertices = len(fileGrid.dimensions["nVertices"])

    latCell   = fileGrid.variables["latCell"][:]
    lonCell   = fileGrid.variables["lonCell"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    lonVertex = fileGrid.variables["lonVertex"][:]

    nEdgesOnCell   = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]

    fileGrid.close()

    return nCells, maxEdges, nVertices, latCell, lonCell, latVertex, lonVertex, nEdgesOnCell, verticesOnCell

#-------------------------------------------------------------------------------

def create_scrip_file_MPAS(filenameMPASGrid, filenameScrip):

    nCells, maxEdges, nVertices, centerLat, centerLon, latVertex, lonVertex, nEdgesOnCell, verticesOnCell = get_mpas_grid_info(filenameMPASGrid)

    cornerLat = np.zeros((nCells,maxEdges))
    cornerLon = np.zeros((nCells,maxEdges))

    # create the corner arrays
    for iCell in range(0,nCells):

       # fill vertices we have
       for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

          iVertex = verticesOnCell[iCell,iVertexOnCell] - 1

          cornerLat[iCell,iVertexOnCell] = latVertex[iVertex]
          cornerLon[iCell,iVertexOnCell] = lonVertex[iVertex]

       # fill up maxEdges corners with repeated vertex
       for iVertexOnCell in range(nEdgesOnCell[iCell], maxEdges):

          cornerLat[iCell,iVertexOnCell] = cornerLat[iCell,nEdgesOnCell[iCell]-1]
          cornerLon[iCell,iVertexOnCell] = cornerLon[iCell,nEdgesOnCell[iCell]-1]

    # create scrip file
    gridDims = np.array([nCells])
    gridImask = np.ones(nCells,dtype="i")

    create_scrip_grid_file(filenameScrip, nCells, maxEdges, 1, gridDims, centerLat, centerLon, gridImask, cornerLat, cornerLon, "MPAS")

    return nCells

#-------------------------------------------------------------------------------

def write_scrip_in_file(srcTitle):

    scripFile = open("scrip_in","w")

    scripFile.write("&remap_inputs\n")
    scripFile.write("    num_maps = 1\n")
    scripFile.write("    grid1_file = 'remap_grid_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    grid2_file = 'remap_grid_MPAS_tmp.nc'\n")
    scripFile.write("    interp_file1 = 'remap_%s_to_MPAS_tmp.nc'\n" %(srcTitle))
    scripFile.write("    interp_file2 = 'remap_MPAS_to_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    map1_name = '%s to MPAS bilinear mapping'\n" %(srcTitle))
    scripFile.write("    map2_name = 'MPAS to %s bilinear mapping'\n" %(srcTitle))
    scripFile.write("    map_method = 'bilinear'\n")
    scripFile.write("    normalize_opt = 'frac'\n")
    scripFile.write("    output_opt = 'scrip'\n")
    scripFile.write("    restrict_type = 'latitude'\n")
    scripFile.write("    num_srch_bins = 90 \n")
    scripFile.write("    luse_grid1_area = .false.\n")
    scripFile.write("    luse_grid2_area = .false.\n")
    scripFile.write("/\n")

    scripFile.close()

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

def get_remapping_data(filenameRemapping, srcGridSize, dstGridSize):

    fileRemap = Dataset(filenameRemapping, "r")

    n_s = len(fileRemap.dimensions["n_s"])
    col = fileRemap.variables["col"][:]
    row = fileRemap.variables["row"][:]
    S = fileRemap.variables["S"][:]

    fileRemap.close()

    # fortran indices
    col[:] = col[:] - 1
    row[:] = row[:] - 1

    # covert to python sparse arrays
    remapMatrixSparse = sparse.coo_matrix((S, (row, col)), shape=(dstGridSize, srcGridSize))
    remapMatrixSparse = remapMatrixSparse.tocsr()

    return remapMatrixSparse

#-------------------------------------------------------------------------------
