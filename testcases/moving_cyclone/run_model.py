import os

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    gridTypes = ["quad","hex"]
    resolutions = ["4km"]

    for gridType in gridTypes:

        print("gridType: ", gridType)

        for resolution in resolutions:

            print("  Resolution: ", resolution)

            os.system("rm -rf output_%s_%s" %(gridType, resolution))

            os.system("rm grid.nc")
            os.system("rm forcing.nc")
            os.system("rm ic.nc")

            os.system("ln -s grid_moving_cyclone_%s_%s.nc grid.nc" %(gridType,resolution))
            os.system("ln -s forcing_moving_cyclone_%s_%s.nc forcing.nc" %(gridType,resolution))
            os.system("ln -s ic_moving_cyclone_%s_%s.nc ic.nc" %(gridType,resolution))

            os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

            os.system("mv output output_%s_%s" %(gridType,resolution))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
