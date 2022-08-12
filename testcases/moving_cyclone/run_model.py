import os

#-------------------------------------------------------------------------------

def run_model():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    resolutions = ["4km"]

    for resolution in resolutions:

        print("Resolution: ", resolution)

        os.system("rm -rf output_%s" %(resolution))

        os.system("rm grid.nc")
        os.system("rm forcing.nc")
        os.system("rm ic.nc")

        os.system("ln -s grid_moving_cyclone_%s.nc grid.nc" %(resolution))
        os.system("ln -s forcing_moving_cyclone_%s.nc forcing.nc" %(resolution))
        os.system("ln -s ic_moving_cyclone_%s.nc ic.nc" %(resolution))

        os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

        os.system("mv output output_%s" %(resolution))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_model()
