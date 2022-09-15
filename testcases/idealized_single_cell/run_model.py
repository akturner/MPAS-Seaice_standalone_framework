import os
import argparse

#-------------------------------------------------------------------------------

def run_model(process=None):

    if (process is None):
        processes = ["congelation",
                     "surface_ice_melt",
                     "frazil"]
    else:
        processes = [process]

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""

    for process in processes:

        print("Process: ", process)

        if (not os.path.isdir("output")):
            os.mkdir("output")

        os.system("rm -rf namelist.seaice streams.seaice output_%s" %(process))
        os.system("ln -s namelist.seaice.%s namelist.seaice" %(process))
        os.system("ln -s streams.seaice.%s streams.seaice" %(process))

        os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

        os.system("mv output output_%s" %(process))

        os.system("mv log.seaice.0000.out log.seaice.0000.out_%s" %(process))

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-p', dest='process', default=None, help='')

    args = parser.parse_args()

    run_model(args.process)
