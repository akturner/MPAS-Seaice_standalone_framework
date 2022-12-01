import os
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_testcase():

    MPAS_SEAICE_EXECUTABLE = os.environ.get('MPAS_SEAICE_EXECUTABLE')
    if (MPAS_SEAICE_EXECUTABLE is None):
        raise Exception("MPAS_SEAICE_EXECUTABLE must be set")
    MPAS_SEAICE_TESTCASES_RUN_COMMAND = os.environ.get('MPAS_SEAICE_TESTCASES_RUN_COMMAND')
    if (MPAS_SEAICE_TESTCASES_RUN_COMMAND is None):
        MPAS_SEAICE_TESTCASES_RUN_COMMAND = ""
    MPAS_SEAICE_DOMAINS_DIR = os.environ.get('MPAS_SEAICE_DOMAINS_DIR')
    if (MPAS_SEAICE_DOMAINS_DIR is None):
        raise Exception("MPAS_SEAICE_DOMAINS_DIR must be set")

    # forcing
    os.system("python %s/domain_sc_71.35_-156.5/get_domain.py" %(MPAS_SEAICE_DOMAINS_DIR))

    os.system("python create_forcing.py -o forcing/ridging_forcing.nc")

    os.system("python create_input_itd.py -o forcing/itd.nc")

    os.system("rm namelist.seaice streams.seaice")
    os.system("ln -s namelist.seaice.ridging namelist.seaice")
    os.system("ln -s streams.seaice.ridging streams.seaice")

    # run MPAS-Seaice
    os.system("%s %s" %(MPAS_SEAICE_TESTCASES_RUN_COMMAND, MPAS_SEAICE_EXECUTABLE))

    # plot output
    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
