from create_mesh import create_mesh
from create_forcing import create_forcing
from run_model import run_model
from plot_testcase import plot_testcase

#-------------------------------------------------------------------------------

def run_testcase():

    # create mesh
    print("Create mesh")
    create_mesh

    # create forcing
    print("Create forcing")
    create_forcing

    # run the model
    print("Run models")
    run_model()

    # plot sub tests
    print("Plot subtests")
    plot_testcase()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    run_testcase()
