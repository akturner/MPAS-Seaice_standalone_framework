from make_testcase import make_testcase
from run_model import run_model
from plot_testcase import plot_testcase

make_testcase()

run_model()

plot_testcase("quad","4km")
plot_testcase("hex","4km")
