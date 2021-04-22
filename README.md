# k-Minimum Spanning Tree

This repo contains the programming exercise for the course Mathematical Programming.
The goal is to develop a MILP formulation to solve the k-Minimum Spanning Tree.

## Setup
1. Install CPLEX Studio
2. Add the path to your CPLEX installation in the `Makefile`. In the `Makefile`, locate the `STEP 2: ...`comment for instructions.
3. Run `make clean`
4. Run `make`
5. Run `./kmst`, which will start the solver using the default solver `mtz` and the default instance `g01.dat`
6. Run `./kmst --help` for usage examples
7. Run `./test.sh` to run all test instances using all 3 models
8. Run `python test.py` to run all test instances, to collect the output in the log files and to see the results in the `out.csv` table.


