import re
import subprocess
import csv
from subprocess import PIPE, run
from pathlib import Path

MODELS = ["mtz", "scf", "mcf"]
APP = "./kmst"
DATA = "data/"
LOG = "log/"
Path(LOG).mkdir(parents=True, exist_ok=True)

INSTANCES = [["g01.dat", 2, 46],    ["g01.dat", 5, 477],
             ["g02.dat", 4, 373],   ["g02.dat", 10, 1390],
             ["g03.dat", 10, 725],  ["g03.dat", 25, 3074],
             ["g04.dat", 14, 909],  ["g04.dat", 35, 3292],
             ["g05.dat", 20, 1235], ["g05.dat", 50, 4898],
             ["g06.dat", 40, 2068], ["g06.dat", 100, 6705],
             ["g07.dat", 60, 1335], ["g07.dat", 150, 4534],
             ["g08.dat", 80, 1620], ["g08.dat", 200, 5787],
             ["g09.dat", 200, 2289], ["g09.dat", 500, 7595],
             ["g10.dat", 400, 4182], ["g10.dat", 1000, 14991]]

reg_objective = re.compile("Objective value:\s*(\d+)")
reg_duration = re.compile("CPU time:\s*(\d+(?:\.\d+)?)")
reg_bb = re.compile("Branch-and-Bound nodes: (\d+)")
reg_status = re.compile("CPLEX status:\s*([a-zA-Z]{3})")
reg_gap = re.compile("(\d+\.\d\d\%)")

with open('out.csv', 'w', newline='') as csvfile:
    result_writer = csv.writer(csvfile, delimiter="\t")

    for i in INSTANCES:
        for m in MODELS:
            command = [APP, "-f", DATA + i[0], "-m", m, "-k", str(i[1])]
            result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
            output = result.stdout

            with open(f"{LOG}{i[0]}-{i[1]}-{m}.txt", "w") as log_file:
                log_file.write(LOG + output)

            try:
                status = reg_status.search(output).group(1)
            except:
                print(output)
                raise Exception("ERROR could not check status!")

            if status == "Opt" or status == "Fea":
                duration = float(reg_duration.search(output).group(1))
                bbnodes = int(reg_bb.search(output).group(1))
                value = int(reg_objective.search(output).group(1))

                try:
                    gap = reg_gap.findall(output)[-1]
                except:
                    gap = "?"

                if status == "Opt":
                    if value != i[2]:
                        raise Exception(f"ERROR FOR INSTANCE {i[0]}: {i[2]} expected but output is {value}!")
                    print(f"{i[0]} \t{m} optimal in \t{duration} s.\t{bbnodes} bb nodes")
                elif status == "Fea":
                    print(f"{i[0]} \t{m} feasible in \t{duration} s.\t{bbnodes} bb nodes, \t{gap} gap (opt: {i[2]} val: {value})")
                result_writer.writerow([i[0], i[1], m, i[2], value, status, duration, bbnodes, gap])
            elif status == "Unk":
                print(f"{i[0]} \t{m} unsolvable")
                result_writer.writerow([i[0], i[1], m, i[2], 0, status, 0, 0, 0])
