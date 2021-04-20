import re
import subprocess
import csv

MODELS = ["mcf"]
APP = "./kmst"
DATA = "data/"

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

with open('out.csv', 'w', newline='') as csvfile:
    result_writer = csv.writer(csvfile, delimiter="\t")

    for m in MODELS:
        for i in INSTANCES:
            output = subprocess.check_output([APP, "-f", DATA + i[0],
                                                "-m", m, "-k", str(i[1])])
            output = output.decode('utf-8')

            match = reg_objective.search(output)
            if not match:
                raise Exception("INVALID OUTPUT FOR ", i)

            duration = float(reg_duration.search(output).group(1))
            bbnodes = int(reg_bb.search(output).group(1))
            value = int(match.group(1))
            
            if value != i[2]:
                raise Exception(f"ERROR FOR INSTANCE {i[0]}: {i[2]} expected but output is {value}!")

            print(f"{i[0]} solved in \t{duration} s. with \t{m} and \t{bbnodes} bb nodes")
            
            result_writer.writerow([i[0], i[1], m, value, duration, bbnodes])
