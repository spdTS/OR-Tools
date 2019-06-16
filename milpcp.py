from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model
from pathlib import Path


class RPQ:
    def __init__(self, r, p, q):
        self.R = r
        self.P = p
        self.Q = q


def Milp(jobs, instanceName):

    variablesMaxValue = 0
    for i in range(len(jobs)):
        variablesMaxValue += (jobs[i].R + jobs[i].P + jobs[i].Q)

    solver = pywraplp.Solver('simple_mip_program',
    pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)

    # variables :
    alfasMatrix = {}
    for i in range(len(jobs)):
        for j in range(len(jobs)):
            alfasMatrix[i, j] = solver.IntVar(0, 1, "alfa" + str(i) + "_" + str(j))
    starts = []
    for i in range(len(jobs)):
        starts.append(solver.IntVar(0, variablesMaxValue, "starts" + str(i)))
    cmax = solver.IntVar(0, variablesMaxValue, "cmax")

    # constraints:
    for i in range(len(jobs)):
        solver.Add(starts[i] >= jobs[i].R)
        solver.Add(cmax >= starts[i] + jobs[i].P + jobs[i].Q)

    for i in range(len(jobs)):
        for j in range(i + 1, len(jobs)):
            solver.Add(starts[i] + jobs[i].P <= starts[j]
                + alfasMatrix[i, j] * variablesMaxValue)
            solver.Add(starts[j] + jobs[j].P <= starts[i]
                + alfasMatrix[j, i] * variablesMaxValue)
            solver.Add(alfasMatrix[i, j] + alfasMatrix[j, i] == 1)

    # solver:
    solver.Minimize(cmax)
    
    status = solver.Solve()
    if status is not pywraplp.Solver.OPTIMAL:
        print("Not optimal!")
    print(instanceName, "Cmax:", solver.Objective().Value())
    pi = []
    for i in range(len(starts)):
        pi.append((i+1, starts[i].solution_value()))
    pi.sort(key=lambda x: x[1])
    print(pi)


def cp(jobs, instance_name):
    
    model = cp_model.CpModel()

    variables_max_value = 0
    for i in range(len(jobs)):
        variables_max_value += (jobs[i].R + jobs[i].P + jobs[i].Q)

    # variables :
    alfa_matrix = {}
    for i in range(len(jobs)):
        for j in range(len(jobs)):
            alfa_matrix[i, j] = model.NewIntVar(0, 1, "alfa" + str(i) + "_" + str(j))
    starts = []
    for i in range(len(jobs)):
        starts.append(model.NewIntVar(0, variables_max_value, "starts" + str(i)))
    cmax = model.NewIntVar(0, variables_max_value, "cmax")

    # constraints:
    for i in range(len(jobs)):
        model.Add(starts[i] >= jobs[i].R)
        model.Add(cmax >= starts[i] + jobs[i].P + jobs[i].Q)

    for i in range(len(jobs)):
        for j in range(i + 1, len(jobs)):
            model.Add(starts[i] + jobs[i].P <= starts[j]
                + alfa_matrix[i, j] * variables_max_value)
            model.Add(starts[j] + jobs[j].P <= starts[i]
                + alfa_matrix[j, i] * variables_max_value)
            model.Add(alfa_matrix[i, j] + alfa_matrix[j, i] == 1)

    # solver:
    model.Minimize(cmax)

    solver = cp_model.CpSolver()
    #solver.parameters.max_time_in_seconds = 120.0
    solver.Solve(model)
    print(instance_name, "Cmax:", solver.ObjectiveValue())
    pi = []
    for i in range(len(starts)):
        pi.append((i+1, solver.Value(starts[i])))
    pi.sort(key=lambda x: x[1])
    show = [x[0] for x in pi]
    print(show)


def GetRPQsFromFile(pathToFile):
    
    fullTextFromFile = Path(pathToFile).read_text()
    words=fullTextFromFile.replace("\n"," ").split(" ")
    words_cleaned=list(filter(None,words))
    numbers=list(map(int,words_cleaned))
    
    numberOfJobs=numbers[0]
    numbers.pop(0)
    numbers.pop(0)
    
    jobs=[]
    for i in range(numberOfJobs):
        jobs.append(RPQ(numbers[0],numbers[1],numbers[2]))
        numbers.pop(0)
        numbers.pop(0)
        numbers.pop(0)
    return jobs


if __name__ == '__main__':
    
    file_paths=["data000.txt"]
    
    for i in range(len(file_paths)):
        jobs=GetRPQsFromFile(file_paths[i])
        #Milp(jobs,file_paths[i])
        cp(jobs,file_paths[i])
