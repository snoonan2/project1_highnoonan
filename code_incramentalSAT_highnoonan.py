#Sophia Noonan
#snoonan2
#team name: highnoonan

#rewrite of the DumbSAT.py Brute Force solver using the incramental method of going
#through every possible solution

#need the imports of time (for execution times), random for the wffs setting, and sys for the input file of test cases
import time
import random
import sys
import matplotlib.pyplot as plt
import numpy as np

#this is the main function that changes the way we go through every solution compared to DumbSAT
def increment_assignment(assignment):
    #for every possible setting of the variable to 0 or 1 check and continue on to the next combination
    for i in range(len(assignment)):
        if assignment[i] == 0:
            assignment[i] = 1
            #returns true and would stop iterating before the DumbSAT
            return True
        assignment[i] = 0
    return False

def check(Wff, Nvars, Nclauses, Assignment):
    #check if there exists an assignment that satisfies the wff
    while True:
        #evaluate all clauses with the current assignment
        Satisfiable = True
        #check each claude
        for i in range(Nclauses): 
            Clause = Wff[i]
            ClauseSatisfiable = False
            
            #check each literal
            for Literal in Clause: 
                #find the literal's value
                VarValue = Assignment[abs(Literal)]  
                #calculate the expected value
                LitValue = 1 if Literal > 0 else 0  
                if LitValue == VarValue:
                    ClauseSatisfiable = True
                    #the clause is SAT and break
                    break  
            
            if not ClauseSatisfiable:
                Satisfiable = False
                #the clause is UNSAT break
                break  
        
        if Satisfiable:
            #SAT assignment
            return True  

        if not increment_assignment(Assignment):
            #there are no more assingments to try
            break  

    #no SAT assignment found
    return False 

#my build is the same as the given code in DumbSAT.py
def build_wff(Nvars, Nclauses, LitsPerClause):
    wff = []
    for i in range(Nclauses):
        clause = []
        for j in range(LitsPerClause):
            var = random.randint(1, Nvars)
            if random.randint(0, 1) == 0:
                var = -var
            clause.append(var)
        wff.append(clause)
    return wff

#my test_wff is also the same as the DumbSAT.py
def test_wff(wff, Nvars, Nclauses):
    Assignment = list((0 for _ in range(Nvars + 2)))
    start = time.time()  #start timer
    SatFlag = check(wff, Nvars, Nclauses, Assignment)
    end = time.time()  #end timer
    exec_time = int((end - start) * 1e6)
    return [wff, Assignment, SatFlag, exec_time]

#function to run the cases from the input or the defualt test cases
#this again is structured from the DumbSAT.py given to use from the canvas files to run and solve the cases but now with my incramental version
def run_cases(TestCases, ProbNum, resultsfile, tracefile, cnffile):
    f1 = open(resultsfile + ".csv", 'w')
    f2 = open(tracefile + ".csv", 'w')
    f3 = open(cnffile + ".cnf", "w")
    
    Nwffs = 0
    Nsat = 0
    Nunsat = 0

    for i in range(len(TestCases)):
        TestCase = TestCases[i]
        Nvars = TestCase[0]
        NClauses = TestCase[1]
        LitsPerClause = TestCase[2]
        Ntrials = TestCase[3]
        
        #times for the lines for my visual plot
        Scount = Ucount = 0
        AveStime = AveUtime = 0
        MaxStime = MaxUtime = 0

        for j in range(Ntrials):
            Nwffs += 1
            random.seed(ProbNum)
            wff = build_wff(Nvars, NClauses, LitsPerClause)
            results = test_wff(wff, Nvars, NClauses)
            Assignment = results[1]
            Exec_Time = results[3]

            #if satisfied
            if results[2]:
                y = 'S'
                Scount += 1
                AveStime += Exec_Time
                MaxStime = max(MaxStime, Exec_Time)
                Nsat += 1
            #if unsatisfied
            else:
                y = 'U'
                Ucount += 1
                AveUtime += Exec_Time
                MaxUtime = max(MaxUtime, Exec_Time)
                Nunsat += 1
            
            x = f"{ProbNum},{Nvars},{NClauses},{LitsPerClause},{NClauses * LitsPerClause},{y},1,{Exec_Time}"
            if results[2]:
                x += "," + ",".join(str(Assignment[k]) for k in range(1, Nvars + 1))
            print(x)
            #write to the files
            f1.write(x + '\n')
            f2.write(x + '\n')

            #add wff to cnf file
            x = f"c {ProbNum} {LitsPerClause} {y}\n"
            f3.write(x)
            x = f"p cnf {Nvars} {NClauses}\n"
            f3.write(x)
            for clause in wff:
                x = " ".join(str(literal) for literal in clause) + " 0\n"
                f3.write(x)
            ProbNum += 1

        #description for the tracefile as well as terminal output
        counts = f'# Satisfied = {Scount}. # Unsatisfied = {Ucount}'
        maxs = f'Max Sat Time = {MaxStime}. Max Unsat Time = {MaxUtime}'
        aves = f'Ave Sat Time = {AveStime / Scount if Scount else 0}. Ave UnSat Time = {AveUtime / Ucount if Ucount else 0}'
        print(counts)
        print(maxs)
        print(aves)
        #write to the files
        f2.write(counts + '\n')
        f2.write(maxs + '\n')
        f2.write(aves + '\n')

    x = f"{cnffile},TheBoss,{Nwffs},{Nsat},{Nunsat},{Nwffs},{Nwffs}\n"
    f1.write(x)
    f1.close()
    f2.close()
    f3.close()

#function to parse input from my input file with the test cases
def parse_input_file(filename):
    TestCases = []
    with open(filename, 'r') as f:
        for line in f:
            #skip empty lines and comments
            if line.strip() == '' or line.strip().startswith('#'):
                continue
            #parse the line into a list of integers
            case = list(map(int, line.strip().split()))
            if len(case) != 4:
                print(f"Warning: Skipping invalid line in input file: {line}")
                continue
            TestCases.append(case)
    return TestCases

#my function to run the cases with timing for the configuration of my plot
def run_cases_with_timing(TestCases):
    #collects the problems that are sat and unsat as well as the timings for my max and average line on my plot
    all_timings_sat = []
    all_timings_unsat = []
    all_problem_sizes_sat = []
    all_problem_sizes_unsat = []
    max_timings = []
    avg_timings = []
    unique_problem_sizes = []

    #do on every case in the test cases
    for case in TestCases:
        Nvars, NClauses, LitsPerClause, Ntrials = case
        problem_size = Nvars * NClauses * LitsPerClause

        #want to save every case's time
        case_timings = []
        for _ in range(Ntrials):
            wff = build_wff(Nvars, NClauses, LitsPerClause)
            start = time.time()
            is_sat = check(wff, Nvars, NClauses, [0] * (Nvars + 2))
            end = time.time()
            exec_time = (end - start) * 1000  #i converted to milliseconds
            case_timings.append(exec_time)
            
            #if the wff was satisfiable
            if is_sat:
                all_timings_sat.append(exec_time)
                all_problem_sizes_sat.append(problem_size)
            #if the wff was unsatisfiable
            else:
                all_timings_unsat.append(exec_time)
                all_problem_sizes_unsat.append(problem_size)

        #add to the lists created above for the creation of my lines
        max_timings.append(max(case_timings))
        avg_timings.append(sum(case_timings) / len(case_timings))
        unique_problem_sizes.append(problem_size)

    return (all_problem_sizes_sat, all_timings_sat,
            all_problem_sizes_unsat, all_timings_unsat,
            unique_problem_sizes, max_timings, avg_timings)

#function to plot the evecution time of each problem to the size of the problem
def plot_timings(all_problem_sizes_sat, all_timings_sat,
                 all_problem_sizes_unsat, all_timings_unsat,
                 unique_problem_sizes, max_timings, avg_timings):
    plt.figure(figsize=(12, 8))

    #plot SAT problems
    plt.scatter(all_problem_sizes_sat, all_timings_sat, alpha=0.3, color='blue', marker='^', label='SAT Problems')

    #plot UNSAT problems
    plt.scatter(all_problem_sizes_unsat, all_timings_unsat, alpha=0.3, color='red', marker='s', label='UNSAT Problems')

    #plot max timings for each unique problem size
    plt.scatter(unique_problem_sizes, max_timings, color='purple', s=100, label='Max Time per Problem Size')
    plt.plot(unique_problem_sizes, max_timings, color='purple', linestyle='--')

    #plot average timings for each unique problem size
    plt.scatter(unique_problem_sizes, avg_timings, color='green', s=100, label='Average Time per Problem Size')
    plt.plot(unique_problem_sizes, avg_timings, color='green', linestyle='--')

    #add exponential growth curve to show that as execution time increases with the problem size increase
    x = np.linspace(min(unique_problem_sizes), max(unique_problem_sizes), 100)
    y = np.exp(x / max(unique_problem_sizes) * np.log(max(max_timings)))
    plt.plot(x, y, color='orange', linestyle=':', label='Exponential Growth')

    plt.xlabel('Problem Size (Variables * Clauses * Literals per Clause)')
    plt.ylabel('Execution Time (ms)')
    plt.title('Incramental SAT Solver Execution Time vs Problem Size')

    plt.legend()
    plt.grid(True)
    plt.yscale('log') 
    plt.xscale('log') 
    plt.show()
    plt.savefig('plot_rewriteSATgraph_highnoonan.png')

# default test cases if there is no data input file for test cases
default_test_cases = [
    [4, 10, 2, 20],
    [8, 16, 2, 20],
    [12, 24, 2, 20],

]


if __name__ == "__main__":
    #check if an input file is provided as an argument from the user
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        try:
            TestCases = parse_input_file(input_file)
        except FileNotFoundError:
            print(f"Error: Input file '{input_file}' not found. Using default test cases.")
            TestCases = default_test_cases
    else:
        print("No input file provided. Using default test cases.")
        TestCases = default_test_cases
    
    ProbNum = 3

    #creating of the files names
    resultsfile = r'output_resultsOfRewriteSAT_highnoonan'
    tracefile = r'output_traceOfProgram_highnoonan'
    cnffile = r'output_cnffile_highnoonan'

#actually runs the Incramental SAT solver on the test cases
run_cases(TestCases, ProbNum, resultsfile, tracefile, cnffile)
(all_problem_sizes_sat, all_timings_sat,
 all_problem_sizes_unsat, all_timings_unsat,
 unique_problem_sizes, max_timings, avg_timings) = run_cases_with_timing(TestCases)

plot_timings(all_problem_sizes_sat, all_timings_sat,
             all_problem_sizes_unsat, all_timings_unsat,
             unique_problem_sizes, max_timings, avg_timings)