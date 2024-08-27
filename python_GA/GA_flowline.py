#import matplotlib as mpl
#from matplotlib import pyplot as plt
import numpy as np
import types
import pygad
import random

np.random.seed(42)
random.seed(42)

penalty = 1000

counter_evaluations = 0
counter_problems = 0


def genetic_algorithm(line, TP_target, do_output=False, do_plots=False):
    global counter_evaluations
    global counter_problems

    costs = np.concatenate([line.costs_buffers, line.costs_spares])

    ### fitness
    def fitness_func(ga_instance, solution, solution_idx):
        global counter_evaluations
        global counter_problems

        line.C = solution[:line.machines - 1]
        line.Q = solution[line.machines - 1:]
        costs_total_local = np.sum(solution * costs)
        line.solve()
        counter_evaluations += 1
        if not line.characteristics.terminated_normally:
            counter_problems += 1
        TP = line.characteristics.TP

        if TP >= TP_target:
            fitness = -costs_total_local
        else:
            fitness = -costs_total_local - penalty * ((TP_target - TP) / TP_target)
        return fitness


    ### config
    fitness_function = fitness_func

    num_generations = 100
    num_parents_mating = 2

    sol_per_pop = 32
    num_genes = len(costs)

    # parent_selection_type="sss": The parent selection type. Supported types are sss (for steady-state selection),
    # rws (for roulette wheel selection), sus (for stochastic universal selection), rank (for rank selection),
    # random (for random selection), and tournament (for tournament selection).
    parent_selection_type = "sss"

    keep_parents = 1

    crossover_type = "single_point"

    mutation_type = "random"
    # mutation_percent_genes = 10
    mutation_num_genes = 1

    gene_space = []
    for i in range(0, line.machines - 1):
        gene_space.append(range(line.C_min, line.C_max + 1))
    for i in range(0, line.machines):
        gene_space.append(range(line.Q_min, line.Q_max + 1))


    ### GA
    ga_instance = pygad.GA(save_solutions=True,
                           save_best_solutions=True,
                           num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=fitness_function,
                           sol_per_pop=sol_per_pop,
                           num_genes=num_genes,
                           gene_space=gene_space,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           crossover_type=crossover_type,
                           mutation_type=mutation_type,
                           # mutation_percent_genes=mutation_percent_genes,
                           mutation_num_genes=mutation_num_genes,
                           gene_type=int,
                           random_seed=42)

    # run
    ga_instance.run()

    # best solution
    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    line.C = solution[:line.machines - 1]
    line.Q = solution[line.machines - 1:]
    costs_total = np.sum(solution * costs)
    line.solve()
    TP = line.characteristics.TP

    # prepare result
    result = types.SimpleNamespace()
    result.C = line.C
    result.Q = line.Q
    result.TP = TP
    result.costs = costs_total
    result.counter_iterations = ga_instance.best_solution_generation
    result.counter_evaluations = counter_evaluations
    result.counter_problems = counter_problems


    ### output
    if do_output:
        print("Initial Population")
        print(ga_instance.initial_population)

        print("Final Population")
        print(ga_instance.population)

        print("Parameters of the best solution: {solution}".format(solution=solution))
        print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))

        print("Costs of the best solution: {prediction}".format(prediction=costs_total))
        print("TP of the best solution: {prediction}".format(prediction=TP))

    if do_plots:
        pass
        # ga_instance.plot_fitness(save_dir='1.png')
        # ga_instance.plot_new_solution_rate(save_dir='2.png')
        # ga_instance.plot_genes(solutions='all', save_dir='3.png')
        # ga_instance.plot_genes(solutions='best', save_dir='4.png')

    return result
