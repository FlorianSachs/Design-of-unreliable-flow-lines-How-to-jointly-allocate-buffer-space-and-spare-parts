import importlib
import os
import sys

import numpy as np
import pygad

os.chdir('python_GA')
sys.path.append(os.getcwd())

import flowline
import GA_flowline

importlib.reload(flowline)
importlib.reload(GA_flowline)

### machine setup
machines = 2

mu_normal = 1
p_normal = 0.01
gamma_normal = 0.1
costs_per_buffer = 1
costs_per_spare = 1

C_min = 1
C_max = 10
Q_min = 1
Q_max = 5

target_throughput = 0.90


### init
line = flowline.flowline()
line.fill(machines, C_max, Q_max, mu_normal, p_normal, gamma_normal, costs_per_buffer, costs_per_spare)
line.set_limits(C_min, C_max, Q_min, Q_max)

### calc
result = GA_flowline.genetic_algorithm(line, target_throughput, False)


### output
print(f'C = {result.C}')
print(f'Q = {result.Q}')
print(f'TP = {result.TP:.4f}')
print(f'costs = {result.costs:.4f}')
