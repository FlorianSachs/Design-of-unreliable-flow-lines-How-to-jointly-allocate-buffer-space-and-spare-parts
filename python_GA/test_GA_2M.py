import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pygad
import random
from flowline import flowline
from GA_flowline import genetic_algorithm


np.random.seed(42)
random.seed(42)


### machine setup
machines = 2

mu_normal = 1
p_normal = 0.005
gamma_normal = 0.05
costs_per_buffer = 1
costs_per_spare = 1

C_min = 1
C_max = 30
Q_min = 1
Q_max = 5

target_throughput = 0.75


### init
line = flowline()
line.fill(machines, C_max, Q_max, mu_normal, p_normal, gamma_normal, costs_per_buffer, costs_per_spare)
# line.mu = np.array([1.05, 0.96], dtype=float)
# line.p = np.array([0.0055, 0.0045], dtype=float)
# line.gamma = np.array([0.049, 0.051], dtype=float)
line.mu = np.array([0.987454011884736, 1.01981617140197], dtype=float)
line.p = np.array([0.00468513292883862, 0.00501908178513622], dtype=float)
line.gamma = np.array([0.0476170568373591, 0.053275189221063], dtype=float)
line.set_limits(C_min, C_max, Q_min, Q_max)


### check solution
line.C = np.array([2])
line.Q = np.array([2, 2])
line.solve()
result = line.characteristics
print(result.TP)


### calc
# plt.ioff()
# mpl.use('Agg')
result = genetic_algorithm(line, target_throughput, False, True)
# plt.ion()
# mpl.use('TkAgg')


### output
print(f'C = {result.C}')
print(f'Q = {result.Q}')
print(f'TP = {result.TP:.4f}')
print(f'costs = {result.costs:.4f}')
print(f'counter_iterations = {result.counter_iterations}')
print(f'counter_evaluations = {result.counter_evaluations}')
print(f'counter_problems = {result.counter_problems}')


### plots
# img1 = plt.imread('1.png')
# img2 = plt.imread('2.png')
# img3 = plt.imread('3.png')
# img4 = plt.imread('4.png')
#
# fig, axs = plt.subplots(2, 2)
# axs[0, 0].imshow(img1)
# axs[1, 0].imshow(img2)
# axs[0, 1].imshow(img3)
# axs[1, 1].imshow(img4)
# plt.axis('off')
# plt.show()
