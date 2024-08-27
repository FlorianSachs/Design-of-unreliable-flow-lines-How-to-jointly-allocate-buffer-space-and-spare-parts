import numpy as np
from flowline import flowline


### input
line = flowline()
line.fill(2, 2, 1, 1, .01, .1, 1, 1)
line.fill(2, 50, 5, 1, .005, .05, 1, 1)

line.mu = np.array([0.987454011884736, 1.01981617140197], dtype=float)
line.p = np.array([0.00468513292883862, 0.00501908178513622], dtype=float)
line.gamma = np.array([0.0476170568373591, 0.053275189221063], dtype=float)


### calc
line.solve()
result = line.characteristics


### output
print(result.TP)
