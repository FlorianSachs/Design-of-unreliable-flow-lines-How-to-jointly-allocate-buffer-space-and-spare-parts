import numpy as np
from flowline import flowline


### input
line = flowline()
line.fill(3, 2, 1, 1, .01, .1, 1, 1)


### calc
line.solve()
result = line.characteristics


### output
print(result.TP)
