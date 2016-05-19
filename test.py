from MatrixSolver import *
import numpy as np;
import time


SYSTEM = 'Systems/t004'
TARGET = 3
RTOLERANCE = 1e-6
TOLERANCE = 1e-3

LMAX = 10
NUMWALKS = 100000

G, z = FileReader.read_system(SYSTEM);
print (G.shape[0]);

# z = np.zeros(G.shape[0]);
# z[0] = 1;

problem = System(G, z, TARGET);
solver = Solver(problem, LMAX);

start = time.time()
print 'reverse only: '+str(solver.reverse_only(RTOLERANCE))
end = time.time();
print 'time: '+str(end - start);

solver = Solver(problem, LMAX);
start = time.time()
print 'reverse component: ' +str(solver.reverse_only(TOLERANCE))
print 'final: '+str(solver.forward_back(TOLERANCE, NUMWALKS));
end = time.time();
print 'time: '+ str((end - start))

'''
test e_t^TG^ke_s
do full summation from k to L 

actually perform the full summation
only expectation is over z 

'''
