from MatrixSolver import *
import numpy as np;
import time


SYSTEM = 'Systems/t002'
TARGET = 3
RTOLERANCE = 1e-6
TOLERANCE = 1e-1

LMAX = 200
NUMWALKS = 10000

G, z = FileReader.read_system(SYSTEM);

z = np.zeros(G.shape[0]);
z[target] = 1;

problem = System(G, z, TARGET);
solver = Solver(problem, LMAX);

start = time.time()
print solver.reverse_only(RTOLERANCE);
end = time.time();
print (end - start);

solver = Solver(problem, LMAX);
start = time.time()
print solver.reverse_only(TOLERANCE);
print solver.forward_back(TOLERANCE, NUMWALKS);
end = time.time();
print (end - start);

'''
test e_t^TG^ke_s
'''