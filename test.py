from MatrixSolver import *
import numpy as np;
import time


SYSTEM = 'Systems/t003'
TARGET = 2
RTOLERANCE = 1e-4
TOLERANCE = 1e-1

LMAX = 100
NUMWALKS = 10000

G, z = FileReader.read_system(SYSTEM);

problem = System(G, z, TARGET);
solver = Solver(problem, LMAX);

start = time.time()
print solver.reverse_only(TOLERANCE);
end = time.time();
print (end - start);

start = time.time()
print solver.forward_back(TOLERANCE, NUMWALKS);
end = time.time();
print (end - start);