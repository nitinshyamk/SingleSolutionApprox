# MatrixSolver.py
# Nitin Shyamkumar, Siddhartha Banerjee
# This python file contains classes for testing the bidirectional randomized
# algorithm for approximating a single component of the solution x
# to the matrix equation x = Gx + z, assuming that G is substochastic 
# and that the maximum rowsum of abs(G) < 1

import numpy as np
import numpy.random as rd


'''System is just a representation of the problem '''
class System:
    #G should be an n x n numpy array representing the matrix
    #z should be a length n vector
    #i is the solution component you want to solve for
    def __init__(self, G, z, target):
        self.G = G;
        self.z = z;
        self.target = target;


import csv
class FileReader:
    '''
    File Reader reads the system x = Gx + z from a directory
    containing the following files file in the following format:
        
        G.dat: n x n matrix G (CSV format)
        z.dat: 1 x n matrix z (vector with n entries)
        info.dat: 1 x 2 matrix (n, target)
        
        n: size of matrix/vector
        target: entry for which to solve for
        G: substochastic matrix 
        z: the solution to the matrix equation
'''
    ''' read_matrix takes in a file name argument '''
    @staticmethod
    def read_system(directory):
        matrix = [];
        with open(directory +'/G.dat', 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter = ',');
            for row in reader:
                row = [float(e) for e in row]
                matrix.append(row);
        G = np.array(matrix);

        z = [];
        with open(directory + '/z.dat', 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter = ',');
            for row in reader: 
                z = [float(e) for e in row];
        z = np.transpose(np.array(z));
        return (G, z);


class ReverseWork(System):
    '''
    Fields: 
        inherits same fields from System, additionally:

        lmax: maximum length to perform reverse push operation
        residuals: matrix of residuals - dimension: lmax x states
        qest: the estimate matrix - dimension: lmax x states
            contains the estimates from the resulting reverse push operation
        qest_v: the vector generated from the qest matrix by summing across lmax
    '''
    def __init__(self, G, z, target, lmax):
        System.__init__(self, G, z, target);
        self.lmax = lmax;
        self.residuals = np.zeros((lmax+1, self.z.shape[0]), dtype = np.float);
        self._residuals_abs = np.zeros((lmax+1, self.z.shape[0]), dtype = np.float);
        self.qest = np.zeros(self.residuals.shape, dtype = np.float);
        self.qest_v = np.zeros(self.z.shape)
        self.residuals[0, self.target] = 1.0;
        self._residuals_abs[0, self.target] = 1.0;


    '''
    reverse_push() performs a single reverse push operation maintaining 
    residual and estimate invariants
    '''
    def reverse_push(self):
        #find correct indexing
        ind = np.argmax(self._residuals_abs);
        l = ind/self.z.shape[0];
        state = ind - l*self.z.shape[0];

        #update qest vector and residual (individual cell operations)
        self.qest[l, state] = self.qest[l, state] + self.residuals[l, state];
        topush = self.residuals[l, state];

        self.residuals[l, state] = 0;
        self._residuals_abs[l, state] = 0;

        #now perform vector operations
        if (l < self.lmax):
            self.residuals[l+1, :] = self.residuals[l+1, :] + topush * self.G[state, :];
            self._residuals_abs[l+1, :] = np.absolute(self.residuals[l+1, :]);

    '''
    reverse_work performs a sequence of reverse_push operations until
    magnitude of residual is less than tol (maxresidual)
    '''
    def reverse_work(self, tol):
        maxresid = np.max(self._residuals_abs);
        while (maxresid > tol):
            self.reverse_push();
            maxresid = np.max(self._residuals_abs);
        self.qest_v = np.sum(self.qest, axis = 0);


''' modified binary search to find the appropriate index i
    for which v[i-1] <val < v[i]'''
def bin_search(vec, val):
    def bin_search_help(vec, val, i, j):
        index = (i+j)/2;
        if (j-i)==1:
            if val < vec[i]:
                return i;
            else:
                return j
        elif vec[index] < val :
            return bin_search_help(vec, val, index, j);
        else:
            return bin_search_help(vec, val, i, index);

    return bin_search_help(vec, val, 0, vec.shape[0]);



class ForwardWork(System):
    '''
    ForwardWork(G, z, target, residuals, lmax)
    constructs an object to simulate forward random walks from z 
    to the residuals 
    '''
    def __init__(self, G, z, target, residuals, lmax):
        System.__init__(self, G, z, target);
        Gabs = np.absolute(G);
        rowsums = np.sum(Gabs, axis = 1);
        self.rowsums = rowsums;
        self.Gprob = Gabs / rowsums[:, None] #normalize
        self.Gprob = np.cumsum(self.Gprob, axis = 1);
        self.Gsign = np.sign(G);
        self.residuals = residuals;
        self.lmax = lmax;


    '''
    runs a walk of length k and get the residual at length - k
    '''
    def run_k_fwalk(self, ind, k, length):
        score = 1.0;
        for i in range(k):
            next = bin_search(self.Gprob[ind, :], rd.random());
            score *= self.rowsums[ind]*self.Gsign[ind, next];
            ind = next;
        return score*self.residuals[length - k, ind];

    '''
    perform numwalks forward walks drawing
    length ~ [1, lmax]
    and k ~[0, length]
    '''
    def forward_work(self, numwalks):
        normalize = np.sum(np.absolute(self.z));
        sample = np.cumsum(np.absolute(self.z)/normalize);
        accum = 0.0;
        for length in rd.randint(self.lmax+1, size = numwalks):
            start = bin_search(sample, rd.random());
            if length > 0:
                for k in rd.randint(length, size = 1):
                    accum+= self.lmax*np.sign(self.z[start])*normalize*self.run_k_fwalk(start, k, length)
        print accum/float(numwalks);
        return accum/float(numwalks) + self.z;

import time
class Solver:
    def __init__(self, prob, lmax):
        #assumes prob has type System
        self.prob = prob
        self.rev = ReverseWork(prob.G, prob.z, prob.target, lmax);
        self.forward = ForwardWork(prob.G, prob.z, prob.target,
            self.rev.residuals, lmax);

    '''
    solves problem using only
     reverse_push operations until residuals < maxresidual 
    '''
    def reverse_only(self, maxresid):
        self.rev.reverse_work(maxresid);
        return np.sum(self.rev.qest_v*self.prob.z);

    '''
    does reverse push operations until maxresid tolerance is 
    reached, and then does numwalks forward walks
    '''
    def forward_back(self, maxresid, numwalks):
        self.rev.reverse_work(maxresid);
        const = np.sum(self.rev.qest_v*self.prob.z);
        score = self.forward.forward_work(numwalks)
        return score + const;

    '''
    balanced not implemented yet
    '''
    def both_balanced(self):
        start = time.time();
        self.forward_work(10);
        end = time.time();


#time forward walks
# total time is nt 
# keep track of time spent on reverse push
# keep track of time spent on 

