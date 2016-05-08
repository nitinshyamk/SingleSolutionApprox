function write_system(directory, states, epsilon, target)
mkdir(directory);
G = getQ(states, epsilon);
z = normrnd(0, 2, 1, states);

csvwrite([directory, '/G.dat'], G);
csvwrite([directory, '/z.dat'], z);
csvwrite([directory, '/info.dat'], [states, target]);

end

