function Q = getQ(states, epsilon)
Q = normrnd(0, 1, states, states);

specnorm = max(abs(eig(Q)));
if specnorm >= 1
    Q = Q.*(1/specnorm - epsilon);
end
maxrow = max(sum(abs(Q), 2));
if maxrow > 1
    Q = Q.*(1/(maxrow + 0.05));
end

end