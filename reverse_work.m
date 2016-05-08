function [Qest, R] = reverse_work(Q, target, lmax, maxresid)

[states, ~] = size(Q);

Qest = zeros(lmax, states);
R = zeros(lmax, states);
R(1, target) = 1;

while abs(max(max(R)))>maxresid
    [Qest, R] = reverse_push(Qest, R, Q);
end

end

function [Qest, R] = reverse_push(Qest, R, P)

%identify which index to push
[maxcol, colind] = max(abs(R));
[~ ,ind] = max(maxcol);
row = colind(ind);


residual = R(row, ind);
Qest(row, ind) = Qest(row, ind) + residual;
R(row, ind) = 0;

[lmax, ~] = size(Qest);
if (row < lmax)
    R(row+1, :)= R(row+1, :)+residual.*(P(ind, :));
end

end

