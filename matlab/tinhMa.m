function Ma = tinhMa(ins)
n = 4;
a = ins(1:n); q = ins(n+1:2*n);
M = matrixM(q);
Ma = M*a;
end
