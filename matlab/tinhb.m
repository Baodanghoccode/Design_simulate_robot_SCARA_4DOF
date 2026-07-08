function b = tinhb(ins)
global m1 I1 l1 m2 I2 l2 m3 I3 l3 m4 I4 l4 gr
n = 4;
q = ins(1:n); qdot = ins(n+1:2*n);
C = matrixC(q,qdot);
G = -(m3+m4);

b = C*qdot+G;
end