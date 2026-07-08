function C = matrixC(q,qdot)
q1 = q(1);  q2 = q(2); d3 = q(3); q4 = q(4);
q1_dot=qdot(1); q2_dot = qdot(2); d3_dot = qdot(3); q4_dot = qdot(4);
global m1 I1 l1 m2 I2 l2 m3 I3 l3 m4 I4 l4 gr

c11 = -(m2 + 2*m3 + 2*m4)*l1*l2*sin(q2)*q2_dot;
c12 = 0;
c21 = (m2/2 + m3 + m4)*l1*l2*sin(q2)*q1_dot;
c22 = -(m2/2 + m3 + m4)*l1*l2*sin(q2)*q2_dot;

C = [c11 c12   0   0;
     c21 c22   0   0;
     0   0    0   0;
     0   0    0   0];
end