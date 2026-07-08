function M = matrixM(q)
q1 = q(1);  q2 = q(2); d3 = q(3); q4 = q(4);
global m1 I1 l1 m2 I2 l2 m3 I3 l3 m4 I4 l4 gr

m11 = ((1/4*m1 + m2 + m3 + m4)*l1^2 + (1/4*m2 + m3 + m4)*l2^2 + (m2+2*m3 + 2*m4)*l1*l2*cos(q2) + I1 + I2 + I3 + I4);
m12 = ((1/4*m2 + m3 + m4)*l2^2 + (1/2*m2 + m3 + m4)*l1*l2*cos(q2) + I2 + I3 + I4);
m13 =  0;
m14 = -I4;
m21 = m12;
m22 = ((1/4*m2 + m3 + m4)*l2^2 + I2 + I3 + I4);
m23 = 0;
m24 = m14;
m31 = 0;
m32 = 0;
m33 = m3 + m4;
m34 = 0;
m41 = m14;
m42 = m14;
m43 = 0;
m44 = m14;

M = [m11 m12 m13 m14;m21 m22 m23 m24; m31 m32 m33 m34;m41 m42 m43 m44];
end