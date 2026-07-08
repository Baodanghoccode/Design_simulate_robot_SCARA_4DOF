function q = traj_planning_4dof(t)
% 1. Thong so thoi gian di chuyen
    T = 10; 

    % 2. Vi tri ban dau cua 4 khop
    q10 = 2.62;
    q20 = 4;
    r30 = 0.01;
    q40 = 3.3336*pi/4;

    % 3. Vi tri ket thuc cua 4 khop
    q1f = 3;
    q2f = 2.53;
    r3f = 0.06;
    q4f = pi/4;


    % 5. Tinh toan quy dao vi tri cho tung khop theo da thuc bac 3
    qd1 = q10 + ((3/T^2)*(q1f - q10))*t^2 + ((-2/T^3)*(q1f - q10))*t^3;
    qd2 = q20 + ((3/T^2)*(q2f - q20))*t^2 + ((-2/T^3)*(q2f - q20))*t^3;
    qd3 = r30 + ((3/T^2)*(r3f - r30))*t^2 + ((-2/T^3)*(r3f - r30))*t^3;
    qd4 = q40 + ((3/T^2)*(q4f - q40))*t^2 + ((-2/T^3)*(q4f - q40))*t^3;

    % 6. Xuat ra vector cot 4x1 de dua vao Simulink
    q = [qd1; qd2; qd3; qd4];
end

