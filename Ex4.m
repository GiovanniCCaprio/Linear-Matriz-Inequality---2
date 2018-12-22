A = -1;
B = 2;
C = 1;
D = 0;

sys = ss(A,B,C,D);

norm(sys,2)
norm(sys,Inf)