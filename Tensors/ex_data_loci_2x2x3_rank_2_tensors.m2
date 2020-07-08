n=2*2*3;
kk=QQ[x_1..x_n,y_1..y_n];

J=ideal(x_9*x_11-x_8*x_12, x_6*x_11-x_5*x_12, x_3*x_11-x_2*x_12, x_9*x_10-x_7*x_12, x_8*x_10-x_7*x_11,x_6*x_10-x_4*x_12, x_5*x_10-x_4*x_11, x_3*x_10-x_1*x_12, x_2*x_10-x_1*x_11, x_6*x_9-x_3*x_12, x_5*x_9-x_2*x_12,x_4*x_9-x_1*x_12, x_6*x_8-x_2*x_12, x_5*x_8-x_2*x_11, x_4*x_8-x_1*x_11, x_3*x_8-x_2*x_9, x_6*x_7-x_1*x_12,x_5*x_7-x_1*x_11, x_4*x_7-x_1*x_10, x_3*x_7-x_1*x_9, x_2*x_7-x_1*x_8, x_3*x_5-x_2*x_6, x_3*x_4-x_1*x_6,x_2*x_4-x_1*x_5);

I=ideal(x_6*x_8*x_10-x_5*x_9*x_10-x_6*x_7*x_11+x_4*x_9*x_11+x_5*x_7*x_12-x_4*x_8*x_12,x_3*x_8*x_10-x_2*x_9*x_10-x_3*x_7*x_11+x_1*x_9*x_11+x_2*x_7*x_12-x_1*x_8*x_12,x_3*x_5*x_10-x_2*x_6*x_10-x_3*x_4*x_11+x_1*x_6*x_11+x_2*x_4*x_12-x_1*x_5*x_12,x_3*x_5*x_7-x_2*x_6*x_7-x_3*x_4*x_8+x_1*x_6*x_8+x_2*x_4*x_9-x_1*x_5*x_9);

SingX=ideal(x_9*x_11-x_8*x_12, x_6*x_11-x_5*x_12, x_3*x_11-x_2*x_12, x_9*x_10-x_7*x_12, x_8*x_10-x_7*x_11,x_6*x_10-x_4*x_12, x_5*x_10-x_4*x_11, x_3*x_10-x_1*x_12, x_2*x_10-x_1*x_11, x_6*x_8-x_5*x_9, x_3*x_8-x_2*x_9,x_6*x_7-x_4*x_9, x_5*x_7-x_4*x_8, x_3*x_7-x_1*x_9, x_2*x_7-x_1*x_8, x_3*x_5-x_2*x_6, x_3*x_4-x_1*x_6,x_2*x_4-x_1*x_5);

c=codim I;
Y=matrix{{x_1..x_n}}-matrix{{y_1..y_n}};
Jac= jacobian gens I;
S=submatrix(Jac,{0..n-1},{0..numgens(I)-1});
Jbar=S|transpose(Y);
EX = I + minors(c+1,Jbar);
EXreg=saturate(EX,SingX);

DLRk1=eliminate(toList(x_1..x_n),EXreg+J);

DLSingRk2=eliminate(toList(x_1..x_n),EXreg+SingX);--- ==I 