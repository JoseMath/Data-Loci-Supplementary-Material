n=9;
kk=QQ[x_1..x_n,u_1..u_n];
f1=x_2-x_4;
f2=x_3-x_5;
f3=x_3-x_7;
f4=x_6-x_8;
Hankel=ideal(f1,f2,f3,f4);
Help=matrix{{x_1..x_3},{x_4..x_6},{x_7..x_9}}
Rank=minors(2,Help);
X=Rank;
c=codim X;
A=X+Hankel;

Y=matrix{{x_1-u_1,x_2-u_2,x_3-u_3,x_4-u_4,x_5-u_5,x_6-u_6,x_7-u_7,x_8-u_8,x_9-u_9}};
Jac= jacobian gens X;
S=submatrix(Jac,{0..n-1},{0..numgens(X)-1});
Jbar=S|transpose(Y);
projGammaCorr= X + minors(c+1,Jbar);

SingX=ideal(x_1,x_2,x_3,x_4,x_5,x_6,x_7,x_8,x_9);
projGammaCorrRegular=saturate(projGammaCorr,SingX);

PreimDL=projGammaCorrRegular+Hankel;
DLA=eliminate(toList(x_1..x_n),PreimDL);

g1=u_2-u_4;
g2=u_3-u_5;
g3=u_3-u_7;
g4=u_6-u_8;
HankelU=ideal(g1,g2,g3,g4);

DLAH=DLA+HankelU;
