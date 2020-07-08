n=6;
kk=QQ[d12,d13,d14,d23,d24,d34,e12,e13,e14,e23,e24,e34];

M=matrix{{2*d14, d14+d24-d12,d14+d34-d13},{d14+d24-d12, 2*d24, d24+d34-d23},{d14+d34-d13,d24+d34-d23,2*d34}};
I=minors(3,M);
J=ideal(d12-d23,d23-d34,d34-14);
c=codim I;
Y=matrix{{d12,d13,d14,d23,d24,d34}}-matrix{{e12,e13,e14,e23,e24,e34}};
Jac= jacobian gens I;
S=submatrix(Jac,{0..n-1},{0..numgens(I)-1});
Jbar=S|transpose(Y);
EX = I + minors(c+1,Jbar);
SingX=radical I+minors(c,Jac);
--SingX=ideal(d12,d13,d14,d23,d24,d34);
EXreg=saturate(EX,SingX);

DL=eliminate(toList(d12,d13,d14,d23,d24,d34),EXreg+J);
to String oo