clear
clc
% Declare parameter as symbolic objects
syms H I P Q s a h1 h2 e B1 B2 p r1 r2 d t

%Define control model
dH=s*(1-a*h1*I/(H+I))-e*H+B1*H*(P+Q);
dI=s*a*h1*I/(H+I)-e*I+p*I+B1*I*(P+Q);
dP=r1*P*(H+I)+(1-h2)*r2*Q*(H+I)-e*P+B2*P*(H+I)-d*P;
dQ=r2*h2*Q*(H+I)-e*Q+B2*Q*(H+I)-d*Q;

%Determine the steady states of the control model
[Hsol,Isol,Psol,Qsol]=solve(dH==0,dI==0,dP==0,dQ==0,H,I,P,Q)

%Establish these as the non-trivial equilibrium points
Hbar=(s - a*h1*s)/p;
Ibar=(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2));
Pbar=((r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2);
Qbar=(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s)/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2);

%Calculate the Jacobian Matrix of the control model
J=jacobian([dH,dI,dP,dQ],[H,I,P,Q])

%Evaluate the Jacobian matrix at the non-trivial equilibrium points
J1=subs(J,[H,I,P,Q],[Hbar,Ibar,Pbar,Qbar])

%Empty Vectors to store resulst in, i.e. the number of positive steady
%states and the positive steady states that are stable

Pos=[];
PosStab=[];

%creates 10^7 random parameter combinations in a give range
for i=1:10000000
    %Generate random parameter combinations
    a=rand; B1=-10*(rand-0.5); B2=-10*(rand-0.5); d=rand*10; p=rand*10; e=rand*10; r1=rand*10; r2=rand*10; s=rand*10; A=1+rand; theta=rand;
    h1=theta;
    h2=A*(theta-1)/(theta-A);
    
    %Numerical evaluation of the Jacobian matrix evaluated at the
    %non-trivial equilibrium point
    j1=[ B1*((e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s)/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) + ((r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2)) - e + (a*h1*s*(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s))/(p*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2)))^2*(B2 + h2*r2)),                                                                                                                                                                                                                                                                                                                             -s*((a*h1)/((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2))) - (a*h1*(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s))/(p*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2)))^2*(B2 + h2*r2))),                                                                                                                                                                                     (B1*(s - a*h1*s))/p,                                                                                                                                                                                        (B1*(s - a*h1*s))/p];
    j2=[                                                                                                                                                                                                                                                                                                                         -(a*h1*s*(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s))/(p*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2)))^2*(B2 + h2*r2)), p - e + B1*((e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s)/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) + ((r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2)) + (a*h1*s)/((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2))) - (a*h1*s*(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s))/(p*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2)))^2*(B2 + h2*r2)),                                                                                                                           (B1*(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s))/(p*(B2 + h2*r2)),                                                                                                                              (B1*(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s))/(p*(B2 + h2*r2))];
    j3=[                                                            (B2*(r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) - (r2*(h2 - 1)*(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) + (r1*(r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2),                                                                                                                                                                       (B2*(r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) - (r2*(h2 - 1)*(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) + (r1*(r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2), r1*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2))) - e - d + B2*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2))),                                                                                                   -r2*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2)))*(h2 - 1)];
    j4=[                                                                                       (B2*(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) + (h2*r2*(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2),                                                                                                                                                                                                  (B2*(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2) + (h2*r2*(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2),                                                                                                                                                                                                       0, B2*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2))) - e - d + h2*r2*((s - a*h1*s)/p + (d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2)))];
    
    J=[j1;j2;j3;j4];
    
    Hbar=(s - a*h1*s)/p;
    Ibar=(d*p - B2*s + e*p - h2*r2*s + B2*a*h1*s + a*h1*h2*r2*s)/(p*(B2 + h2*r2));
    Pbar=((r2 - h2*r2)*(d*p - d*e + e*p - e^2 + B2*a*h1*s + a*h1*h2*r2*s))/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2);
    Qbar=(e^2*r1 + d*e*r1 - d*p*r1 - e*p*r1 - e^2*h2*r2 - d*e*h2*r2 + d*h2*p*r2 + e*h2*p*r2 + a*h1*h2^2*r2^2*s - B2*a*h1*r1*s + B2*a*h1*h2*r2*s - a*h1*h2*r1*r2*s)/(B1*d*r1 - B1*d*r2 + B1*e*r1 - B1*e*r2);
    %determine if the equilibrium points are all positive
    if Hbar>0 && Ibar>0 && Pbar>0 && Qbar>0 
        Pos=[Pos;1];
        
        %Determine whether the equilibrium points are stable if they are
        %positive
        if sum(sign(real(eig(J))))==-4
            PosStab=[PosStab;1];
        end
    end
        
end
%Calculates the proportion of positive steady states that are stable
length(PosStab)/length(Pos)
