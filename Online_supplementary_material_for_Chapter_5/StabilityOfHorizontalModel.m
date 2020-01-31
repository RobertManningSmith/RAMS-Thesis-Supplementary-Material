clear
clc
 
E=[]; % empty vector to store eigenvalues in
SSCATCH=[]; % empty vector to store the value of the steady states
% for loop to generate random parameter samples within the parameter space
for i=1:10000000
    A=1+rand; % trade-off curve shaping parameter
    theta=rand; % trade-off parameter value
    h1=theta; % trade-off function attached to transmission through the plant
    h2=A*(theta-1)/(theta-A); % trade-off function attached to transmission through the fungus
   
    
    a=rand; % Vertical transmission efficiency between 0 and 1
    B1=10*(-1+rand)*rand; % Fitness effect of interacting with the fungus experienced by the plant bound between -5 and 5
    B2=10*(-1+rand)*rand; % Fitness effect of interacting with the plant experienced by the fungusbound between -5 and 5
    d=rand*10; % Fungus mortality rate bound between 0 and 10
    e=rand*10; % Plant mortality rate bound between 0 and 10
    p=rand*10; % Fitness benefit experienced by the plant from interaction with the virus bound between 0 and 10
    r2=rand*10; % Infected fungus growth rate bound between 0 and 10
    r1=rand*r2; % Uninfected fungus growth rate bound between 0 and 10, with r1<r2
    s=rand*10; % Planting rate bound between 0 and 10
    t=rand; % Rate of loss of virus in infected fungi bound between 0 and 1
    z=rand; % Rate of horizontal transmission of infection in the fungus population
    
    
    % Numerical solutions found for the resisdent system to estimate the non-trivial steady state of the resident system    
    
    % Time span the system is solved for
    tspan=[0,1000]; 
    % Express the system as a function which MATLAB can interpret
    f=@(t,y)[s*(1-a*h1*y(2)/(y(1)+y(2)))-e*y(1)+B1*y(1)*(y(3)+y(4));s*a*h1*y(2)/(y(1)+y(2))-e*y(2)+p*y(2)+B1*y(2)*(y(3)+y(4));r1*y(3)*(y(1)+y(2))-e*y(3)-d*y(3)+B2*y(3)*(y(1)+y(2))+t*y(4)-z*h2*y(3)*y(4);r2*y(4)*(y(1)+y(2))-e*y(4)-d*y(4)+B2*y(4)*(y(1)+y(2))-t*y(4)+z*h2*y(3)*y(4)];
    % Arbitrary initial conditions
    finit=[1,1,1,1];
    % Solve system using ODE solver
    [T,y]=ode15s(f,tspan,finit);
    H=y(:,1); % Solutions for uninfected plant host
    I=y(:,2); % Solutions for infected plant host
    P=y(:,3); % Solutions for uninfected fungus host
    Q=y(:,4); % Solutions for infected fungus host
 
    
    Hbar=H(end); % Estimate of steady state for uninfected plant host
    Ibar=I(end); % Estimate of steady state for infected plant host
    Pbar=P(end); % Estimate of steady state for uninfected fungus host
    Qbar=Q(end); % Estimate of steady state for infected fungus host
    
    % if any of the estimates of the steady states are negative, skip the
    % loop
    if Hbar<0
        continue
    end
    if Ibar<0
        continue
    end
    if Pbar<0
        continue
    end
    if Qbar<0
        continue
    end
    steadystates=[Hbar,Ibar,Pbar,Qbar];
    SSCATCH=[SSCATCH;steadystates];
    TF1=isinf(steadystates);
    % if any of the steady states are valued at infinity, skip the loop
    if sum(TF1)>0
        continue
    end
    % if any of the steady states are valued as NaN, skip the loop
    TF2=isnan(steadystates);
    if sum(TF2)>0
        continue
    end
    % if all the steady states are positive, calculate the jacobian matrix
    if Hbar>0 && Ibar>0 && Pbar>0 && Qbar>0
    
    j11=[ B1*(Pbar + Qbar) - e + (Ibar*a*h1*s)/(Hbar + Ibar)^2,                     -s*((a*h1)/(Hbar + Ibar) - (Ibar*a*h1)/(Hbar + Ibar)^2),B1*Hbar,B1*Hbar];
    j12=[-(Ibar*a*h1*s)/(Hbar + Ibar)^2, p - e + B1*(Pbar + Qbar) + (a*h1*s)/(Hbar + Ibar) - (Ibar*a*h1*s)/(Hbar + Ibar)^2,B1*Ibar,B1*Ibar];
    j13=[Pbar*r1 + B2*Pbar, Pbar*r1 + B2*Pbar, B2*(Hbar + Ibar) - e - d + r1*(Hbar + Ibar) - Qbar*h2*z,t - Pbar*h2*z];
    j14=[Qbar*r2 + B2*Qbar, Qbar*r2 + B2*Qbar, Qbar*h2*z, B2*(Hbar + Ibar) - e - t - d + r2*(Hbar + Ibar) + Pbar*h2*z];
    
    % these if statements calculate if any part of Jacobian matrix is
    % valued at infinity and skips the loop if so
    TF3=isinf(j11);
    if sum(TF3)>0
        continue
    end
    TF4=isinf(j12);
    if sum(TF4)>0
        continue
    end
    TF5=isinf(j13);
    if sum(TF5)>0
        continue
    end
    TF6=isinf(j14);
    if sum(TF6)>0
        continue
    end
    % these if statements calculate if any part of Jacobian matrix is
    % valued as NaN and skips the loop if so
    DF3=isnan(j11);
    if sum(DF3)>0
        continue
    end
    DF4=isnan(j12);
    if sum(DF4)>0
        continue
    end
    DF5=isnan(j13);
    if sum(DF5)>0
        continue
    end
    DF6=isnan(j14);
    if sum(DF6)>0
        continue
    end
    Jnum=[j11;j12;j13;j14];
    % Vector to store the relevant information about the stability of the      non-trivial steady state, specifically if the sign of the real part of all eigenvalues of the Jacobian matrix are negative, then the steady state is stable.
    E=[E;sum(sign(real(eig(Jnum))))];
    end
end
% This calculates the proportion of positive steady states which are stable, counting number of stable equilibrium points and dividing that by the number of positive steady states.
propstab1=sum(E(:)==-4)/length(E)
 
 
 



