clear all
clc
T=100; % Time length of the numerical solution of the resident system in absense of the invader system 
tspan=[0,T];
params=[];


    %Loop determined by m determines how many PIPs are created
for m=1:200
    
    %Loop determined by k parameterises the time-varying parameter model
    
for k=1:100
    
    
    
    %Loop determined by i finds parameters which will likely give a non-trivial equibilibrium
    %point, using the control model as a basis
for i=1:100
    
    alpha=rand; % Vertical transmission efficiency between 0 and 1
    beta1=10*(-1+rand)*rand; % Fitness effect of interacting with the fungus experienced by the plant bound between -5 and 5
    beta2=10*(-1+rand)*rand; % Fitness effect of interacting with the plant experienced by the fungusbound between -5 and 5
    delta=rand*10; % Fungus mortality rate bound between 0 and 10
    eta=rand*10; % Plant mortality rate bound between 0 and 10
    phi=rand*10; % Fitness benefit experienced by the plant from interaction with the virus bound between 0 and 10
    r2=rand*10; % Infected fungus growth rate bound between 0 and 10
    r1=rand*r2; % Uninfected fungus growth rate bound between 0 and 10, with r1<r2
    sigma=rand*10; % Planting rate bound between 0 and 10
    tau=rand; % Rate of loss of virus in infected fungi bound between 0 and 1
    
 
    A=1+rand; % trade-off curve shaping parameter
    theta=rand; % trade-off parameter value
    h1=theta; % trade-off function attached to transmission through the plant
    h2=A*(theta-1)/(theta-A); % trade-off function attached to transmission through the fungus
    
    % Control model equilibrium points
    Hbar=(sigma - alpha*h1*sigma)/phi;
    Ibar=(delta*phi + eta*phi - beta2*sigma + phi*tau - h2*r2*sigma + alpha*beta2*h1*sigma + alpha*h1*h2*r2*sigma)/(phi*(beta2 + h2*r2));
    Pbar=(tau*(beta2 + h2*r2)*(delta*phi - delta*eta + eta*phi - eta*tau + phi*tau - eta^2 + alpha*beta2*h1*sigma + alpha*h1*h2*r2*sigma))/(beta1*r1*tau^2 + beta1*delta^2*r1 + beta1*eta^2*r1 + 2*beta1*delta*r1*tau + 2*beta1*eta*r1*tau - beta1*delta^2*h2*r2 - beta1*eta^2*h2*r2 - beta1*h2*r2*tau^2 + 2*beta1*delta*eta*r1 - 2*beta1*delta*eta*h2*r2 - 2*beta1*delta*h2*r2*tau - 2*beta1*eta*h2*r2*tau);
    Qbar= -(alpha*h1*sigma*beta2^2*tau - beta2*delta*eta*tau - alpha*h1*sigma*beta2*delta*h2*r2 + phi*beta2*delta*tau + alpha*h1*r1*sigma*beta2*delta - beta2*eta^2*tau - alpha*h1*sigma*beta2*eta*h2*r2 - beta2*eta*tau^2 + phi*beta2*eta*tau + alpha*h1*r1*sigma*beta2*eta + alpha*h1*sigma*beta2*h2*r2*tau + phi*beta2*tau^2 + alpha*h1*r1*sigma*beta2*tau + delta^2*eta*h2*r2 - r1*delta^2*eta - phi*delta^2*h2*r2 + phi*r1*delta^2 + 2*delta*eta^2*h2*r2 - 2*r1*delta*eta^2 + delta*eta*h2*r2*tau - 2*phi*delta*eta*h2*r2 - 2*r1*delta*eta*tau + 2*phi*r1*delta*eta - alpha*h1*sigma*delta*h2^2*r2^2 - phi*delta*h2*r2*tau + alpha*h1*r1*sigma*delta*h2*r2 + 2*phi*r1*delta*tau + eta^3*h2*r2 - r1*eta^3 + eta^2*h2*r2*tau - phi*eta^2*h2*r2 - 2*r1*eta^2*tau + phi*r1*eta^2 - alpha*h1*sigma*eta*h2^2*r2^2 - phi*eta*h2*r2*tau + alpha*h1*r1*sigma*eta*h2*r2 - r1*eta*tau^2 + 2*phi*r1*eta*tau + alpha*h1*r1*sigma*h2*r2*tau + phi*r1*tau^2)/(beta1*r1*tau^2 + beta1*delta^2*r1 + beta1*eta^2*r1 + 2*beta1*delta*r1*tau + 2*beta1*eta*r1*tau - beta1*delta^2*h2*r2 - beta1*eta^2*h2*r2 - beta1*h2*r2*tau^2 + 2*beta1*delta*eta*r1 - 2*beta1*delta*eta*h2*r2 - 2*beta1*delta*h2*r2*tau - 2*beta1*eta*h2*r2*tau);
    % determines whether the non-trivial steady states are positive
    if Hbar>0 && Ibar>0 && Pbar>0 && Qbar>0
    
    j11=beta1*(Pbar + Qbar) - eta + (Ibar*alpha*h1*sigma)/(Hbar + Ibar)^2;
    j12=-sigma*((alpha*h1)/(Hbar + Ibar) - (Ibar*alpha*h1)/(Hbar + Ibar)^2);
    j13=Hbar*beta1;
    j14=Hbar*beta1;
    j21=-(Ibar*alpha*h1*sigma)/(Hbar + Ibar)^2;
    j22=phi - eta + beta1*(Pbar + Qbar) + (alpha*h1*sigma)/(Hbar + Ibar) - (Ibar*alpha*h1*sigma)/(Hbar + Ibar)^2;
    j23=Ibar*beta1;
    j24=Ibar*beta1;
    j31=Pbar*beta2 + Pbar*r1;
    j32=Pbar*beta2 + Pbar*r1;
    j33=beta2*(Hbar + Ibar) - eta - delta + r1*(Hbar + Ibar);
    j34=tau;
    j41=Qbar*beta2 + Qbar*h2*r2;
    j42=Qbar*beta2 + Qbar*h2*r2;
    j43=0;
    j44=beta2*(Hbar + Ibar) - eta - tau - delta + h2*r2*(Hbar + Ibar);
    J=[j11,j12,j13,j14;j21,j22,j23,j24;j31,j32,j33,j34;j41,j42,j43,j44];
    % determines whether the non-rtivial steady states are stable
 
    if sum(sign(real(eig(J))))==-4
        break
    end
    end
end
 
%Begin simulating the resident system dynamics, finding a chaotic attractor
r=r1; gamma=rand; omega=10*rand; upsilon=r+r2-r1;
 
 
 
 
tspan=[0,50]; % length of time the time-varying model is solved over
% initial conditions for the time-varying model
finit = [sigma, sigma,r , r];
 
 
 
% express the time-varying model in a form matlab can interpret
f = @(t,y)[sigma*(1-alpha*h1*y(2)/(y(1)+y(2)))-eta*y(1)+beta1*y(1)*(y(3)+y(4)); sigma*(alpha*h1*y(2)/(y(1)+y(2)))-eta*y(2)+beta1*y(2)*(y(3)+y(4))+phi*y(2); (r*(1+gamma*sin(2*pi*t*omega)))*y(3)*(y(1)+y(2))-eta*y(3)+beta2*y(3)*(y(1)+y(2))-delta*y(3)+tau*y(2); (r*(1+gamma*sin(2*pi*t*omega))+upsilon)*h2*y(4)*(y(1)+y(2))-eta*y(4)+beta2*y(4)*(y(1)+y(2))-delta*y(4)-tau*y(4)];
[t,S]=ode15s(f,tspan,finit);
S1=S(:,1); % Solutions of the uninfected plant host
S2=S(:,2); % Solutions of the infected plant host
S3=S(:,3); % Solutions of the uninfected fungus host
S4=S(:,4); % Solutions of the infected fungus host
 
% check that the solutions are non-trivial, if they are, break the loop k
% and matlab will try again until the right parameter values are found
if S1(end)>0.1 && S2(end)>0.1 && S3(end)>0.1 && S4(end)>0.1
    params=[params;sigma, eta, alpha, beta1, beta2, phi, delta, tau, r, gamma, omega, upsilon, A,theta];
    break
   
end
end
 
thetar=linspace(0,1,40); %resident values of theta
thetai=linspace(0,1,40); %invader values of theta
 
 
%% pip
 
 
% time span used to calculate lyapunov exponents
T=100;
tspan=[0,T];
% create a blank figure which we use to plot the PIP
figure
 
 
% initial conditions used to numerically calculate the lyapunov exponents, 
% the method used is stated in the chapter, equations 23-26
    finit=[S1(end),S2(end),S3(end),S4(end),1,0,0,1,0,0];
 
    
    % for loops to create the PIP
    
for i=1:length(thetar)
    for j=1:length(thetai)
% %    for reference
% %    h1r=thetar(i);
% %    h2r=A*(thetar(i)-1)/(thetar(i)-A);
% %    
% %    h1i=thetai(j);
% %    h2i=A*(thetai(j)-1)/(thetai(j)-A);
   
    [t,y]=ode15s(@(t,y) [sigma*(1-alpha*thetar(i)*y(2)/(y(1)+y(2)))-eta*y(1)+beta1*y(1)*(y(3)+y(4));
        sigma*alpha*thetar(i)*y(2)/(y(1)+y(2))-eta*y(2)+phi*y(2)+beta1*y(2)*(y(3)+y(4));
        ((r*(1+gamma*sin(2*pi*t/omega))))*y(3)*(y(1)+y(2))-eta*y(3)+beta2*y(3)*(y(1)+y(2))-delta*y(3)+tau*y(4);
        ((r*(1+gamma*sin(2*pi*t/omega))+upsilon))*(A*(thetar(i)-1)/(thetar(i)-A))*y(4)*(y(1)+y(2))-eta*y(4)+beta2*y(4)*(y(1)+y(2))-delta*y(4)-tau*y(4);
        y(7)*(y(5)*y(7)*(phi-eta+beta1*(y(3)+y(4))+sigma*alpha*thetai(j)/(y(1)+y(2)))+y(6)*y(8)*(beta2*(y(1)+y(2))-eta-tau-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2))));
        y(8)*(y(5)*y(7)*(phi-eta+beta1*(y(3)+y(4))+sigma*alpha*thetai(j)/(y(1)+y(2)))+y(6)*y(8)*(beta2*(y(1)+y(2))-eta-tau-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2))));
        -y(5)*(y(5)*y(7)*(phi-eta+beta1*(y(3)+y(4))+sigma*alpha*thetai(j)/(y(1)+y(2)))+y(6)*y(8)*(beta2*(y(1)+y(2))-eta-tau-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2))));
        -y(6)*(y(5)*y(7)*(phi-eta+beta1*(y(3)+y(4))+sigma*alpha*thetai(j)/(y(1)+y(2)))+y(6)*y(8)*(beta2*(y(1)+y(2))-eta-tau-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2))));
        y(5)*y(5)*((phi-eta+beta1*(y(3)+y(4))+sigma*alpha*thetai(j)/(y(1)+y(2))))+y(6)*y(6)*(beta2*(y(1)+y(2))-eta-tau-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)));
        y(7)*y(7)*((phi-eta+beta1*(y(3)+y(4))+sigma*alpha*thetai(j)/(y(1)+y(2))))+y(8)*y(8)*(beta2*(y(1)+y(2))-eta-tau-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))],tspan,finit);
        
    
        L1=y(:,9); 
        L2=y(:,10);
        lya1=L1(end)/T; %First lyapunov exponent
        lya2=L2(end)/T; %Second lyapunov exponent
        Lmax=max(lya1,lya2); %Find the largest of the two lyapunov exponents
        if Lmax>0 % Invasion criteria
        plot(thetar(i),thetai(j),'black.','MarkerSize',7)
        hold on
        end      
    end
end
%labels for the graphs
xlim([0 1])
ylim([0 1])
xlabel('\theta_{r}')
ylabel('\theta_{i}')
    
end

%The PIPs were checked manually, and in these PIPs, transmission through the fungus or mixed transmission was never observed.
