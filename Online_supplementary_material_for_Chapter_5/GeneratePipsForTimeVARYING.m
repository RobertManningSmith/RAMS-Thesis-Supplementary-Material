clear
clc

thetar=linspace(0,1,45);
thetai=linspace(0,1,45);
T1=50; 
T2=50;
tspan1=[0 T1];
tspan2=[0 T2];


for o=1:5
    figure
    sigma=10*rand; eta=10*rand; alpha=rand; beta1=-10*(rand-0.5); beta2=-10*(rand-0.5); phi=10*rand;  delta=10*rand;
    A=1+rand;
    r=10*rand; gamma=rand; omega=10*rand; upsilon=rand;
    finit=[sigma,sigma,r,r];
for i=1:length(thetar)
    for j=1:length(thetai)
    
        [t,y]=ode45(@(t,y) [sigma*(1-alpha*thetar(i)*y(2)/(y(1)+y(2)))-eta*y(1)+beta1*y(1)*(y(3)+y(4));
                            sigma*alpha*thetar(i)*y(2)/(y(1)+y(2))-eta*y(2)+phi*y(2)+beta1*y(2)*(y(3)+y(4));
                            r*(1+gamma*sin(2*pi*t/omega))*y(3)*(y(1)+y(2))-eta*y(3)+beta2*y(3)*(y(1)+y(2))-delta*y(3)+(r*(1+gamma*sin(2*pi*t/omega))+upsilon)*y(4)*(y(1)+y(2))*(1-(A*(thetar(i)-1)/(thetar(i)-A)));
                            (r*(1+gamma*sin(2*pi*t/omega))+upsilon)*y(4)*(y(1)+y(2))*(A*(thetar(i)-1)/(thetar(i)-A))-eta*y(4)+beta2*y(4)*(y(1)+y(2))-delta*y(4)],tspan1,finit);
        Y1=y(:,1);
        Y2=y(:,2);
        Y3=y(:,3);
        Y4=y(:,4);
        Ginit=[Y1(end),Y2(end),Y3(end),Y4(end),1,0,0,1,0,0];
        
        
        [t,z]=ode45(@(t,y) [sigma*(1-alpha*thetar(i)*y(2)/(y(1)+y(2)))-eta*y(1)+beta1*y(1)*(y(3)+y(4)); 
                            sigma*alpha*thetar(i)*y(2)/(y(1)+y(2))-eta*y(2)+phi*y(2)+beta1*y(2)*(y(3)+y(4));
                            r*(1+gamma*sin(2*pi*t/omega))*y(3)*(y(1)+y(2))-eta*y(3)+beta2*y(3)*(y(1)+y(2))-delta*y(3)+(r*(1+gamma*sin(2*pi*t/omega))+upsilon)*y(4)*(y(1)+y(2))*(1-(A*(thetar(i)-1)/(thetar(i)-A)));
                            (r*(1+gamma*sin(2*pi*t/omega))+upsilon)*y(4)*(y(1)+y(2))*(A*(thetar(i)-1)/(thetar(i)-A))-eta*y(4)+beta2*y(4)*(y(1)+y(2))-delta*y(4);
                            y(7)*((phi-eta+beta1*(y(3)+y(4))+(sigma*alpha*thetai(j))/(y(1)+y(2)))*y(5)*y(7)+(beta2*(y(1)+y(2))-eta-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))*y(6)*y(8));
                            y(8)*((phi-eta+beta1*(y(3)+y(4))+(sigma*alpha*thetai(j))/(y(1)+y(2)))*y(5)*y(7)+(beta2*(y(1)+y(2))-eta-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))*y(6)*y(8));
                            y(5)*((phi-eta+beta1*(y(3)+y(4))+(sigma*alpha*thetai(j))/(y(1)+y(2)))*y(5)*y(7)+(beta2*(y(1)+y(2))-eta-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))*y(6)*y(8));
                            y(6)*((phi-eta+beta1*(y(3)+y(4))+(sigma*alpha*thetai(j))/(y(1)+y(2)))*y(5)*y(7)+(beta2*(y(1)+y(2))-eta-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))*y(6)*y(8));
                            (phi-eta+beta1*(y(3)+y(4))+(sigma*alpha*thetai(j))/(y(1)+y(2)))*y(5)*y(5)+(beta2*(y(1)+y(2))-eta-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))*y(6)*y(6);
                            (phi-eta+beta1*(y(3)+y(4))+(sigma*alpha*thetai(j))/(y(1)+y(2)))*y(7)*y(7)+(beta2*(y(1)+y(2))-eta-delta+(A*(thetai(j)-1)/(thetai(j)-A))*(y(1)+y(2)))*y(8)*y(8)],tspan2,Ginit);
        
        Z1=z(:,9);
        Z2=z(:,10);
        
        
        
        lya1=Z1(end)/T2;
        lya2=Z2(end)/T2;
        Lmax=max(lya1,lya2);
        LL(i,j)=Lmax;
        if Lmax>0
            plot(thetar(i),thetai(j),'black.','MarkerSize',8)
            hold on
        end
        
        
    end
end
LL1=sign(LL);
xlim([0 1])
ylim([0 1])
xlabel('\theta_{r}')
ylabel('\theta_{i}')
end