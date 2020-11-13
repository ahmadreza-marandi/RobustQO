% To run the code, you need to have the Data_30.mat, RobustPorfolio.m, and
% ChancePortfolio.m
%The code uses Yalmip (https://yalmip.github.io/) to pass the optimization
%problem to the solver. 
%As the solver, we use Mosek (https://www.mosek.com/)
%
load('Data_30.mat','z');
n=size(z,2);
sizem=round((n+2)*(n+1)/2)-1;
xbar=mean(z); %1*n matrix
S=cov(z);   %n*n matrix
Schol=chol(S);
 
lambda=3;
 %% Robust version
alpha=[0:0.02:0.98];
z_a30=[];
timeRO=[];
Rphi95=0;
Rrelax95=0;
timeCO=[];
 for ii=1:size(alpha,2)
    [z_a,Rphi,RW,s,solvert_time]=RobustPortfolio(z, alpha(ii), lambda);
    [Rrelaxphi,RelaxW,s_moment,R,solvert_timeCO]=ChancePortfolio(z, alpha(ii), lambda);
    if alpha(ii)==0.94
        Rphi95=Rphi;
        Rrelax95=Rrelaxphi;
    end
    z_a30=[z_a30,z_a];
    timeRO=[timeRO,solvert_time ];
    timeCO=[timeCO,solvert_timeCO ];
    
 end
 %% Chance version
 Rrelaxphi=Rrelax95;

Rrelaxmean=12*Rrelaxphi'*xbar';
Rrelaxrisk=sqrt(12*Rrelaxphi'*S*Rrelaxphi);
%% Nominal problem
m=zeros(sizem,n);
m(1:n,1:n)=eye(n);
count=n;
for i=1:n
    m(count+1:count+n-i+1,i:n)=[ones(n-i+1,1),zeros(n-i+1,n-i)]+eye(n-i+1);
    count=count+n-i+1;
end
Ahelp=ones(1,n);
for i=n+1:sizem
    if size(find(m(i,:)),2)==2
        Ahelp=[Ahelp,sqrt(2)];
    else
        Ahelp=[Ahelp,1];
    end
end
 A=diag(Ahelp);
 
phi=sdpvar(n,1);
sol_nom=optimize([phi>=0,sum(phi)==1],-phi'*xbar'+lambda*(phi'*S*phi),sdpsettings('solver','mosek'));
phi=value(phi);
Nominalphi=phi;
Nominalval=phi'*xbar'-lambda*phi'*S*phi;

%% Finding the worst-case scenario
r=sqrt(chi2inv(0.95,s_moment)/size(z,1)); % chi-square for n value and alpha% confident
mu=sdpvar(n,1);
sigma=sdpvar(n,n);
vec=sdpvar(n,1);
nu=sdpvar(sizem,1);
constraint=[sigma>=0,[mu;svec(sigma)]==A*R*nu+[xbar';svec(S)],norm(nu)<=r];
sol=optimize(constraint,Rphi95'*mu-lambda*(Rphi95'*sigma*Rphi95),sdpsettings('solver','mosek'));
mu95=value(mu); 
sigma95=value(sigma);

display('risk of robust solution for the scenario ')
sqrt(12*(Rphi95'*sigma95*Rphi95))
display('mean return of robust solution for the scenario ')
12*Rphi95'*mu95
display('Risk of nominal solution for the scenario ')
sqrt(12*(Nominalphi'*sigma95*Nominalphi))
display('mean of nominal solution for the scenario ')
12*Nominalphi'*mu95
display('risk of robust solution for niminal case')
sqrt(12*(Rphi95'*S*Rphi95))
display('mean return of robust solution for nominal case')
12*Rphi95'*xbar'
display('Risk of nominal solution for nominal case ')
sqrt(12*(Nominalphi'*S*Nominalphi))
display('mean of nominal solution for nominal case ')
12*Nominalphi'*xbar'
vec=value(vec);
nu=value(nu);
display('relaxation mean and risk')
Rrelaxmean
Rrelaxrisk