clear all
%% input
load('energy.mat');
%date time year-month-day hour:minute:second 
% Appliances, energy use in Wh 
% lights, energy use of light fixtures in the house in Wh 
% T1, Temperature in kitchen area, in Celsius 
% RH_1, Humidity in kitchen area, in % 
% T2, Temperature in living room area, in Celsius 
% RH_2, Humidity in living room area, in % 
% T3, Temperature in laundry room area 
% RH_3, Humidity in laundry room area, in % 
% T4, Temperature in office room, in Celsius 
% RH_4, Humidity in office room, in % 
% T5, Temperature in bathroom, in Celsius 
% RH_5, Humidity in bathroom, in % 
% T6, Temperature outside the building (north side), in Celsius 
% RH_6, Humidity outside the building (north side), in % 
% T7, Temperature in ironing room , in Celsius 
% RH_7, Humidity in ironing room, in % 
% T8, Temperature in teenager room 2, in Celsius 
% RH_8, Humidity in teenager room 2, in % 
% T9, Temperature in parents room, in Celsius 
% RH_9, Humidity in parents room, in % 
% To, Temperature outside (from Chievres weather station), in Celsius 
% Pressure (from Chievres weather station), in mm Hg 
% RH_out, Humidity outside (from Chievres weather station), in % 
% Wind speed (from Chievres weather station), in m/s 
% Visibility (from Chievres weather station), in km 
% Tdewpoint (from Chievres weather station), °C 
% rv1, Random variable 1, nondimensional 
% rv2, Random variable 2, nondimensional

%% seperating learning and test dataset
m=floor(size(data,1)*0.75);
n=size(data,2)-1;
N=size(data,1);
test_N=N-m;
Y=data(1:m,1);
Y_test=data(m+1:end,1);
X=data(1:m,2:end);
X_test=data(m+1:end,2:end);


Xscale=max(abs(X));
Y_scale=max(abs(Y));


%% Nominal problem
alpha=sdpvar(n,1);
beta=sdpvar(1,1);
sol=optimize([],norm([ones(m,1),X,Y]*[beta;alpha;-1]),sdpsettings('solver','mosek'));
x_n=[value(alpha); -1;value(beta)];

%% inner approximation
rho=0.01;
A=[X,Y, ones(m,1)];    
coef=repmat([Xscale(1:n) Y_scale 0],m,1);
B=[ones(m,n+1) zeros(m,1)];
%finding omega
x=sdpvar(n+1,1);
opt_omega=sdpvar(1,1);
optimize([norm(x,2)<=1,opt_omega>=-[Xscale(1:n) Y_scale]*x],opt_omega,sdpsettings('solver','mosek'));
Omega=rho*sqrt(m)*(-value(opt_omega));
%making the matrix regarding the Temprature features only
B1=zeros(m,n+2);
J=2*[1:10];
B1(:,J)=ones(m,10);
rho2=0.001*sum(sum(A.*B1));
W=sdpvar(n+2,n+2);
W_bar=sdpvar(m,n+2);
y=sdpvar(n+2,1);
tau=sdpvar(1,1);
W3=2*A*W-W_bar;
constraint=[ W3(:,n+2)==zeros(m,1),(ones(m,n+2)-B1).*W_bar==0, [W, y; y',1]>=0 ,y(n+1)==-1,tau>=trace((A'*A+Omega^2*eye(n+2))*W)+rho*norm(coef.*W3,1)+rho2*norm(W_bar.*B1,inf)];
sol_I=optimize(constraint,tau,sdpsettings('solver','mosek'));
W=value(W);
y=value(y);
X_I=value(y);
%% checking the solution by LASSO
n_l=20;
tic
[Lasso, Lasso_Info]=lasso(X,Y,'CV',5,'NumLambda',n_l);
lasso_time=toc;
[a,b]=min(Lasso_Info.MSE);
lasso_sol=[Lasso(:,b);-1;Lasso_Info.Intercept(b)];

%% evaluating the solutions
rand_size=100;
Eval_Scenario_robust=[];
Eval_Scenario_nom=[];
Eval_Scenario_lasso=[];
i=1;
number_iter=0;
while i<rand_size
    Scenar=[rand(m,n+1)*2*rho.*coef(:,1:n+1)-rho*coef(:,1:n+1) zeros(m,1)];
    if norm(Scenar.*B1,1)<rho2 %inside the uncertainty set
        Scenar=Scenar+A;
        Eval_Scenario_robust=[Eval_Scenario_robust;abs(Scenar*X_I)];
        Eval_Scenario_lasso=[Eval_Scenario_lasso;abs(Scenar*lasso_sol)];
        Eval_Scenario_nom=[Eval_Scenario_nom;abs(Scenar*x_n)];
        i=i+1
    end
    number_iter=number_iter+1;
end
%% worst-case
K_I=floor(0.1*sum(sum(abs(A)))/max(max(abs(A))))-1;
worst_help=worst_value(rho,B1,zeros(m,n+2),X_I,A*X_I,K_I,m,n+2);
while norm(coef.*worst_help.*B1,1)<=rho2 &&K_I<m*(n+2)
    K_I=floor(2*K_I+1)
    Worst_I=worst_help;
    worst_help=worst_value(rho,B1,zeros(m,n+2),X_I,A*X_I,K_I,m,n+2);
end
K_I=floor(K_I/2);
 worst_help=worst_value(rho,B1,zeros(m,n+2),X_I,A*X_I,K_I,m,n+2);
while norm(coef.*worst_help.*B1,1)<=rho2 &&K_I<m*(n+2)
    K_I=K_I+10
    Worst_I=worst_help;
    worst_help=worst_value(rho,B1,zeros(m,n+2),X_I,A*X_I,K_I,m,n+2);
end
%
K_n=floor(0.1*sum(sum(abs(A)))/max(max(abs(A))))-1;
 worst_help=worst_value(rho,B1,zeros(m,n+2),x_n,A*x_n,K_n,m,n+2);
while norm(coef.*worst_help.*B1,1)<=rho2 &&K_n<m*(n+2)
    K_n=2*K_n+1
    Worst_n=worst_help;
    worst_help=worst_value(rho,B1,zeros(m,n+2),x_n,A*x_n,K_n,m,n+2);
end
K_n=floor(K_n/2);
 worst_help=worst_value(rho,B1,zeros(m,n+2),x_n,A*x_n,K_n,m,n+2);
while norm(coef.*worst_help.*B1,1)<=rho2 &&K_n<m*(n+2)
    K_n=K_n+10
    Worst_n=worst_help;
    worst_help=worst_value(rho,B1,zeros(m,n+2),x_n,A*x_n,K_n,m,n+2);
end
%
K_l=floor(0.1*sum(sum(abs(A)))/max(max(abs(A))))-1;
 worst_help=worst_value(rho,B1,zeros(m,n+2),lasso_sol,A*lasso_sol,K_l,m,n+2);
while norm(coef.*worst_help.*B1,1)<=rho2 &&K_l<m*(n+2)
    K_l=2*K_l+1
    Worst_l=worst_help;
    worst_help=worst_value(rho,B1,zeros(m,n+2),lasso_sol,A*lasso_sol,K_l,m,n+2);
end
K_l=floor(K_l/2);
 worst_help=worst_value(rho,B1,zeros(m,n+2),lasso_sol,A*lasso_sol,K_l,m,n+2);
while norm(coef.*worst_help.*B1,1)<=rho2 &&K_l<m*(n+2)
    K_l=K_l+10
    Worst_l=worst_help;
    worst_help=worst_value(rho,B1,zeros(m,n+2),lasso_sol,A*lasso_sol,K_l,m,n+2);
end
