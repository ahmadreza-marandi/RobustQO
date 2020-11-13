clear all
load("Random_instances.mat")
%% initialization/definition
n_size=100;
ran_size=10;
rho=0:0.0025:0.1;
X_I=cell(ran_size,1,n_size);
X_O=cell(ran_size,1,n_size);
Wors_I=zeros(ran_size,size(rho,2),n_size);
Wors_N=zeros(ran_size,size(rho,2),n_size);
Wors_N_inv=zeros(ran_size,size(rho,2),n_size);
Wors_O=zeros(ran_size,size(rho,2),n_size);
Nom_I=zeros(ran_size,size(rho,2),n_size);
Nom_O=zeros(ran_size,size(rho,2),n_size);
Opval_I=zeros(ran_size,size(rho,2),n_size);
Opval_O=zeros(ran_size,size(rho,2),n_size);
Time_N=zeros(ran_size,1,n_size);
Time_I=zeros(ran_size,size(rho,2),n_size);
Time_O=zeros(ran_size,size(rho,2),n_size);
C_rho_instance=zeros(ran_size,size(rho,2),n_size); 
Opval_N_opt=zeros(ran_size,1,n_size);
Opval_N_inv=zeros(ran_size,1,n_size);
for n=4:100
    for random_construction=1:ran_size
        X_I_help=zeros(n,size(rho,2));
        X_O_help=zeros(n,size(rho,2));
        Wors_I_help=zeros(size(rho,2),1);
        Wors_N_help=zeros(size(rho,2),1);
        Wors_O_help=zeros(size(rho,2),1);
        Nom_I_help=zeros(size(rho,2),1);
        Nom_O_help=zeros(size(rho,2),1);
        Opval_I_help=zeros(size(rho,2),1);
        Opval_O_help=zeros(size(rho,2),1);
        Time_I_help=zeros(size(rho,2),1);
        Time_O_help=zeros(size(rho,2),1);
        %% constructing the parameters of the uncertainty set
        B1=Uncenrtainty{random_construction,n,1};
        B2=Uncenrtainty{random_construction,n,2};
        K=Uncenrtainty{random_construction,n,3};
        %% constructing A and b
        A=matrix{random_construction,n};
        b=right{random_construction,n};
        %% solve using norm minimization
        x=sdpvar(n,1);tau=sdpvar(1,1);
        sol=optimize([tau>=norm(A*x-b,2)],tau,sdpsettings('solver','mosek'));
        x_n=value(x);
        Opval_N_opt(random_construction)=norm(A*x_n-b);
        %% solving using linear algebra
        x_nn=A\b;
        Opval_N_inv(random_construction)=norm(A*x_nn-b);
        Time_N(random_construction,n)=sol.solvertime;
        %% making uncertainty set 
        for i=1:size(rho,2)
            %% finding the worst-case of the nominal solution
            B_N=worst_value(rho(i), B1,B2,x_n,A*x_n-b,K,n,n);
            Wors_N(random_construction,i)=norm((A+B_N)*x_n-b);
            Wors_N_inv(random_construction,i)=norm((A+(rho(i)*repmat(sign(A*x_nn-b),1,n).*repmat(sign(x_nn)',n,1)))*x_nn-b);
            %% inner approximation     
            B=rho(i)*ones(n,n);
            Gamma=norm(B,2);
            W=sdpvar(n,n);
            W1=sdpvar(n,n);
            W2=sdpvar(n,n);
            W3=sdpvar(n,n);
            y=sdpvar(n,1);
            tau=sdpvar(1,1);
            constraint=[ [W, y; y',1]>=0 ,W1+W2+W3==2*W*A'-2*b*y', (ones(n,n)-B1).*W2==0,(ones(n,n)-B2).*W3==0 ,tau>=trace((A'*A+Gamma^2*eye(n))*W)+rho(i)*norm(W1,1)+K*rho(i)*(norm(W2.*B1,inf)+norm(W3.*B2,inf))-2*b'*A*y];
            sdp=sdpsettings;
            sdp.mosek.MSK_DPAR_OPTIMIZER_MAX_TIME=0.01;
            sol=optimize(constraint,tau,sdpsettings('solver','mosek'));
            y=value(y);
            W=value(W);
            tau=value(tau);
            Time_I(random_construction,i,n)=sol.solvertime;
            X_I_help(:,i)=value(y);
            Opval_I(random_construction,i,n)=tau+norm(b,2)^2;
            %% finding the worst-case using the hueristic
            B_I=worst_value(rho(i), B1,B2,y,A*y-b,K,n,n);
            Wors_I(random_construction,i)=norm((A+B_I)*y-b);
            Nom_I(random_construction,i,n)=norm((A)*y-b);
             %% checking assumption C
            y=sdpvar(n,1);
            sol=optimize([A*y>=0,norm(y,1)<=1],y'*A'*A*y-2*rho(i)*ones(1,n)*A*y,sdpsettings('solver','mosek'));
            y=value(y);
            l=y'*A'*A*y-2*rho(i)*ones(1,n)*A*y;
            if l>=-0.05
                C_rho_instance(random_construction,i,n)=1;
            end
            %% outer approximation
            B=rho(i)*ones(n,n);
            W=sdpvar(n,n);
            W1=sdpvar(n,n);
            W2=sdpvar(n,n);
            W3=sdpvar(n,n);
            y=sdpvar(n,1);
            tau=sdpvar(1,1);
            constraint=[ [W, y; y',1]>=0 ,W1+W2+W3==2*W*A'-2*b*y', (ones(n,n)-B1).*W2==0,(ones(n,n)-B2).*W3==0 ,tau>=trace((A'*A)*W)+rho(i)*norm(W1,1)+K*rho(i)*(norm(W2.*B1,inf)+norm(W3.*B2,inf))-2*b'*A*y];
            sol=optimize(constraint,tau,sdpsettings('solver','mosek'));
            y=value(y);
            W=value(W);
            tau=value(tau);
            Time_O(random_construction,i,n)=sol.solvertime;
            X_O_help(:,i)=value(y);
            Opval_O(random_construction,i,n)=tau+norm(b)^2;
            %% finding the worst-case using the hueristic
            B_O=worst_value(rho(i), B1,B2,y,A*y-b,K,n,n);
            Wors_O(random_construction,i)=norm((A+B_O)*y-b);
            Nom_O(random_construction,i,n)=norm((A)*y-b);
         end
        X_I(random_construction,1,n)={X_I_help};
        X_O(random_construction,1,n)={X_O_help};
    end
end
save('solution_scalable')