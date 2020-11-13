function [z_a,Rphi,W,s_moment,solver_time]=RobustPortfolio(z, alpha, lambda)
% -z is the data set. Each row corresponds to a time slot and each column an
%  asset
% -alpha is the confidence level for which we want to consider the
%  uncertainty set
% -lambda is the risk-aversion coefficient
%%%
% -Rphi and W is the solutions of the robust counterpart
% -z_a=Prob(mu'*Rphi-lambda*Rphi'*S*Rphi>=z*), where z* is the robust
% optimal value
% -s_moment is the size of the moment matrix V in the paper
%
% Written by Ahmadreza Marandi and Bertrand Melenberg based on the paper:
% Extending the scope of robust quadratic optimization
% All rights are preserved.

n=size(z,2);
sizem=round((n+2)*(n+1)/2)-1;
m=zeros(sizem,n);
m(1:n,1:n)=eye(n);
count=n;
%% 
for i=1:n
    m(count+1:count+n-i+1,i:n)=[ones(n-i+1,1),zeros(n-i+1,n-i)]+eye(n-i+1);
    count=count+n-i+1;
end

xbar=mean(z); %1*n matrix
S=cov(z);   %n*n matrix
Schol=chol(S);
matxbar=repmat(xbar,[size(z,1),1]);
v=z-matxbar;
sizev=size(z,1);
H=[];
for i=1:sizem
    y=[];
    for j=1:sizev
         y=[y;prod(v(j,:).^m(i,:))];
    end
   
    H=[H,y];
    y=[];
end
v=cov(H);

jacobian=eye(sizem);

helpm=m(n+1:end,:);
mmatrix=diag(xbar);
helpm=-helpm*mmatrix;
[rh,ch]=size(helpm);
for i=1:rh
    dumm=0;
    for j=1:ch
        if m(n+i,j)==1 && dumm==0
            aa=helpm(i,j);
            ja=j;
            dumm=1;
        elseif m(n+i,j)==1 && dumm==1
            bb=helpm(i,j);
            jb=j;
        end
    end
    if dumm==1
        helpm(i,ja)=bb;
        helpm(i,jb)=aa;
    end
end

jacobian(n+1:end,1:n)=helpm(:,:);
v=jacobian*v*jacobian';

s_moment=size(v,2);
R=chol(v);
Ahelp=ones(1,n);
for i=n+1:sizem
    if size(find(m(i,:)),2)==2
        Ahelp=[Ahelp,sqrt(2)];
    else
        Ahelp=[Ahelp,1];
    end
end
 A=diag(Ahelp);
  
 %% Robust version
robustInf=0;

r=sqrt(chi2inv(alpha,size(v,2))/720); % chi-square for n value and alpha% confident
phi=sdpvar(n,1);
W=sdpvar(n,n);
constraint=[sum(phi)==1,[1, phi';phi, W] >=0, phi>=0];
sol=optimize(constraint,-phi'*xbar'+lambda*(trace(S*W))+(norm((A*R)'*[-phi;lambda*svec(W)])*r),sdpsettings('solver','mosek'));
solver_time=sol.solvertime;
phi=value(phi);
W=value(W);
Rphi=phi;
Rval=phi'*xbar'-lambda*(trace(S*W))-(norm((A*R)'*[-phi;lambda*svec(W)])*r);

b=[phi;-lambda*svec(phi*phi')];
 nor=(b'*A*v*A*b)^0.5;   
z_a=normcdf((b'*[xbar';svec(S)]-Rval)*sqrt(size(z,1))/nor);

end