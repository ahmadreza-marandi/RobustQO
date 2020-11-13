function s = svec( W)
%svec returning the svec vector of symmetric matrix W

n=size(W,1);

S=W-diag(diag(W));
S=sqrt(2)*S;
S=S+diag(diag(W));
s=[];
for i=1:n
    for j=i:n
        s=[s;S(i,j)];
    end
end

end

