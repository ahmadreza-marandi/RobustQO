
function Delta=worst_value(rho, B1,B2,y,x_value,K,m,n)
%x_value A*y-b
% m rows, n column or variables

Delta=zeros(m,n);
for i=1:m
    for j=1:n
        Delta(i,j)=rho*sign(y(j))*sign(x_value(i));
    end
end

[RowB1,ColB1]=find(B1);
[RowB2,ColB2]=find(B2);
numberOFones1=size(RowB1,1);
numberOFones2=size(RowB2,1);

BothCol=intersect(ColB1,ColB2);
[~,SortBoth]=sort(abs(y(BothCol)),'descend');


count1=0;%counter for the 
count2=0;%
redundant1=0;
redundant2=0;
if numberOFones1<=K % if the number of ones is fewer than K, this means the constraint is redundant
    redundant1=1;
else% else we need to assing rho to the highest values only, so we start from zero matrix
    for i=1:numberOFones1
        Delta(RowB1(i),ColB1(i))=0;
    end
end
if numberOFones2<=K 
   redundant2=1;
else
     for i=1:numberOFones2
        Delta(RowB2(i),ColB2(i))=0;
    end
end

j=1;
while count1<K &&redundant1==0
    if j>size(SortBoth,1)  % this means K is higher than the number of ones in the columns of B1 among the BothCol
           % so we need to go further than BothCol. Since the columns of B2 that have ones outside BothCol are not 
           % among the columns of B1 that have ones  (otherwise it will be
           % in BothCol), we are sure that the allocations here won't
           % conflict with B2 ones. So, we go freely by increasing count
        NotInIntersect1=setdiff(ColB1,BothCol);
        [~,SortB1]=sort(abs(y(NotInIntersect1)),'descend');
        j_secondlayer=1;
        while j_secondlayer<size(SortB1,1)
            i=1;
            while i<=m && count1<K
                 if B1(i,SortB1(j_secondlayer))==1
                    Delta(i,SortB1(j_secondlayer))=rho*sign(y(SortB1(j_secondlayer)))*sign(x_value(i));
                    count1=count1+1;
                 end
                i=i+1;
            end
             j_secondlayer=j_secondlayer+1;
        end
        
    else
        i=1;
        while i<=m && count1<K% start from the column with highest value of y and go through the 
                  %rows one by one and increasing count1
            if B1(i,SortBoth(j))==1 && B2(i,SortBoth(j))==1
                if count2<K
                    Delta(i,SortBoth(j))=rho*sign(y(SortBoth(j)))*sign(x_value(i));
                    count1=count1+1;
                    count2=count2+1;
                end
            elseif B1(i,SortBoth(j))==0 && B2(i,SortBoth(j))==0
                Delta(i,SortBoth(j))=rho*sign(y(SortBoth(j)))*sign(x_value(i));
            elseif B1(i,SortBoth(j))==1
                Delta(i,SortBoth(j))=rho*sign(y(SortBoth(j)))*sign(x_value(i));
                count1=count1+1;
            else
                if count2<K
                     Delta(i,SortBoth(j))=rho*sign(y(SortBoth(j)))*sign(x_value(i));
                    count2=count2+1;
                end
            end
            i=i+1;
        end
        j=j+1;
    end
end
%same story for B2
if j~=1
    j=j-1;  % this is because the last iteration in j may not be complete,
end   % so we need to check it for B2

while count2<K&& redundant2==0
    if j>size(SortBoth,1) 
        NotInIntersect2=setdiff(ColB2,BothCol);
        [~,SortB2]=sort(abs(y(NotInIntersect2)),'descend');
        j_secondlayer=1;
        while j_secondlayer<size(SortB2,1)
            i=1;
            while i<=m && count2<K
                 if B2(i,SortB2(j_secondlayer))==1
                    Delta(i,SortB2(j_secondlayer))=rho*sign(y(SortB2(j_secondlayer)))*sign(x_value(i));
                    count2=count2+1;
                 end
                i=i+1;
            end
             j_secondlayer=j_secondlayer+1;
        end
        
    else
        i=1;
        while i<=m && count2<K% start from the column with highest value of y and go through the 
                              %rows one by one and increasing count1
            if B1(i,SortBoth(j))==0 && B2(i,SortBoth(j))==0
                Delta(i,SortBoth(j))=rho*sign(y(SortBoth(j)))*sign(x_value(i));
            elseif B2(i,SortBoth(j))==1&& B1(i,SortBoth(j))==0
                Delta(i,SortBoth(j))=rho*sign(y(SortBoth(j)))*sign(x_value(i));
                count2=count2+1;
            end
            i=i+1;
        end
        j=j+1;
    end
end

end