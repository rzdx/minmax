function [ newX ] = xovr( X,muX,CR )
% bin
newX=X*0;
NP=size(X,2);
D=size(X,1);
rdj=randi(D);
for i=1:NP;
    for j=1:D
        if rand<=CR||rdj==j
            newX(j,i)=muX(j,i);
        else
            newX(j,i)=X(j,i);
        end
    end
end
end