function [ newX ] = muta( X,F )
% DE/rand/1
NP=length(X);
newX=X*NaN;
for i=1:NP
    r1=fn.dfrandi(1,NP);
    r2=fn.dfrandi(1,NP,r1);
    r3=fn.dfrandi(1,NP,r1,r2);
    rd=[r1,r2,r3];
    newX(:,i)=X(:,rd(1))+F*(X(:,rd(2))-X(:,rd(3)));
end
end