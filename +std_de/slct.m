function [ newX ] = slct( X,oldX,evfnnm,evn,evX )
newX=X*0;
NP=size(X,2);
for i=1:NP
    if feval(evfnnm,X(i),evn)<=evX(i)
        newX(:,i)=X(:,i);
    else
        newX(:,i)=oldX(:,i);
    end
end
end