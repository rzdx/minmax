function [minv,minpara,gminv,gmeanv,gstd] = dfev(Gmax,NP,D,ibv,evfnnm,evn,F,CR)
Xmax=zeros(D,1)+100;
Xmin=zeros(D,1)-100;
X=zeros(D,NP);
newX=X;
tX=zeros(D,1);
tU=tX;
for n=1:NP
    for d=1:D
        tX(d)=Xmin(d)+rand*(Xmax(d)-Xmin(d));
    end
    X(:,n)=tX;
end

evX=evalX(X,evfnnm,evn);
for G=1:Gmax
    [bestXv,bestXidx]=min(evX);
    for n=1:NP
        r1=fn.dfrandi(1,NP);
        r2=fn.dfrandi(1,NP,r1);
        r3=fn.dfrandi(1,NP,r1,r2);
        rd=[r1,r2,r3];
        switch ibv
            case 1 %rand
                tV=X(:,rd(1))+F*(X(:,rd(2))-X(:,rd(3)));
            case 2 %best
                tV=X(:,bestXidx)+F*(X(:,rd(2))-X(:,rd(3)));
            case 3 %target-to-best
                tV=X(:,n)+F*(X(:,bestXidx)-X(:,n))+F*(X(:,rd(2))-X(:,rd(3)));
            otherwise
                error('type_not_match');
        end
        
        rdj=randi(D);
        for j=1:D
            if rand<=CR||rdj==j
                tU(j)=tV(j);
            else
                tU(j)=X(j,n);
            end
        end
        
        if feval(evfnnm,tU,evn)<=evX(n)
            newX(:,n)=tU;
        else
            newX(:,n)=X(:,n);
        end
    end
    X=newX;
    evX=evalX(X,evfnnm,evn);
    [minv,minidx]=min(evX);
    minpara=X(:,minidx);
    gminv(G)=minv;
    gmeanv(G)=mean(evX);
    gstd{G}=std(X,[],2);
end
end

function [evX]=evalX(X,evfnnm,evn)
evX=zeros(size(X,2),1);
for i=1:size(X,2)
    evX(i)=feval(evfnnm,X(:,i),evn);
end
end