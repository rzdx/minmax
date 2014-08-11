function [minv,minpara] = NPbestfev(MAXFEs,D,ibv,fhd,fi,F,CR,pr,pt)
NPmax=5000;
NPmin=D;
NP=5*D;
Xmax=zeros(D,1)+100;
Xmin=zeros(D,1)-100;
X=zeros(D,NPmax);
newX=X;
tX=zeros(D,1);
tU=tX;
for n=1:NPmax
    for d=1:D
        tX(d)=Xmin(d)+rand*(Xmax(d)-Xmin(d));
    end
    X(:,n)=tX;
end
G=0;
FEs=0;
evX=evalX(X,NP,fhd,fi);
while FEs<=MAXFEs
    G=G+1;
    [bestXv,bestXidx]=min(evX);
    for n=1:NP
        r1=fn.dfrandi(1,NP);
        r2=fn.dfrandi(1,NP,r1);
        r3=fn.dfrandi(1,NP,r1,r2);
        rd=[r1,r2,r3];
        switch ibv
            case 1 %rand
                tV=X(:,rd(1))+(X(:,rd(2))-X(:,rd(3)));
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
        tUv=feval(fhd,tU,fi);
        if tUv<=evX(n)
            newX(:,n)=tU;
            evX(n)=tUv;
        else
            newX(:,n)=X(:,n);
        end
    end
    oldX=X;
    X=newX;
    
    if mod(G,10)==0
        s=sg(X,NP,oldX,NP);
        m=(log(pt)-log(stdd(X)))/(s+realmin);
        NPe=(MAXFEs-FEs)/m;
        NP=floor((1-pr)*NP+pr*NPe);
        if NP>NPmax
            NP=NPmax;
        end
        if NP<NPmin
            NP=NPmin;
        end
    end
    FEs=FEs+NP;
end
[minv,minidx]=min(evX);
minpara=X(:,minidx);
end

function [O]=sg(X,NP1,oldX,NP2)
O=mean(std(X(:,1:NP1),[],2)./(std(oldX(:,1:NP2),[],2)+realmin));
end

function [stdx]=stdd(X)
stdx=mean(std(X,[],2));
end

function [evX]=evalX(X,NP,evfnnm,evn)
evX=zeros(NP,1);
for i=1:NP
    evX(i)=feval(evfnnm,X(:,i),evn);
end
end