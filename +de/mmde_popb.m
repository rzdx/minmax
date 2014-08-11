function [ bestind,nv,vnew,tst,neval ] = mmde_popb( xnew,bx,XL,XU,YL,YU,NP,D,CR,mutp,evn,minmaxb,itermax,maxStagnation )
neval=0;
X=std_de.init(YL,YU,NP,D);
W=std_de.init(0,1,NP,1);
if ~isempty(bx)
    for i=1:size(bx,2)
        X(:,i)=bx(:,i);
    end
end
X(:,end)=std_de.init(YL,YL+1e-6*rand,1,D);
X(:,end-1)=std_de.init(YU-1e-6*rand,YU,1,D);
[nvX,vnewX,tstX,nfiteval]=testf.mmde_getfit(evn,X,xnew,XL,XU,YL,YU,1);
neval=neval+nfiteval;
[bestidx]=testf.mmde_evafit(nvX,vnewX,tstX,minmaxb);
bestind=X(:,bestidx);

lastBestChange=0;
it=0;
while it<itermax & it-lastBestChange<maxStagnation
    X(:,end)=std_de.init(YL,YL+1e-6*rand,1,D);
    X(:,end-1)=std_de.init(YU-1e-6*rand,YU,1,D);
    for i=1:NP
        F=fn.expPDF(0,0.5);
        % mutation
        if mutp==0 % DE/rand/1/bin
            r1=fn.dfrandi(1,NP);
            r2=fn.dfrandi(1,NP,r1);
            r3=fn.dfrandi(1,NP,r1,r2);
            tV=X(:,r1)+F*(X(:,r2)-X(:,r3));
        else % DE/DEGL-SAW/1/bin
            degl_k=round((1/3)*NP);
            nborset=(i-degl_k:i+degl_k);
            for ii=1:length(nborset)
                if nborset(ii)<1
                    nborset(ii)=nborset(ii)+NP;
                elseif nborset(ii)>NP
                    nborset(ii)=nborset(ii)-NP;
                end
            end
            
            [nnvX,nvnewX,ntstX,nfiteval]=testf.mmde_getfit(evn,X(:,nborset),xnew,XL,XU,YL,YU,1);
            neval=neval+nfiteval;
            [nbestidx]=testf.mmde_evafit(nnvX,nvnewX,ntstX,minmaxb);
            nbestind=X(:,nborset(nbestidx));
            nr1=fn.dfrandi(1,length(nborset),i);
            nr2=fn.dfrandi(1,length(nborset),i,nr1);
            r1=fn.dfrandi(1,NP,i);
            r2=fn.dfrandi(1,NP,i,r1);
            tVL=X(:,i)+F*(nbestind-X(:,i))+F*(X(:,nborset(nr1))-X(:,nborset(nr2)));
            tVG=X(:,i)+F*(bestind-X(:,i))+F*(X(:,r1)-X(:,r2));
            
            tW=W(i)+F*(W(bestidx)-W(i))+F*(W(r1)-W(r2));
            if tW>0.95
                tW=0.95;
            end
            if tW<0.05
                tW=0.05;
            end
            tV=tW*tVG+(1-tW)*tVL;
        end
        % crossover
        tU=tV*0;
        rdj=ceil(rand*D);
        for j=1:D
            if rand<=CR||rdj==j
                tU(j)=tV(j);
            else
                tU(j)=X(j,i);
            end
        end
        ynew=tU;
        [ynvX,yvnewX,ytstX,nfiteval]=testf.mmde_getfit(evn,[X(:,i),ynew],xnew,XL,XU,YL,YU,1);
        neval=neval+nfiteval;
        [ybestidx]=testf.mmde_evafit(ynvX,yvnewX,ytstX,minmaxb);
        nvX(i)=ynvX(ybestidx);
        vnewX(i)=yvnewX(ybestidx);
        tstX(i)=ytstX(ybestidx);
        if ybestidx==2
            if mutp==1
                W(i)=tW;
            end
            X(:,i)=ynew;
            [fbestidx]=testf.mmde_evafit(nvX,vnewX,tstX,minmaxb);
            fnewbestind=X(:,fbestidx);
            if ~isequal(bestind,fnewbestind)
                lastBestChange=it;
                bestidx=fbestidx;
                bestind=fnewbestind;
            end
        end
    end
    it=it+1;
end
[fbestidx]=testf.mmde_evafit(nvX,vnewX,tstX,minmaxb);
bestind=X(:,fbestidx);
nv=nvX(fbestidx);
vnew=vnewX(fbestidx);
tst=tstX(fbestidx);
end