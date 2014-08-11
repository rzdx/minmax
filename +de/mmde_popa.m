function [ fit,bestX,bestY,neval ] = mmde_popa( itermax,maxStagnation,bfraction,NP,D,mutp,evn,CR,XL,XU,YL,YU,minmaxb )
neval=0;
X=std_de.init(XL,XU,NP,D);
W=std_de.init(0,1,NP,1);
popc=zeros(D,NP);
nvX=zeros(NP,1);
vnewX=nvX;
tstX=nvX;
for i=1:NP
    [popbind,popbnv,popbvnew,popbtst,popbneval]=de.mmde_popb(X(:,i),[],XL,XU,YL,YU,NP,D,mutp,CR,evn,1-minmaxb,itermax,maxStagnation);
    popc(:,i)=popbind;
    nvX(i)=popbnv;
    vnewX(i)=popbvnew;
    tstX(i)=popbtst;
    neval=neval+popbneval;
end
[bestidx]=testf.mmde_evafit(nvX,vnewX,tstX,minmaxb);
bestind=X(:,bestidx);

lastBestChange=0;
it=0;
while it<itermax & it-lastBestChange<maxStagnation
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
            [nnvX,nvnewX,ntstX,nfiteval]=testf.mmde_getfit(evn,X(:,nborset),popc(:,nborset),XL,XU,YL,YU,0);
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
        rdj=ceil(rand*D);
        tU=tV*0;
        for j=1:D
            if rand<=CR||rdj==j
                tU(j)=tV(j);
            else
                tU(j)=X(j,i);
            end
        end
        xnew=tU;
        dst=zeros(NP,1);
        for j=1:NP
            dst(j)=norm(xnew-X(:,j));
        end
        [~,sdstidx]=sort(dst);
        nbest=round(bfraction*NP);
        [ypopbind,ypopbnv,ypopbvnew,ypopbtst,popbneval]=de.mmde_popb(xnew,popc(:,sdstidx(1:nbest)),XL,XU,YL,YU,NP,D,mutp,CR,evn,1-minmaxb,itermax,maxStagnation);
        popcnew=ypopbind;
        neval=neval+popbneval;
        
        xnvX(1)=nvX(i);
        xvnewX(1)=vnewX(i);
        xtstX(1)=tstX(i);
        xnvX(2)=ypopbnv;
        xvnewX(2)=ypopbvnew;
        xtstX(2)=ypopbtst;
        [xbestidx]=testf.mmde_evafit(xnvX,xvnewX,xtstX,minmaxb);
        if xbestidx==2
            if mutp==1
                W(i)=tW;
            end
            X(:,i)=xnew;
            popc(:,i)=popcnew;
            nvX(i)=xnvX(xbestidx);
            vnewX(i)=xvnewX(xbestidx);
            tstX(i)=xtstX(xbestidx);
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
bestX=X(:,fbestidx);
bestY=popc(:,fbestidx);
fit=testf.mmde_tst(evn,bestX,bestY);
end