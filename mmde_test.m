mpool
clear
% clc
fldn='result';
sfldn={'meanstd'};
fio.nfolds(fldn,sfldn);

runmax=10;
itermax=200;

funmin=1;
funmax=3;
funlen=funmax-funmin+1;

maxStagnation=20;
NP=30;
CR=0.9;
mutps=[0,1];
bfractions=[0,1/5];
mutns={'rand','DEGL'};
bfns={'without copy','with copy'};
Ds=[1 1 1];
XLs=[-3.14 0 -5];
XUs=[3.14 6 5];
YLs=[-3.14 2 -5];
YUs=[3.14 8 5];
minmaxbs=[0 0 0];
optsol={[-0.437082 -0.4370820;-2.553833 -3.14],...
    [4.143 4.143;4.807 6.907],...
    [-0.88734;5]};
optsolv=[0.0085865 1.10255 3.22169];
Tct=0;
tic
for mi=1:length(mutps)
    disp(mutns{mi});
    for bi=1:length(bfractions)
        disp(bfns{bi});
        
        Tct=Tct+1;
        fitrst=zeros(funlen,runmax);
        nevalrst=fitrst;
        dstrst=fitrst;
        Srst=fitrst;
        bestXrst=fitrst;
        bestYrst=bestXrst;
        parfor j=1:runmax
            disp(['run: ',num2str(j)]);
            for evn=funmin:funmax
                disp(['P',num2str(evn)]);
                
                D=Ds(evn);
                XL=XLs(evn);
                XU=XUs(evn);
                YL=YLs(evn);
                YU=YUs(evn);
                minmaxb=minmaxbs(evn);
                bfraction=bfractions(bi);
                mutp=mutps(mi);
                [fit,bestX,bestY,neval]=de.mmde_popa(itermax,maxStagnation,bfraction,NP,D,mutp,evn,CR,XL,XU,YL,YU,minmaxb);
                %----------------------------result
                fitrst(evn,j)=fit;
                nevalrst(evn,j)=neval;
                optdst=zeros(size(optsol{evn},2),1);
                for di=1:size(optsol{evn},2)
                    optdst(di)=norm(optsol{evn}(:,di)-[bestX;bestY]);
                end
                [dstrst(evn,j),mindstidx]=min(optdst);
                sX=fn.err_within(optsol{evn}(1,mindstidx),0.02,bestX);
                sY=fn.err_within(optsol{evn}(2,mindstidx),0.02,bestY);
                if sX==1 & sY==1
                    S=1;
                else
                    S=0;
                end
                Srst(evn,j)=S;
                bestXrst(evn,j)=bestX;
                bestYrst(evn,j)=bestY;
            end
        end
        rownm=cell(funlen,1);
        colnm={'mean','std','success_rate'};
        datav=zeros(length(rownm),length(colnm));
        for dr=1:funlen
            datav(dr,1)=mean(dstrst(dr,:));
            datav(dr,2)=std(dstrst(dr,:),1);
            datav(dr,3)=sum(Srst(dr,:))/runmax;
            rownm{dr}=[mutns{mi},'-',bfns{bi},'_P',num2str(dr)];
        end
        T={datav,rownm,colnm};
        shT{Tct}=T;
        save([fio.addslash(1,fldn,sfldn{1}),mutns{mi},'-',bfns{bi},'-','T.mat'], 'T');
        save([fio.addslash(1,fldn,sfldn{1}),mutns{mi},'-',bfns{bi},'-','fitrst.mat'], 'fitrst');
        save([fio.addslash(1,fldn,sfldn{1}),mutns{mi},'-',bfns{bi},'-','dstrst.mat'], 'dstrst');
        save([fio.addslash(1,fldn,sfldn{1}),mutns{mi},'-',bfns{bi},'-','Srst.mat'], 'Srst');
        save([fio.addslash(1,fldn,sfldn{1}),mutns{mi},'-',bfns{bi},'-','bestXrst.mat'], 'bestXrst');
        save([fio.addslash(1,fldn,sfldn{1}),mutns{mi},'-',bfns{bi},'-','bestYrst.mat'], 'bestYrst');
    end
end
toc

XT={[fio.addslash(1,fldn,sfldn{1}),'mmde_rst.xls'],0,shT};
tio.xlswt(XT);
disp('---OVER---');