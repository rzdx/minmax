function [ nvX,vnewX,tstX,nfiteval ] = mmde_getfit( id,X,Y,XL,XU,YL,YU,popab )
% popab: popa=0,popb=1
nfiteval=0;
nvX=zeros(size(X,2),1);
vnewX=nvX;
tstX=nvX;
for n=1:size(X,2)
    if popab==0 % popa
        x=X(:,n);
        y=Y(:,n);
    else % popb
        x=Y;
        y=X(:,n);
    end
    
    switch id
        case 1
            g(1)=y-x*(x+6.28);
            g(2)=y-x*(x-6.28);
        case 2
            g(1)=-(x-5)^2-(y-3)^2+4;
            g(2)=(x-5)^2+(y-3)^2-16;
        case 3
            g(1)=-x^2-y^2+25;
    end
    nv=0;
    vnew=0;
    for i=1:length(g)
        if g(i)>0
            nv=nv+1;
            vnew=vnew+g(i);
        end
    end
    nvX(n)=nv;
    vnewX(n)=vnew;
    tstX(n)=testf.mmde_tst(id,x,y);
    nfiteval=nfiteval+1;
    
    for d=1:length(x)
        if x(d)<XL | x(d)>XU | y(d)<YL | y(d)>YU
            nvX(n)=inf;
        end
    end
end
end