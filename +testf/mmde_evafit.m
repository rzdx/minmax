function [ bestidx ] = mmde_evafit( nvX,vnewX,tstX,minmaxb )
minv=min(nvX);
nvidx=find(nvX==minv);
if length(nvidx)==1
    bestidx=nvidx;
    return;
end
minv=min(vnewX(nvidx));
vnewXidx=find(vnewX==minv);
if length(vnewXidx)==1
    bestidx=vnewXidx;
    return;
end
if minmaxb==0
    minv=min(tstX(vnewXidx));
else
    minv=max(tstX(vnewXidx));
end
tstXidx=find(tstX==minv);
if length(tstXidx)==1
    bestidx=tstXidx;
    return;
else
    bestidx=tstXidx(1);
end
end