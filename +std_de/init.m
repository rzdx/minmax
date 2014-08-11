function [X]=init_X(l,u,NP,D)
X=l+rand(D,NP)*(u-l);
end