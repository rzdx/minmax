function [ x ] = expPDF( a,b )
u1=rand;
u2=rand;
if u1>0.5
    x=a+b*log(u2);
else
    x=a-b*log(u2);
end
x=abs(x);
end