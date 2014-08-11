function [ f ] = mmde_tst( id,x,y )
switch id
    case 1
        f=(cos(y)+cos(2*y+x))^2;
    case 2
        f=x^2+y^2+2*x*y-20*x-20*y+100;
    case 3
        f=sin(x)^2-x*cos(y)+2*sin(x)-cos(y)^2+y-1;
end
end