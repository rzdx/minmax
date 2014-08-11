function [ b,sol_range ] = err_within( value,range,x )
vr=abs(value*range);
sol_range(1)=value-vr;
sol_range(2)=value+vr;
if x>=sol_range(1) & x<=sol_range(2)
    b=1;
else
    b=0;
end
end

% b: 1=true 0=false