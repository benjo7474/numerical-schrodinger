
function [total] = nonUniformTrap(x,f)
total=0;
for i=1:size(x)-1
    total=total+abs(x(i+1)-x(i))*(f(i)+f(i+1))*0.5;
end;