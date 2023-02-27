function [R] = Rcgs(param,theta)
%constants
k=-4;
f=2;
p=2*pi;
A = exp(-((f.*theta+k*p - param(3)).^2)/(2*param(4)));
for k=-4:4
A = A + exp(-((f.*theta+k*p - param(3)).^2)/(2*param(4)));
end
R= param(1) + param(2) * A ;

end