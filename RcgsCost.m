function f = RcgsCost(b,x,y)

y_predicted = Rcgs(b,x);    % get prediction
f = sum((y-y_predicted).^2);    % return the sum-squared-error
