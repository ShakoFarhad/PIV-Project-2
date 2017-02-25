% Newtons Method
% 
% f    is the function to find root for
% df   is the derivative of f
% tol  is the accepted tolerance of the error
% kmax is the max iterations
% x    is the root
% ex   is the error of the root
% k    is the number of iterations
% 
% 
% Example of running:
% [x, ex, k] = NewtonsMethod('9.81*x*tanh(x*0.72) - (2*pi*1.425)^2', '9.81*tanh(0.72*x)+9.81*0.72*x*(sech(0.72*x))^2', 4, 0.00001, 10000);
% 


function [x, ex, k] = NewtonsMethod(f, df, x0, tol, kmax) %// Change

    if nargin == 3
        tol = 1e-4;
        kmax = 1e1;
    elseif nargin == 4
        kmax = 1e1;
    elseif nargin ~= 5
        error('newton: invalid input parameters');
    end
    
    f = inline(f);
    df = inline(df);
    x(1) = x0 - (f(x0)/df(x0));
    ex(1) = abs(x(1)-x0);
    k = 2;
    while (ex(k-1) >= tol) && (k <= kmax)
        x(k) = x(k-1) - (f(x(k-1))/df(x(k-1)));
        ex(k) = abs(x(k)-x(k-1));
        k = k+1;
    end
end