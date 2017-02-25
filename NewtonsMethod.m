%Newtons Method

function [x, ex, k] = NewtonsMethod(f, df, x0, tol, nmax) %// Change

    if nargin == 3
        tol = 1e-4;
        nmax = 1e1;
    elseif nargin == 4
        nmax = 1e1;
    elseif nargin ~= 5
        error('newton: invalid input parameters');
    end
    
    f = inline(f);
    df = inline(df);
    x(1) = x0 - (f(x0)/df(x0));
    ex(1) = abs(x(1)-x0);
    k = 2;
    while (ex(k-1) >= tol) && (k <= nmax)
        x(k) = x(k-1) - (f(x(k-1))/df(x(k-1)));
        ex(k) = abs(x(k)-x(k-1));
        k = k+1;
    end
end