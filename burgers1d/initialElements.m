function [elementsU, elementsV] = initialElements(X, args)
% initial values and parameters
% input: 
%       X: struct, X.start, X.steps, X.end 
%       T: struct, t_start, t_steps, t_end 
% output: 
%       elementsV: Ev, Av, Fv
%       elementsU: Eu, Au, Fu

    epsilon = args.epsilon;
    x_scale = args.x_scale;

    %% initial u for WCNS scheme
    dx = (X.end - X.start) /x_scale/(X.steps - 1);  % signed distance function
    x_span = (X.start:dx:X.end);  
    n = length(x_span);

    Eu = eye(n);

    Au = -30*eye(n) - diag(ones(n-2, 1), 2) - diag(ones(n-2, 1), -2) ...
        + 16*diag(ones(n-1,1), 1) + 16*diag(ones(n-1,1), -1);

    Au = epsilon * Au / (12*dx^2);    
%     Au = sparse(Au);

    Fu = @(u) -DFWCNS(u, epsilon, dx);

    elementsU.E = Eu;
    elementsU.A = Au;
    elementsU.F = Fu;
    %% initial v
    dx = (X.end - X.start) / (X.steps - 1);  % signed distance function
    x_span = (X.start:dx:X.end);  
    n = length(x_span);

    Ev = eye(n);

    Av = -30*eye(n) - diag(ones(n-2, 1), 2) - diag(ones(n-2, 1), -2) ...
        + 16*diag(ones(n-1,1), 1) + 16*diag(ones(n-1,1), -1);
    Av(1, 1) = -1;  Av(1, 2) = 2;
    Av(n, n-1) = 2; Av(n, n) = -1;
    Av(2, 1) = 14; Av(2, 2) = -29;
    Av(n-1,n) = 14; Av(n-1, n-1) = -29;
    Av = epsilon * Av / (12*dx^2);
%     Av = sparse(Av);

    Fv = @(v) -0.5*((1-tanh(v/(4*epsilon))) ...
                      .*([v(3)-v(1);v(3:n)-v(1:n-2);v(n)-v(n-2)]/2/dx) ...
                   +tanh(v/(4*epsilon)) ...
                      .*(([v(3)-v(1);v(3:n)-v(1:n-2);v(n)-v(n-2)]/2/dx).^2));

    elementsV.E = Ev;
    elementsV.A = Av;
    elementsV.F = Fv;
end