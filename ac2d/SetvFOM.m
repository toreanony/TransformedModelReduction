% initialize the FOM of slow variable v for different kinds of boundary condition. 

if bdCase == 1
    bdModifier = [29 -14; -2 1];   % parital^2 v / parital n^2 = 0.
elseif bdCase == 2
    bdModifier = [15 0; -1 0];  % partial v / partial n = 1.
else 
    disp("please choose the boundary condition.");
end

A0 = -30*eye(n)+16*diag(ones(n-1,1), 1)-diag(ones(n-2,1), 2);
A0(1:2, 1:2) = A0(1:2,1:2) + bdModifier; 
A0 = A0 + rot90(A0, 2);  % the middle 30 of another axis is added here. 
A0 = sparse(A0); 
Ad = repmat({A0}, n, 1);
Av = sparse(1:N-n, n+1:N, 16*ones(N-n, 1), N, N) - sparse(1:N-2*n, 2*n+1:N, ones(N-2*n, 1), N, N);
Av(1:2*n,1:2*n)=Av(1:2*n,1:2*n)+sparse(kron(bdModifier, eye(n)));
Av = Av + rot90(Av, 2) + blkdiag(Ad{:});  
Av = Av/12/h^2;
clear Ad A0 

if bdCase == 1
    dyA = diag(ones(n-1, 1), 1); dyA(1, 1:2) = [-2 2];
    dyA = sparse(dyA - rot90(dyA, 2))/2/h;
    Fv = @(v) sqrt(2)/epsilon*tanh(v/sqrt(2)/epsilon).*(1 - ([2*v(n+1:2*n) - 2*v(1:n); ...
        v(2*n+1:N) - v(1:N-2*n); 2*v(N-n+1:N) - 2*v(N-2*n+1:N-n)]/2/h).^2 - ...
        reshape(dyA*reshape(v, n, n), N, 1).^2);
    clear dyA
elseif bdCase == 2
    dyA = diag(ones(n-1, 1), 1); dyA(1, 2) = 0;
    dyA = sparse(dyA - rot90(dyA, 2))/2/h;
    bA = zeros(n, 1); bA(1:2) = [7/6/h;-1/12/h]; bA = sparse(bA + flip(bA));
    Fv = @(v) sqrt(2)/epsilon*tanh(v/sqrt(2)/epsilon).*(1 - ([-2*h*ones(n, 1); ...
        v(2*n+1:N)-v(1:N-2*n);2*h*ones(n, 1)]/2/h).^2 - ...
        (reshape(dyA*reshape(v,n,n), N, 1)+repmat([-1;zeros(n-2,1);1], n, 1)).^2) + ...
        repmat(bA, n, 1) + kron(bA, ones(n, 1));  % from laplace term. 
    clear dyA bA
else 
    disp("please choose the boundary condition.");
end

Rv = @(v) Av * v + Fv(v);
