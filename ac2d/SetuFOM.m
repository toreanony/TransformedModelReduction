% initialize the FOM of slow variable v for different kinds of boundary condition. 

if bdCase == 1
    bdModifier = [29 -14; -2 1];   % parital^2 u / parital n^2 = 0.
elseif bdCase == 2
    bdModifier = [15 0; -1 0];  % partial u / partial n = 0.
else 
    disp("please choose the boundary condition.");
end

A0 = -30*eye(n)+16*diag(ones(n-1,1), 1)-diag(ones(n-2,1), 2);
A0(1:2, 1:2) = A0(1:2,1:2) + bdModifier; 
A0 = A0 + rot90(A0, 2);  % the middle 30 of another axis is added here. 
A0 = sparse(A0); 
Ad = repmat({A0}, n, 1);
Au = sparse(1:N-n, n+1:N, 16*ones(N-n, 1), N, N) - sparse(1:N-2*n, 2*n+1:N, ones(N-2*n, 1), N, N);
Au(1:2*n,1:2*n)=Au(1:2*n,1:2*n)+sparse(kron(bdModifier, eye(n)));
Au = Au + rot90(Au, 2) + blkdiag(Ad{:});  
Au = Au/12/h^2;
clear Ad A0 

Fu = @(u) (u - u.^3) / epsilon^2;

Ru = @(u) Au * u + Fu(u);