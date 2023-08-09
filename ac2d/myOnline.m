function approSol = myOnline(ROM, reduInitFunc, T)
% This function computes the approximate solution. 
% input:
%     ROM: struct of reduced order model. (coef matrix A, nonlinear function F, rhs function R)
%     reduInitFunc: reduced initial function. 
%     T:   struct of parameters related to time. (tStart, tStop, dt, tMiddleStep, tSpan)
%     args: struct of parameters in need. (epsilon, n, N, h, shotsNum, bdCase, deim_on)    
% output:
%     approSol: approximate solution by POD or POD-qDEIM, has the same size to full order solution. 

%     approSol = TTY_RK(ROM.R, reduInitFunc, T);

    rhs = @(t, y) ROM.R(y);
    [~, approSol] = ode15s(rhs, T.tSpan, reduInitFunc);
    approSol = approSol';

    approSol = ROM.PODbasis * approSol;  % lift back to full order. 
end