function [fullSol, reducedElements] = offline(elements, v0, t_span, args)
% this is the offline procedure
% input:
%       elements:
%           E, A, F
%       v0
%       t_span 
%       args: 
%           baseChose: base num
% output: 
%       fullSol: full order solution
%       redusedElements: V, Er, Ar, Fr, sigV, sigF

    E = elements.E;
    A = elements.A;
    F = elements.F;

    baseChosePOD = args.baseChosePOD;
    baseChoseDEIM = args.baseChoseDEIM;
    DEIM_on = args.DEIM_on;  % chose if DEIM

    %% solver
    R = @(t, v) A*v + F(v);
    opt = odeset('Mass', E);
    [~, fullSol] = ode15s(R, t_span, v0, opt);  % row solution

    %% POD
    [V, sigV] = POD(fullSol, baseChosePOD);

    %% reduced elements
    reducedElements.V  = V;
    reducedElements.Er = V' * E * V;
    reducedElements.Ar = V' * A * V;
    % DEIM
    if DEIM_on
        % NOTICE: need to deliver column solution to DEIM.
        [~, UM, P, sigF] = DEIM(F, fullSol', baseChoseDEIM);
        VtU = V' * UM;
        reducedElements.sigF = sigF;
        reducedElements.Fr = @(v) rF(VtU, P, F, V, v);  % POD-DEIM
    else
        reducedElements.Fr = @(v) V' * F( V * v);  % only POD
    end
    reducedElements.sigV = sigV;

    %% rough nonlinear function reduce by using DEIM
    function rf = rF(leftMat, P, F, V, u)
        rf = F(V*u);
        rf = rf(P);
        rf = leftMat * rf;
    end    
end