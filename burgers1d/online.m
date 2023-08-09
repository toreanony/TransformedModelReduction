function reduSol = online(reduceElements, v0r, t_span)
% this is the online procedure 
% input: 
%       reduceElements: V, Er, Ar, Fr, (only these are needed)
%       v0r
%       t_span
% output: 
%       reduSol: reduced order solution

    Er = reduceElements.Er;
    Ar = reduceElements.Ar;
    Fr = reduceElements.Fr;
    V  = reduceElements.V;

    Rr = @(t, v) Ar * v + Fr(v);
    opt = odeset('Mass', Er);
    [~, reduSol] = ode15s(Rr, t_span, v0r, opt);  % call solver
    reduSol = reduSol*V';% 1D row solutions. 
end