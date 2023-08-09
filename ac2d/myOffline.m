function [snapshots, PODsig,PODbasis, DEIMsig, DEIMbasis] = myOffline(FOM, initFunc, T, args)
% This function compute snapshots, singular values,  some POD basis and DEIM basis. In the online process only use part of the basis.
% input:
%     FOM: struct of full order model. (coef matrix A, nonlinear function F, rhs function R)
%     initFunc: initial function. 
%     T:   struct of parameters related to time. (tStart, tStop, dt, tMiddleStep, tSpan)
%     args: struct of parameters in need. (epsilon, n, N, h, shotsNum, bdCase, deim_on)
% output:
%     snapshots: a matrix, each column is a snapshot corresponds to a time step.
%     PODsig: column vector of POD singular values.
%     DEIMsig: column vector of DEIM singular values.
%     PODbasis: a matrix, each column is a POD mode.
%     DEIMbasis: a matrix, each column is a DEIM mode. 

    if isrow(initFunc), initFunc = initFunc'; end  % make sure it is column vector. 
    snapshots = TTY_RK(FOM.R, initFunc, T);
    [PODbasis, PODsig] = POD(snapshots', args.n);  % no need to use more than n basis. 
    if args.deim_on
        DEIMsnapshots = cell2mat(cellfun(FOM.F, num2cell(snapshots, 1), "uni", false));
        DEIMsig = svds(DEIMsnapshots, args.shotsNum);
        [DEIMbasis,~,~] = svds(DEIMsnapshots, args.n);
    else
        DEIMsig = 0;
        DEIMbasis = 0;
    end
end 