% main program.

clear; clc; 
close all;

initial  % modify inside parameters and initial functions. 

% offline of slow variable v. 
[vFull, sigV, V_all, sigFv, V_alldeim] = myOffline(vFOM, v0, T, args);
uvFull = v2uFun(vFull);  % U^v in the thesis. 

% offline of origional variable u.
[uFull, sigU, U_all, sigFu, U_alldeim] = myOffline(uFOM, u0, T, args);

% % save data. 
% save Udeim.mat U_alldeim sigFu
% save Vdeim.mat V_alldeim sigFv
% 
% save ufull.mat uFull
% save vfull.mat vFull
% 
% save U_all.mat U_all sigU
% save V_all.mat V_all sigV
% 
% save parameters.mat args X Y T 
%% choose POD and DEIM basis number by tolerance. 
% tol = 1e-10;
% tolPOD = cumsum(sigV.^2)/sum(sigV.^2);
% PODnum = find(1- tolPOD<tol, 1);
% tolDEIM = cumsum(sigFv.^2)/sum(sigFv.^2);
% DEIMnum = find(1 - tolDEIM<tol, 1);

PODnum = 18;
DEIMnum = 30;

%% compute the reduced solution.
SetvROM
% online of slow variable v. 
vRedu = myOnline(vROM, v0r, T);
uvRedu = v2uFun(vRedu);  % U^v_{appr} in the thesis.


SetuROM
% online of origional variable u.
uRedu = myOnline(uROM, u0r, T);

% MakeFigures

