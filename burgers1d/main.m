% main program in offline-online strategy for 1D Buger's equation
% author: t.y.Tang
% date: 2022/12/11

clear; clc;
close all;

%% initial settings. 
X.start = 0; X.end = 1; X.steps = 1000;
T.start = 0; T.end = 1; T.steps = 1000;

x0 = 0.5; 

args.baseChosePOD = 10;
args.baseChoseDEIM = 20;
args.DEIM_on = 1;
args.epsilon = 1e-4;
args.x_scale = 1;
args.t_scale = 1;

u0Fun = @(x) 0.5*(1-tanh((x - x0)/4/args.epsilon));
v0Fun = @(x) (x - x0);
uvFun = @(v) 0.5*(1-tanh(v/4/args.epsilon));

dxU = (X.end-X.start)/(X.steps-1)/args.x_scale;
dtU = (T.end-T.start)/(T.steps-1)/args.t_scale;
x_spanU = (X.start:dxU:X.end);
t_spanU = (T.start:dtU:T.end);
x_spanV = linspace(X.start, X.end, X.steps);
t_spanV = linspace(T.start, T.end, T.steps);
[elementsU, elementsV] = initialElements(X, args);
u0 = u0Fun(x_spanU);
v0 = v0Fun(x_spanV);
if isrow(u0), u0 = u0'; end
if isrow(v0), v0 = v0'; end
%% exact solution
exac = zeros(length(t_spanV), length(x_spanV));
exac(1,:) = u0Fun(x_spanV);
for i = 2:length(t_spanV)
    exac(i,:) = 0.5*(1-tanh((x_spanV-x0-0.5*t_spanV(i))/4/args.epsilon));
end
%% WCNS (ode15s)
[u_full, reducedElementsU] = offline(elementsU, u0, t_spanU, args);
u0r = pinv(reducedElementsU.V) * u0;
u_redu = online(reducedElementsU, u0r, t_spanU);
u_full = u_full(1:args.t_scale:end, 1:args.x_scale:end);
u_redu = u_redu(1:args.t_scale:end, 1:args.x_scale:end);

%% U^v, U^v_appr
[v_full, reducedElementsV] = offline(elementsV, v0, t_spanV, args);
v0r = pinv(reducedElementsV.V) * v0;
v_redu = online(reducedElementsV, v0r, t_spanV);
uv_full = v2u(uvFun, v_full);
uv_redu = v2u(uvFun, v_redu);

%% test new u0
nx0 = 0.3;
nu0Fun = @(x) 0.5*(1-tanh((x - nx0)/4/args.epsilon));
nv0Fun = @(x) (x - nx0);

nexac = zeros(length(t_spanV), length(x_spanV));
nexac(1,:) = nu0Fun(x_spanV);
for i = 2:length(t_spanV)
    nexac(i,:) = 0.5*(1-tanh((x_spanV-nx0-0.5*t_spanV(i))/4/args.epsilon));
end

nv0 = nv0Fun(x_spanV); if isrow(nv0), nv0 = nv0'; end
nv0r = pinv(reducedElementsV.V) * nv0;
nv_redu = online(reducedElementsV, nv0r, t_spanV);
nuv_redu = v2u(uvFun, nv_redu);
%%
if ~exist("./data", "dir"),    mkdir("./data");    end
if ~exist("./figures", "dir"), mkdir("./figures"); end
dataName = strcat("./data/", ...
    "d=", num2str(args.DEIM_on),...
    "_x0=",num2str(x0),...
    "_epsilon=" , num2str(args.epsilon),...
    ".mat");

save(dataName);
