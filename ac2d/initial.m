% this file initialize time region and spatial domain and all the variables generated from given parameters.

%%
tStart= 0;
tStop = 0.13;
xStart = -1; 
xStop = 1;
yStart = xStart;
yStop = xStop;


epsilon = 0.01;
dt = 5e-7;
n = 500;            % n = ceil(4/epsilon). each axis discretize number.
shotsNum = 1000;    % the number of needed snapshots.

deim_on = 1;  % 0: POD (no DEIM); 1: POD-DEIM (qDEIM).

% modify the boundary condition. 
% bdCase = 1: parital^2 v / parital n^2 = 0; parital^2 u / parital n^2 = 0.
% bdCase = 2: partial v / partial n = 1;     partial u / partial n = 0. 

bdCase = 2;  

% generate from parameters above. 
N = n^2;            % total size of FOM.
h = (xStop - xStart) / (n-1);
xSpan = (xStart:h:xStop); 
ySpan = (yStart:h:yStop);

tMiddleStep = round((tStop-tStart)/dt/shotsNum);  % middle steps. 
if tMiddleStep <= 1, tMiddleStep = 1; dt = (tStop - tStart)/shotsNum; end
tSpan = linspace(tStart, tStop, shotsNum);

uFileName = strcat("./data/u_n=",num2str(n),"_",num2str(epsilon),"_",num2str(tStop),".mat");
vFileName = strcat("./data/v_n=",num2str(n),"_",num2str(epsilon),"_",num2str(tStop),".mat");

%% initial functions and full order models. 

v2uFun = @(v) tanh(v/sqrt(2)/epsilon);  % transfer equation of AC.  

[xGrid, yGrid] = meshgrid(xSpan, ySpan);
v0 = (1.5*xGrid.^2 + yGrid.^2) - 0.3;  % the initial function of slow variable.
% v0 = (2*xGrid.^2 + yGrid.^2) - 0.3;  % the initial function of slow variable.
% v0 = (xGrid.^2/3*4 + yGrid.^2)-0.45;  % 可能是之前的19, 80 还可以的结果. 0.01, 200, 0.13(?)   

clear xGrid yGrid
v0 = reshape(v0, N, 1);
u0 = v2uFun(v0);  % origional variable u0. given by transform v0. 

SetvFOM 

SetuFOM 

%% deliver to structs, easier to use function. 
T.dt = dt; 
T.tStart= tStart;
T.tStop = tStop;
T.tMiddleStep = tMiddleStep;
T.tSpan = tSpan;

X.xStart = xStart;
X.xStop = xStop;
X.xSpan = xSpan;

Y.yStart = yStart;
Y.yStop = yStop;
Y.ySpan = ySpan;

args.epsilon = epsilon;
args.n = n;
args.N = N;
args.h = h;
args.shotsNum = shotsNum;
args.bdCase = bdCase;
args.deim_on = deim_on;

vFOM.A = Av; 
vFOM.F = Fv;
vFOM.R = Rv;

uFOM.A = Au;
uFOM.F = Fu;
uFOM.R = Ru;