% 1: excat solution
figure; pcolor(x_spanV, t_spanV, exac); 
colorbar; % title("exact solution");
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex'); shading flat;
set(gca, "FontSize", 17); savefig("./figures/burgers_exact.fig");

% 2: v u full solution
% 2.1
figure; 
pcolor(x_spanV, t_spanV, u_full);
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex');
shading flat; colorbar; % u full
set(gca, "FontSize", 17); savefig("./figures/burgers_original.fig");
% 2.2
figure;
pcolor(x_spanV, t_spanV, v_full);
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex');
shading flat; colorbar; % v full
set(gca, "FontSize", 17); savefig("./figures/burgers_transformed.fig");
% 2.3
figure;
pcolor(x_spanV, t_spanV, v2u(uvFun, v_full));
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex');
shading flat; colorbar; % u full
set(gca, "FontSize", 17); savefig("./figures/burgers_phi(v).fig");

% 3: residual wcns POD-qDEIM
figure;
pcolor(x_spanV, t_spanV, exac - u_redu);
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex');
shading flat; colorbar; % residual of wcns scheme
set(gca, "FontSize", 17); savefig("./figures/burgers_original_residual.fig");

% 4: residual rescaled POD-qDEIM
figure;
pcolor(x_spanV, t_spanV, exac - uv_redu);
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex');
shading flat; colorbar; % residual of our method
set(gca, "FontSize", 17); savefig("./figures/burgers_transformed_residual.fig");

%% 5: sigular value comparison
figure;
semilogy(reducedElementsU.sigV(1:500), 'b*-');
hold on
semilogy(reducedElementsV.sigV(1:500), 'r+-');
h = legend("PODsigU", "PODsigV", FontSize=13);
% set(h, "Box", "off");
xlabel("index",'Interpreter','LaTex'); ylabel("singular value",'Interpreter','LaTex');
set(gca, "FontSize", 17); savefig("./figures/burgers_sig_compare.fig");

%% 6: singular values of DEIM snapshots
figure;
semilogy(reducedElementsU.sigF(1:500), 'b*-');
hold on
semilogy(reducedElementsV.sigF(1:500), 'r+-');
h = legend("DEIMsigU", "DEIMsigV",'location', "best", FontSize=13);
% set(h, "Box", "off");
xlabel("index",'Interpreter','LaTex'); ylabel("singular value",'Interpreter','LaTex');
set(gca, "FontSize", 17); savefig("./figures/burgers_sig_compare_deim.fig");

%% 7: new v0 with different x0
figure;
pcolor(x_spanV, t_spanV, nexac - nuv_redu);
xlabel("$x$",'Interpreter','LaTex'); ylabel("$t$",'Interpreter','LaTex');
shading flat; colorbar; % residual of wcns scheme
set(gca, "FontSize", 17); savefig("./figures/burgers_newInitialFunc_residual.fig");

% 8: error over time original
errorOverTime1 = zeros(length(t_spanV), 1);
for i = 1:length(t_spanV)
    errorOverTime1(i) = computeError(exac(i,:), v_redu(i,:));
end

% 9: error over time vu redu
errorOverTime = zeros(length(t_spanV), 1);
for i = 1:length(t_spanV)
    errorOverTime(i) = computeError(exac(i,:), uv_redu(i,:));
end
figure;
semilogy(t_spanV, errorOverTime1, "b*-");
hold on
semilogy(t_spanV, errorOverTime, "r+-");
legend("$U_{appr}$", "$U_{appr}^v$", "Interpreter","latex", FontSize=13);
xlabel("$t$", "Interpreter","latex");
ylabel("relative error", "Interpreter","latex");
set(gca, "FontSize", 17); savefig("./figures/burgers_error_over_t.fig");
