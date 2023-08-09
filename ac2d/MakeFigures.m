idx = [1, 500, 900];
time = T.dt*T.tMiddleStep*(0:args.shotsNum-1);
close all

if ~exist(".\figures", "dir"), mkdir(".\figures"); end

%%
% POD
figure;
semilogy(sigU(1:500), "b*-");
hold on
semilogy(sigV(1:500), "r+-");
legend("PODsigU", "PODsigV");
xlabel("index",'Interpreter','LaTex');
ylabel("singular value",'Interpreter','LaTex');
set(gca, "FontSize", 17); 
savefig("./figures/oval_ac2_sig_compare.fig");

% DEIM
figure;
semilogy(sigFu(1:500), "b*-");
hold on
semilogy(sigFv(1:500), "r+-");
legend("DEIMsigU", "DEIMsigV");
xlabel("index",'Interpreter','LaTex');
ylabel("singular value",'Interpreter','LaTex');
set(gca, "FontSize", 17);
savefig("./figures/oval_ac2_sig_compare_deim.fig");

%%
close all
for j = 1:length(idx)
    i = idx(j);
    figure;
    pcolor(xSpan, ySpan, reshape(uFull(:,i)-uvRedu(:,i), n, n));
%     title("residual of $U_{appr}$", "Interpreter","latex");
    xlabel("$x$", "Interpreter","latex");
    ylabel("$y$", "Interpreter","latex");
%     zlim([-1.1 1.1])
    shading flat
    colorbar
    axis square
    set(gca, "FontSize", 17);
    savefig(strcat("./figures/oval_ac2_uvres_",...
    num2str(round(time(i), 2)), ".fig"));
    pause(0.5)
end
%%
close all
for j = 1:length(idx)
    i = idx(j);
    figure;
    pcolor(xSpan, ySpan, reshape(uFull(:,i)-uRedu(:,i), n, n));
        xlabel("$x$", "Interpreter","latex");
    ylabel("$y$", "Interpreter","latex");
    shading flat
    colorbar
    axis square
    set(gca, "FontSize", 17);
    savefig(strcat("./figures/oval_ac2_ures_",...
    num2str(round(time(i), 2)), ".fig"));
    pause(0.5)
end

%%
close all
for j = 1:length(idx)
    i = idx(j);
    figure;
    surf(xSpan, ySpan, reshape(uFull(:,i), n, n));
        xlabel("$x$", "Interpreter","latex");
    ylabel("$y$", "Interpreter","latex");
    zlim([-1.1 1.1])
    shading flat
    set(gca, "FontSize", 17);
    savefig(strcat("./figures/oval_ac2_u_",...
    num2str(round(time(i), 2)), ".fig"));
    pause(0.5)
end


