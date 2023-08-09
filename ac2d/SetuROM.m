% generate the ROM of u. 
U = U_all(:, 1:PODnum);  % POD basis of u. 
u0r = pinv(U) * u0;
Aur = U' * Au * U;
if deim_on 
    D = U_alldeim(:, 1:DEIMnum);
    [~, R, RP] = qr(D', 'vector');
    m = DEIMnum;
    P = RP(1:m); 
    D = [eye(m) ; (R(:,1:m)\R(:,m+1:N))'];
    Pinverse(RP) = 1 : N ; D = D(Pinverse,:);  % this is D*(P^TD)^{-1}.
    coeFu = U' * D; 
    uDEIMrows = U(P, :);  % for u, only need to take out the P rows.
    Fur = @(ur) coeFu * (uDEIMrows*ur - (uDEIMrows*ur).^3)/epsilon^2;
else
    Fur = @(ur) U' * Fu(U * ur);
end
Rur = @(ur) Aur * ur + Fur(ur);
uROM.PODbasis = U;
uROM.A = Aur;
uROM.F = Fur;
uROM.R = Rur;