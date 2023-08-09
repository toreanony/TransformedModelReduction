% generate the ROM of v. 

V = V_all(:, 1:PODnum);  % POD basis of v.
v0r = pinv(V) * v0;
Avr = V' * Av * V;
if deim_on 
    D = V_alldeim(:, 1:DEIMnum);
    [~, R, RP] = qr(D', 'vector');
    m = DEIMnum;
    P = RP(1:m); 
    D = [eye(m) ; (R(:,1:m)\R(:,m+1:N))'];
    Pinverse(RP) = 1 : N ; D = D(Pinverse,:);  % this is D*(P^TD)^{-1}.
    coeFv = V' * D;  
    
    if bdCase == 1
        Fvr = @(vr) coeFv * Fp_bd1(V, vr, P, args);
    elseif bdCase == 2
        Fvr = @(vr) coeFv * Fp_bd2(V, vr, P, args);
    else
        Fvr = @(vr) coeFv * Fp_bad(Fv, V, vr, P);
    end
else
    Fvr = @(vr) V'*Fv(V*vr);
end
Rvr = @(vr) Avr * vr + Fvr(vr);
vROM.PODbasis = V;
vROM.A = Avr; 
vROM.F = Fvr;
vROM.R = Rvr;
