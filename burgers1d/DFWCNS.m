function df = DFWCNS(u, epsilon, dx)
% generate derivative of f(u). using WCNS scheme. 
    beta = 64/45; eta = dx^2;
    f = 0.5 * u.^2;
    n = length(u);
    fedgel = 0.5 * u(1)^2;
    fedger = 0.5 * u(n)^2;
    
    S0p = zeros(n+1, 1); S1p = zeros(n+1, 1); S2p = zeros(n+1, 1);
    S0p(4:n+1) = 1/8 * (3*f(1:n-2) -10*f(2:n-1) + 15*f(3:n));
    S0p(2) = 1/8 * (-7*fedgel + 15*f(1));
    S0p(3) = 1/8 * (3*fedgel -10*f(1) + 15*f(2));
    S0p(1) = fedgel;
    S1p(3:n) = 1/8 * (-f(1:n-2) + 6*f(2:n-1) + 3*f(3:n));
    S1p(2) = 1/8 * (-fedgel + 6*f(1) + 3*f(2));
    S1p(n+1) = 1/8 * (-f(n-1) + 6*f(n) + 3*fedger);
    S1p(1) = (5*fedgel + 3*f(1))/8;
    S2p(2:n-1) = 1/8 * (3*f(1:n-2) + 6*f(2:n-1) - f(3:n));
    S2p(n) = 1/8 * (3*f(n-1) + 6*f(n) - fedger);
    S2p(n+1) = 1/8 * (3*f(n) + 5*fedger);
    S2p(1) = (3*fedgel + 6*f(1) - f(2))/8;

    S0n = zeros(n+1, 1); S1n = S2p; S2n = S1p;
    S0n(2:n-2) = 1/8 *(15*f(2:n-2) - 10*f(3:n-1) + 3*f(4:n));
    S0n(n-1) = 1/8 * (15*f(n-1) - 10*f(n) + 3*fedger);
    S0n(n) = 1/8 * (15*f(n) - 7*fedger);
    S0n(n+1) = fedger;
    S0n(1) = (15*f(1)-10*f(2)+3*f(3))/8;

    elementVec = zeros(n, 1);
    elementVec(2:n-1) = 13*(f(1:n-2)-2*f(2:n-1)+f(3:n)).^2/12;
    elementVec(1) = 13*(fedgel-2*f(1)+f(2))^2/12;
    elementVec(n) = 13*(f(n-1)-2*f(n)+fedger)^2/12;

    elementVec1 = zeros(n, 1);
    elementVec1(2:n-1) = 1/4*(f(1:n-2)-4*f(2:n-1)+3*f(3:n)).^2;
    elementVec1(1) = 1/4*(fedgel - 4*f(1) + 3*f(2))^2;
    elementVec1(n) = 1/4*(f(n-1) - 4*f(n) + 3*fedger)^2;

    elementVec2 = zeros(n, 1);
    elementVec2(2:n-1) = 1/4*(3*f(1:n-2)-4*f(2:n-1)+f(3:n)).^2;
    elementVec2(1) = 1/4*(3*fedgel - 4*f(1) + f(2))^2;
    elementVec2(n) = 1/4*(3*f(n-1) - 4*f(n) + fedger)^2;

    elementVec3 = zeros(n, 1);
    elementVec3(2:n-1) = 1/4*(f(1:n-2) - f(3:n)).^2;
    elementVec3(1) = 1/4*(fedgel - f(2))^2;
    elementVec3(n) = 1/4*(f(n-1) - fedger)^2;

    IS0p = zeros(n+1, 1); IS1p = zeros(n+1, 1); IS2p = zeros(n+1, 1);
    IS0p(3:n+1) = elementVec1(1:n-1) + elementVec(1:n-1);
    IS0p(2) = 10*(-fedgel+f(1))^2/3;
    IS1p(2:n+1) = elementVec3 + elementVec;
    IS1p(1) = 4*(fedgel-f(1))^2/3;
    IS2p(1:n) = elementVec2 +elementVec;
    IS2p(n+1) = 10*(f(n)-fedger)^2/3;

    IS0n = zeros(n+1, 1); IS1n = zeros(n+1, 1); IS2n = zeros(n+1, 1);
    IS0n(1:n) = IS2p(2:n+1); 
    IS1n(1:n) = IS1p(2:n+1); 
    IS1n(n+1) = 1/4*(fedger-f(n))^2 + 13*(f(n)-fedger)^2/12;
    IS2n(1:n) = IS0p(2:n+1);
    IS2n(n+1) = elementVec1(n) + elementVec(n);

    coe0p = 1 ./(eta + IS0p).^2;
    coe1p = 10./(eta + IS1p).^2;
    coe2p = 5 ./(eta + IS2p).^2;
    sump = coe0p + coe1p + coe2p;
    coe0n = 1 ./(eta + IS0n).^2;
    coe1n = 10./(eta + IS1n).^2;
    coe2n = 5 ./(eta + IS2n).^2;
    sumn = coe0n + coe1n + coe2n;
    coe0p = coe0p ./ sump;
    coe1p = coe1p ./ sump;
    coe2p = coe2p ./ sump;
    coe0n = coe0n ./ sumn;
    coe1n = coe1n ./ sumn;
    coe2n = coe2n ./ sumn;

    % flux
    flux = (coe0p .* S0p + coe1p .* S1p + coe2p .* S2p) ...
         + (coe0n .* S0n + coe1n .* S1n + coe2n .* S2n);
    df1 = beta*(flux(2:n+1)-flux(1:n))/dx;
    
    df2 = zeros(n, 1); df2(2:n-1) = f(3:n) - f(1:n-2);
    df2(1) = f(2) - fedgel; df2(n) = fedger - f(n-1);
    df2 = (192-175*beta)*df2/(256*dx);

    df3 = zeros(n, 1); df3(3:n-2) = f(5:n) - f(1:n-4);
    df3(1) = f(3) - fedgel; df3(n) = fedger - f(n-2);
    df3(2) = f(4) - fedgel; df3(n-1) = fedger - f(n-3);
    df3 = (35*beta-48)*df3/(320*dx);

    df4 = zeros(n, 1); df4(4:n-3) = f(7:n) - f(1:n-6);
    df4(1) = f(4) - fedgel; df4(n) = fedger - f(n-3);
    df4(2) = f(5) - fedgel; df4(n-1) = fedger - f(n-4);
    df4(3) = f(6) - fedgel; df4(n-2) = fedger - f(n-5);
    df4 = (64 - 45*beta)*df4/(3840*dx);

    b = zeros(n, 1);
    b(1) = 5*epsilon*u(1)/4/dx^2;
    b(n-1) = -epsilon*u(n)/12/dx^2;
    b(2) = -epsilon*u(1)/12/dx^2;
    b(n) = 5*epsilon*u(n)/4/dx^2;

    df = df1/2 + df2 + df3 + df4 - b;
end