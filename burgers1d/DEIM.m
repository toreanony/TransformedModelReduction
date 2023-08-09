function [U, UM, P, SigF] = DEIM(F, D, arg)
% Q-DEIM process
% input:  
%   F: conduct column vector 
%   D: column solution
% output: 
%   U: DEIM basis
%   UM: the U*inv(U(P,:))
%   P: selected rows of function F
%   SigF: singular values of F(D)
    [k, num] = size(D);
    % F(D) basis
    U = cell2mat(cellfun(F, num2cell(D, 1),'uni',false));
    % DEIM basis
    SigF = svds(U, num);
    total = sum(SigF);
    % 如果arg>=1, 则直接用于指定取基个数. 否则作为tolerance确定基数
    if arg >= 1
        [U, ~, ~] = svds(U, arg);
        m = arg;
    else
        m = 1;
        while sum(SigF(1:m))/total < 1 - arg
            m = m + 1;
        end
        [U, ~, ~] = svds(U, m);
    end


%     % DEIM main process
%     P = zeros(m, 1);
%     u1 = U(:, 1);
%     [rho, P(1)] = max(abs(u1));
%     M_inv = 1/rho;
%     
%     for i = 2:m
%         ui = U(:, i);
%         r = abs(ui - U(:, 1:i-1) * M_inv * ui(P(1:i-1)));     
%         [rho, P(i)] = max(r);
%         while find(P(1:i-1) == P(i))
%             r(P(i)) = -1;
%             [rho, P(i)] = max(r);
%         end
%         M_inv = getNewM_inv(U, M_inv, rho, P(1:i), i);
%     end
%     UM = U * M_inv;

    % q_deim from "A New Selection Operator 
    % for the Discrete Empirical Interpolation Method---
    % Improved A Priori Error Bound and Extensions"
    [~, R, RP] = qr(U', 'vector');
    P = RP(1:m);
    UM = [eye(m) ; (R(:,1:m)\R(:,m+1:k))'];
    Pinverse(RP) = 1 : k ; UM = UM(Pinverse,:);

    function M_inv = getNewM_inv(U, M_inv, rho, P, i)
    A = U(P(i), 1:i-1);
    C = M_inv * U(P(1:i-1), i);
%     % method1
%     M_inv = [eye(i-1), -C; zeros(1, i-1), 1] *...
%         [M_inv, zeros(i-1, 1); -A*M_inv/rho, 1/rho];
    
    % method2
    term1 = eye(i) + [C;-1] * [A, -1] / rho; 
    term1(i, i) = term1(i, i) - 1;
    term2 = eye(i);
    term2(1:i-1, 1:i-1) = M_inv;
    M_inv = term1 * term2;
end
end