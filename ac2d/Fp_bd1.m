function Fdeim = Fp_bd1(V, vr, P, args)
% DEIM for transformed variable equation
    h = args.h; n = args.n;
    Fdeim = zeros(size(P)); 
    if isrow(Fdeim), Fdeim = Fdeim'; end
    for k = 1:length(P)
        id = P(k) - 1;  % in-vector index
        j = mod(id, n);  % column index
        i = (id - j) / n;  % row index
        id = id + 1; j = j + 1; i = i + 1; 
        v_id = V(id, :) * vr;
        if i == 1
            v_x = (V(id+n,:)*vr - v_id)/h;
        elseif i == n
            v_x = (v_id - V(id-n,:)*vr)/h;
        else
            v_x = (V(id+n,:)-V(id-n,:))*vr/2/h;
        end
        if j == 1
            v_y = (V(id+1,:)*vr - v_id)/h;
        elseif j == n
            v_y = (v_id - V(id-1,:)*vr)/h;
        else
            v_y = (V(id+1,:)-V(id-1,:))*vr/2/h;
        end
        Fdeim(k) = tanh(v_id/sqrt(2)/args.epsilon)*(1-v_x^2-v_y^2);
    end
    Fdeim = sqrt(2)/args.epsilon*Fdeim;
end
