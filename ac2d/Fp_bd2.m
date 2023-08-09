function Fdeim = Fp_bd2(V, vr, P, args)
% DEIM for transformed variable equation
% for pupn=1 condition
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
            v_x = -1;  % 左边界
            Fdeim(k) = Fdeim(k) + 7/6/h;
        elseif i == n
            v_x = 1;   % 右边界
            Fdeim(k) = Fdeim(k) + 7/6/h;
        elseif i == 2
            v_x = (V(id+n,:)-V(id-n,:))*vr/2/h;
            Fdeim(k) = Fdeim(k) - 1/12/h;
        elseif i == n-1
            v_x = (V(id+n,:)-V(id-n,:))*vr/2/h;
            Fdeim(k) = Fdeim(k) - 1/12/h;
        else
            v_x = (V(id+n,:)-V(id-n,:))*vr/2/h;
        end
        if j == 1
            v_y = -1;  % 下边界
            Fdeim(k) = Fdeim(k) + 7/6/h;
        elseif j == n
            v_y = 1;   % 上边界
            Fdeim(k) = Fdeim(k) + 7/6/h;
        elseif j == 2
            v_y = (V(id+1,:)-V(id-1,:))*vr/2/h;
            Fdeim(k) = Fdeim(k) - 1/12/h;
        elseif j == n-1
            v_y = (V(id+1,:)-V(id-1,:))*vr/2/h;
            Fdeim(k) = Fdeim(k) - 1/12/h;
        else
            v_y = (V(id+1,:)-V(id-1,:))*vr/2/h;
        end
        Fdeim(k) = Fdeim(k) + sqrt(2)/args.epsilon*tanh(v_id/sqrt(2)/args.epsilon)*(1-v_x^2-v_y^2);
    end
end