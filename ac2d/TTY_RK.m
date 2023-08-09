function v = TTY_RK(R, v0, T)
    % RK4, RK24
    % input: 
    %   R: right hand side function
    %   v0: column initial condition
    %   T: time struct (tStart, tStop, dt, tMiddleStep, tSpan)
    % output:
    %   v: column solution, each column corresponds to a time step in `tSpan`.
    
    N = length(v0);
    stepNum = length(T.tSpan);
    v = zeros(N, stepNum);
    v(:, 1) = v0;
    
    % RK4
    for i = 2:stepNum
        temp = v(:, i-1);
        for j = 2:T.tMiddleStep
            K0 = temp;
            K1 = R(K0);
            K2 = R(K0 + T.dt*K1/2);
            K3 = R(K0 + T.dt*K2/2);
            K4 = R(K0 + T.dt*K3);
            temp = K0 + T.dt*(K1+2*K2+2*K3+K4)/6;
        end
        v(:, i) = temp;
%         if mod(i, 10) == 0, i, end
        if any(isnan(temp)), warning("NaN failed!"), break; end
    %     TestScript
    end
    
end