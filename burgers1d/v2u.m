function v = u2v(vu, u)
% this is transform from u to v
% input: 
%       vu: v(u) expression 
%       u: u row solution 
% output: 
%       v: v row solution

    [n, m] = size(u);
    v = zeros(n , m);
    
    for i = 1:n
        v(i, :) = vu(u(i, :));
    end 
end

