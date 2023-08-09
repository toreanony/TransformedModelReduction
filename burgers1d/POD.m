function [DV, SigV] = POD(snapshots, arg)
    % Input: 
    %   snapshots: each row is a snapshot
    %   arg: if arg >= 1, it indicate reduced order; otherwise the
    %   tolerance.
    % Output:
    %   DV: POD base matrix, each column is a POD base.
    %   SigV: singular values of POD snapshots matrix.

    SigV = svds(snapshots', size(snapshots, 1));
  
    if arg >= 1
        [DV, ~, ~] = svds(snapshots', arg);
    else
        kr = 1;
        while sum(SigV(1:kr)) / total < 1 - arg
            kr = kr + 1;
        end 
        DV = svds(snapshots', kr);
    end
end