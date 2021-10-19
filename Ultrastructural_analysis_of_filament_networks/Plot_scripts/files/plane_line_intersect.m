function [I,check]=plane_line_intersect(n,V0,P0,P1)
    % Check if multitoolbox is available
    persistent useMulti
    if isempty(useMulti)
        if exist('multiprod', 'file')
            useMulti = true;
        end
    end

    % Make sure n and V0 are column vectors
    if ~(size(n,1) == 3 && size(n,2) == 1)
        error('Normal must be a 3x1 vector')
    end
    if ~(size(V0,1) == 3 && size(V0,2) == 1)
        error('Point on the plane must be a 3x1 vector')
    end


    % Do some calculation
    nPoints = size(P0,2);
    u = P1-P0;
    w = P0 - repmat(V0, [1 nPoints]);
    if useMulti
        D = multiprod(n', u);
        N = -multiprod(n',w);
    else
        D = nan(1,nPoints);
        N = nan(1,nPoints);
        for i = 1:size(u,2)
            D(i) = dot(n,u(:,i));
            N(i) = -dot(n,w(:,i));
        end
    end
    
    % Check for parallel cases
    idxParallel = abs(D) < 10^-7;
    idxOnPlane = idxParallel & N==0;
    idxDisjoint = idxParallel & N~=0;

    % Ccompute the intersection parameter
    sI = nan(1,nPoints);
    sI(~idxParallel) = N(~idxParallel) ./ D(~idxParallel);
    I = P0+ repmat(sI, [3,1]).*u;

    % Check if intersection is between the points or outside them
    idxIntersectBetween = (sI(1,:) > 0 & sI(1,:) < 1);
    idxIntersectOutside = ~idxIntersectBetween & ~idxParallel;
    
    % Prepare check vector
    check = zeros(1, nPoints);
    check(idxOnPlane) = 2;
    check(idxDisjoint) = 0;
    check(idxIntersectBetween) = 1;
    check(idxIntersectOutside) = 3;

end
