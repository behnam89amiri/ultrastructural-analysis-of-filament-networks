function P = sliceDelaunay(d, varargin)
    % Make sure type is okay
    if ~isa(d, 'delaunayTriangulation')
        error('First argument must be of delaunayTriangulation type');
    end
    
    % Get the triangles of delaunay structure
    t = d.convexHull;
    p = d.Points;
    
    % Find all possible lines of the mesh (between each triangle points 1,
    % 2 et 3), starting points are on 1st, 3rd and 5th columns while end
    % points are on 2nd, 4th and 6th (This will be usefull later)
    t_tp = t(:,[1 2 1 3 2 3]); 

    % Find what is the cutting plane
    if ischar(varargin{1})
        useStandardizedPlane = true;
        if strcmp(varargin{1}, 'x')
            plane = [1 0 0; varargin{2} 0 0];
        elseif strcmp(varargin{1}, 'y')
            plane = [0 1 0; 0 varargin{2} 0];
        elseif strcmp(varargin{1}, 'z')
            plane = [0 0 1; 0 0 varargin{2}];
        else
            errordlg('Plan non reconnu');
            error('Plan non reconnu');
        end
    else % If normal/point formalism is used
        useStandardizedPlane = false;
        % Test for dimensions
        if ~(numel(varargin{1})==3 && length(varargin{1}) == 3 && numel(varargin{2})==3 && length(varargin{2}) == 3)
            error('Normal and point must be a 3x1 or 1x3 vector')
        end
        if size(varargin{1},1) == 3
            varargin{1} = varargin{1}';
        end
        if size(varargin{2},1) == 3
            varargin{2} = varargin{2}';
        end
        
        % set-up the plane
        plane = [varargin{1}; varargin{2}];
    end
        
    % Compute intersection between the plane and each lines
    [p_inter, check] = plane_line_intersect(plane(1,:)', plane(2,:)', p(t_tp(:,1:2:end),:)', p(t_tp(:,2:2:end),:)');
    % Remove every outside crossing 
    p_inter = p_inter(:,check==1);

    % Compute the convex hull in the desired plane
    if size(p_inter,2) <= 3
        % If there is 3 points or less, the convex hull is all of them
        % independetly of ther order
        p2 = p_inter; 
    else
        if useStandardizedPlane
            % If it's a standard plane (align with global reference frame, we
            % can just remove the axes of that plane to compute the convex hull
            p2 = p_inter(~plane(1,:),:);
        else
            % If its a plane not aligned on global axis, we must align point on a global reference frame in
            % order to compute the convex hull

            % Find a matrix of rotation
            Z = plane(1,:)'; % Z is the normal of the plane
            % Create any point which is not collinear to Z
            X = [1 0 0]';
            if dot(X,Z)/(norm(X)*norm(Z)) >0.95
                X = double(~X);
            end
            % Calculate, recalculate and normalize axes
            Y = cross(Z,X);
            X = cross(Y,Z);
            Z = Z/norm(Z); 
            X = X/norm(X);
            Y = Y/norm(Y);
            
            % Create the rotation matrix
            R = [X, Y, Z];
            
            % Rotate each element
            p2 = R' * p_inter;
            
            % Remove the z axis
            p2 = p2(1:2, :);
        end
    end
    
    % Compute convexe hull
    if isempty(p2)
        P = double.empty(3,0);
        return;
    else
        idx = convhull(p2(1,:), p2(2,:));
    end

    % Returns points
    P = p_inter(:,idx);
    
end