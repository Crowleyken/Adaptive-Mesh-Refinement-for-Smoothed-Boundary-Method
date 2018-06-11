function c_out = Ficks(c_in,OP,direc,dx,dy)

% C_in - 2D concentration matrix
% OP - 2D Order Parameter matrix
% direc - r, l, u or d. direction of finite difference scheme
% dx - grid size i direction
% dy - grid size j direction


    % Computes right
    if direc == 0
        % Average of order parameter right
        OP = (circshift(OP,[0 1])+OP)/2;
        c_out = ((circshift(c_in,[0 1]) - c_in)/dx).*OP;
        
    % Compute left
    elseif direc == 1
        % Average of order parameter left
        OP = (circshift(OP,[0 -1])+OP)/2;
        c_out = ((c_in - circshift(c_in,[0 -1]))/dx).*OP;
        
    % Compute up
    elseif direc == 2
        % Average of order parameter up
        OP = (circshift(OP,[-1 0])+OP)/2;
        c_out = ((circshift(c_in,[-1 0]) - c_in)/dy).*OP;
        
    % Compute down
    elseif direc == 3
        % Average of order parameter down
        OP = (circshift(OP,[1 0])+OP)/2;
        c_out = ((c_in - circshift(c_in,[1 0]))/dy).*OP;
    
    else
            disp('must input d,l,u or d')
    end        
end