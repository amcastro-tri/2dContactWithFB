%close all;

h = tt(2)-tt(1);

skip = floor(0.01/h);
scale = 0.5;
%wsz = [-0.6 -0.2 0 0.4];
wsz = [-0.1 0.6 -0.002 0.4];
start_step = 1;
end_step = size(tt);
%end_step = 1;
%skip=1;

p_BoC = params.geometry();
X       = p_BoC(1, :);
Y       = p_BoC(2, :);
nc_max = size(p_BoC, 2);
 
 hSquare = fill(X,Y,'r');

 axis(wsz)
 axis equal;
 
 V0 = get(hSquare,'Vertices')';
 
 nsteps = length(xx);
 
 for istep=start_step:skip:end_step  
     x = xx(istep, 1);
     z = xx(istep, 2);
     theta = xx(istep, 3);
     
     % Rotation matrix:
     c = cos(theta);
     s = sin(theta);
     R_WB = [c, -s; s, c];
     
     C       = repmat([x z], nc_max, 1)';
     
     V = R_WB * V0 + C;             % do the rotation relative to the centre of the square
     
     fill(V(1,:),V(2,:),'b', 'LineWidth', 4, 'EdgeColor','red');
     hold on     
     
     quiver(V(1,:), V(2,:), ft(istep,:), fn(istep,:), 'AutoScale','off')   

      axis equal;
      axis(wsz)
      hold off;
     
     pause(0.1);
 end