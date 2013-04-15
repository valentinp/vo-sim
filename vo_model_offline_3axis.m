function [RMS_error] = vo_model_offline_3axis(map_dim, map_type, theta, fov, dz, noiseFactor, pathType, inputLandmarks)
%Revised VO model

%Define the path
    function [ro_orient] = circle_path (step)
            x_offset = -2;
            z_offset = 0;
            r = 2;
            phi = step*dz/r;
            
            ro_orient = [
            x_offset + r*cos(phi);
            0;
            z_offset + r*sin(phi);
            0;
            -phi;
            0
            ];
    end

    function [ro_orient] = straight_path (step)
            
            ro_orient = [
            0;
            0;
            dz*step;
            0;
            0;
            0
            ];
    end
    function [ro_orient] = straight_path_diag (step)
            
            ro_orient = [
             dz*step/sqrt(2);
            0;
            dz*step/sqrt(2);
            0;
            0;
            0
            ];
    end
if pathType == 's'
    path_fn = @straight_path;
elseif pathType == 'c'
    path_fn = @circle_path;
elseif pathType == 'd'
    path_fn = @straight_path_diag;
end
    
%% Create Map
if map_type == 0
    landmarks = map3D(map_dim, map_type, 100);
else
    if (nargin == 8)
        landmarks = inputLandmarks;
    else
        error('If map_type is set to 1, landmarks input expected as the 8th parameter.');
    end
end
% Rover Position and Orientation: [X Y Z Theta-X Theta-Y Theta-Z]
rover_pos = path_fn(0);



%% Move the Rover
no_steps = 1;


%Initialize
rover(landmarks, rover_pos, fov, true, noiseFactor, theta);

%Main Loop
R = [];
T = [];
Global_R = eye(3);
rover_est_pos = path_fn(0);
rover_est_pos = rover_est_pos(1:3);

rover_est_pos_track=zeros(3,no_steps);
rover_est_pos_error=zeros(no_steps,3);

    for step=1:no_steps
        rover_est_pos_track(1:3, step) = rover_est_pos;
        rover_pos = path_fn(step);
       
        [R, T] = rover(landmarks, rover_pos,fov, false, noiseFactor, theta);
        %Propagate the motion estimation
       if norm(Global_R'*T) < map_dim 
            rover_est_pos = Global_R'*T + rover_est_pos;
            Global_R = R*Global_R;
        end
  
        rover_est_pos_error(step, 1:3) = rover_est_pos - rover_pos(1:3);
    end
    
    RMS_error = [sqrt(sum(rover_est_pos_error(:,1).^2)/length(rover_est_pos_error(:,1)));
        sqrt(sum(rover_est_pos_error(:,2).^2)/length(rover_est_pos_error(:,2)));
        sqrt(sum(rover_est_pos_error(:,3).^2)/length(rover_est_pos_error(:,3)))
        ]; 
end

