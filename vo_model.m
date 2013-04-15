function [] = vo_model(map_dim, map_type, theta, fov, dz, noiseFactor, pathType)
%Revised VO model

%Helper functions
%----------------------------
%Rotation Cosine Matrix 

Ry = @(y) [cos(y) 0 sin(y); 0 1 0;  -sin(y) 0  cos(y)];

%% Create Map
landmarks = map3D(map_dim, map_type, 100);

%% Draw Figures
close all;
figTrue = figure('Position',[50 50 500 500]);
figEst = figure('Position',[550 50 500 500]);

figure(figTrue)
clf(figTrue)
scatter(landmarks(:,1), landmarks(:,3))
axis([-map_dim-2 map_dim+2 -map_dim-2 map_dim+2])
hold on;

%Define the path
    function [ro_orient] = circle_path (step)
            x_offset = 0;
            z_offset = 0;
            r = 3;
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

if pathType == 's'
    path_fn = @straight_path;
elseif pathType == 'c'
    path_fn = @circle_path;
end

% Rover Position and Orientation: [X Y Z Theta-X Theta-Y Theta-Z]
rover_pos = path_fn(0);


%Determine the right and left vectors that define the field of view
fov_vec_0 = Ry(fov/2)*Ry(-theta)*Ry(rover_pos(5))*[0 0 3]';
fov_vec_1 = Ry(-fov)*fov_vec_0;



scatter(rover_pos(1), rover_pos(3), 'g','filled')

figure(figEst)
clf(figEst)
scatter(landmarks(:,1), landmarks(:,3))
axis([-map_dim-2 map_dim+2 -map_dim-2 map_dim+2])
hold on;
scatter(rover_pos(1), rover_pos(3), 'b','filled')


%% Move the Rover
no_steps = 1;


%Initialize
rover(landmarks, rover_pos, fov, true, noiseFactor,theta);

%Main Loop
R = [];
T = [];
Global_R = eye(3);
rover_est_pos = rover_pos(1:3);
%rover_est_pos_track=zeros(no_steps,2);
%rover_est_pos_error=zeros(no_steps,1);

    for step=1:no_steps
        rover_pos = path_fn(step);
        [R, T] = rover(landmarks, rover_pos,fov, false, noiseFactor, theta);
        %Propagate the motion estimation
       if norm(Global_R'*T) < map_dim 
            rover_est_pos = Global_R'*T + rover_est_pos;
            Global_R = R*Global_R;
        end
        
        set(0,'CurrentFigure',figTrue)
        scatter(rover_pos(1), rover_pos(3), 'g','filled')
        % Plot FOV lines
        h = findobj(figTrue, 'type','line');
        delete(h);
        fov_vec_0 = Ry(fov/2)*Ry(-theta)*Ry(rover_pos(5))*[0 0 3]';
        fov_vec_1 = Ry(-fov)*fov_vec_0;
        
        line([0,fov_vec_0(1)] + rover_pos(1),[0,fov_vec_0(3)] +rover_pos(3),'linewidth',1,'color',[1,0,0]);
        line([0,fov_vec_1(1)] + rover_pos(1),[0,fov_vec_1(3)] +rover_pos(3),'linewidth',1,'color',[1,0,0]);
    
        set(0,'CurrentFigure',figEst)
        scatter(rover_est_pos(1), rover_est_pos(3), 'b','filled')
        %rover_est_pos_track(step,:) = rover_est_pos';
        %rover_est_pos_error(step) =  norm(rover_pos(1:2) - rover_est_pos);
        pause(0.1);
    end
    %RMS_error = mean(rover_est_pos_error);
    %display(RMS_error);
end

