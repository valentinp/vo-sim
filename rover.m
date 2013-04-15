function [opt_R, opt_T] = rover(landmarks, true_rover_pos,fov, initialize, noiseFactor, theta)
%rover.m: A model of a stereo camera equipped rover.
%% Set up Parameters
b = 0.25; % Baseline [m]
f = 200;
%NOTE: Pixel dimension is not necessary if f is in units of pixels.
%pixel_dim = 2*f*tan(fov/2)/512; %Dimension of Each Pixel assuming 512px wide camera.
c_u = 0;
c_v = 0;

persistent previous_landmarks;

%Helper functions
%----------------------------
%Rotation Cosine Matrix
%Negative because the y-axis points "toward the ground", thus to preserve
%expected rotations, the angle is negated
Rx = @(x) [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Ry = @(y) [cos(y) 0 sin(y); 0 1 0;  -sin(y) 0  cos(y)];
Rz = @(z) [cos(z) -sin(z) 0;  sin(z) cos(z) 0; 0 0 1];


%Triangulation Functions
%Assumed structure: u = [u_l v_l u_r v_r]';
tri_X = @(u) (b/2)*(u(1) + u(3))/(u(1) - u(3));
tri_Y = @(u) (b/2)*(u(2) + u(4))/(u(1) - u(3));
tri_Z = @(u) (b*f)*1/(u(1) - u(3));
%-----------


%Determine the right and left vectors that define the field of view
%In the camera frame
fov_vec_0 = Ry(fov/2)*[0 0 1]';
fov_vec_1 = Ry(-fov)*fov_vec_0;

%The orientation of the camera with respect to the rover frame and then the
%inertial frame
C_cr = Ry(theta); 
orient = true_rover_pos(4:6);
C_ri = Ry(orient(2)); 

C_ci = C_cr*C_ri;

%Position of the rover in the inertial frame
r_ci = true_rover_pos(1:3);

no_landmarks = length(landmarks(:,1));

%% Track Landmarks
%Projection Matrix (3D -> 2D in 2 Cameras)
%Assumes camera origin is taken to be in the middle of the two cameras
M = [f 0 c_u f*(b/2);
    0 f c_v 0;
    f 0 c_u -f*(b/2);
    0 f c_v 0;
];

current_landmarks = [];

for i=1:no_landmarks
    %Calculate the position vector of the landmark in the camera frame
    r_li = landmarks(i,:)';
    r_lc = C_ci*(r_li- r_ci);
    
    %Account for both FOVs
    r_lc_right = C_ci*(r_li- (r_ci + C_ci*[b/2 0 0]'));
    r_lc_left = C_ci*(r_li- (r_ci + C_ci*[-b/2 0 0]'));
    
    %Is this landmark within both of our cameras' field of view?
    %Recall that the y-axis points downward (hence the negative sign)
    
    if ([0 -1 0]*cross(fov_vec_0,r_lc_right) > 0 && [0 -1 0]*cross(r_lc_right, fov_vec_1) > 0) && ([0 -1 0]*cross(fov_vec_0,r_lc_left) > 0 && [0 -1 0]*cross(r_lc_left, fov_vec_1) > 0)
        %Project the landmark onto the left and right camera 
       u = 1/r_lc(3)*M*[r_lc; 1]; 
       % Add some noise
       noiseVec = noiseFactor.*randn(4,1);
       u_new = u + [noiseVec(1) 0 noiseVec(3) 0]';
       
       if (u_new(1) - u_new(3) < 4)
           %disp('error');
           %disp(u);
           %disp(u_new);
           continue;
       end
       
       
       
       %Add projection and the landmark index to the list of tracked landmarks
       current_landmarks = [current_landmarks; i u_new';];
    end

end

%% Triangulate and Determine Estimated Motion
if initialize
    %Retain information about the landmarks
    %Do not estimate motion
    previous_landmarks = current_landmarks;
    opt_R = 0;
    opt_T = 0;
else
    %FEATURE TRACKING
     %Determine all landmarks that are found in both the current and
     %previous tracked sets.
    C_index = current_landmarks(:,1);
    P_index = previous_landmarks(:,1);
    ci_intersect = intersect(C_index, P_index);
  
        %Get 1D Positions
        p_l = previous_landmarks(ismember(previous_landmarks(:,1),ci_intersect),:);
        c_l = current_landmarks(ismember(current_landmarks(:,1),ci_intersect),:);
       

        %Initialize 3D position arrays
        p_3d_l = zeros(3,length(p_l));
        c_3d_l = zeros(3,length(p_l));
        
        %Initialize covariance arrays
        sigma_c = zeros(3,3,length(p_l));
        sigma_p = zeros(3,3,length(p_l));

        %Triangulates all recorded landmarks into their 3D locations in the
        %previous and current frames. Then calculate the covariance matrix
        %of each.
        
        tracked_landmark_amount = length(c_l(:,1));
        
        %disp(tracked_landmark_amount);
        
       for j=1:tracked_landmark_amount
            p_3d_l(1:3,j) = [tri_X(p_l(j,2:5)) tri_Y(p_l(j,2:5)) tri_Z(p_l(j,2:5)) ]';
            c_3d_l(1:3,j) = [tri_X(c_l(j,2:5)) tri_Y(c_l(j,2:5)) tri_Z(c_l(j,2:5)) ]';
        
        Jc = jacob3D(c_l(j,2:5), b, f);
        sigma_c(1:3,1:3,j) = Jc*Jc';
        Jp = jacob3D(p_l(j,2:5), b, f);
        sigma_p(1:3,1:3,j) = Jp*Jp';

       end

    %Get a motion estimate in the old reference frame
    %If less than 3 landmarks exist, return the trivial transformation
   % if tracked_landmark_amount > 3
        [opt_R,opt_T] = motion_est_matrix(p_3d_l, c_3d_l,sigma_p,sigma_c);
        opt_T = C_cr'*opt_T;
        %opt_R = C_cr'*opt_R;
   % else
   %     opt_R = eye(3);
   %     opt_T = zeros(3,1);
   % end
     
end
    

    %Save these landmarks for the next iteration
    %Note these are all the current landmarks, even those that were not
    %present in the prior set.
    previous_landmarks = current_landmarks;
end

