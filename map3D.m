function [landmarks] = map3D(map_dim, map_type, max_landmarks)
%Create the global map
%OUTPUT: nx2 array of feature coordinates in the Inertial Frame

%0 : Uniform
%1 : Random
if (map_type == 0)
    PTS_X = [];
    PTS_Z = [];
    for i = -map_dim:map_dim
        PTS_Z = [PTS_Z [-map_dim:map_dim]];
        PTS_X = [PTS_X ones(1,2*map_dim+1)*i];   
    end
    landmarks = [PTS_X; zeros(1, length(PTS_X)); PTS_Z]';
else
    landmarks = -map_dim*ones(max_landmarks,2) + (2*map_dim).*rand(max_landmarks,2);
    landmarks = [landmarks(:,1) zeros(max_landmarks,1) landmarks(:,2)];
end

end

