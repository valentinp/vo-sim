%% Calculate error as a function of theta

%Set parameters
map_dim = 20;
max_theta = 2*pi;
theta_num = 40;
pathType = 'c'; %c (circular) or s (straight) or d (diagonal)
mapType = 'r'; %r (random) or u (uniform)
runs = 500;
ss = [0.5 1 1.5];
nf = [1.5];
fv = [90]*pi/180;


%Determine the theta list
theta = 0:max_theta/theta_num:max_theta;

for k = 1:length(ss)
    for m = 1:length(nf)
        for n = 1:length(fv)
    step_size = ss(k);
    noiseFactor = nf(m);
    fov = fv(n);
    
    if mapType == 'u'
        map_type = 0;
    else
        map_type = 1;
        landmarks = map3D(15,1,900);
    end
    %Keep track of the mean square root error and its standard deviation
    rms_error_x = zeros(1, length(theta));
    rms_error_z = zeros(1, length(theta));
    
    rms_error_x_std = zeros(1, length(theta));
    rms_error_z_std = zeros(1, length(theta));
    
    for i = 1:length(theta)
        disp(theta(i))
        errors = zeros(runs, 3);

        for j = 1:runs
            if map_type == 0
                errors(j, 1:3) = vo_model_offline_3axis(map_dim,map_type,theta(i),fov, step_size, noiseFactor, pathType);
            else
                errors(j, 1:3) = vo_model_offline_3axis(map_dim,map_type,theta(i),fov, step_size, noiseFactor, pathType, landmarks);
            end
        end
        rms_error_x(i) = mean(errors(:,1));
        rms_error_x_std(i) = std(errors(:,1));
        rms_error_z(i) = mean(errors(:,3));
        rms_error_z_std(i) = std(errors(:,3));
    end
    
    %Plot and add aesthetics 
    close all;
    f = figure();
    hold on;
    errorbar(theta*(180/pi), rms_error_x, rms_error_x_std, '--bs','LineWidth',1,'MarkerFaceColor','r');
    errorbar(theta*(180/pi), rms_error_z, rms_error_z_std, '--ks','LineWidth',1,'MarkerFaceColor','g');

    legend('X-axis','Z-axis');
    xlabel('Theta (deg)','FontSize', 12);
    ylabel('Average Error Norm (m)','FontSize', 12);
    chart_title = sprintf('Visual Odometry Error vs. Offset Angle \n Map Type: %i  Map Dim: %i Theta Incr: %i Runs Per Theta: %i \n FOV: %.2f Step Size: %.2f Noise Factor: %.2f', ...
        map_type, map_dim, length(theta), runs, round(fov*180/pi), step_size, noiseFactor);
    title(chart_title, 'FontSize', 12);
    xlim([-10 max_theta*180/pi+10]);
    set(gca,'FontSize',14);
    set(gca,'box','on');

    file_name = sprintf('Figures/VO-E-%s-LM-Theta-FOV-%.2f-SS-%.2f-NF-%.2f-RPT-%d-%d',pathType, fov, step_size, noiseFactor,runs, randi(100));
    saveas(f, strcat(file_name,'.png'));
    saveas(f, strcat(file_name,'.fig'));
    hold off;
        end
       
    end
end



 