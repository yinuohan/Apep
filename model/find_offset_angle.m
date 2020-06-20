% Given the observed (projected) position angle and orbital parameters, 
% find the unprojected offset angle from periastron

function true_anomaly = find_offset_angle(obs_angle, big_omega, inclination, little_omega)

% big_omega = 0;
% little_omega = 0;
% inclination = 0;

% ALL ANGLES ARE TRUE ANOMALY

% Rotate periastron pa and calculated projected angle
%peri_point = [1;0;0];
%peri_obs_point = rotate_z(big_omega/180*pi) * (rotate_x(inclination/180*pi) * (rotate_z(little_omega/180*pi) * peri_point));
%peri_obs_angle = atan2d(peri_point(2),peri_point(1)); % True anomaly

% Rotate test pa and calculated projected angle
% Compare to observed
%obs_angle = 278 - 180; % True anomaly

% Get approximate range
test_true_anomalies = 0:0.1:360; % deg
tolerance = 0.2; % deg

found_point = 0;
if inclination >= 85
    disp("High inclination -- harder to find true anomaly!")
end

for test = test_true_anomalies
    test_point = [cosd(test); sind(test); 0];
    test_obs_point = rotate_z(big_omega/180*pi) * (rotate_x(inclination/180*pi) * (rotate_z(little_omega/180*pi) * test_point));
    test_obs_angle = atan2d(test_obs_point(2),test_obs_point(1));
    
    if abs(test_obs_angle - obs_angle) < tolerance
        true_anomaly = test;
        found_point = found_point + 1;
        break
    end
end

% Get finer estimate
%delta = 180 * inclination*3 / 90*3; % High inclination needs larger range
delta = 1;
test_true_anomalies = true_anomaly-delta:1e-4:true_anomaly+delta; % deg
tolerance = 2e-4; % deg

for test = test_true_anomalies
    test_point = [cosd(test); sind(test); 0];
    test_obs_point = rotate_z(big_omega/180*pi) * (rotate_x(inclination/180*pi) * (rotate_z(little_omega/180*pi) * test_point));
    test_obs_angle = atan2d(test_obs_point(2),test_obs_point(1));
    
    if abs(test_obs_angle - obs_angle) < tolerance
        true_anomaly = test;
        found_point = found_point + 1;
        break
    end
end

if found_point ~= 2
    disp("Didn't find the right offset angle!")
    disp(found_point)
end

end
