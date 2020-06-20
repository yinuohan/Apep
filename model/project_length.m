% Given projection angles of orbit, 
% caluculate the apparent length of a line in orbital plane

% Determine length and direction of the line
angle = offset_angle;

% Make it into a point without projection
point = [cosd(offset_angle); sind(offset_angle); 0];

% Rotate
new_point = rotate_z(big_omega/180*pi) * (rotate_x(inclination/180*pi) * (rotate_z(little_omega/180*pi) * point));

% Projected angle
new_angle = atan2d(new_point(2),new_point(1)) % Projected true anomaly
length_factor = sqrt(new_point(1)^2 + new_point(2)^2) % Projected length

% Caluculate true length
projected_length = 42
orig_length = projected_length / length_factor

semimajor_axis = orig_length / true_to_radius(offset_angle, 1, eccentricity)

