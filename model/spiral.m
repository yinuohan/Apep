function [im, theta] = spiral(skeleton, gif, windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim)

% Plotting parameters
n_p = 5000;                      % Number of points per circle
n_c = 5000;                     % Number of circles per period

im_siz = 360;        %360       % Image size (pixels) (square)
im_res = 45;                    % Image resolution (mas/pix)

% Physics parameters
w     = windspeed / 365;        % windspeed [mas/year]->[mas/day]
P     = period * 365;           % period [year]->[day]
omega = deg2rad(little_omega);  % omega
inc   = deg2rad(inclination);   % inclination  
Ohm   = deg2rad(big_omega);     % Omega 
ecc   = eccentricity;           % Eccentricity
pa    = periastron * 365;       % Periastron date (Julian days)
cone  = deg2rad(cone_angle);    % cone angle (deg)
offst = offset * 365;           % time offset [year]->[day]
%offst = t_obs - pa;            % time offset
rnuc = turn_off;
lim = theta_lim/180*pi;         % limit for theta to produce dust
use_density = 0;

% Don't use in general
if omega_lock
    adjust = kepler_solve(pa-245, P, ecc);
    omega = rectify(omega + adjust);
end

% Time vector
t = (0:n_circ*n_c)/n_c*P + offst; 

% Angles from Keplers laws as a function of time
theta = kepler_solve(t, P, ecc); 

% Radius of each dust circle
r2 = w * (n_circ*P - (0:n_circ*n_c)/n_c*P);
% First ring in time has biggest radius

% 3D spiral plume
% NOTE: Spiral is always generated to align with first point
% on the x-axis (theta = 0). The result is then rotated by 
% the specified angles using direction cosine matrices 
% (note: order of rotations is important & must be preserved
% - omega, inc, Omega)

% Generate the coordinates of the initial unit circle
% which is copied around the spiral
chi = (0:n_p-1)/n_p*pi*2; % Angle
ones2 = ones(1,length(chi)); % Just ones

% The circle is parallel to the y-z plane
circ = [(cos(cone/2)*ones2); % x - becomes North
        (sin(cone/2)*cos(chi)); % y - becomes East
        (sin(cone/2)*sin(chi))]; % z - becomes (-) line of sight

% With anisotropic winds, can try and ellipse
% elongation_factor = 1;
% circ = [(cos(cone/2)*ones2); % x
%         (sin(cone/2)*cos(chi)); % y
%         (sin(cone/2)*sin(chi)) * elongation_factor]; % z

% Initialise full array to store coordinates of all points
circfull = zeros(3, (n_p+1)*length(theta));
gen = 1:n_p; % indices of one circle

% Calculate coordinates of each circle on the spiral
for j = 1:length(theta)
    if r2(j) >= rnuc 
        if rectify(theta(j)) >= lim(1)
            if rectify(theta(j)) <= lim(2)
                circj = rotate_z(theta(j)) * circ * r2(j);
                circfull(:, (j-1)*n_p+gen) = circj;
            end
        end
    end
end

% TEST orbit
% Input
% im_res = 1
% e = ecc
% a = 100
% c = e * a;
% b = sqrt(a^2 - c^2);
% 
% circfull = zeros(3, (n_p+1)*length(theta));
% for theta = 1:1:360
%     circfull(1,floor(theta)) = a*cosd(theta);
%     circfull(2,floor(theta)) = b*sind(theta);
% end
% 
% circfull(1,361) = c;
% circfull(2,361) = 0;
% 
% for i = 362:1:362+a-1
%     circfull(1,floor(i)) = i-362;
%     circfull(2,floor(i)) = 0;
% end

% Rotate points by specifed rotations -------------------
circfull = rotate_z(Ohm) * (rotate_x(inc) * (rotate_z(omega) * circfull));
%circfull = circfull' * rotate_z(omega) * rotate_x(inc) * rotate_z(Ohm) ;
%circfull = circfull';

% Density function
if use_density
    dmax_theta = 270/180 * pi; %10;
    density = cos(chi - dmax_theta) + 2;
end

% Generate image
n_total = n_p*length(theta);
im = zeros(im_siz, im_siz);

% TEST line
% circfull = zeros(3, (n_p+1)*length(theta));
% for i = 1:100
%     circfull(1,i) = i;
%     circfull(2,i) = i/2;
% end


% Project 3D spiral to pixel values by looping over all points
for i = 1:n_total
    
    % View along z axis
    imy = fix(circfull(1, i)/im_res + im_siz/2);
    imx = im_siz - fix(circfull(2, i)/im_res + im_siz/2);
    
    %imy = fix(circfull(1, i)/im_res + im_siz/2);
    %imx = fix(circfull(2, i)/im_res + im_siz/2);
    
    % Image y = up = North = coordinate x
    % Image -x = left = East = coordinate y
    
    % Add density
    if (imx > 0) && (imx <= im_siz) && (imy > 0) && (imy <= im_siz)
        if use_density
            dens = density(mod(i-1,n_p)+1);
            im(imy, imx) = im(imy, imx) + dens;
        else
            im(imy, imx) = im(imy, imx) + 1;
        end
    % Coordinate 1 = image y
    % Coordinate 2 = image x
    end
end

% Get rid of centre
im(im_siz/2, im_siz/2) = 0;

% Normalise
if max(max(im)) > 0 && ~gif
    im = im/max(max(im));
end


% ---------- PLOTTING ----------
scale = im_res/1e3;
dim = im_siz/2;
x = linspace(-dim*scale, dim*scale, dim*2);

if ~gif && 1
    
    % DATA
    figure
    imagesc(-x,x,-skeleton)
    axis image
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    set(gca,'YDir','normal')
    
    % MODEL OUTLINE
    hold on
    high_pass_model = imgaussfilt(im - imgaussfilt(im,2),1);
    model_threshold = 0.003;
    high_pass_model(high_pass_model > model_threshold) = 1;
    high_pass_model(high_pass_model <= model_threshold) = 0;
    high_pass_model = imgaussfilt(high_pass_model,1);
    high_pass_model = high_pass_model/max(max(high_pass_model));
    
    h2 = imagesc(-x,x,-high_pass_model);
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    axis image
    colormap gray
    set(gca,'YDir','normal')
    set(gca,'XDir','reverse')
    alpha = (1 - (high_pass_model<0.01)) * 0.9;
    set(h2, 'AlphaData', alpha)
end

% MODEL
if ~gif && 1
    figure
    %h = imagesc(x,x,-min(im,0.15));
    imagesc(-x, x, -min(im,0.15));
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    colormap gray
    axis image
    set(gca,'YDir','normal')
    set(gca,'XDir','reverse')
end


