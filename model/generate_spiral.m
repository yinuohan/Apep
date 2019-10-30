function high_pass_model = spiral(windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ)

n_p = 500;                     % Number of points per circle
n_c = 500;                     % Number of circles per period

im_siz = 360;        %360      % Image size (pixels) (square)
im_res = 45;                   % Image resolution (mas/pix)

% Parameters
w     = windspeed / 365;    % windspeed [mas/year]->[mas/day]
P     = period * 365;       % period [year]->[day]
omega = deg2rad(little_omega);      % omega
inc   = deg2rad(inclination); % inclination  
Ohm   = deg2rad(big_omega);   % Omega 
ecc   = eccentricity;       % Eccentricity
pa    = periastron;              % Periastron date (Julian days)
cone  = deg2rad(cone_angle);    % cone angle (deg)
offst = offset * 365;            % time offset [year]->[day]
%offst = t_obs - pa;        % time offset
rnuc = turn_off * 365;      % Turn off [year]->[day]

if omega_lock
    adjust = kepler_solve(pa-245, P, ecc);
    omega = rectify(omega + adjust);
end

t = (0:n_circ*n_c-1)/n_c*P - offst; % Time vector

% Angle solutions from Keplers laws as a function of time
theta = kepler_solve(t, P, ecc); 
r2 = w*(1:n_circ*n_c+1)/n_c*P;
orb_theta = theta(1);

% 3D spiral plume
% NOTE: Spiral is always generated to align with first point
% on the x-axis (theta = 0). The result is then rotated by 
% the specified angles using direction cosine matrices 
% (note: order of rotations is important & must be preserved
% - omega, inc, Omega)

% Generate initial unit circle which is copied around the spiral

chi = (0:n_p)/n_p*pi*2;
ones2 = ones(1,length(chi));

circ = [(cos(cone/2.)*ones2)', (sin(cone/2.)*cos(chi))', (sin(cone/2.)*sin(chi))'];

% Initialise full array to store all points (decreases process time)
circfull = zeros((n_p+1)*length(theta), 3);
gen = 0:n_p;

nr = 0;

% Calculate each circle around the spiral

for j = 1:length(theta)
    circj = circ*r2(j)' * rotate_z(theta(j));
    if r2(j) <= rnuc
        continue
    else
        circfull(j*(n_p+1)+gen, :) = circj;
        nr = nr + 1;
    end
end


% Rotate points by specifed rotations -------------------
circfull = circfull * rotate_z(-omega) * rotate_y(-inc) * rotate_z(-Ohm);

% Density function

dmax_theta = 10; %*pi/180
q = 3;
rin = rnuc;

density = cos(chi - dmax_theta); 
%dmax_theta * (rin/r2) ** q%(chi - dmax_theta) ** r

% Generate image
n_total = (n_p+1)*length(theta);
im = zeros(im_siz, im_siz);

imxs = zeros(1,n_total);
imys = zeros(1,n_total);


% Convert 3D coords to pixel coords
for i = 1:n_total
    imx = floor(fix(circfull(i,1)/im_res + im_siz/2));
    imy = floor(fix(circfull(i,2)/im_res + im_siz/2));
    
    imxs(i) = imx;
    imys(i) = imy;
    
    dens = density(mod(i,n_p+1)+1);
    
    if (imx >= 0) && (imx <= im_siz) && (imy >= 0) && (imy <= im_siz)
        im(imx, imy) = im(imx, imy) + 1; % Add density value at that location
    end
end

im(floor(im_siz/2), floor(im_siz/2)) = 0;

if max(max(im)) > 0
    im = im/max(max(im));
end

im = flip(im);

scale = 45/1e3;
dim = 180;

x = linspace(-dim*scale, dim*scale, dim*2);

% DATA
if 0
    figure()
    path = 'your path';
    image = fitsread([path  'IMAGE_141_J8.9.fits']);
    %image = flip(image);
    image = image/max(max(image));
    image = min(max(image,0),0.08);
    image = image/max(max(image));
    imagesc(x,x,image);
    axis image



    % MODEL OUTLINE
    high_pass_model = imgaussfilt(im - imgaussfilt(im,2),1);
    model_threshold = 0.005;
    high_pass_model(high_pass_model > model_threshold) = 1;
    high_pass_model(high_pass_model <= model_threshold) = 0;
    high_pass_model = imgaussfilt(high_pass_model,1);
    high_pass_model = high_pass_model/max(max(high_pass_model));
    

    hold on
    h2 = imagesc(x,x,high_pass_model);
    xlabel("Relative RA ('')")
    ylabel("Relative Dec ('')")
    axis image
    set(gca,'YDir','normal')
    alpha = (1 - (high_pass_model<0.01)) * 0.9;
    set(h2, 'AlphaData', alpha)
end

% MODEL
figure
h = imagesc(x,x,min(im,0.12));
xlabel("Relative RA ('')")
ylabel("Relative Dec ('')")
colormap hot
set(h, 'AlphaData', 1)
axis image
set(gca,'YDir','normal')

