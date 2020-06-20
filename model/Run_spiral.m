% Run this file to generate a spiral.
% Calls spiral.m which does the computation. 

% ---------- Data image ----------
% Read image
path = 'Where is the data?';
image = fitsread([path  'APEP_VISIR_2018_J8.9_large.fits']);

% Centre and cut out
image = image(71:430,1:360);

% Bad pixels
image(320:325, 46:51) = 0;
image(59:69, 2:10) = 0;

% Clip and normalise
image = max(image,0);
image = image/max(max(image));

% Log transform
image = log(image+0.00001);

% Clip again
range = [-6, -3];
image = min(max(image,range(1)),range(2));

% Shift min value to 0
image = image - min(min(image));

% Normalise again
image = image/max(max(image));

% Blur a bit
image = imgaussfilt(image,0.5);

% Plot original image
if 0
    figure()
    imagesc(-image)
    colormap gray
    set(gca,'YDir','normal')
    axis image
end

% Plot filtered image
if 0
    figure()
    imagesc(skeleton)
    set(gca,'YDir','normal')
    axis image
end

% ---------- Model parameters ----------
% Fixed
omega_lock = false; % fixed by program design
periastron = 0; % fixed by construction
windspeed = 80; % fixed by proper motion

% Shape
eccentricity = 0.7
cone_angle = 125; % fixed by im geometry
period = 125; % fixed by im size

% On and off
delay = 18
% Delay determines the angle past periastron
% of the dust production centre
theta_lim = [-132 + delay, 132 + delay]; % fixed by im angles
%theta_lim = [-180, 180];
turn_off = 0; % fixed by program design
n_circ = 2; % fixed by program design

% Angles
delta = 0
% Delta determines the inclination axis
little_omega = 0 + delta; % fixed by im orientation
inclination = 25
big_omega = -70 - delay - delta; % fixed by im orientation

% Calculate offset, fixed by position angle
% Rotations are centred at focus point (use true anomaly)
%offset = period/2
obs_angle = 278 - 180;
offset_angle = find_offset_angle(obs_angle, big_omega, inclination, little_omega); % True anomaly
offset = eccentric_to_time(true_to_eccentric(offset_angle, eccentricity), period, eccentricity); % Offset time

% But this is the current offset
% Need to convert to beginning of simultion

offset = offset - n_circ*period

% big_omega = 0
% little_omega = 0
% inclination = 0
% eccentricity = 0

% ---------- Model image ----------
make_gif = 0;
save_gif = 0;

% Plot model image
if ~make_gif
    %figure
    [im, theta] = spiral(image,make_gif,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
end

% Plot angle of the rings
if ~make_gif && 0
    figure
    plot(theta)
end

% Animation
if make_gif
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    filename = 'Where do you want to save this?';
    
    scale = 45/1e3;
    dim = 360/2;
    x = linspace(-dim*scale, dim*scale, dim*2);
    
    offset_initial = 0
    increment = 1
    offset_final = ceil(period)

    for offset = offset_initial:increment:offset_final
        model = spiral(image,make_gif,windspeed,period,inclination,big_omega,turn_off,eccentricity,omega_lock,little_omega,periastron,cone_angle,offset,n_circ,theta_lim);
        
        hold off
        limit = 3000;
        imagesc(-x,x,-min(model,limit));
        title("Offset: " + string(offset) + " yrs")
        xlabel("Relative RA ('')")
        ylabel("Relative Dec ('')")
        colormap gray
        axis image
        set(gca,'YDir','normal')
        set(gca,'XDir','reverse')
        set(gcf,'Position',[20 20 580 580])

        hold on
        plot([6,4], [-7.3,-7.3], '-k', 'LineWidth',1)

        text(6,-6.8,"2 arcsec")
        text(6,-6.1,"Offset (yr): " + string(offset))
        pause(0.01)

        if save_gif
            if offset == offset_initial
                %gif(filename, 'nodither')
                gif(filename)
            else
                %gif('nodither', 'frame',gca)
                gif
            end
        end

    end
end


% ARCHIVE

% The following is wrong
peri_point = [1;0;0];
peri_point = rotate_z(big_omega/180*pi) * (rotate_x(inclination/180*pi) * (rotate_z(little_omega/180*pi) * peri_point));
peri_angle = atan2d(peri_point(2),peri_point(1)); % Projected true anomaly
obs_angle = 278 - 180; % Projected true anomaly
E = obs_angle - peri_angle;
offset2 = eccentric_to_time(E, period, eccentricity);
