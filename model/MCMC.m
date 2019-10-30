% An MCMC to fit a geometric model to the ridges of Apep. 

% Set up data
data = flip(sketch2); % After running Make_sketch.m

% Search parameters
steps = 50000;

windspeed = zeros(1,steps);
period = zeros(1,steps);
inclination = zeros(1,steps);
big_omega = zeros(1,steps);
turn_off = zeros(1,steps);
eccentricity = zeros(1,steps);

start_windspeed = 43;
start_period = 100;
start_inclination = 0;%-30
start_big_omega = 20;
start_turn_off = 9.5;
start_eccentricity = 0.4;

s_windspeed = 2;
s_period = 5;
s_inclination = 5;
s_big_omega = 20;
s_turn_off = 0.1;
s_eccentricity = 0.1;

windspeed(1) = start_windspeed;
period(1) = start_period;
inclination(1) = start_inclination;
big_omega(1) = start_big_omega;
turn_off(1) = start_turn_off;
eccentricity(1) = start_eccentricity;

b1_windspeed = start_windspeed - 10*s_windspeed;
b1_period = start_period - 10*s_period;
b1_inclination = start_inclination - 10*s_inclination;
b1_big_omega = start_big_omega - 10*s_big_omega;
b1_turn_off = start_turn_off - 10*s_turn_off;
b1_eccentricity = start_eccentricity - 4*s_eccentricity;

b2_windspeed = 45;
b2_period = 110;
b2_inclination = start_inclination + 10*s_inclination;
b2_big_omega = start_big_omega + 10*s_big_omega;
b2_turn_off = start_turn_off + 10*s_turn_off;
b2_eccentricity = start_eccentricity + 4*s_eccentricity;

model = spiral(start_windspeed,start_period,start_inclination,start_big_omega,start_turn_off,start_eccentricity);
old_chi2 = sum(sum((model-data)^2));

% Initialise
chi2 = zeros(1,steps);
chi2(1) = old_chi2;
acc = zeros(1,steps);
acc(1) = 1;


for i = 2:steps+1
    % New step
    new_windspeed = windspeed(i-1)+normrnd(0,s_windspeed);
    new_period = period(i-1)+normrnd(0,s_period);
    new_inclination = inclination(i-1)+normrnd(0,s_inclination);
    new_big_omega = big_omega(i-1)+normrnd(0,s_big_omega);
    new_turn_off = turn_off(i-1)+normrnd(0,s_turn_off);
    new_eccentricity = eccentricity(i-1)+normrnd(0,s_eccentricity);
    
    % Evaluate r
    if ~(b1_windspeed<new_windspeed && new_windspeed<b2_windspeed) || ~(b1_period<new_period && new_period<b2_period) || ~(b1_inclination<new_inclination && new_inclination<b2_inclination) || ~(b1_big_omega<new_big_omega && new_big_omega<b2_big_omega) || ~(b1_turn_off<new_turn_off && new_turn_off<b2_turn_off) || ~(b1_eccentricity<new_eccentricity && new_eccentricity<b2_eccentricity)
        r = 0;
    else
        model = spiral(new_windspeed,new_period,new_inclination,new_big_omega,new_turn_off,new_eccentricity);
        new_chi2 = sum(sum((model-data).^2));
        r = old_chi2 / new_chi2;
    end
    
    % Obtain U
    U = rand();
    
    % Compare U with r
    if r > U
        windspeed(i) = new_windspeed;
        period(i) = new_period;
        inclination(i) = new_inclination;
        big_omega(i) = new_big_omega;
        turn_off(i) = new_turn_off;
        eccentricity(i) = new_eccentricity;
        
        old_chi2 = new_chi2;
        accept = 1;
    else
        windspeed(i) = windspeed(i-1);
        period(i) = period(i-1);
        inclination(i) = inclination(i-1);
        big_omega(i) = big_omega(i-1);
        turn_off(i) = turn_off(i-1);
        eccentricity(i) = eccentricity(i-1);
        
        accept = 0;
    end

    acc(i) = (acc(i-1)*(i-1) + accept)/i;
    chi2(i) = old_chi2;
end


% Plot
subplot(2,4,1)
histogram(windspeed)
title('windspeed')
subplot(2,4,2)
histogram(period)
title('period')
subplot(2,4,3)
histogram(inclination)
title('inclination')
subplot(2,4,4)
histogram(big_omega)
title('big omega')
subplot(2,4,5)
histogram(turn_off)
title('turn off')
subplot(2,4,6)
histogram(eccentricity)
title('eccentricity')
subplot(2,4,7)
plot(acc)
subplot(2,4,8)
plot(chi2)

% Best solution
disp(windspeed(chi2==min(chi2)))
disp(period(chi2==min(chi2)))
disp(inclination(chi2==min(chi2)))
disp(big_omega(chi2==min(chi2)))
disp(turn_off(chi2==min(chi2)))
disp(eccentricity(chi2==min(chi2)))



