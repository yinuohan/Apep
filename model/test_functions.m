% A test to make sure that the functions 
% in this directory are consistent

period = 10;
eccentricity = 0.5;
time = 0:0.1:period;

true_anomaly = kepler_solve(time, period, eccentricity)/pi*180
eccentric_anomaly = true_to_eccentric(true_anomaly, eccentricity)
time2 = eccentric_to_time(eccentric_anomaly, period, eccentricity)

figure
plot(time, time2)

figure
plot(time, true_anomaly)
