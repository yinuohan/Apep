% Converts eccentric anomaly to true anomaly

function theta = eccentric_to_true(E, eccentricity)

theta = 2*atand(sqrt((1+eccentricity)/(1-eccentricity))*tand(E/2));

end
