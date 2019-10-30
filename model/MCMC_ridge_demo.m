% Demo for feasibility of MCMC on blurred ridge with one free parameter.

% Set up data
dim = 360;
a = zeros(dim,dim);
a(floor(dim/2),:) = 1;
a = imgaussfilt(a,10);

% Search parameters
mu = 0;
sigma = 10;
steps = 10000;
start_x = ceil(rand()*360);
blur_model = 5;

% Store x
x = zeros(1,steps+1);
x(1) = start_x;
acc = zeros(1,steps+1);
acc(1) = 1;
chi2 = zeros(1,steps+1);

% Initialise
b = zeros(dim,dim);
b(start_x,:) = 1;
b = imgaussfilt(b,blur_model);
old_chi2 = sum(sum((a-b).^2));
chi2(1) = old_chi2;

for i = 2:steps+1
    % New step
    random_step = normrnd(mu,sigma);
    new_x = x(i-1) + floor(abs(random_step))*sign(random_step);
    
    % Evaluate r
    if new_x < 1 || new_x > dim
        r = 0;
    else
        b = zeros(dim,dim);
        b(new_x,:) = 1;
        b = imgaussfilt(b,blur_model);
        new_chi2 = sum(sum((a-b).^2));
        r = old_chi2 / new_chi2;
    end
    
    % Obtain U
    U = rand();
    
    % Compare U with r
    if r > U
        x(i) = new_x;
        old_chi2 = new_chi2;
        accept = 1;
    else
        x(i) = x(i-1);
        accept = 0;
    end

    acc(i) = (acc(i-1)*(i-1) + accept)/i;
    chi2(i) = old_chi2;
end


close all

figure()
imagesc(a)

figure()
histogram(x,100)

figure()
plot(x)

figure()
scatter(x,chi2)


