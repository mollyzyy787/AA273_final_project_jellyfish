close all; clear all; clc;
%%
load filtered_volume_history.mat;

%%
[numRows, numCols]=size(img);
calibration_factor = 0.2/numRows; %meter/pixel
Cd  = 0.30;
dt = 1/30;

Volume_history = zeros(1,length(height_history));
h_hist = zeros(1,length(height_history));
r_hist = zeros(1,length(height_history));
% jelly_dx = zeros(1,length(height_history)-1);
% jelly_dy = zeros(1,length(height_history)-1);
u_hist = zeros(1,length(height_history));
px = zeros(1,length(height_history));
py = zeros(1,length(height_history));

for i = 1:length(height_history)
    h = height_history(i)*calibration_factor;
    r = width_history(i)/2*calibration_factor;
    h_hist(i) = h;
    r_hist(i) = r;
    Volume_history(i) = 2/3*pi*h*r^2;
    px(i) = location(1,i)*calibration_factor;
    py(i) = location(2,i)*calibration_factor;
end


for i = 2:length(height_history)
    u_hist(i) = (filtered_volume(i)-filtered_volume(i-1))/dt;
end

%%
figure;
plot(filtered_volume)

figure;
plot(u_hist)

%% EKF
f = @(x,u,h,r) [x(1) + dt*x(3);
                x(2) + dt*x(4);
                x(3) + dt*(2*u^2/(pi*r^2)-0.5*pi*r^2*Cd*x(3)^2)/((1+(h/r)^1.4)*x(5));
                x(4) - dt*(0.5*1.5*h*r*Cd*x(4)^2)/((1+(h/r)^1.4)*x(5));
                x(5) + dt*u];
x0 = [px(1);py(1);0;0.01;Volume_history(1)];
x_pred = zeros(5,length(height_history));
x_pred(:,1) = x0;
x_update = x_pred;
cov_pred = zeros(5,5,length(height_history));
cov_pred(:,:,1) = eye(5);
cov_update = cov_pred;
Q = 0.01*eye(5);
R = eye(2);
for i = 2:length(height_history)
    alpha = (h_hist(i)/r_hist(i))^1.4;
    Sx = pi*r_hist(i)^2;
    Sy = 1.5*h_hist(i)*r_hist(i);
    %predict
    x_pred(:,i) = f(x_update(:,i-1),u_hist(i-1),h_hist(i-1),r_hist(i-1));
    %Jacobian
    A = [1,0,dt,0,0;
         0,1,0,dt,0;
         0,0,1-dt*Sx*Cd*x_pred(3,i)/((1+alpha)*x_pred(5,i)),0, -dt*(2*u_hist(i)^2/Sx-0.5*Sx*Cd*x_pred(3,i)^2)/((1+alpha)*x_pred(5,i)^2);
         0,0,0,1-dt*Sy*Cd*x_pred(4,i)/((1+alpha)*x_pred(5,i)), dt*0.5*Sy*Cd*x_pred(4,i)^2/((1+alpha)*x_pred(5,i)^2);
         0,0,0,0,1];
    C = [1,0,0,0,0;
         0,1,0,0,0];
    cov_pred(:,:,i) = A*cov_update(:,:,i-1)*A' + Q;
    %update
    innov = [px(i);py(i)]-[x_pred(1,i);x_pred(2,i)];
    Kt = cov_pred(:,:,i)*C'*inv(C*cov_pred(:,:,i)*C' + R);
    x_update(:,i) = x_pred(:,i) + Kt*innov;
    cov_update(:,:,i) = (eye(5) - Kt*C)*cov_pred(:,:,i);   
end

%% predictor
% initialize the predicted values with the EKF updates
x_predictor5 = x_update;
% we are going to predict n_pred number of steps
n_pred = 5;
% we go from the initial step  to N steps before the end
for i = 1:length(height_history) - n_pred
    step = x_update(:,i);
    for j = 1:n_pred
        step = f(step,u_hist(i+n_pred-1),h_hist(i+n_pred-1),r_hist(i+n_pred-1));
    end
    x_predictor5(:,i+n_pred) = step;
end
%%
figure()
plot(px,-py,x_update(1,:),-x_update(2,:),'--')
legend('measurement','EKF estimate')
xlabel('x position')
ylabel('y position')

%%
figure()
plot(px,-py,x_update(1,:),-x_update(2,:),x_predictor5(1,:),-x_predictor5(2,:),'--')
xlabel('x position')
ylabel('y position')
legend('measurement','EKF estimate','5 step predictor')



  