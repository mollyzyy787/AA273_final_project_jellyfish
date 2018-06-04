vel = 1:1:100;
vel = vel';

%% filter in place
figure;
plot(1:100,vel)
title('before running average')

for i = 10:length(vel)
    vel(i) = mean(vel(i-9:i,1));
end
figure;
plot(1:100,vel)
title('after running average')

%% post process filter
figure;
plot(1:100,vel)
title('before running average')

filtered_vel = vel;
for i = 10:length(vel)
    filtered_vel(i) = mean(vel(i-9:i,1));
end
figure;
plot(1:100,vel)
title('after running average')