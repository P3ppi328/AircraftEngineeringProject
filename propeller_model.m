clear; close all;clc;

rpm = [2285,2486,2693,2956,3284,3561,3872,4075,4291];
thrust = [6.03,7.22,8.32,10.57,12.86,15.27,17.77,20.49,23.27];
Torque = [2.26,2.71,3.13,3.95,4.84,5.81,6.74,7.85,8.94];
max_rpm = 5200;

% rpm = [1750,2474,3148,3745,4367,4878,5328,5718,6088];
% thrust = [0.70,1.46,2.40,3.41,4.67,5.87,7.16,8.25,9.35];
% Torque = [0.2,0.39,0.62,0.88,1.21,1.52,1.84,2.14,2.42];
% max_rpm = 6700;
% P = T*rpm/9.5488



c_thrust = polyfit(rpm,thrust,2);
c_torque = polyfit(rpm,Torque,2);


rpm_range= 1000:max_rpm;
thrust_model = zeros(1,length(rpm_range));
Torque_model = zeros(1,length(rpm_range));
for i = 1:length(rpm_range)
    thrust_model(i)= c_thrust(1)*rpm_range(i)^2 + c_thrust(2)*rpm_range(i) + c_thrust(3) ;
    Torque_model(i) = c_torque(1)*rpm_range(i)^2 + c_torque(2)*rpm_range(i) + c_torque(3) ;
    
end


%% max thrust 
c = 9.5488;
Thrust_max =c_thrust(1)*max_rpm^2 + c_thrust(2)*max_rpm + c_thrust(3);
Torque_max =c_torque(1)*max_rpm^2 + c_torque(2)*max_rpm + c_torque(3);
P_max = max_rpm * Torque_max / c; 

Thrust_max_N = Thrust_max * 9.81;
%% thrust cruise
thrust_required_N = 80; % N
thrust_required = thrust_required_N / 9.81;

a = thrust_model-thrust_required;
minimum = min(abs(a));
index = find(abs(a)==minimum);
rpm_cruise = rpm_range(index);
T = Torque_model(index);
Torque_cruise = c_torque(1)*rpm_cruise^2 + c_torque(2)*rpm_cruise + c_torque(3);
P_cruise = rpm_cruise * Torque_cruise / c;

%% figure

figure(1)
hold on
plot(rpm,thrust,'o','Color', 'b')
plot(rpm,Torque,'o','Color','r')
plot(rpm_range,thrust_model,'b')
plot(rpm_range,Torque_model,'r')
%plot(rpm_max,Thrust_max,'x','color','b')
%plot(rpm_max,Torque_max,'x','color','r')
xline(max_rpm)
%yline(thrust_required)
title('Thrust and torque over RPM Folding Propeller 32x10.6 CCW 2B MC - puller')
ylabel('Thrust [kgf] - Torque [Nm]')
xlabel('RPM')
legend('Thrust measured','Torque measured','Thrust model','Torque model','Max rpm', 'location','northwest')
hold off


