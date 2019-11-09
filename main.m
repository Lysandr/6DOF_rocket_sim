%% Rocket Powered Vehicle Simulation Environment
% Padraig Lysandrou March 8th, 2019

clc; close all; clear all;
m_dry = 90.7 / 2;
m_wet = 90.7 + m_dry;
Thrust_mag = 4448.2216;

statedim = 13;
dt = 0.1;
T = 200;
t = 0:dt:T;
npoints = length(t);

veh_r = 0.3;            % vehicle radius
veh_h = 1;              % vehicle height

g0 = 9.81;
p.g = [0 0 -g0].';
Ic = diag([(0.5*(m_wet*veh_h*veh_h) + 0.25*(m_wet*veh_r*veh_r)) ...
    (0.5*(m_wet*veh_h*veh_h) + 0.25*(m_wet*veh_r*veh_r)) 0.5*m_wet*veh_r*veh_r]);
p.Ic = Ic;
Isp = 205;
p.mdot = Thrust_mag/(g0 * Isp);
p.T_B = [0 0 Thrust_mag].';
p.r_E_COM = [0 0 -veh_h/2];
p.Cd = 0.75;
p.Aref = pi*(0.1588^2); % 12.5 in diameter

f_dot = @(t_in,state_in,param) dynamics(t_in,state_in,param);
vehicle_state = zeros(13,npoints);

% sigma_0 = [(sind(1/2)/(1+cosd(1/2))) 0 0].';
sigma_0 = [0 0 0].';
vehicle_state(1,1) = m_wet;
vehicle_state(8:10,1) = sigma_0;

tic
for i = 1:npoints-1
    % Distribute the states from the last cycle for computation
    m =         vehicle_state(1,i);
    x_N =       vehicle_state(2:4,i);
    v_N =       vehicle_state(5:7,i);
    sigma_BN =  vehicle_state(8:10,i);
    omega_BN =  vehicle_state(11:13,i);
    
    if m <= m_dry
       p.mdot = 0;
       p.T_B = [0 0 0].';
    end
    
    % Naive, but update the inertia matrix...
    Ic = diag([(0.5*(m*veh_h*veh_h) + 0.25*(m*veh_r*veh_r)) ...
        (0.5*(m*veh_h*veh_h) + 0.25*(m*veh_r*veh_r)) 0.5*m*veh_r*veh_r]);
    p.Ic = Ic;
    
    % Update exogenous torques
    L =  zeros(3,1);  p.L = L;
    
    % RK4 step for the spacecraft dynamics
    k_1 = f_dot(t(i), vehicle_state(:,i), p);
    k_2 = f_dot(t(i)+0.5*dt, vehicle_state(:,i)+0.5*dt*k_1, p);
    k_3 = f_dot((t(i)+0.5*dt),(vehicle_state(:,i)+0.5*dt*k_2), p);
    k_4 = f_dot((t(i)+dt),(vehicle_state(:,i)+k_3*dt), p);
    vehicle_state(:,i+1) = vehicle_state(:,i) + (1/6)*(k_1+(2*k_2)+(2*k_3)+k_4)*dt;
    
    % Perform the nonsingular MRP propagation attitude check
    s = norm(vehicle_state(8:10,i+1));
    if s > 1
        vehicle_state(8:10,i+1) = -(vehicle_state(8:10,i+1) ./(s^2));
    end
end
toc

m =         vehicle_state(1,:);
x_N =       vehicle_state(2:4,:);
v_N =       vehicle_state(5:7,:);
sigma_BN =  vehicle_state(8:10,:);
omega_BN =  vehicle_state(11:13,:);

% figure; plot(t,m); grid on; hold off; title('Mass over time')
% xlabel('Time seconds'); ylabel('Mass, kg')
% 
% figure; plot3(x_N(1,:),x_N(2,:),x_N(3,:));
% grid on; hold off; %axis equal;
% xlabel('x'); ylabel('y'); zlabel('z')
% title('Position vs Time');
% 
% figure;
% subplot(3,1,1); plot(t,x_N(1,:)); ylabel('x pos'); xlabel('time, s');
% subplot(3,1,2); plot(t,x_N(2,:)); ylabel('y pos'); xlabel('time, s');
% subplot(3,1,3); plot(t,x_N(3,:)); ylabel('z pos'); xlabel('time, s');
% title('Position vs Time'); hold off;
% 
% figure;
% subplot(3,1,1); plot(t,v_N(1,:)); ylabel('x vel'); xlabel('time, s');
% subplot(3,1,2); plot(t,v_N(2,:)); ylabel('y vel'); xlabel('time, s');
% subplot(3,1,3); plot(t,v_N(3,:)); ylabel('z vel'); xlabel('time, s');
% title('Velocity vs Time'); hold off;
% 
% figure;
% subplot(3,1,1); plot(t,sigma_BN(1,:)); ylabel('\sigma_1'); xlabel('time, s');
% subplot(3,1,2); plot(t,sigma_BN(2,:)); ylabel('\sigma_2'); xlabel('time, s');
% subplot(3,1,3); plot(t,sigma_BN(3,:)); ylabel('\sigma_3'); xlabel('time, s');
% title('Sigma MRP vs Time'); hold off;
% 
% figure;
% subplot(3,1,1); plot(t,omega_BN(1,:)); ylabel('\omega_1'); xlabel('time, s');
% subplot(3,1,2); plot(t,omega_BN(2,:)); ylabel('\omega_2'); xlabel('time, s');
% subplot(3,1,3); plot(t,omega_BN(3,:)); ylabel('\omega_3'); xlabel('time, s');
% title('Omega vs Time'); hold off;


figure;
subplot(3,1,1); plot(t,x_N(3,:)); ylabel('z pos'); xlabel('time, s');
title('Z position (m) vs Time'); grid on;
subplot(3,1,2); plot(t,v_N(3,:)); ylabel('z vel'); xlabel('time, s');
title('Z velocity (m/s) vs Time'); grid on;
subplot(3,1,3); plot(t,m);
title('Mass vs Time (kg)'); grid on;



apogee_km = max(x_N(3,:))/1000
max_v_mph = max(v_N(3,:))*2.23694
