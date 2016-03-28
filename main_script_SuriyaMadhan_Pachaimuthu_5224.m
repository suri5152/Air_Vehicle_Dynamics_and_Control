clc; clear; close all;

% Parameters that could be changed 

u = 10;
v = 0;
w = 0;
Wx = 1;
Wy = 0;
Wz = 0;
p = 2;
q = 0;
r = 0;
theta = 1.24;
phi = 0.45;
si = 0.24;

Wind = 2;    % To determine air speed
Vg = 30;

%% Rigid Body Model
[uvw pqr uvwdot phithetasi] = SystemModel(u,v,w,Wx,Wy,Wz,p,q,r,theta,phi,si);

%% Sensor Model
[accel gyro Pressure GPS] = Sensors(uvw, pqr, uvwdot, phithetasi, Vg, Wind);

%% Estimation(kalman)
EstaccelX = Estimation(accel(:,1));
EstaccelY = Estimation(accel(:,2));
EstaccelZ = Estimation(accel(:,3));

EstPn = Estimation(GPS(:,1));
EstPe = Estimation(GPS(:,2));

%% Follow Trajectory

Dubin;


