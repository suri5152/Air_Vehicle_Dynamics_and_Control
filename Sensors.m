function [accel gyro Pressure GPS ] = Sensors(uvw, pqr, uvwdot, phithetasi, Vg, Wind)

u = uvw(1,1);
v = uvw(1,2);
w = uvw(1,3);

p = pqr(1,1);
q = pqr(1,2);
r = pqr(1,3);

udot = uvwdot(1,1);
vdot = uvwdot(1,2);
wdot = uvwdot(1,3);

phi = phithetasi(:,1);
theta = phithetasi(:,2);
si = phithetasi(:,3);

%% Constants
rho = 1.2682;
h = 1000;
g = 9.8;
Va = Vg - Wind;

%% Random Noise
Xnacc = 0.1;
Ynacc = 0.01;
Znacc = 0.002;
Xngyro = 0.001;
Yngyro = 0.004;
Zngyro = 0.006;
Psnoise = 0.01;
Pdnoise = 0.008;
Pnnoise = 0.001;
Penoise = 0.009;

%% Accelerometer
Xaccel = udot + q*w - r*v + g*sin(theta) + Xnacc;
Yaccel = vdot + r*u - p*w - g*cos(theta).*sin(phi) + Ynacc;
Zaccel = wdot + p*v - q*u - g*cos(theta).*cos(phi) + Znacc;

accel = [Xaccel Yaccel Zaccel];

%% Gyroscope
Xgyro = p + Xngyro;
Ygyro = q + Yngyro;
Zgyro = r + Zngyro;

gyro = [Xgyro; Ygyro; Zgyro];

%% Pressure Sensor
Ps = rho*h*g + Psnoise;
Pd = 1/2*rho*Va^2 + Pdnoise;

Pressure = [Ps; Pd];

%% GPS
Pn = Vg*sin(si) + Pnnoise;
Pe = Vg*cos(si) + Penoise;

GPS = [Pn Pe];

