function [ uvw pqr uvwdot phithetasi ] = SystemModel(u,v,w,Wx,Wy,Wz,p,q,r,theta,phi,si)


%% Constants
m = 13.5;
g = 9.8;
Ix = 0.8244;
Iy = 1.135;
Iz = 1.759;
Izx = 0.1204;
S = 0.55;
b = 2.8956;
c = 0.18994;
Sprop = 0.2027;
rho = 1.2682;
kmotor = 80;
ktp = 0;
komega = 0;
e = 0.9;
Cl0 = 0.28;
Cd0 = 0.03;
Cm0 = -0.02338;
Clalpha = 3.45;
Cdalpha = 0.30;
Cmalpha = -0.38;
Clq = 0;
Cdq = 0;
Cmq = -3.6;
Cldele = -0.36;
Cddele = 0;
Cmdele = -0.5;
Cprop = 1.0;
M = 50;
alpha0 = 0.4712;
Cdp = 0.0437;
Cndelr = -0.032;
Cy0 = 0;
Cl0 = 0;
Cn0 = 0;
Cybeta = -0.98;
Clbeta = -0.12;
Cnbeta = 0.25;
Cyp = 0;
Clp = -0.26;
Cnp = 0.022;
Cyr = 0;
Clr = 0.14;
Cnr = -0.35;
Cydela = 0;
Cldela = 0.08;
Cndela = 0.06;
Cydelr = -0.17;
Cldelr = 0.105;

Clalphadot = 0;
Cw0 = 0;
cdash = c;
theta0 = 0;
u0 = 0;
wdot = 1;


%% Rigid Body Dynamic Equations

Lv = 1/2*rho*u0*b*S*Clbeta;
Lp = 1/4*rho*u0*b^2*S*Clp;
Lr = 1/4*rho*u0*b^2*S*Clr; 
Nv = 1/2*rho*u0*b*S*Cnbeta;
Np = 1/4*rho*u0*b^2*S*Cnp;
Nr = 1/4*rho*u0*b^2*S*Cnr;

ue = u + Wx;
ve = v + Wy;
we = w + Wz;

Rie = [cos(theta)*cos(phi)                              cos(theta)*sin(phi)                                 -sin(theta);
       -cos(phi)*sin(si)+sin(phi)*sin(theta)*cos(si)    cos(phi)*cos(si)+sin(phi)*sin(theta)*sin(si)     sin(phi)*cos(theta);
        sin(phi)*sin(si)+cos(phi)*sin(theta)*cos(si)    -sin(phi)*cos(si)+cos(phi)*sin(theta)*sin(si)     cos(phi)*cos(theta)];

invRie = inv(Rie);
u = invRie(1,1)*ue + invRie(1,2)*ve + invRie(1,3)*we;
v = invRie(2,1)*ue + invRie(2,2)*ve + invRie(2,3)*we;
w = invRie(3,1)*ue + invRie(3,2)*ve + invRie(3,3)*we;

phidot = p + (q*sin(phi) + r*cos(phi))*tan(theta);
thetadot = q*cos(phi) - r*sin(phi);
sidot = (q*sin(phi) + r*cos(phi))*sec(theta);

p = phidot - sidot*sin(theta);
q = thetadot*cos(phi) + sidot*cos(theta)*sin(phi);
r = sidot*cos(theta)*cos(phi) - thetadot*sin(phi);

L = Lv*v + Lp*p + Lr*r;
N = Nv*v + Np*p + Nr*r;

X0 = m*g*sin(theta0);
Y0 = 0;
Z0 = -m*g*cos(theta0);

Czalphadot = -Clalphadot;
Czq = -Clq;
Czalpha = -(Clalpha + Cd0);
Czu = -2*Cl0;
Cxu = -2*Cd0;
Cxalpha = -(Cdalpha - Cl0);

Xu = rho*u0*S*Cw0*sin(theta0) + 1/2*rho*u0*S*Cxu;
Xw = 1/2*rho*u0*S*Cxalpha;
Yv = 1/2*rho*u0*S*Cybeta;
Yp = 1/4*rho*u0*b*S*Cyp;
Yr = 1/4*rho*u0*b*S*Cyr;
Zu = -rho*u0*S*Cw0*cos(theta0) + 1/2*rho*u0*S*Czu;
Zw = 1/2*rho*u0*S*Czalpha;
Zq = 1/4*rho*u0*cdash*S*Czq;
Zwdot = 1/4*rho*cdash*S*Czalphadot;

DX = Xu*u + Xw*w;
DY = Yv*v + Yp*p + Yr*r;
DZ = Zu*u + Zw*w + Zwdot*wdot + Zq*q;

X = X0 + DX;
Y = Y0 + DY;
Z = Z0 + DZ;


uedot = (X - m*g*sin(theta) - m*q*we + m*r*ve)/m;
vedot = (Y + m*g*cos(theta)*sin(phi) - m*r*ue + m*p*we)/m;
wedot = (Z + m*g*cos(theta)*cos(phi) - m*p*ve + m*q*ue)/m;

udot = invRie(1,1)*uedot + invRie(1,2)*vedot + invRie(1,3)*wedot;
vdot = invRie(2,1)*uedot + invRie(2,2)*vedot + invRie(2,3)*wedot;
wdot = invRie(3,1)*uedot + invRie(3,2)*vedot + invRie(3,3)*wedot;



pdot = (-Iz/(Izx^2 - Ix*Iz))*(L-(Iz-Iy)*q*r + Izx*p*q) + (-Izx/(Izx^2 - Ix*Iz))*(N-(Iy-Ix)*p*q - Izx*q*r);

qdot = (M - r*p*(Ix - Iz) - Izx*(p^2 - r^2))/Iy;

rdot = (-Izx/(Izx^2 - Ix*Iz))*(L-(Iz-Iy)*q*r + Izx*p*q) + (-Ix/(Izx^2 - Ix*Iz))*(N-(Iy-Ix)*p*q - Izx*q*r);

% p = phidot - sidot*sin(theta);
% q = thetadot*cos(phi) + sidot*cos(theta)*sin(phi);
% r = sidot*cos(theta)*cos(phi) - thetadot*sin(phi);

xedot = ue*cos(theta)*cos(si) + ve*(sin(phi)*sin(theta)*cos(si) - cos(phi)*sin(si)) + we*(cos(phi)*sin(theta)*cos(si) + sin(phi)*sin(si));
yedot = ue*cos(theta)*sin(si) + ve*(sin(phi)*sin(theta)*sin(si) + cos(phi)*cos(si)) + we*(cos(phi)*sin(theta)*sin(si) - sin(phi)*cos(si));
zedot = -ue*sin(theta) + ve*sin(phi)*cos(theta) + we*cos(phi)*cos(theta);

[t phi] = ode45(@(t,pd) phidot, [0 1],[phidot]);
[t theta] = ode45(@(t,qd) thetadot, [0 1],[thetadot]);
[t si] = ode45(@(t,rd) sidot, [0 1],[sidot]);

uvw = [u v w];
pqr = [p q r];
uvwdot = [udot vdot wdot];

phithetasi = [phi theta si];



end


