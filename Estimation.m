function [Est] = Estimation(Data)
dt = 1;

t  = 1:dt:41;

Nsamples = length(Data);

Xsaved = zeros(Nsamples, 1);
Zsaved = zeros(Nsamples, 1);

for k=1:Nsamples
  
  z = Data(k,1);
  volt = Kalman(z);
  
  Xsaved(k) = volt;
  Zsaved(k) = z;
end

Est = Xsaved;

end


function volt = Kalman(z)


persistent A H Q R 
persistent x P
persistent firstRun


if isempty(firstRun)
  A = 1;
  H = 1;
  
  Q = 0;
  R = 4;

  x = 14;
  P =  6;
  
  firstRun = 1;  
end

  
xp = A*x;
Pp = A*P*A' + Q;

K = Pp*H'*inv(H*Pp*H' + R);

x = xp + K*(z - H*xp);
P = Pp - K*H*Pp;


volt = x;
end