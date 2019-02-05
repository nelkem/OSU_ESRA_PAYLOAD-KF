% extended kalman filter - jan 20
clear;
clc;
clf;
close all;

ME = 57; %uT earths magnetic field strength
n = 1000; %steps
dt = 0.125; %seconds - sampling period
ts = linspace(0,n*dt,n); %time for 100 samples

%initial values
mx = ME;
my = 0;
mxdot = 0;
mydot = 0;
theta = 0;
thetadot = 0;
Fx = zeros(1, n);
Fy = zeros(1, n);
w = zeros(1, n);


Fx(1) = ME;
Fy(1) = 0;
w(1) = 30; %degrees per second initial
% sensor errors
Sm1 = 10; %uT error
Sm2 = 1; %deg/s error

I = [   1, 0, 0, 0, 0, 0;...
        0, 1, 0, 0, 0, 0;...
        0, 0, 1, 0, 0, 0;...
        0, 0, 0, 1, 0, 0;...
        0, 0, 0, 0, 1, 0;...
        0, 0, 0, 0, 0, 1];

    
R = [   Sm1^2, 0, 0, 0, 0, 0;...
        0, Sm1^2, 0, 0, 0, 0;...
        0, 0, Sm1/3, 0, 0, 0;...
        0, 0, 0, 20, 0, 0;...
        0, 0, 0, 0, 20, 0;...
        0, 0, 0, 0, 0, Sm2^2];

q = [.6;.6; .1; .7;.7;.4];
Q = (q*q');
Q = I .* Q;
    
%{
Q = [  1, 0, 0, 0, 0, 0;...
        0, 1, 0, 0, 0, 0;...
        0, 0, 1, 0, 0, 0;...
        0, 0, 0, 10, 0, 0;...
        0, 0, 0, 0, 10, 0;...
        0, 0, 0, 0, 0, 5];
%}
alpha = zeros(1,n);
thetas = zeros(1,n);
%create function to follow
for k = 2:n
    alpha(k) = 30*sind(15*ts(k))*randn(1);
    w(k) = w(k-1) + alpha(k)*dt; %deg/s process variance
    thetas(k) = thetas(k-1) + w(k)*dt + .5*alpha(k)*dt^2;
    Fx(k) = ME*cosd(thetas(k));
    Fy(k) = ME*sind(thetas(k));
end

for k = 1:n
    Xm(1,k) = Fx(k) + randn(1) * Sm1;
    Xm(2,k) = Fy(k) + randn(1) * Sm1;
    
    if( Xm(1,k)< 0 && Xm(2,k) > 0 )
        Xm(3,k) = atan(Xm(2,k)/Xm(1,k))*180/pi + 180;
    else if (Xm(1,k)<0 && Xm(2,k) < 0)
            Xm(3,k) = atan(Xm(2,k)/Xm(1,k))*180/pi - 180;
        else
            Xm(3,k) = atan(Xm(2,k) / Xm(1,k))*180/pi;
        end
    end
    
    Xm(6,k) = w(k) + randn(1) * Sm2;
    Xm(4,k) = -sqrt(Xm(1,k)^2 + Xm(2,k)^2) * sind(Xm(3,k))*Xm(6,k)*pi/180;
    Xm(5,k) = sqrt(Xm(1,k)^2 + Xm(2,k)^2) * cosd(Xm(3,k))*Xm(6,k)*pi/180;
end

%initial state vector
X(:,1) = Xm(:,1);

%{
physical equations

rho = sqrt(mx^2 + my^2);
theta = arctan(my / mx);
rhodot = (x*xdot + y*ydot) / (x^2 + y^2);

%}
%still assuming constant speed

%setup:
A = [ 1, 0, 0, dt, 0, 0;...
      0, 1, 0, 0, dt, 0;...
      0, 0, 1, 0, 0, dt;...
      0, 0, 0, 1, 0, 0;...
      0, 0, 0, 0, 1, 0;...
      0, 0, 0, 0, 0, 1];
  
px = 10; %uT initial error
py = 10; %uT initial error
ptheta = 10; %degrees initial error
pxdot = 10; %uT/s initial error
pydot = 10; %uT/s initial error
pthetadot = 5; %deg/s initial error

P0 = [ px; py; ptheta; pxdot; pydot; pthetadot];
 
P = P0 * P0';

L = n;
hs = zeros(6,20,L);
Xp = zeros(6,1);
for k = 2:L
    % setup
    % X [ x, y, theta, vx, vy, vtheta ]
    % linear prediction
    xp = X(1, k-1) + X(4, k-1) * dt;
    yp = X(2, k-1) + X(5, k-1) * dt;
    thetap = X(3, k-1) + X(6, k-1) * dt;
    
    vxp = X(4, k-1);
    vyp = X(5, k-1);
    vthetap = X(6, k-1);

    Xp(:) = [xp; yp; thetap; vxp; vyp; vthetap];

    
    Pp = A*P*A' + Q;
    Ppthetas(k) = sum(Pp(3,:));
    
    % set up h, Hj
    c = 0;
    c = double(( xp^2 + yp^2 )^(1/2));
    cs(1,k) = c;
    d = (xp*vxp + yp*vyp);
    its = 1;
	dc = 1.0;
    h = zeros(6,1);
    dvx = 1;
    dvy = 1;
    hs(:,1,k) = Xp(:);
    while( max(dvx,dvy) > 0.00001 && its < 20)
        cs(its,k) = c;
        its = its + 1;
        h(:) = [   c * cosd(thetap);...
                c * sind(thetap);...
                thetap;...
                (xp*vxp + yp*vyp)*cosd(thetap)/c - c*sind(thetap)*vthetap*pi/180;...
                (xp*vxp + yp*vyp)*sind(thetap)/c + c*cosd(thetap)*vthetap*pi/180;...
                vthetap];

        dvx = abs(vxp - h(4));
        dvy = abs(vyp - h(5));
        xp = h(1);
        yp = h(2);
        vxp = h(4);
        vyp = h(5);
        hs(:,its,k) = h;
        dc = abs( sqrt( h(1)^2 + h(2)^2) - c );
        c = sqrt( xp^2 + yp^2 );
        cs(its,k) = c;
        d = (xp*vxp + yp*vyp);
    end
    
    
    Hj(4,:) = [ (vxp/c - vxp*xp^2/c^3 - xp*yp*vyp/c^3)*cosd(thetap) - xp * sind(thetap) * pi * vthetap / (c*180),...
                (vyp/c - vyp*yp^2/c^3 - yp*xp*vxp/c^3)*cosd(thetap) - yp * sind(thetap) * pi * vthetap / (c*180),...
                -d*sind(thetap)*pi/(c*180) - c*cosd(thetap)*vthetap*pi^2/180^2,...
                xp/c * cosd(thetap),...
                yp/c * cosd(thetap),...
                -c*sind(thetap)*pi/180];
            
    Hj(5,:) = [ (vxp/c - vxp*xp^2/c^3 - xp*yp*vyp/c^3)*sind(thetap) + xp * cosd(thetap) * pi * vthetap / (c*180),...
                (vyp/c - d*yp/c^3)*sind(thetap) + yp * cosd(thetap) * pi * vthetap / (c*180),...
                d*cosd(thetap)*pi/(c*180) - c*sind(thetap)*vthetap*pi^2/180^2,...
                xp/c * sind(thetap),...
                yp/c * sind(thetap),...
                c*cosd(thetap)*pi/180];
        
        
    Hj(1:2,:) = [  xp*cosd(thetap)/c, yp*cosd(thetap)/c, -c*sind(thetap)*pi/180, 0, 0, 0;...
                   xp*sind(thetap)/c, yp*sind(thetap)/c, c*cosd(thetap)*pi/180, 0, 0, 0];
    Hj(3,:) = [0, 0, 1, 0, 0, 0];
    Hj(6,:) = [0, 0, 0, 0, 0, 1];
    
    Xm_ = Xm(:,k);
    if( Xm(3,k) * h(3) < 0 && abs(h(3)) > 90 )
        Xm_(3) = Xm(3,k) + (h(3)/abs(h(3)))*360;
    end
    Y = Xm_ - h;
    
    % calculate
 
    P = Hj*(Pp*Hj') + R;
    K = (Pp*Hj') * P^-1;
                
    % update set

    X(:,k) = h + K*Y;
    P = (I - K*Hj)*Pp;
    Ps(:,:,k) = P;
    
    if abs(X(3,k)) > 180
        X(3,k) = X(3,k) - X(3,k) * 360 / abs(X(3,k));
    end
    
    Ks(:,:,k) = K;
    Xps(:,k) = Xp;
    Hjs(:,:,k) = Hj;
end




figure();
plot(ts(1:L), X(1,1:L), 'g', 'LineWidth', 2);
hold on;
plot (ts(1:L), Fx(1:L), 'LineWidth',2);
plot( ts(1:L), Xm(1, L));

diffs = Xm(:,1:L) - X;

figure();
title('Xy vs My vs Fy');
plot( ts, Fy, 'b', 'LineWidth',2);
hold on;
plot( ts, X(2,:), 'r','LineWidth',2);
plot( ts, Xm(2,:), 'y', 'LineWidth',2);
plot( ts, movmean(Xm(2,:),10), 'g', 'LineWidth',2);

angles = zeros(1,L);
angles(1) = 0;
for k = 2:L
    angles(k) = atan(Fy(k)/Fx(k))*180/pi;
    if Fx(k) < 0 && Fy(k) > 0
        angles(k) = 180 + angles(k);
    else if Fx(k) < 0 && Fy(k) < 0
        angles(k) = angles(k) - 180;
        end
    end
end

for k = 1:n
    if( thetas(k) > 180 )
        thetas(k:end) = thetas(k:end) - 360;
    else if thetas(k) < -180
            thetas(k:end) = thetas(k:end) + 360;
        end
    end
end


figure();
title('angle');
plot( ts, angles,'b','LineWidth',2);
hold on;
plot( ts, X(3,:),'g','LineWidth',2 );
plot( ts, Xm(3,:),'r', 'LineWidth',1);
ts_ = [ts, flip(ts)];
Pps(:,:) = Ps(3,3,:);
diffs = [X(3,:),flip(thetas)];
fill1 = fill(ts_,diffs,'g');
set(fill1,'facealpha',0.2);
title('Actual angle, Estimated Angle, Measured Angle');
xlabel('time (s, 0.125s sample freq)');
ylabel('Angle (degrees)');
legend('Actual', 'Estimation', 'Measured');

avg4 = movmean(abs( Xm(3,:) ),4);
avg6 = movmean(abs( Xm(3,:) ),6);
avg10 = movmean(abs( Xm(3,:) ),10);

stdX = 0.0;
stdM = 0.0;
std4 = 0.0;
std6 = 0.0;
std10 = 0.0;

for k = 10:n
    stdX = stdX + (abs( X(3,k) ) - abs( thetas(k) ))^2;
    stdM = stdM + (abs( Xm(3,k) ) - abs( thetas(k) ))^2;
    std4 = std4 + (abs( avg4(k)) - abs( thetas(k) ))^2;
    std6 = std6 + (abs( avg6(k)) - abs( thetas(k) ))^2;
    std10 = std10 + (abs( avg10(k)) - abs( thetas(k) ))^2;
end

stdX = sqrt( stdX / (n-10) );
stdM = sqrt( stdM / (n-10) );
std4 = sqrt( std4 / (n-10) );
std6 = sqrt( std6 / (n-10) );
std10 = sqrt( std10 / (n-10) );

fprintf('Estimation STD: %3.3f deg\n',stdX);
fprintf('Measurement STD: %3.3f deg\n',stdM);
fprintf('Moving average 4 STD: %3.3f deg\n',std4);
fprintf('Moving average 6 STD: %3.3f deg\n',std6);
fprintf('Moving average 10 STD: %3.3f deg\n',std10);
plot( ts, avg4);
plot( ts, avg6);
plot( ts, avg10);

figure();
plot( ts, abs( avg10 ), 'r' );
hold on;
plot( ts, abs( thetas ), 'b' );
plot( ts, abs( X(3,:) ), 'g' );

errorsX = abs( X(3,:) ) - abs( thetas );
errorsM = abs( Xm(3,:) ) - abs( thetas );

figure();
subplot(2,1,1);
histogram( errorsX, 40 );
title('Estimation');
subplot(2,1,2);
histogram( errorsM, 40 );
title('Measurement');

