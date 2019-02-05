clear;
clc;
clf;
close all;

n = 300;
dt = .125; % seconds
N = 1:n;
ts = linspace(0, (n-1)*dt,n);
%function to estimate
F(1) = 0;
w(1) = 40; % deg / s

V1 = 20; %deg
V2 = 10; %deg / s
Vm = [V1^2, V1*V2; V1*V2, V2^2];

for k = 2:n
    w(k) = w(k-1) + randn(1)*5;
    F(k) = F(k-1) + w(k)*dt;
end

plot(ts, F, 'LineWidth',2);
hold on;

Xm = [ F + randn(1,n)*V1; w + randn(1,n)*V2];

plot(ts, Xm(1,:), 'r');
title('Angle');

figure();
plot(ts, w, 'LineWidth',2);
hold on;
plot(ts, Xm(2,:), 'r');
title('Angular Velocity');

S1 = 20; %degrees initial error
S2 = 10; %degrees per second initial error

Q = [0, 0; 0, 14]; % process noise
A = [1, dt; 0, 1];
%initial predicted state is the actual state +/- noise
X(:,1) = [F(1) + rand(1)*S1; w(1) + rand(1)*S2];

P = [S1^2, S1*S2; S1*S2, S2^2];
Ps(1) = P(1);
for k = 2:n
    %predict

    Xp = A* X(:,k-1);
    Pp = A * P * A' + Q;
    %calculate
    K = Pp ./ (Pp + Vm);
    X(:,k) = Xp + K*(Xm(:,k) - Xp);

    %update
    P = ([1,0;0,1] - K) * Pp;
    Ps(k) = P(1);
end

M1 = movmean(Xm(1,:),5);
M2 = movmean(Xm(1,:),10);

figure();
plot(ts,X(1,:),'r','LineWidth',1.5);
hold on;
plot(ts,F(1,:), 'LineWidth',2);
plot(ts,Xm(1,:),'b');
title('Linear KF - Linear Actual');
xlabel('time - s');
ylabel('angle, degrees');

plot(ts,M1);
plot(ts,M2);
legend('Filter','Actual','Measured','MM5','MM10');

figure();

for k = 1:n
    Xk_lims(:,k) = [X(1,k)+Ps(k), X(1,k) - Ps(k)];
end
Xlims = [Xk_lims(1,:), flip(Xk_lims(2,:))];
ts_ = [ts, flip(ts)];
fill1 = fill( ts_, Xlims, 'y');
set(fill1, 'facealpha',0.25);
hold on;
plot(ts, F, 'b','LineWidth',2);
mLims = [Xm(1,:)+Vm(1), flip(Xm(1,:)-Vm(1))];
fill2 = fill( ts_, mLims, 'r');
set(fill2, 'facealpha',0.1);


figure();
diff_est = [F, flip(X(1,:))];
diff_mean10 = [F, flip(M2)];
fill1 = fill(ts_, diff_est, 'r');
hold on;
fill2 = fill(ts_, diff_mean10, 'b');
plot(ts, F, 'b', 'LineWidth', 2);
set(fill1, 'facealpha',0.25);
set(fill2, 'facealpha',0.25);
legend('KF', 'Measured Mean 10');