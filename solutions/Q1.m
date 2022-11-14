%% 
% CEE 261C; Homework 7 solution
%
%% Question 1
rho_a = 1.225;
cp_a = 1005;
V_a = 2880;

rho_f = 2300;   cp_f = 750;     k_f = 0.8;      A_f = 90;   T_f = 293.15;
dx = 0.01;      th = 0.1;       dt = 15;
h = 4;

alpha_f = k_f /(cp_f * rho_f);
lambda = alpha_f * dt / dx^2;

timing = 0:dt:dt*4*60*24;

x = 0:dx:th;
t = 0:dt:dt*4*60*24*4;
T = zeros(length(x),length(t));
n = length(x);

T_int = zeros(1,n);

T_int(1) = T_f;
T(:,1) = T_f;

E_v = zeros(1,length(t));
E_v((t/60/60) >= 7 & (t/60/60) <= 19) = -500;
E_v((t/60/60) >= 24+7 & (t/60/60) <= 24+19) = -500;
E_v((t/60/60) >= 48+7 & (t/60/60) <= 48+19) = -500;
E_v((t/60/60) >= 72+7 & (t/60/60) <= 72+19) = -500;

E_int = 250;

% build matrix and boundary condition
beta = h *dx / k_f;

A = zeros(n);
for i=2:n-1
    A(i,i-1) = lambda;
    A(i,i) = 1-2*lambda;
    A(i,i+1) = lambda;
end
A(1,1) = 1-lambda-lambda*beta;  A(1,2) = lambda;
A(n,n) = 1-lambda-lambda*beta;  A(n,n-1) = lambda;

B = zeros(n,1);
% T(t+1) = A * T(t) + B
for i=1:length(t)-1
    E_f(i) = h*6*A_f*(T(1,i)-T_int(i));     % two side, 3 Af
// %     T_int(i+1) = T_int(i) + dt/(rho_a*cp_a*V_a)*(h*2*3*A_f*(T(1,i)-T_int(i))+E_int+E_v(i));
    T_int(i+1) = T_int(i) + dt/(rho_a*cp_a*V_a)*(E_f(i)+E_int+E_v(i));
    B(1) = lambda*beta*T_int(i);
    B(end) = lambda*beta*T_int(i);
    T(:,i+1) = A*T(:,i)+ B;
end

%%
figure();
subplot(1,2,1); hold on
plot(t/3600,T_int,'r','linewidth',2);
plot(t/3600,T(1,:),'b','linewidth',2)
legend('T_{int}','T_{fabric}')
xlabel('Time [hr]');
ylabel('Temperature [K]');
xlim([0 96]);   xticks(0:24:96);
grid on

subplot(1,2,2); hold on
plot(t/3600,[0, E_f],'r','linewidth',2);
plot(t/3600,E_v,'b','linewidth',2)
plot(t/3600,E_int*ones(length(t),1),'g','linewidth',2)
legend('E_f','E_v','E_{int}');
xlabel('Time [hr]');
ylabel('Heat exchange [W]');
xlim([0 96]);   ylim([-600 600]);
xticks(0:24:96);%    xticklabels({'Noon','Midnight','Noon'});
grid on


%%
figure();hold on
plot(x,T(:,1))
plot(x,T(:,1 + 4*60*4)) % 15s * 4 * 60* 4 
plot(x,T(:,1 + 4*60*8)) % 15s * 4 * 60* 8 
plot(x,T(:,1 + 4*60*12)) % 15s * 4 * 60* 12 
plot(x,T(:,1 + 4*60*16)) % 15s * 4 * 60* 16
plot(x,T(:,1 + 4*60*20)) % 15s * 4 * 60* 20
legend('12 pm', '  4 pm', '12 am', '  4 am', '  8 am');

xlabel('x');
ylabel('Temperature [K]');
xlim([x(1) x(end)])
title('Temperature Profiles throughout the fabric');








