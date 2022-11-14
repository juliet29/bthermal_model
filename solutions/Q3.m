clear;
load('JulyWeek.mat');
time_hr = time/3600;
TintMeas=TintRec;

g=9.81;

w1.alpha = 42;
w1.area = 1.61;
w1.l = 11.34;
w1.cd = cd_pivot(w1.alpha);

w2.alpha = 42;
w2.area = 1.755;
w2.l = 6.62;
w2.cd = cd_pivot(w2.alpha);

w3.alpha = 42;
w3.area = 1.755;
w3.l = 2.07;
w3.cd = cd_pivot(w3.alpha);


idx = zeros(1,length(time));
idx(time_hr <=7 ) = 1;
for i=1:7
    idx(24*(i-1) + 19 <= time_hr & time_hr <= 24*i+7)=1;
end

V1 = zeros(1,length(time));
V2 = zeros(1,length(time));
V3 = zeros(1,length(time));
for i=1:length(idx)
   if(idx(i) ~= 0)
       V1(i) = w1.cd * w1.area * sqrt(2*g*w1.l *(TintMeas(i)-Tout(i))/(Tout(i)+273));
       V2(i) = w2.cd * w2.area * sqrt(2*g*w2.l *(TintMeas(i)-Tout(i))/(Tout(i)+273));
       V3(i) = w3.cd * w3.area * sqrt(2*g*w3.l *(TintMeas(i)-Tout(i))/(Tout(i)+273));
   end
end

figure();
hold on
plot(time_hr/24, V1 ,'r');
plot(time_hr/24, V2 ,'g');
plot(time_hr/24, V3 ,'b');
plot(time_hr/24, V1+V2+V3 ,'k','linewidth',1.5);
legend('V_1','V_2','V_3','V_{total}');
xlabel('Time [day]');
ylabel('Ventilation rate [m^3/s]');
title('Volume flow rate');


