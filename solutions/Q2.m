W = 1;
H = 1;

angle = 0:90;
Cd0 = 0.611;

for i=1:length(angle)
    alpha = angle(i);
    W_pivot = @(z) (1/W^2 + 1./(2*(H-z)*tand(alpha)+sind(alpha)*W).^2).^(-1/2);
    h = H * (1- cosd(alpha));
    A_eff(i) = W*h+ integral(W_pivot,h,H);
%     Cd(i) = cd_pivot(angle(i));
end
Cd = A_eff / (H*W) * Cd0;

figure();
subplot(1,2,1);
plot(angle, A_eff,'linewidth',2);
xlabel('Angle, \alpha [deg]');
ylabel('A_{eff}');
xlim([min(angle), max(angle)]);
grid on

subplot(1,2,2);
plot(angle, Cd,'linewidth',2);
xlabel('Angle, \alpha [deg]');
ylabel('C_d (\alpha)');
xlim([min(angle), max(angle)]);
grid on
