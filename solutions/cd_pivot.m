function Cd = cd_pivot(alpha)
W=1;
H=1;

W_pivot = @(z) (1/W^2 + 1./(2*(H-z)*tand(alpha)+sind(alpha)*W).^2).^(-1/2);
h = H * (1- cosd(alpha));
A_eff = W*h+ integral(W_pivot,h,H);

Cd0 = 0.611;
Cd = A_eff / (H*W) * Cd0;

end
