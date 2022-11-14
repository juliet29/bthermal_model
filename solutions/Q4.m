clear; close all; clc

load('Q3data');

TintRec = TintMeas;
// Data
va      = 2880;                   // Volume of air inside the building [m^3]
rhoa    = 1.225;                  // Density of air [kg/m^3]
ca      = 1005;                   // Specific heat of air [J/KgK]
th      = 0.1;                    // Thickness of thermal mass floors[m]
A_f     = 2*90;                  // Area of the 3 floors [m^2]
rho_f   = 2300;                   // Density of Thermal Mass concrete [kg/m^3]
k_f     = 0.8;                    // Conductivity of Thermal Mass concrete [W/mK]
c_f     = 750;                    // Specific heat of Thermal Mass concrete [J/KgK]
hconv   = 4;                      // Convective heat coefficient [W/m2.K]
Eint    = 250;                   // Internal loads due to occupancy, lighting and equipment [W]

// Implement your own version of Ev as you did in question 1

// Discretization Variables
dt    = 15;
dx = th/10;        
x  = 0:dx:th;           
Nx = length(x);

// Matrix for conduction
alpha  = k_f./(rho_f.*c_f);
lambda = alpha.*dt./dx.^2;
A      = conductionMatrix(lambda, dt, dx, Nx);

// First loop for thermal mass initial temperature 
// Initialize Night Flush Conditions   
Tint    = 273+TintRec(1); 
T0      = 273+TintRec(1);
TNx     = 273+TintRec(1);
Tf      = T0 + (TNx-T0)*x''/th; // Initial Temperature of Thermal Mass -> x transpose...! 

Natural ventilation parameters
g = 9.81;

alpha1 = 42;
alpha2 = 42;
alpha3 = 42;

Cd1 = cd_pivot(alpha1);
Cd2 = cd_pivot(alpha2);
Cd3 = cd_pivot(alpha3);

A1 = 1.61;
A2 = 1.755;
A3 = 1.755;

L1 = 14.41-3.07;
L2 = 14.41-7.79;
L3 = 14.41-12.34;

//// Convertion of data in Kelvin, and cut
time_in_seconds = time*3600;
time_in_hours   = time;
Tout        = 273+Tout;
TintReal    = 273+TintRec;
//The time starts at midnight!!!

//// Recursive Time Loop

for t = 1:length(time)
    // Equivalent hour
    h=time_in_hours(t);
    
    // Natural ventilation
    //Compute Env at each time step
   
    // Q = \dot{V}
    Q1(t) = Cd1 * A1 * sqrt(2*g*L1*(Tint(t)-Tout(t))/Tout(t));
    Q2(t) = Cd2 * A2 * sqrt(2*g*L2*(Tint(t)-Tout(t))/Tout(t));
    Q3(t) = Cd3 * A3 * sqrt(2*g*L3*(Tint(t)-Tout(t))/Tout(t));
    Q_tot(t) = Q1(t)+Q2(t)+Q3(t);
    Env(t) = rhoa*ca*Q_tot(t)*(Tout(t) - Tint(t));
    
    
    //// Convective heat flux
    qconv  = hconv*(Tf(1,t) - Tint(t)); // first floor
    qconvst(t) = qconv;
    //// Walls boundary conditions
    // At x = 0
    b(1)  = -lambda*dx/k_f*(qconv);
    // At x = th
    b(Nx) = -lambda*dx/k_f*(qconv);

    //// Thermal mass temperature
    Tf(:,t+1) = A*Tf(:,t) + b''; // transpose

    //// Air temperature
    Tint(t+1) = Tint(t) + (qconv*6*A_f + Eint +  Env(t))*dt/(va*rhoa*ca); 
    Eist(t) = Eint;
    
end    


//// Plots - NO NEED TO MODIFY THIS
Vnvtot = Q_tot;
figure;
plot(time_in_hours,Vnvtot)
xlabel('Time from midnight (h)');
ylabel('Mass flow rate (m^3.^{-1})');

figure 
plot(time_in_hours, Tint(1:end-1))
hold
plot(time_in_hours, Tf(end,1:end-1))
plot(time_in_hours, Tout)
plot(time_in_hours, TintReal)
legend('T_{int}(t)','T_f(x=th,t)','T_{out}(t)','Real T_{int}(t)');
xlabel('Time from midnight (h)');
ylabel('Temperature (K)');

