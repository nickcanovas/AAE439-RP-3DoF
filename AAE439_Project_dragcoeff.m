function CD0_FB = AAE439_Project_dragcoeff(v, rho)
% v = vehicle velocity (m/s)
% rho = air density (kg/m^3)

in_to_m = 0.0254; % converts inches to meters
l_b = 55*in_to_m; % body length (m)
l_n = 11.5*in_to_m; % nose cone length (m)
d_m = 3.1*in_to_m; % max body diameter (m)
d_b = 3.1*in_to_m; % base diameter (m)
t = 0.0025*in_to_m; % fin thickness (m)
c_tip = 1*in_to_m; % tip chord length for primary fins (m)
c_root = 6*in_to_m; % root chord length for primary fins (m)
c_tip_aft = 3.49*in_to_m; % tip chord length for aft fins (m)
c_root_aft = 4.53*in_to_m; % root chord length for aft fins (m)
h_fin = 5*in_to_m; % primary fin height (m)
h_fin_aft = 5*in_to_m; % aft fin height (m)
num_finpairs = 4; % number of primary-aft fin pairs

Re_cr = 500000; % critical reynolds number
mu = 3.178E-5; % dynamic viscosity of air at sea level (kg/m*s)

c = (c_tip + c_root)/2; % average fin chord for primary fins (m)
c_aft = (c_tip_aft + c_root_aft)/2; % average fin chord for aft fins (m)

Re_c = rho*v*c/mu;

S_F = num_finpairs*(c*h_fin + c_aft*h_fin_aft); % total planform area for all fins (m^2)
S_M = pi*d_m^2; % max body cross_sectional area (m^2)

Re_l = rho*v*l_b/mu; % reynolds number based on body length

l_s = l_b - l_n; % body length not including nose or fin (m)
SS_SM = 2.67*(l_n/d_m) + 4*(l_s/d_m); % ratio of forebody wetted area to maximum body cross-sectional area

% corrections of skin friction coefficient to account for cylinder
delCf_turbulent = (1.6E-3*(l_s/d_m))/(Re_l).^0.4;
delCf_laminar = 2*(l_s/d_m)/Re_l;

% skin friction
Cf_turbulent = 0.074*(Re_l).^(-0.2);
Cf_laminar = 1.328./sqrt(Re_l);
Cf_F = 1.328./sqrt(Re_c);

CD0_F = 2*Cf_F*(1 + 2*t/c); % zero-lift drag coefficient for fins

if Re_l<Re_cr
    Cf_b = 0.074./(Re_l^0.2) - Re_cr*(Cf_turbulent - Cf_laminar)/Re_l + delCf_laminar; % skin friction coefficient for body
else
    Cf_b = 0.074./(Re_l^0.2) - Re_cr*(Cf_turbulent - Cf_laminar)/Re_l + delCf_turbulent; % skin friction coefficient for body
end

if Cf_b < 0
    Cf_b = 0;
end

CD_f_B = Cf_b*(1 + 60/(l_b/d_m)^3 + 0.0025*(l_b/d_m))*SS_SM; % forebody drag coefficient
CD_b = (0.029*(d_b/d_m)^3)/sqrt(CD_f_B); % base drag coefficient
CD0_B = CD_f_B + CD_b; % zero-lift drag coefficient for body

CD0_FB = CD0_F*(S_F/S_M) + CD0_B; % overall zero-lift drag coefficient
