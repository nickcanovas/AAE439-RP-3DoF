%% Initialization
clear
clc

% Define Important Constants
const = struct;
const.g = 9.81; % gravity [m/s^2]
const.rho = 1.225; % density of air [kg/m^3]
const.vw = 3; % Wind Velocity [m/s]
const.dt = .0005; % Time step [s]
dt = const.dt;
const.mp_i = 0.0393; % Initial mass of the propellant, [kg]

% Define Aerodynamic Constants
in_to_m = 0.0254; % converts inches to meters
const.l_b = 55*in_to_m; % body length (m)
const.l_n = 11.5*in_to_m; % nose cone length (m)
const.d_m = 3.1*in_to_m; % max body diameter (m)
const.d_b = 3.1*in_to_m; % base diameter (m)
const.t = 0.0025*in_to_m; % fin thickness (m)
const.c_tip = 1*in_to_m; % tip chord length for primary fins (m)
const.c_root = 6*in_to_m; % root chord length for primary fins (m)
const.c_tip_aft = 3.49*in_to_m; % tip chord length for aft fins (m)
const.c_root_aft = 4.53*in_to_m; % root chord length for aft fins (m)
const.h_fin = 5*in_to_m; % primary fin height (m)
const.h_fin_aft = 5*in_to_m; % aft fin height (m)
const.num_finpairs = 4; % number of primary-aft fin pairs

const.Re_cr = 500000; % critical reynolds number
const.mu = 1.789E-5; % dynamic viscosity of air at sea level (kg/m*s)

const.c = (const.c_tip + const.c_root)/2; % average fin chord for primary fins (m)
c_aft = (const.c_tip_aft + const.c_root_aft)/2; % average fin chord for aft fins (m)
const.S_F = const.num_finpairs*(const.c*const.h_fin + c_aft*const.h_fin_aft); % total planform area for all fins (m^2)

%m_initial = 1.1; % Initial mass of the rocket, [kg]
m_initial = 1.24738 + 0.0876 + const.mp_i;
m_prev = m_initial;


% Here is where thrust v. time is input
F = readmatrix('Thrust.txt');
F = F *4.44822; %[N]
t_b = readmatrix('Time.txt');
mdot = const.mp_i/t_b(end); % ASSUMING constant m_dot (!!)

% Define the states
x = zeros(1, length(F)); % X-direction [m]
z = zeros(1, length(F)); % Altitude [m]
v = zeros(1, length(F)); % Velocity [m/s]
t = zeros(1, length(F)); % time [s]
v(1) = F(1)/mdot;
z(1) = 1;
t(1) = 0;
theta = zeros(1, length(F)); % Flight Path Angle [degrees]
psi = zeros(1, length(F)); % Pitch Angle [degrees]
theta(1) = 40.00; % launch angle [deg]

% Define the forces
%F = 10*ones(1,300); % Thrust [N]
L = zeros(1, length(F)); % Lift [N]
D = zeros(1, length(F)); % Drag [N]

% Define other changing values
m = zeros(1, length(F)); % mass, [kg] 
v_rel = zeros(1, length(F));

% Rod Characteristics
L_rod = 1.8288; [m]
CD_pc = 1.3;        % FIX THIS

%% Calculations
i = 1;
while (z(i) > 0)
    t(i+1) = t(i) + const.dt;
    if (sqrt(z(i)^2+x(i)^2)) < L_rod
        disp('on rod.');
        v_rel(i) = sqrt((v(i)*sind(theta(i)))^2+(const.vw + v(i)*cosd(theta(i)))^2);
        psi(i) = acosd((v(i)*sind(theta(i)))/v_rel(i));

        L(i) = lift_function(const, v_rel(i));
        [CD0_FB, CD_b, CD0_F] = AAE439_Project_dragcoeff(v_rel(i), const);
        D(i) = drag_function(const, v_rel(i),CD0_FB);
        % L(i) = lift_function(const, v_rel(i), CD_b, CD0_F);

        if i > 2711
            F(i) = 0;
            %t_b(i) = 0;
            m(i) = m_prev;
        else
            m(i) = m_prev - mdot * dt;
        end

        % fprintf('The cd0 is: %.2f\n', CD0_FB);
        % fprintf('The lift is: %.2f\n', L(i));
        % fprintf('The Drag is: %.2f\n', D(i));
        % fprintf('The mass is: %.2f\n', m(i));
        % fprintf('The Thrust is: %.2f\n', F(i));
        % fprintf('The velocity is: %.2f\n', v(i));

        % [H(i),G(i), H_star(i), G_star(i)] = h_g_function(const,F(i),L(i),D(i),CD0_FB,CD_b,CD0_F, v(i),theta(i),psi(i),z(i), m(i));
        [H(i),G(i), H_star(i), G_star(i)] = h_g_function(const,F(i),L(i),D(i),CD0_FB, v(i),theta(i),psi(i),z(i), m(i));

        % fprintf('The H is: %.2f\n', H(i));
        % fprintf('The G is: %.2f\n', G(i));
        % fprintf('The H_star is: %.2f\n', H_star(i));
        % fprintf('The G_star is: %.2f\n', G_star(i));

        v(i+1) = v(i) + dt/2 * (G(i) + G_star(i));
        theta(i+1) = theta(i);
        x(i+1) = x(i) + dt/2 * (v(i)*cosd(theta(i)) + v(i+1)*cosd(theta(i+1)));
        z(i+1) = z(i) + dt/2 * (v(i)*sind(theta(i)) + v(i+1) * sind(theta(i+1)));
    elseif(t(i) > (t_b(end) + 6))
        disp('parachute deployed.');
        v_rel(i) = sqrt((v(i)*sind(theta(i)))^2+(const.vw + v(i)*cosd(theta(i)))^2);
        psi(i) = acosd((v(i)*sind(theta(i)))/v_rel(i));

        L(i) = 0;
        %[CD0_FB, CD_b, CD0_F] = AAE439_Project_dragcoeff(v_rel(i), const);
        % Insert parachute CD
        D(i) = drag_function(const, v_rel(i),CD_pc);

        if i > 2711
            F(i) = 0;
            %t_b(i) = 0;
            m(i) = m_prev;
        else
            m(i) = m_prev - mdot * dt;
        end
        
        % fprintf('The cd0 is: %.2f\n', CD0_FB);
        % fprintf('The lift is: %.2f\n', L(i));
        % fprintf('The Drag is: %.2f\n', D(i));
        % fprintf('The mass is: %.2f\n', m(i));
        % fprintf('The Thrust is: %.2f\n', F(i));
        % fprintf('The velocity is: %.2f\n', v(i));

        [H(i),G(i), H_star(i), G_star(i)] = h_g_function(const,F(i),L(i),D(i),CD0_FB, v(i),theta(i),psi(i),z(i), m(i));

        % fprintf('The H is: %.2f\n', H(i));
        % fprintf('The G is: %.2f\n', G(i));
        % fprintf('The H_star is: %.2f\n', H_star(i));
        % fprintf('The G_star is: %.2f\n', G_star(i));

        v(i+1) = v(i) + dt/2 * (G(i) + G_star(i));
        theta(i+1) = theta(i) + dt/2 * (H(i) + H_star(i));
        x(i+1) = x(i) + dt/2 * (v(i)*cosd(theta(i)) + v(i+1)*cosd(theta(i+1)));
        z(i+1) = z(i) + dt/2 * (v(i)*sind(theta(i)) + v(i+1) * sind(theta(i+1)));
    else
        v_rel(i) = sqrt((v(i)*sind(theta(i)))^2+(const.vw + v(i)*cosd(theta(i)))^2);
        psi(i) = acosd((v(i)*sind(theta(i)))/v_rel(i));

        L(i) = lift_function(const, v_rel(i));
        [CD0_FB, CD_b, CD0_F] = AAE439_Project_dragcoeff(v_rel(i), const);
        D(i) = drag_function(const, v_rel(i),CD0_FB);
        % L(i) = lift_function(const, v_rel(i), CD_b, CD0_F);

        disp('else');

        if i > 2711
            F(i) = 0;
            %t_b(i) = 0;
            m(i) = m_prev;
        else
            m(i) = m_prev - mdot * dt;
        end

        % fprintf('The cd0 is: %.2f\n', CD0_FB);
        % fprintf('The lift is: %.2f\n', L(i));
        % fprintf('The Drag is: %.2f\n', D(i));
        % fprintf('The mass is: %.2f\n', m(i));
        % fprintf('The Thrust is: %.2f\n', F(i));
        % fprintf('The velocity is: %.2f\n', v(i));

        % [H(i),G(i), H_star(i), G_star(i)] = h_g_function(const,F(i),L(i),D(i),CD0_FB,CD_b,CD0_F, v(i),theta(i),psi(i),z(i), m(i));
        [H(i),G(i), H_star(i), G_star(i)] = h_g_function(const,F(i),L(i),D(i),CD0_FB, v(i),theta(i),psi(i),z(i), m(i));

        % fprintf('The H is: %.2f\n', H(i));
        % fprintf('The G is: %.2f\n', G(i));
        % fprintf('The H_star is: %.2f\n', H_star(i));
        % fprintf('The G_star is: %.2f\n', G_star(i));

        v(i+1) = v(i) + dt/2 * (G(i) + G_star(i));
        theta(i+1) = theta(i) + dt/2 * (H(i) + H_star(i));
        x(i+1) = x(i) + dt/2 * (v(i)*cosd(theta(i)) + v(i+1)*cosd(theta(i+1)));
        z(i+1) = z(i) + dt/2 * (v(i)*sind(theta(i)) + v(i+1) * sind(theta(i+1)));
    end

    m_prev = m(i);
    i = i + 1
end
%%
figure(1)
plot(x,z)
axis("equal")
figure(2)
plot(t,z)
figure(3)
plot(t,v)
[CD0_FB, ~, ~] = AAE439_Project_dragcoeff(v_rel(10730), const)

%% Functions
function[L] = lift_function(const, v)
    rho = const.rho;
    A = const.S_F;

    c_l = 0.3; % cl(rho,v)
	L = (1/2)*rho*A*v^2*c_l;
end

function[D] = drag_function(const, v, cd)
    rho = const.rho;
    A = const.S_F;

	D = (1/2)*rho*A*v^2*cd;
end

function[H,G,H_star, G_star] = h_g_function(const,F,L,D,cd, v,theta,psi,z, m)
% function[H,G,H_star, G_star] = h_g_function(const,F,L,D,cd,CD_b, CD0_F, v,theta,psi,z, m)
	g  = const.g;
    dt = const.dt;
    vw = const.vw;

	G = ((F-D)*cosd(theta-psi) + L*sind(theta-psi) - m*g*cosd(90-theta))/m;
	H = ((F-D)*sind(psi-theta)+L*cosd(psi-theta)-m*g*cosd(theta))/(m*v);

    % Approximations
	v_star = v + G * dt;
    theta_star = theta + H*dt;
    psi_star = atand((v_star*sind(theta_star))/(vw+v_star*cosd(theta_star)));
    z_star = z + v*sind(theta)*dt;
    % L_star = lift_function(const, v_star,CD_b, CD0_F);
    L_star = lift_function(const, v_star);
    D_star = drag_function(const, v_star,cd);
    
    G_star = ((F-D_star)*cosd(theta_star-psi_star) + L_star * sind(theta_star-psi_star) - m*g*cosd(90-theta))/m;
    H_star = ((F-D_star)*sind(psi_star-theta)+L*cosd(psi_star-theta)-m*g*cosd(theta))/(m*v);

end

function [CD0_FB, CD_b, CD0_F] = AAE439_Project_dragcoeff(v, const)
    v = abs(v);% v = vehicle velocity (m/s)
    rho = const.rho; % rho = air density (kg/m^3)

    l_b = const.l_b; % body length (m)
    l_n = const.l_n; % nose cone length (m)
    d_m = const.d_m; % max body diameter (m)
    d_b = const.d_b; % base diameter (m)
    t = const.t; % fin thickness (m)
    c_tip = const.c_tip; % tip chord length for primary fins (m)
    c_root = const.c_root; % root chord length for primary fins (m)
    c_tip_aft = const.c_tip_aft; % tip chord length for aft fins (m)
    c_root_aft = const.c_root_aft; % root chord length for aft fins (m)
    h_fin = const.h_fin; % primary fin height (m)
    h_fin_aft = const.h_fin_aft; % aft fin height (m)
    num_finpairs = const.num_finpairs; % number of primary-aft fin pairs
    
    Re_cr = const.Re_cr; % critical reynolds number
    mu = const.mu; % dynamic viscosity of air at sea level (kg/m*s)

    c = (c_tip + c_root)/2; % average fin chord for primary fins (m)
    c_aft = (c_tip_aft + c_root_aft)/2; % average fin chord for aft fins (m)
    
    Re_c = rho*v*c/mu;
    
    S_F = num_finpairs*(c*h_fin + c_aft*h_fin_aft); % total planform area for all fins (m^2)
    S_M = pi*d_m^2; % max body cross_sectional area (m^2)
    
    Re_l = rho*v*l_b/mu; % reynolds number based on body length
    
    l_s = l_b - l_n; % body length not including nose or fin (m)
    SS_SM = 2.67*(l_n/d_m) + 4*(l_s/d_m); % ratio of forebody wetted area to maximum body cross-sectional area
    
    % corrections of skin friction coefficient to account for cylinder
    delCf_turbulent = (1.6E-3*(l_s/d_m))/(Re_l)^0.4;
    delCf_laminar = 2*(l_s/d_m)/Re_l;
    
    % skin friction
    Cf_turbulent = 0.074*(Re_l)^(-0.2);
    Cf_laminar = 1.328/sqrt(Re_l);
    Cf_F = 1.328/sqrt(Re_c);
    
    CD0_F = 2*Cf_F*(1 + 2*t/c); % zero-lift drag coefficient for fins
    
    if Re_l<Re_cr
        Cf_b = 0.074/(Re_l^0.2) - Re_cr*(Cf_turbulent - Cf_laminar)/Re_l + delCf_laminar; % skin friction coefficient for body
    else
        Cf_b = 0.074/(Re_l^0.2) - Re_cr*(Cf_turbulent - Cf_laminar)/Re_l + delCf_turbulent; % skin friction coefficient for body
    end
    if v < 3.1
        Cf_b = abs(Cf_b);
    end
    CD_f_B = Cf_b*(1 + 60/(l_b/d_m)^3 + 0.0025*(l_b/d_m))*SS_SM; % forebody drag coefficient
    CD_b = (0.029*(d_b/d_m)^3)/sqrt(CD_f_B); % base drag coefficient
    CD0_B = CD_f_B + CD_b; % zero-lift drag coefficient for body
    
    CD0_FB = CD0_F*(S_F/S_M) + CD0_B; % overall zero-lift drag coefficient
end

% function [L] = lift_function(const, v, CD_b, CD0_F) 
%     h_fin = const.h_fin; % Fin length in m
%     c = const.c; % Average cord length in m
%     rho = const.rho; % Air density in kg/m^3
%     S_F = const.S_F; % Fin area in m^2;
% 
%     b = 2*h_fin; % Wingspan in m
%     AR = b/c; % Aspect ratio
%     e = 0.8; % Assumed Oswald Efficiency Factor for fins
%     k = 1/(pi*e*AR);
%     c_l = sqrt((CD_b - CD0_F) / k); 
%     disp(c_l)
% 	L = (1/2)*rho*S_F*v^2*c_l;
% end