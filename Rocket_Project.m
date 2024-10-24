%% Initialization
clear
clc

% Define Constants
const = struct;
const.g = 9.81; % gravity [m/s^2]
const.rho = 1.225; % density of air [kg/m^3]
const.A = 1; % Surface Area of the Rocket [m^2]
const.vw = 10; % Wind Velocity [m/s]
const.dt = .01; % Time step [s]
const.mp_i = 0.0393; % kg
m_prev = m_initial;


% Here is where thrust v. time is input
F = readmatrix('Thrust.txt');
t_b = readmatrix('Time.txt');
m_dot = const.mp_i/t_b(end); % ASSUMING constant m_dot (!!)

% Define the states
x = zeros(300); % X-direction [m]
z = zeros(300); % Altitude [m]
v = zeros(300); % Velocity [m/s]
theta = zeros(300); % Flight Path Angle [degrees]
psi = zeros(300); % Pitch Angle [degrees]
psi(1) = pi/4;
theta(1) = pi/4;

% Define the forces
%F = 10*ones(1,300); % Thrust [N]
L = zeros(1,300); % Lift [N]
D = zeros(1,300); % Drag [N]

% Define other changing values
m = zeros(1,300); % mass, [lbm]

%% Calculations
i = 1;
while (z(i) > 0) || (i == 1)
    v_rel(i) = sqrt((v(i)*sind(theta(i)))^2+(const.vw + v(i)*cosd(theta(i)))^2);

    L(i) = lift_function(const, v_rel(i));
    D(i) = drag_function(const, v_rel(i));

    m(i) = m_prev - mdot * dt;

    [H(i),G(i), H_star(i), G_star(i)] = h_g_function(const,F(i),L(i),D(i),v(i),theta(i),psi(i),z(i), m(i));

    v(i+1) = v(i) + dt/2 * (G(i) + G_star(i));
    theta(i+1) = theta(i) + dt/2 * (H(i) + H_star(i));
    x(i+1) = x(i) + dt/2 * (v(i)*cosd(theta(i)) + v(i+1)*cosd(theta(i+1)));
    z(i+1) = z(i) + dt/2 * (v(i)*sind(theta(i)) + v(i+1) * sind(theta(i+1)));

    i = i + 1;
    m_prev = m(i);
end


%% Functions
function[L] = lift_function(const, v)
    rho = const.rho;
    A = const.A;

    c_l = 1; % cl(rho,v)
	L = (1/2)*rho*A*v^2*c_l;
end

function[D] = drag_function(const, v)
    rho = const.rho;
    A = const.A;

    c_d = 1; % cd(rho,v)
	D = (1/2)*rho*A*v^2*c_d;
end

function[H,G,H_star, G_star] = h_g_function(const, F,L,D,v,theta,psi,z, m)
	g  = const.g;
    dt = const.dt;
    vw = const.vw;

	G = ((F-D)*cosd(theta-psi) + L*sind(theta-psi) - m*g*cosd(pi/2-theta))/m;
	H = ((F-D)*sind(psi-theta)+L*cosd(psi-theta)-m*g*cosd(theta))/(m*v);

    % Approximations
	v_star = v + G * dt;
    theta_star = theta + H*dt;
    psi_star = atand((v_star*sind(theta_star))/(vw+v_star*cosd(theta_star)));
    z_star = z + v*sind(theta)*dt;
    L_star = lift_function(const, v_star);
    D_star = drag_function(const, v_star);
    
    G_star = ((F-D_star)*cosd(theta_star-psi_star) + L_star * sind(theta_star-psi_star) - m*g*cosd(pi/2-theta))/m;
    H_star = ((F-D_star)*sind(psi_star-theta)+L*cosd(psi_star-theta)-m*g*cosd(theta))/(m*v);

end
