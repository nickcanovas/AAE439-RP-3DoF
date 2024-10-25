function[L] = lift_function(const, v, CD_b, CD0_F) 
    h_fin = cosnt.h_fin; % Fin length in m
    c = const.c; % Average cord length in m
    rho = const.rho; % Air density in kg/m^3
    S_F = const.S_F; % Fin area in m^2;

    b = 2*h_fin; % Wingspan in m
    AR = b/c; % Aspect ratio
    e = 0.8; % Assumed Oswald Efficiency Factor for fins
    k = 1/(pi*e*AR);
    c_l = sqrt((CD_b - CD0_F) / k); 
	L = (1/2)*rho*S_F*v^2*c_l;
end