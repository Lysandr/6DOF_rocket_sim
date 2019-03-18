function f_x = dynamics(~,state_in,p)
%%This function defines all of the system dynamics
% TODO: add in proper gravitational dynamics,
% add in rho variation
% add in attitude varied drag function

    % pull out struct constants
    Ic = p.Ic;
    Thrust_B = p.T_B;

    % distribute the states for use! these are row vectors.
    m =         state_in(1);
    x_N =       state_in(2:4);
    v_N =       state_in(5:7);
    sigma_BN =  state_in(8:10);
    omega_BN =  state_in(11:13);


    % Mass depletion dynamics
    m_dot = -p.mdot;

    % Velocity -- Inertial
    x_dot = v_N;

    % Acceleration -- Inertial
    % (p.mu/r^3)*r_v
    v_dot = ((MRP2C(sigma_BN).'*Thrust_B)/m) + p.g; % - ((0.5*v_N*v_N.'*rho*p.Cd*p.A)/m)

    % MRP integration, remember to perform norm check on this badboi: 3X1
    sigma_dot = (0.25.*((1 -(sigma_BN.'*sigma_BN))*eye(3) + 2*skew(sigma_BN.') ...
        + (2*sigma_BN*(sigma_BN.'))))*omega_BN;
    
    % Angular Rate Integration
    omega_dot = (-skew(omega_BN)*(Ic*omega_BN)) + Ic\(skew(p.r_E_COM)*Thrust_B) ...
        + Ic\p.L;
    
    % Send out the derivative
    f_x = [m_dot x_dot.' v_dot.' sigma_dot.' omega_dot.'].'; 
end












