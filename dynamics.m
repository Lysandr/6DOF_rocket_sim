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
    
    %% Acceleration
    omega_E = [0 0 7.29211505392569e-05].'; % equatorial rotation of earth
    H = 7250;                       % m
    rho_naught = 1.225;             % kg/m^3
    A = p.Aref;                     % m^2
    % make this dependant upon the MRP at some point 
    r = norm(x_N) + 6371000;
    Re = 6371000;
    vatm = v_N - skew((Re/r)*omega_E)*x_N;
    % taking care if singular cases
    if r-Re >= 0
        rho = rho_naught*exp(-(r-Re)/H);
    else
        rho = 0;
    end
    Cd = p.Cd;
    v = norm(vatm) + 1E-12;
    v_hat = vatm./v;
    q = (rho*v*v)/2;
    a_drag = -(q*Cd*A*v_hat)/m;
    v_dot = ((MRP2C(sigma_BN).'*Thrust_B)/m) + p.g + a_drag;
    
    %% MRP integration, remember to perform norm check on this badboi: 3X1
    sigma_dot = (0.25.*((1 -(sigma_BN.'*sigma_BN))*eye(3) + 2*skew(sigma_BN) ...
        + (2*sigma_BN*(sigma_BN.'))))*omega_BN;
    % Angular Rate Integration
    omega_dot = Ic\((-skew(omega_BN)*(Ic*omega_BN)) + (skew(p.r_E_COM)*Thrust_B) ...
        + p.L);
    % Send out the derivative
    f_x = [m_dot x_dot.' v_dot.' sigma_dot.' omega_dot.'].'; 
    
    
end












