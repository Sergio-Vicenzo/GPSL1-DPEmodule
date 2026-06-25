function [satPositions,satVelocity,satClkCorr] = ...
    satPosVel(transmitTime,prnList,eph,ephRecNum,activeChnList)

%SATPOS and Velocity Calculation of X,Y,Z satellites coordinates at TRANSMITTIME 
%for given ephemeris EPH. Coordinates are calculated for each satellite in the
%list PRNLIST.
%[satPositions,satVelocity,satClkCorr] = satPosVel(transmitTime, prnList,eph)
%
%   Inputs:
%       transmitTime  - transmission time for all satellites
%       prnList       - list of PRN-s to be processed
%       eph           - ephemeridies of satellites
%       ephRecNum     - Number of the ephemeris records
%       activeChnList - Active channel list
%
%   Outputs:
%       satPositions  - positions of satellites (in ECEF system [X; Y; Z;])
%       satVelocity   - velocity of satellites (in ECEF system [VX; VY; VZ;])
%       satClkCorr    - correction of satellites clocks

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre 04-09-96
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%Modified by Xiaofan Li at University of Colorado at Boulder
% CVS record:
% $Id: satpos.m,v 1.1.2.15 2006/08/22 13:45:59 dpl Exp $

%% Initialize constants ===================================================

numOfSatellites = size(transmitTime, 2);

% GPS constatns

gpsPi          = 3.1415926535898;  % Pi used in the GPS coordinate
% system

%--- Constants for satellite position calculation -------------------------
Omegae_dot     = 7.2921151467e-5;  % Earth rotation rate, [rad/s]
GM             = 3.986005e14;      % Earth's universal
% gravitational parameter,
% [m^3/s^2]
F              = -4.442807633e-10; % Constant, [sec/(meter)^(1/2)]

%% Initialize results =====================================================
satClkCorr   = zeros(1, numOfSatellites);
satPositions = zeros(3, numOfSatellites);
satVelocity  = zeros(3, numOfSatellites);

%% Process each satellite =================================================

for satNr = activeChnList
    
    prn  = prnList(satNr);
    index= 1;
    
    %% Find initial satellite clock correction --------------------------------

    %--- Find time difference ---------------------------------------------
    dt = check_t(transmitTime(satNr) - eph(index,prn).t_oc);

    %--- Calculate clock correction ---------------------------------------
    satClkCorr(satNr) = (eph(index,prn).a_f2 * dt + eph(index,prn).a_f1) * dt + ...
        eph(index,prn).a_f0 - eph(index,prn).T_GD;

    time = transmitTime(satNr) - satClkCorr(satNr);
    
    %% Find satellite's position ----------------------------------------------

    %Restore semi-major axis
    a   = eph(index,prn).sqrtA * eph(index,prn).sqrtA;

    %Time correction
    tk  = check_t(time - eph(index,prn).t_oe);

    %Initial mean motion
    n0  = sqrt(GM / a^3);
    %Mean motion
    n   = n0 + eph(index,prn).deltan;
    
    %Mean anomaly
    M   = eph(index,prn).M_0 + n * tk;
    %Reduce mean anomaly to between 0 and 360 deg
    M   = rem(M + 2*gpsPi, 2*gpsPi);
    
    %New Line
    M_dot = n;
    

    
    %Initial guess of eccentric anomaly
    E   = M;

    %--- Iteratively compute eccentric anomaly ----------------------------
    for ii = 1:10
        E_old   = E;
        E       = M + eph(index,prn).e * sin(E);
        dE      = rem(E - E_old, 2*gpsPi);

        if abs(dE) < 1.e-12
            % Necessary precision is reached, exit from the loop
            break;
        end
    end

    %Reduce eccentric anomaly to between 0 and 360 deg
    E   = rem(E + 2*gpsPi, 2*gpsPi);
    
    % New Line
    E_dot = M_dot/(1.0 - eph(index,prn).e*cos(E));
    
    %Compute relativistic correction term
    dtr = F * eph(index,prn).e * eph(index,prn).sqrtA * sin(E);

    %Calculate the true anomaly
    nu   = atan2(sqrt(1 - eph(index,prn).e^2) * sin(E), cos(E)-eph(index,prn).e);
    
    %New Line
    nu_dot=sin(E)*E_dot*(1.0+eph(index,prn).e*cos(nu))/(sin(nu)*(1.0-eph(index,prn).e*cos(E)));
    

    
    
    %Compute angle phi
    phi = nu + eph(index,prn).omega;
    %Reduce phi to between 0 and 360 deg
    phi = rem(phi, 2*gpsPi);

    %Correct argument of latitude
    u = phi + ...
        eph(index,prn).C_uc * cos(2*phi) + ...
        eph(index,prn).C_us * sin(2*phi);
    %Correct radius
    r = a * (1 - eph(index,prn).e*cos(E)) + ...
        eph(index,prn).C_rc * cos(2*phi) + ...
        eph(index,prn).C_rs * sin(2*phi);
    %Correct inclination
    i = eph(index,prn).i_0 + eph(index,prn).iDot * tk + ...
        eph(index,prn).C_ic * cos(2*phi) + ...
        eph(index,prn).C_is * sin(2*phi);
    
    % New Lines
    u_dot=nu_dot+2.0*(eph(index,prn).C_us*cos(2.0*u)-eph(index,prn).C_uc*sin(2.0*u))*nu_dot;
    r_dot = a*eph(index,prn).e*sin(E)*n/(1.0-eph(index,prn).e*cos(E)) + ...
        2.0*(eph(prn).C_rs*cos(2.0*u)-eph(index,prn).C_rc*sin(2.0*u))*nu_dot;
    i_dot = eph(index,prn).iDot + (eph(index,prn).C_is*cos(2.0*u)-eph(index,prn).C_ic*sin(2.0*u))*2.0*nu_dot;
    
 
    xp=r*cos(u);
    yp=r*sin(u);
    
    xp_dot=r_dot*cos(u)-yp*u_dot;
    yp_dot=r_dot*sin(u)+xp*u_dot;
    
    %Compute the angle between the ascending node and the Greenwich meridian
    Omega = eph(index,prn).omega_0 + (eph(index,prn).omegaDot - Omegae_dot)*tk - ...
        Omegae_dot * eph(index,prn).t_oe;
    %Reduce to between 0 and 360 deg
    Omega = rem(Omega + 2*gpsPi, 2*gpsPi);
    
    % New Line
    Omega_dot = (eph(index,prn).omegaDot-Omegae_dot);
    
    %--- Compute satellite coordinates ------------------------------------
    satPositions(1, satNr) = xp * cos(Omega) - yp * cos(i)*sin(Omega);
    satPositions(2, satNr) = xp * sin(Omega) + yp * cos(i)*cos(Omega);
    satPositions(3, satNr) = yp * sin(i);
    
    
    % New Lines-- Compute satellite veloticy------------------------------
    satVelocity(1, satNr)= (xp_dot-yp*cos(i)*Omega_dot)*cos(Omega)-...
        (xp*Omega_dot+yp_dot*cos(i)-yp*sin(i)*i_dot)*sin(Omega);
    satVelocity(2, satNr)= (xp_dot-yp*cos(i)*Omega_dot)*sin(Omega)+...
        (xp*Omega_dot+yp_dot*cos(i)-yp*sin(i)*i_dot)*cos(Omega);
    satVelocity(3, satNr)= yp_dot*sin(i)+yp*cos(i)*i_dot;
   

    
    %% Include relativistic correction in clock correction ----------------
    satClkCorr(satNr) = (eph(index,prn).a_f2 * dt + eph(index,prn).a_f1) * dt + ...
        eph(index,prn).a_f0 - ...
        eph(index,prn).T_GD + dtr;
    
end % for satNr = 1 : numOfSatellites