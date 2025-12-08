function [pos,el, az, dop,satpos] = leastSquarePos(satpos,obs,settings)
%Function calculates the Least Square Solution.
%Edited to also output satellite position - Sergio Vicenzo
% 14 Feb 2024
%[pos, el, az, dop] = leastSquarePos(satpos, obs, settings);
%
%   Inputs:
%       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)
%       obs         - Observations - the pseudorange measurements to each
%                   satellite corrected by SV clock error
%                   (e.g. [20000000 21000000 .... .... .... .... ....]) 
%       settings    - receiver settings
%
%   Outputs:
%       pos         - receiver position and receiver clock error 
%                   (in ECEF system: [X, Y, Z, dt]) 
%       el          - Satellites elevation angles (degrees)
%       az          - Satellites azimuth angles (degrees)
%       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])
%       satpos      - Satellites positions after Least Squares for DPE
%                   (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: leastSquarePos.m,v 1.1.2.12 2006/08/22 13:45:59 dpl Exp $
%==========================================================================

%=== Initialization =======================================================
nmbOfIterations = 10;

dtr     = pi/180;
pos     = zeros(4, 1);   % center of earth
X       = satpos;
nmbOfSatellites = size(satpos, 2);

A       = zeros(nmbOfSatellites, 4);
omc     = zeros(nmbOfSatellites, 1);
az      = zeros(1, nmbOfSatellites);
el      = az;

% To save the updated satellite positions for DPE
% Added by Sergio Vicenzo - 16 Feb 2024
Rot_X_up = zeros(3, nmbOfSatellites);

%=== Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations

    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);
            trop = 2;
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                   (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c ;

            %--- Correct satellite position (do to earth rotation) --------
            % Convert SV position at signal transmitting time to position 
            % at signal receiving time. ECEF always changes with time as 
            % earth rotates.
            Rot_X = e_r_corr(traveltime, X(:, i));
            Rot_X_up(:,i) = Rot_X;
            
            %--- Find the elevation angle of the satellite ----------------
            [az(i), el(i), ~] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));

            if (settings.useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end
            
        end % if iter == 1 ... ... else 

        %--- Apply the corrections ----------------------------------------
        omc(i) = ( obs(i) - norm(Rot_X - pos(1:3), 'fro') - pos(4) - trop ); 

        %--- Construct the A matrix ---------------------------------------
        A(i, :) =  [ (-(Rot_X(1) - pos(1))) / norm(Rot_X - pos(1:3), 'fro') ...
                     (-(Rot_X(2) - pos(2))) / norm(Rot_X - pos(1:3), 'fro') ...
                     (-(Rot_X(3) - pos(3))) / norm(Rot_X - pos(1:3), 'fro') ...
                     1 ];
    end % for i = 1:nmbOfSatellites

    % These lines allow the code to exit gracefully in case of any errors
    if rank(A) ~= 4
        pos     = zeros(1, 4);
        dop     = inf(1, 5);
        fprintf('Cannot get a converged solution! \n');
        return
    end

    %--- Find position update (in the least squares sense)-----------------
    x   = A \ omc;
    
    %--- Apply position update --------------------------------------------
    pos = pos + x;
    
end % for iter = 1:nmbOfIterations

%--- Fixing result --------------------------------------------------------
pos = pos';
% Saves satellite position for more accurate satellite position for DPE 
% Added by Sergio Vicenzo - 16 Feb 2024
satpos=Rot_X_up; 

%=== Calculate Dilution Of Precision ======================================

% Commented to always calculate DOP
% Sergio Vicenzo - 16 Feb 2024

% if nargout  == 4
    %--- Initialize output ------------------------------------------------
    dop     = zeros(1, 5);
    
    %--- Calculate DOP ----------------------------------------------------
    Q       = inv(A'*A);
    
    dop(1)  = sqrt(trace(Q));                       % GDOP    
    dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
    dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
    dop(4)  = sqrt(Q(3,3));                         % VDOP
    dop(5)  = sqrt(Q(4,4));                         % TDOP
% end  % if nargout  == 4
