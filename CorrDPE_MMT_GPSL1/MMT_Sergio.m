function [MMT_cost_func3,store_a,store_a3,store_b,store_b3,store_c,store_c3,...
    store_d,store_d3,LOS_codeDelay,NLOS_codeDelay,codeDelay_LOS_indx,...
    codeDelay_NLOS_indx] = MMT_Sergio(remCodePhase,codePhaseStep,blksize,...
                            caCode,iBasebandSignal,qBasebandSignal,MMT_const)

% MMT_Sergio performs Multipath Mitigation Technology (MMT) to produce code
% delay estimates in tracking, along with the complex amplitudes for
% MMT-DPE positioning.

% Only a single reflected path is assumed in the following computation
% i.e., a LOS path and a single reflected path

% [MMT_cost_func3,store_a,store_a3,store_b,store_b3,store_c,store_c3,...
%     store_d,store_d3,LOS_codeDelay,NLOS_codeDelay,codeDelay_LOS_indx,...
%     codeDelay_NLOS_indx] = MMT_Sergio(remCodePhase,codePhaseStep,blksize,...
%                             caCode,iBasebandSignal,qBasebandSignal,MMT_const)

%   Inputs:
%       remCodePhase    - Current tracking epoch's residual code phase
%       codePhaseStep   - Current phasestep based on code freq (variable) and
%                           sampling frequency (fixed)
%       blksize         - Current size of a "block" or code period in whole 
%                           samples
%       caCode          - C/A code
%       iBasebandSignal - Inphase (Real) component of baseband signal
%       qBasebandSignal - Quadrature pahse (Imaginary) component of baseband signal       
%       MMT_const       - MMT amplitude constraint 

%
%   Outputs:
%       MMT_cost_func3      - MMT cost function for the grid of LOS and
%                               reflected code delays
%       store_a             - Grid of real part of complex amplitude of LOS path
%       store_b             - Grid of real part of complex amplitude of reflected path
%       store_c             - Grid of Imag part of complex amplitude of LOS path
%       store_d             - Grid of Imag part of complex amplitude of reflected path
%       store_a3            - Grid of real part of complex amplitude of LOS path,
%                               constrained with Lagrange multiplier
%       store_b3            - Grid of real part of complex amplitude of reflected path,
%                               constrained with Lagrange multiplier
%       store_c3            - Grid of Imag part of complex amplitude of LOS path,
%                               constrained with Lagrange multiplier
%       store_d3            - Grid of Imag part of complex amplitude of reflected path,
%                               constrained with Lagrange multiplier
%       LOS_codeDelay       - LOS path code delay estimate
%       NLOS_codeDelay      - Reflected path code delay estimate
%       codeDelay_LOS_indx  - Index for global maxima for grid of LOS path 
%                               complex amplitudes 
%       codeDelay_NLOS_indx - Index for global maxima for grid of reflected 
%                               path complex amplitudes 

% -------------------------------------------------------------------------
% ------------------------ MMT_Sergio v1.0 --------------------------------
%   Copyright (C) 2025 Sergio Vicenzo
%   Written by Sergio Vicenzo
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, write to the Free Software Foundation, Inc.,
%   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% Last Updated: 13 Jan 2025

%% === Generate multi-correlator values ===================================
chip_spacings = [-(flip(0.1:0.1:0.4)),0,0.1:0.1:0.4];

count = 1;
precalc_correlations    = ...
    zeros(length(chip_spacings),1);

for spacing=(chip_spacings)

    % Define index into the code vector 
    delay_index = ...
        remCodePhase-spacing : ...
            codePhaseStep : ...
            ((blksize-1)*codePhaseStep + ...
           remCodePhase-spacing);

    caCode1 = caCode;
    tcodee = ceil(delay_index)+3;

    s = caCode1(tcodee);

    I = sum(s  .* iBasebandSignal);
    Q = sum(s  .* qBasebandSignal);

    % Store the correlations
    precalc_correlations(count,1) = ...
        sqrt(I.^2 + Q.^2);
    count=count+1;
end % End for getting multicorr values
  
% Get maximum corr value to act a first estimates of LOS and reflected
% path code delay
[~,codedelay_index] = ...
        max(precalc_correlations(:,1),[],'all','linear');

%% == Initialize LOS and reflected path code delays ======================
LOS_codeDelay = (chip_spacings(codedelay_index));
NLOS_codeDelay=LOS_codeDelay;
spacing=0.1;
iterations=1:5;

delay_index_0 = ...
        remCodePhase : ...
        codePhaseStep : ...
        ((blksize-1)*codePhaseStep + ...
        remCodePhase);

caCode1 = caCode;
tcodee_0 = ceil(delay_index_0)+3;
s_0 = caCode1(tcodee_0);

 %% === Multipath Mitiugation Technology ==================================
 % Compute MMT starting with coarse to fine resolution
for iter=iterations

    % Initialize MMT variables
    % Generate grid of LOS and reflected path delay
    if iter==1
        % Grid of reflected path code delays
        % Assume that the max relative delay is 1 chip
        codeDelay_NLOS = ...
            NLOS_codeDelay-0.1:spacing:NLOS_codeDelay+1;
        % Grid of LOS code delays
        codeDelay_LOS  = ...
            LOS_codeDelay-0.1:spacing:LOS_codeDelay+0.1;
    else
        % Grid of LOS code delays
        codeDelay_LOS = ...
            LOS_codeDelay-spacing:spacing:LOS_codeDelay+spacing;
        % Grid of reflected path code delays
        codeDelay_NLOS  = ...
            NLOS_codeDelay-spacing:spacing:NLOS_codeDelay+spacing;
    end
    MMT_cost_func3 = nan(length(codeDelay_LOS),...
                        length(codeDelay_NLOS));

    NLOS_length =length(codeDelay_NLOS);
    LOS_length  =length(codeDelay_LOS);

    store_a = nan(length(codeDelay_LOS),length(codeDelay_NLOS));
    store_b = nan(length(codeDelay_LOS),length(codeDelay_NLOS)); 
    store_c = nan(length(codeDelay_LOS),length(codeDelay_NLOS));
    store_d = nan(length(codeDelay_LOS),length(codeDelay_NLOS));

    store_a3 = nan(length(codeDelay_LOS),length(codeDelay_NLOS));
    store_b3 = nan(length(codeDelay_LOS),length(codeDelay_NLOS)); 
    store_c3 = nan(length(codeDelay_LOS),length(codeDelay_NLOS));
    store_d3 = nan(length(codeDelay_LOS),length(codeDelay_NLOS));

    for delay1_index=1:LOS_length
        delay1=codeDelay_LOS(delay1_index);

        for delay2_index=1:NLOS_length
            delay2=codeDelay_NLOS(delay2_index);

            % Constraint to only compute MMT when the reflected
            % path is larger than LOS path
            if delay2 >= delay1 

                delay_index_1 = ...
                    delay_index_0-delay1 ;

                delay_index_2 = ...
                    delay_index_0-delay2;

                delay_index_12 = ...
                    delay_index_0-(delay2-delay1) ;



                tcodee_1 = ceil(delay_index_1)+3;
                tcodee_2 = ceil(delay_index_2)+3;
                tcodee_12 = ceil(delay_index_12)+3;


                s_1 = caCode1(tcodee_1);
                s_2 = caCode1(tcodee_2);
                s_12 = caCode1(tcodee_12);
                R_xm_1 = sum(iBasebandSignal.*s_1);
                R_ym_1 = sum(qBasebandSignal.*s_1);

                R_xm_2 = sum(iBasebandSignal.*s_2);
                R_ym_2 = sum(qBasebandSignal.*s_2);

                R_mm_12 = sum(s_0.*s_12);
                R_mm_0  = sum(s_0.*s_0); 

                R_1 = (sum(iBasebandSignal.*s_1) ...
                + (sum(qBasebandSignal.*s_1)).*1i);

                R_2 = (sum(iBasebandSignal.*s_2) ...
                + (sum(qBasebandSignal.*s_2)).*1i);


                % Apply amplitude constraint through Lagrange
                % multiplier
                X = (MMT_const^4)*(norm(R_2,1)^2) ...
                    - (MMT_const^2)*(norm(R_1,1)^2);
                    
                Y = 2*(MMT_const^2)*...
                 (R_mm_0*((norm(R_1,1)^2)+(norm(R_2,1)^2)) ...
                 - 2*R_mm_12*real(R_1*R_2));

                Z = norm(R_mm_0*R_2 - R_mm_12*R_1,1)^2 - ...
                 (MMT_const^2)*norm(R_mm_0*R_1 - ...
                 R_mm_12 - R_mm_12*R_2,1)^2;

     
                gamma_1 = (-Y + sqrt(Y^2 - 4*X*Z))/(2*X);
                gamma_12 = sqrt(Y^2 - 4*X*Z);

                a = ((R_mm_0 * R_xm_1) - (R_xm_2*R_mm_12))/...
                        (R_mm_0^2 - R_mm_12^2);
    
                b = ((R_mm_0 * R_xm_2) - (R_xm_1*R_mm_12))/...
                     (R_mm_0^2 - R_mm_12^2);

                c = ((R_mm_0 * R_ym_1) - (R_ym_2*R_mm_12))/...
                    (R_mm_0^2 - R_mm_12^2);
                
                d = ((R_mm_0 * R_ym_2) - (R_ym_1*R_mm_12))/...
                     (R_mm_0^2 - R_mm_12^2);

                % Don't apply amplitude constraint if gamma_1 is
                % imaginary
                if ~isreal(gamma_12)==1
                
                        a3 = a;
        
                        b3 = b;
        
                        c3 = c;
                        
                        d3 = d;

                else

                    a3 = (((R_mm_0 - gamma_1)* R_xm_1) - (R_xm_2*R_mm_12))/...
                            (R_mm_0^2 - gamma_1*R_mm_0*(1-MMT_const^2)...
                            - (MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                    b3 = (((R_mm_0 + ((MMT_const^2) * gamma_1))...
                        * R_xm_2) - (R_xm_1*R_mm_12))/...
                        (R_mm_0^2 - gamma_1*R_mm_0*(1-MMT_const^2)...
                        - (MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                    c3 = (((R_mm_0 - gamma_1)* R_ym_1) - (R_ym_2*R_mm_12))/...
                        (R_mm_0^2 - gamma_1*R_mm_0*(1-MMT_const^2)...
                        - (MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                    d3 = (((R_mm_0 + ((MMT_const^2) * gamma_1))...
                        * R_ym_2) - (R_ym_1*R_mm_12))/...
                        (R_mm_0^2 - gamma_1*R_mm_0*(1-MMT_const^2)...
                        - (MMT_const^2)*(gamma_1^2) - R_mm_12^2);

                end

                MMT_cost_func3(delay1_index,delay2_index) = ...
                    abs(2*real((a3-1i*c3)*...
                    (R_xm_1 + 1i*R_ym_1))...
                    +2*real((b3-1i*d3)*...
                    (R_xm_2 + 1i*R_ym_2)) - ...
                    2*(a3*b3 + c3*d3)*R_mm_12 - ...
                    (a3^2 + b3^2 + c3^2 +d3^2)*R_mm_0);

                store_a3(delay1_index,delay2_index)=...
                    a3 ;
                store_b3(delay1_index,delay2_index)=...
                    b3 ;
                store_c3(delay1_index,delay2_index)=...
                    c3 ;
                store_d3(delay1_index,delay2_index)=...
                    d3 ;

                store_a(delay1_index,delay2_index)=...
                    a ;
                store_b(delay1_index,delay2_index)=...
                    b ;
                store_c(delay1_index,delay2_index)=...
                    c ;
                store_d(delay1_index,delay2_index)=...
                    d ;
            end 

         end % End for delay2

    end % End for delay1

     % Maximise MMT cost function
    [MaxCorrValue,~] = max(MMT_cost_func3,[],'all','linear');

    % Get LOS and reflected path code delays
    [codeDelay_LOS_indx,codeDelay_NLOS_indx,~] = ...
        find(MMT_cost_func3==MaxCorrValue); 

    LOS_codeDelay = mean(codeDelay_LOS(codeDelay_LOS_indx));
    NLOS_codeDelay = mean(codeDelay_NLOS(codeDelay_NLOS_indx));

    % Narrow the grid spacing (resolution) by half
    spacing=spacing/2;

end % end for iter

end
