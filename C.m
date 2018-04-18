function [value] = C(x) 
%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: C_Q2.m
%--------------------------------------------------------------------------
% Description: Flow set
% Return 0 if outside of C, and 1 if inside C
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

    global T_t
    
    % H_a states
    tauN = x(5);
    taud = x(6);
    q    = x(7);
    
    % H_b states
    tauT = x(14);
    tauM = x(15);
    tauS = x(16);
    qP = x(17);  
    
    if ((tauN > 0) && q == 0) || ((taud > 0) && q == 1) || ((tauM > 0) && (qP == 1)) || ((tauS > 0)&& (qP == 2)) || (tauT <= T_t)
        value = 1;
    else
        value = 0;
    end
    
end