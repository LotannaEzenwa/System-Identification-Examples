function [Y,X] = dlpv_response(sys,a,x0,u,feedback)
%LPV_RESPONSE Calculates the response of Discrete-Time LPV System
%   Has a Discrete Linear Parameter Varying System and calculates the response over
%   a varying parameter from an initial state and varying input.

sys.initial_state = x0;
sys.parameter = a;
sys.input = u;
n = size(a,1);


end

