%script to run fwd and inv cqts

clear all;
addpath ./scatt_functions

Npad = 2^15;
T = 2048;
scparam.N = Npad;
scparam.T = T;
scparam.Q = 32;
%scparam.dse = 1;
filts = create_scattfilters( scparam );

input = randn(Npad, 1);

cosa = fwdcqt( input, scparam, filts{1});
aviam = invcqt( cosa, scparam, filts{1}) ;

