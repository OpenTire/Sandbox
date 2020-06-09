%% This script is used to auto define the C++ variable types. 

% Define inputs
parameterSource = mfeval.readTIR('TNO_car205_60R15.tir');
Fz = 3000;
kappa = 0;
alpha = 0;
gamma = 0;
phit = 0;
Vx = 10;
useMode = 111;

% Call the wrapper
[outMF] = mfeval_wrapper(parameterSource, Fz, kappa, alpha, gamma, phit, Vx, useMode);

