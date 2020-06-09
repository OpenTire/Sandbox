function [outMF] = mfeval_wrapper(parameterSource, Fz, kappa, alpha, gamma, phit, Vx, useMode)


% Pack inputs
inputsMF = [Fz kappa alpha gamma phit Vx];

% Call mfeval
outMF = mfeval(parameterSource, inputsMF, useMode);

end


