function sampledCode = sampledCAcode(PRN,sampFreq,numSamples)

% Function to sample the C/A code at the specified sampling frequency
%
% Returns 1 for binary 0 and -1 for binary 1
%
% Input
% PRN         - sat id number (PRN = 1-37, WAAS = 120-158)
% sampFreq      - sampling frequency
% numSamples    - number of samples to create
%
% Output
% sampledCode   - vector containing the sampled code sequence


code = generateCAcode(PRN);

codeLength  = 1023;
codeRate    = 1.023e6;

stepSize = codeRate / sampFreq;        % step size (in chips) for sampling

% sample vector at stepSize increments
samples = 0:stepSize:(numSamples*stepSize)-stepSize;

samples = floor(samples);                   % convert to integer values
samples = rem(samples,codeLength);         % handle samples > codeLength
samples = samples + 1;                      % convert from 0-based to 1-based
sampledCode = code(samples);               % grab values from original ranging code