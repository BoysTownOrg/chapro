iirfb                                  % Generate IIR filter coefficients
e=system('tst_iffb test/carrots.wav'); % Apply AFC processing
if (e==0)                              % Check for error
    tst_iffb                           % Show results
end
