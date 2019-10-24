% Test NI based interfacing with the Cortech USB/Parallel interface for
% TTL communication with the Biosemi USB receiver

s = daq.createSession('ni') % identify the interface
ch = addDigitalChannel(s,'Dev1', 'Port2/Line0:7', 'OutputOnly') % specify the channels to use
    % (port 2 is connected with the USB receiver)
outputSingleScan(s, dec2binvec(64,8)) % output 64 in 8 binary values (1 per channel)
outputSingleScan(s, [0 0 0 0 0 0 0 0]) % immediately reset to 0

% Clear variables when done with session
clear all
close all
