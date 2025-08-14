%SAMPLEWIFI - Sample data over WiFi interface
%
%   Example that shows how to sample data from a device over WiFi 
%   interface, after setting the device configuration and channel 
%   configuration. The example assumes that the data recorder interface is
%   configured for WiFi measurements (run ChangeDataRecorderInterface.m
%   example first). The example shows how to enable repair logging and how 
%   to perform a repair after samples have been missed. 

% If the script has been terminated manually, it can occur that the device
% keeps sampling. To be sure, a lib.cleanUp() is performed prior to running
% the example.
if exist('lib', 'var')
    disp('[N/A] Ensure library is cleaned...')
    lib.cleanUp();
end

% Initialize the library
lib = TMSiSAGA.Library();

% Desired device configuration
deviceConfig = struct('ImpedanceMode', false, ...
                      'BaseSampleRate', 4000, ...
                      'Dividers', {{'uni', 0; 'bip', 0; 'dig', 0}}, ...
                      'RepairLogging', true, ...
                      'Triggers', false);  

% Desired channel configuration
channelConfig=struct('uni', [1:6], 'bip', [1 2], 'acc',1);

% Number of retries the device may do to find the first available device
num_retries = 10;

% Code within the try-catch to ensure that all devices are stopped and
% closed properly in case of a failure.
try    
    % Get a single device from the connected devices
    % possible interfaces 'usb'/'network' and 'electrical'/'optical/'wifi'
    device = lib.getFirstAvailableDevice('usb', 'wifi', num_retries);
    
    % Open a connection to the device
    device.connect();
    
    % Set device configuration
    device.setDeviceConfig(deviceConfig);
    
    % Set channel configuration
    device.setChannelConfig(channelConfig);    
    
    % Create a real time plot
    rPlot = TMSiSAGA.RealTimePlot('Plot', device.sample_rate, device.getActiveChannels());
    rPlot.show();
    
    % Create a timestamp for the measurement (saving needed for Repair)
    timestamp = datestr(datetime('now'),'dd-mm-yyyy_HH.MM.SS');
    
    % Create a file storage
    poly5 = TMSiSAGA.Poly5(['./SampleWiFi_' timestamp '.poly5'], ...
        'Plot', device.sample_rate, device.getActiveChannels());
    
    % Start sampling on the device
    device.start();
    
    % As long as we do not press the X or 'q' keep on sampling from the
    % device.
    while rPlot.is_visible
        % Sample from device
        [samples, num_sets, type] = device.sample();
        
        % Append samples to the plot and redraw
        if num_sets > 0
            rPlot.append(samples);
            poly5.append(samples);
            rPlot.draw();
        end
    end
    
    % Stop sampling on the device
    device.stop();
    
    % Close file
    poly5.close();
    
    % Get repair data
    [repair_data, num_sets] = device.getMissingSamples();
    
    % Repair data and save to SampleWiFi.poly5.repaired.poly5
    TMSiSAGA.Repair.repairPoly5(['SampleWiFi_' timestamp '.poly5'], ...
        repair_data, numel(device.getActiveChannels()));
    
    % Disconnect from device
    device.disconnect();
    
    % clean up and unload the library
    lib.cleanUp();
    
catch e
    % In case of an error close all still active devices and clean up
    % library itself
    lib.cleanUp();
    
    % Rethrow error to ensure you get a message in console
    rethrow(e)
end