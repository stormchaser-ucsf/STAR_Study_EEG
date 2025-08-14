%SAMPLEELECTRICAL - Sample data over electrical interface
%
%   Example that shows how to sample data from a device, after setting the
%   device configuration and channel configuration. Furthermore, the
%   example shows how to save data to a Poly5 file and stream data to a
%   RealTimePlot object.

% If the script has been terminated manually, it can occur that the device
% keeps sampling. To be sure, a lib.cleanUp() is performed prior to running
% the example.
if exist('lib', 'var')
    disp('[N/A] Ensure library is cleaned...')
    lib.cleanUp();
end

% Initialize the library
lib = TMSiSAGA.Library();

% Specify desired device configuration
deviceConfig = struct('ImpedanceMode', false,...
                          'BaseSampleRate', 4000, ...
                          'Dividers', {{'uni', 0; 'bip', 0; 'aux', 0; 'dig', 0}});  

% Specify desired channel configuration
channelConfig=struct('uni', [1:32], 'bip', [1:4], 'aux', [1:3], 'acc',1, 'dig', 0);
                      
% Code within the try-catch to ensure that all devices are stopped and
% closed properly in case of a failure.
try
    
    % Get a single device from the connected devices
    % possible interfaces 'usb'/'network' and 'electrical'/'optical/'wifi'
    device = lib.getFirstAvailableDevice('usb', 'electrical');
    
    % Open a connection to the device
    device.connect();
    
    % Set device configuration
    device.setDeviceConfig(deviceConfig);
    
    % Set channel configuration
    device.setChannelConfig(channelConfig);    
    
    % Create a real time plot
    rPlot = TMSiSAGA.RealTimePlot('Plot', device.sample_rate, device.getActiveChannels());
    rPlot.show();
    
    % Create a file storage
    poly5 = TMSiSAGA.Poly5(['./SampleElectrical_' datestr(datetime('now'),'dd-mm-yyyy_HH.MM.SS') '.poly5'], ...
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