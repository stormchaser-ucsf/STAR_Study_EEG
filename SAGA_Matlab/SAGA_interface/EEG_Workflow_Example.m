%EEG_WORKFLOW_EXAMPLE - Example that shows a workflow of an EEG measurement. 
%
%   The workflow starts with an Impedance check. This impedance check is 
%   stopped when the figure window is closed. Then it continues with a 
%   measurement that is plotted using a RealTimePlot object. Closing the 
%   window terminates the measurement.

% If the script has been terminated manually, it can occur that the device
% keeps sampling. To be sure, a lib.cleanUp() is performed prior to running
% the example. 
if exist('lib', 'var')
    disp('[N/A] Ensure library is cleaned...')
    lib.cleanUp();
end

% Initialize the library
lib = TMSiSAGA.Library();

% Impedance configuration
config_impedance = struct('ImpedanceMode', true, ... 
                          'ReferenceMethod', 'average');

% Measurement configuration
config_measurement = struct('ImpedanceMode', false,...
                          'BaseSampleRate', 4000, ...
                          'Dividers', {{'uni' 3}});

% Channel configuration                      
channel_config = struct('uni', 1:32);

% Code within the try-catch to ensure that all devices are stopped and
% closed properly in case of a failure.
try
    % Get a single device from the connected devices
    device = lib.getFirstAvailableDevice('usb', 'electrical');
    
    % Open a connection to the device
    device.connect();
    
    % Set whether a SAGA32+ or SAGA64+ by getting the total number of
    % hardware channels and substracting 15 (for other channel types).
    SAGA_type = device.num_hw_channels - 15;
    
    % Load the locationsfile for the head lay-out impedance plot
    load(['EEGChannels' num2str(SAGA_type) 'TMSi.mat']);
    
    % Update the channel names to the EEG name convention
    for ii = 1:numel(channel_config.uni)
        device.channels{channel_config.uni(ii)+1}.setAlternativeName(ChanLocs(channel_config.uni(ii)).labels);
    end      
        
    % Update the device configuration
	device.setDeviceConfig(config_impedance);
    
    % Update the channel configuration
    device.setChannelConfig(channel_config);
       
    % Create a list of channel names that can be used for printing the
    % impedance values next to the figure
    for i=1:length(device.getActiveChannels())
        channel_names{i}=device.getActiveChannels{i}.alternative_name;
    end
    
    % Initialise the ImpedancePlot object
    iPlot = TMSiSAGA.ImpedancePlot('Topographic plot Impedance', channel_config.uni, channel_names, SAGA_type);
    iPlot.show();
       
    % Start sampling on the device
    device.start();
    
    % As long as we do not press the X or 'q' keep on sampling from the
    % device.
    while iPlot.is_visible
        % Sample from device
        [samples, num_sets, type] = device.sample();
        
        % Append samples to the plot and redraw
        % need to divide by 10^6.
        if num_sets > 0
            s=samples ./ 10^6;
            
            % Make a topographic plot of the impedance
            iPlot.head_layout(s, ChanLocs); 
        end
    end
    
    % Stop sampling on the device
    device.stop();
    
    % Update the device configuration
	device.setDeviceConfig(config_measurement);    
       
    % Create a real time plot
    rPlot = TMSiSAGA.RealTimePlot('Plot of the measured EEG channels', device.sample_rate, device.getActiveChannels());
    rPlot.show();
    
    % Create a poly5 file storage
    poly5 = TMSiSAGA.Poly5(['./EEG_measurement_' datestr(datetime('now'),'dd-mm-yyyy_HH.MM.SS') '.poly5'], ...
    'Plot', device.sample_rate, device.getActiveChannels());
        
    % Start sampling on the device
    device.start();
    
    while rPlot.is_visible
        % Sample from device
        [samples, num_sets, type] = device.sample();
        
        % Only enter this part of the code if the samples variable is not
        % empty
        if num_sets > 0
            % Append the data to the Poly5 file
            poly5.append(samples);
                        
            % Append samples to the plot and redraw
            rPlot.append(samples);
            rPlot.draw();
        end
    end
    
    % Cleanup internal data storage buffer
    poly5.close();
    
    % Stop sampling on the device
    device.stop();
    
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