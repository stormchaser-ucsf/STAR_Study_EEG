%EMG_WORKFLOW_EXAMPLE - This example shows a workflow of an EMG measurement. 
%
%   The workflow starts with an impedance check. This impedance check is 
%   stopped when the figure window is closed. Next, a heat map of the 
%   muscle activation is shown.  Channels that are out of range are marked red, as this 
%   could indicate that the connection between electrode and skin are not 
%   good.
%
%   Visualisation functions work for 32 channels or 64 channels. 
%
%   By default, a window size of 0.5 seconds is chosen to plot the heat 
%   map. The window size should not be chosen to small, as creation of the
%   heat map requires processing time on your pc during which no calls are
%   made to the device.

% If the script has been terminated manually, it can occur that the device
% keeps sampling. To be sure, a lib.cleanUp() is performed prior to running
% the example. 
if exist('lib', 'var')
    disp('[N/A] Ensure library is cleaned...')
    lib.cleanUp();
end

if ~license('test', 'Signal_Toolbox')
    error('This example requires the Signal Processing Toolbox to be installed')    
end

% Initialize the library
lib = TMSiSAGA.Library();

% Normalisation initialisation
normalisation = 0;

% Direction the connector points to seen from an anterior view of the
% frontal plane ('left', 'down', 'right' and 'up')
connector_orientation = 'up';

% Impedance configuration
config_impedance = struct('ImpedanceMode', true, ... 
                          'ReferenceMethod', 'common', ...
                          'Triggers', false);

% Measurement configuration
config_measurement = struct('ImpedanceMode', false, ... 
                          'ReferenceMethod', 'common', ...
                          'AutoReferenceMethod', true, ...
                          'BaseSampleRate', 4000, ...
                          'Dividers', {{'uni', 0;}}, ...
                          'RepairLogging', false);
                      
% Channel configuration
channel_config = struct('uni',1:64);
                      
% Filter configuration
order = 2;
Fc = 10; % Hz

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
          
    % Update configuration for impedance measurement
    device.setDeviceConfig(config_impedance);
    
    % Update channel configuration
    device.setChannelConfig(channel_config);
    
    % Create a list of channel names that can be used for printing the
    % impedance values next to the figure
    for i=1:length(device.getActiveChannels())
        channel_names{i}=device.getActiveChannels{i}.alternative_name;
    end
    
    % Create an ImpedancePlot object
    iPlot = TMSiSAGA.ImpedancePlot('Impedance of the HD EMG grid', channel_config.uni, channel_names, SAGA_type);
    iPlot.show();

    % Start sampling on the device
    device.start();
    
    % Remain in impedance mode until the figure is closed
    while iPlot.is_visible
        % Sample from device
        [samples, num_sets, type] = device.sample();
        
        % Append samples to the plot and redraw
        % need to divide by 10^6.
        if num_sets > 0
            s=samples ./ 10^6;
            iPlot.grid_layout(s);
        end        
    end
    
    % Stop sampling on the device
    device.stop();   
    
    % Update device configuration
    device.setDeviceConfig(config_measurement);
    
    % Initialise the sample buffer and the window_size
    sample_buffer = zeros(numel(device.getActiveChannels()), 0);
    window_seconds = 0.5;
    window_samples = round(window_seconds * device.sample_rate);
    
    % Initialise the objects to be used in the workflow
    vPlot = TMSiSAGA.Visualisation(device.sample_rate, device.getActiveChannels(), window_samples);
    
    % Create a Poly5 file storage
    poly5 = TMSiSAGA.Poly5(['./HD_EMG_measurement_' datestr(datetime('now'), 'dd-mm-yyyy_HH.MM.SS') '.poly5'], ...
        'Plot', device.sample_rate, device.getActiveChannels());
   
    % Start sampling from the device again
    device.start();
    
    while vPlot.is_visible
        % Sample from device
        [samples, num_sets, type] = device.sample();    
               
        % Append samples to the plot and redraw
        if num_sets > 0
            poly5.append(samples);
            
            % Append samples to a buffer, so that there is always
            % a minimum of window_samples to process the data
            sample_buffer(:, size(sample_buffer, 2) + size(samples, 2)) = 0;
            sample_buffer(:, end-size(samples, 2) + 1:end) = samples;
            
            % As long as we have enough samples calculate RMS values.
            while size(sample_buffer, 2) >= window_samples
                % Find the RMS value of the window
                data_plot = vPlot.RMSEnvelope(sample_buffer(:,1:window_samples), Fc, order);
                
                % Call the visualisation function
                vPlot.EMG_Visualisation(sample_buffer, data_plot, normalisation,...
                    'Heat map of the muscle activation', ...
                    connector_orientation)
                
                % Check if there are no missed samples due to processing
                % time.
                if any(diff(sample_buffer(end,:)) > 1)
                    % Cleanup internal data storage buffer
                    poly5.close();
                    error('Sample counter error. Samples are lost due to slow processing of the HD EMG visualisation.')
                end
                
                % Clear the processed window from the sample buffer
                sample_buffer = sample_buffer(:, window_samples + 1:end);
            end
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