%SAMPLEIMPEDANCEVISUAL - Example that shows a visualisation for the impedance measurement.
%
%   Example that shows how to make a graphical visualisation of the 
%   skin-electrode impedance for a HD EMG or EEG measurement. The impedance  
%   can be plotted in a head or grid layout. The example calls the 
%   .mat-file EEGChannels32/64TMSi.mat when the head layout is chosen.

% If the script has been terminated manually, it can occur that the device
% keeps sampling. To be sure, a lib.cleanUp() is performed prior to running
% the example. 
if exist('lib', 'var')
    disp('[N/A] Ensure library is cleaned...')
    lib.cleanUp();
end

% Specify the desired configuration
device_config = struct('ImpedanceMode', true, ... 
                'ReferenceMethod', 'common');
 
% Set the channel configuration
channel_config = struct('uni', 1:32);

% Initialize the library
lib = TMSiSAGA.Library();

% Lay-out of the plot ('grid' or 'head')
lay_out = 'head';

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
    
    % Load the EEG location file if the head lay-out is used for plotting
    if strcmp(lay_out,'head')
        load(['EEGChannels' num2str(SAGA_type) 'TMSi.mat']);
    end
    
    % Set the device configuration
    device.setDeviceConfig(device_config);
    
    % Set the channel configuration
    device.setChannelConfig(channel_config);
    
    % Create a list of channel names that can be used for printing the
    % impedance values next to the figure
    for i=1:length(device.getActiveChannels())
        channel_names{i}=device.getActiveChannels{i}.alternative_name;
    end
    
    % Create an ImpedancePlot object that is used for plotting
    iPlot = TMSiSAGA.ImpedancePlot([lay_out ' plot Impedance'], channel_config.uni, channel_names, SAGA_type);
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
                   
            % Set the type of layout (head_layout or grid_layout) that is executed
            if strcmp(lay_out, 'head')
                iPlot.head_layout(s, ChanLocs);
            elseif strcmp(lay_out, 'grid')
                iPlot.grid_layout(s);
            end
        end
    end
    
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