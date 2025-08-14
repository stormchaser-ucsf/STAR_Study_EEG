%CHANGEDATARECORDERINTERFACE - Change data recorder interface
%   Shows how to change the interface being used by the data recorder.

% If the script has been terminated manually, it can occur that the device
% keeps sampling. To be sure, a lib.cleanUp() is performed prior to running
% the example. 
if exist('lib', 'var')
    disp('[N/A] Ensure library is cleaned...')
    lib.cleanUp();
end

% Initialize the library
lib = TMSiSAGA.Library();

% Code within the try-catch to ensure that all devices are stopped and 
% closed properly in case of a failure.
try
    % Get a single device from the connected devices
    device = lib.getFirstAvailableDevice('usb', 'electrical');
    
    % Open a connection to the device
    device.connect();
	   
    % Change interface for data recorder to wifi
    device.changeDataRecorderInterfaceTo('wifi');
       
    % Wait for the device to be ready
    disp('Wait for the device to be connected to wifi.');
    
    % Get device by wifi
    device = lib.getFirstAvailableDevice('usb', 'wifi', 10);
       
    % Open a connection to the device
    device.connect();
       
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

