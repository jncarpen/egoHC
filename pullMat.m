s% Once all data is preprocessed and organized into .mat files, pull the
% files out one-by-one for analysis.

% Adopted from:
% https://matlab.fandom.com/wiki/FAQ#How_can_I_process_a_sequence_of_files.3F

% Specify the folder where the .mat files live
myFolder = 'C:\Users\yourUserName\Documents\My Pictures';
% Check to make sure that folder actually exists.  Throw error if not.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end

% Get a list of all .mat files
filePattern = fullfile(myFolder, '*.mat');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    imageArray = imread(fullFileName);
    imshow(imageArray);  % Display image.
    drawnow; % Force display to update immediately.
end