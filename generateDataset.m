%% GENERATE DATASETS

% This script will generate and store dataset (.mat file) for a single
% animal and store it in specified file location.
% jCarpenter


%%  Make sure everything is on the path
addpath(genpath("C:\Users\17145\Downloads\MClust-4.4")) % MClust for LoadSpikes.m
addpath(genpath("C:\Users\17145\Documents\github_local\MATLAB\moser_matlab\OVC\bnt-20190903T101355Z-001")); %BNT
addpath(genpath("C:\Users\17145\Documents\github_local\egoHC")); % project code

% specify storage location (should be the same for all animals)
storeLoc = "D:\Data\Dataset";

%% Run sessionDirectory.m function
% Fill this section out with relevant information for the animal of interest 
recFolderPath = 'D:\Data\24116\Axona recordings'; % change this for each animal (or make it into a function)
[animal_trackfilelist, animal_fmEvents, animal_behaviourlist, animal_posfilelist, animal_tfilelist, animal_EEGfilelist] = sessionDirectory(recFolderPath);
LabNotes_flnm = "D:\Data\Labnotes\24116.xlsx";

%% Run animalMAT.m script

% get length of sessions
totalSessions = length(animal_tfilelist); % grab # of sessions
totalEEGChannels = length(animal_EEGfilelist{1,1}); % grab # EEG channels

% create cell arrays, where each cell is 1 session
spikeTimes = cell(1,totalSessions);
unitID = cell(1,totalSessions);
eeg_and_time = cell(1,totalEEGChannels);
eeg = cell(1,totalSessions);
pos = cell(1,totalSessions);
speed = cell(1,totalSessions);
accel = cell(1,totalSessions);
hd = cell(1,totalSessions);
sessInfo = cell(1, totalSessions);
behaviour = cell(1, totalSessions);
fmEvents = cell(1,totalSessions);

% loop through the filelists for each session
for sessionNum = 1:totalSessions
    [spikeTimes{1, sessionNum}, unitID{1,sessionNum}] = LoadSpikes(animal_tfilelist{1,sessionNum});
    eegSession = animal_EEGfilelist{1,sessionNum};
    pullPos = io.axona.getPos(animal_posfilelist{1,sessionNum}{1,1}); % from BNT* /interpolate NaNs?
    [xInt1, yInt1] = general.interpolatePositions(pullPos(:,1), [pullPos(:,2), pullPos(:,3)]); % interpolate positions for LED1
    [xInt2, yInt2] = general.interpolatePositions(pullPos(:,1), [pullPos(:,4), pullPos(:,5)]); % interp positions for LED2
    intPos = [pullPos(:,1), xInt1, yInt1, xInt2, yInt2]; % squash it back together
    fixedLEDpos = general.fixLedAssignment(intPos); % Fix point assignment to LEDs in tracked positions (computationally heavy)
    
    smoothPos(1:length(fixedLEDpos),1) = fixedLEDpos(:,1); % we don't want to smooth time vector
    for col = 2:5
        smoothPos(1:length(fixedLEDpos),col) = general.smoothGauss(fixedLEDpos(:,col), 5); % smooth positions 
    end
    pos{1,session} = smoothPos; % save position data for the session
    
    for eegChan = 1:totalEEGChannels
       eeg_and_time{1,eegChan} = read_eeg_file(eegSession{1,eegChan}); % read in all (4) EEG channels
    end
    
   EEG{1,sessionNum} = eeg_and_time; %% [eeg time]; units? --> match with timestamps *
   
  %% speed + accel + hd
  numLeds = 2; 
  t = smoothPos(:,1); % in seconds
  x = smoothPos(:,2); 
  x2 = smoothPos(:,4);
  y = smoothPos(:,3); 
  y2 = smoothPos(:,5);
  v = zeros(size(smoothPos,1), numLeds); % velocity
  a = zeros(size(smoothPos,1), numLeds); % acceleration 
 
  for i = 2:size(smoothPos,1)-1
    v(i, 1) = sqrt((x(i+1) - x(i-1))^2 + (y(i+1) - y(i-1))^2) / (t(i+1) - t(i-1));
    v(i, 2) = sqrt((x2(i+1) - x2(i-1))^2 + (y2(i+1) - y2(i-1))^2) / (t(i+1) - t(i-1));
  end
  
  % pad the vector
  v(1,1) = v(2,1);
  v(1,2) = v(2,2);
  v(end, 2) = v(end-1, 2);
  v(end, 1) = v(end-1, 1);

  v1 = v(:,1);
  v2 = v(:,2);
  for i = 2:size(v,1)-1
      a(i,1) = (v1(i+1) - v1(i-1))/(t(i+1) - t(i-1));
      a(i,2) = (v2(i+1) - v2(i-1))/(t(i+1) - t(i-1));
  end
         
  speed{1,sessionNum} = v; 
  accel{1,sessionNum} = a;
  hd{1,sessionNum} = rem(atan2d(y2-y, x2-x) + 180, 360); % HD in DEGREES
  
  %% Session Information
  % figure out how you want to deal with setfile info- is there an easier
  % way to read this all in?
  
  info = cell(2,17);
  [header, ~] = read_binary_data(animal_posfilelist{1,sessionNum}{1,1},'pos'); % grab position header
  info{1,1} = "EEGsamplesPerPos";
  info{2,1} = header.KeyValue("EEG_samples_per_position");
  info{1,2} = "sampleRate";
  info{2,2} = header.KeyValue("sample_rate");
  info{1,3} = "bytesPerTS";
  info{2,3} = header.KeyValue("bytes_per_timestamp");
  info{1,4} = "bytes_per_coord";
  info{2,4} = header.KeyValue("bytes_per_coord");
  info{1,5} = "pixels_per_metre";
  info{2,5} = header.KeyValue("pixels_per_metre");
  info{1,6} = "num_pos_samples";
  info{2,6} = header.KeyValue("num_pos_samples");
  info{1,7} = "duration";
  info{2,7} = header.KeyValue("duration");
  info{1,8} = "trial_date";
  info{2,8} = header.KeyValue("trial_date");
  info{1,9} = "trial_time";
  info{2,9} = header.KeyValue("trial_time");
  info{1,10} = "window_min_x";
  info{2,10} = header.KeyValue("window_min_x");
  info{1,11} = "window_max_x";
  info{2,11} = header.KeyValue("window_max_x");
  info{1,12} = "window_min_y";
  info{2,12} = header.KeyValue("window_min_y");
  info{1,13} = "window_max_y";
  info{2,13} = header.KeyValue("window_max_y");
  info{1,14} = "min_x";
  info{2,14} = header.KeyValue("min_x");
  info{1,15} = "max_x";
  info{2,15} = header.KeyValue("max_x");
  info{1,16} = "min_y";
  info{2,16} = header.KeyValue("min_y");
  info{1,17} = "max_y";
  info{2,17} = header.KeyValue("max_y");
  sessInfo{1,sessionNum} = info;
  
  %% For stuff JS already generated 
  
  % for fmEvents.mat
  if animal_fmEvents{1,session}{1,1} == "OFS"
        fmEvents{1,session} = 'OFS';
  else
      eventStruct = load(animal_fmEvents{1,session}{1,1});
      fmEvents{1,session} = eventStruct;
  end
  
  % for behaviour.mat
  if isempty(animal_behaviourlist{1,session}{1,1})
      behaviour{1,session} = 'Empty';
  else
      behaviourMat = load(animal_behaviourlist{1,session}{1,1});
      behaviour{1,session} = behaviourMat;
  end
  
  
end % end of sessions loop

  [unitInfo] = pullLabNotes(LabNotes_flnm, 4);
  [uniqueID, neuronID, sessType] = getUniqueID(unitInfo, t64List);

    %% Bin stuff
    spikeTrain = cell(1, totalSessions);
    for session = 1:totalSessions
        if ~isempty(pos{1,session}) && ~isempty(spikeTimes{1,session}) % do we have everything we need?
            spikeTrain{1,session} = binSpikes(pos{1,session}(:,1), spikeTimes{1,session}); 
        else
            spikeTrain{1,session} = []; % no tracker file/spikeTimes for that session...
        end
    end

%% Save all variables to a single .matfile

% unitInfo, uniqueID, neuronID, sessType, spikeTrain, 
