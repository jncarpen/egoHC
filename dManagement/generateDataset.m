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
filename = '25398_v3.mat';

%% Run sessionDirectory.m function
% Fill this section out with relevant information for the animal of interest 
recFolderPath = 'D:\Data\25398'; % change this for each animal (or make it into a function)
[animal_trackfilelist, animal_fmEvents, animal_behaviourlist, animal_posfilelist, animal_tfilelist, animal_EEGfilelist] = sessionDirectory(recFolderPath);
LabNotes_flnm = "D:\Data\Labnotes\25398.xlsx";

%% Run animalMAT.m script

% get length of sessions
totalSessions = length(animal_tfilelist); % grab # of sessions
totalEEGChannels = length(animal_EEGfilelist{1,1}); % grab # EEG channels

% create cell arrays, where each cell is 1 session
spikeTimes = cell(1,totalSessions); % sampling rate of spikeTimes?
unitID = cell(1,totalSessions);
eeg_and_time = cell(1,totalEEGChannels); % [time EEG]
rawEEG = cell(1,totalSessions);
pos = cell(1,totalSessions);
speed = cell(1,totalSessions);
accel = cell(1,totalSessions);
hd = cell(1,totalSessions);
sessInfo = cell(1, totalSessions);
behaviour = cell(1, totalSessions);
fmEvents = cell(1,totalSessions);

% loop through the filelists for each session
for sessionNum = 1:totalSessions
    
    % grab session information (from .pos file) --> getSessInfo.m
     posFilePath = animal_posfilelist{1,sessionNum}{1,1};
     [info] = getSessInfo(posFilePath);
     sessInfo{1,sessionNum} = info;
    
    % read in files using path locations specified in the animal_xlists
    [spikeTimes{1, sessionNum}, unitID{1,sessionNum}] = LoadSpikes(animal_tfilelist{1,sessionNum});
    spikeTimes{1, sessionNum} = LoadSpikes(animal_tfilelist{1,sessionNum});
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
    pos{1,sessionNum} = smoothPos; % save position data for the session
    
    for eegChan = 1:totalEEGChannels
       eeg_and_time{1,eegChan} = read_eeg_file(eegSession{1,eegChan}); % read in all (4) EEG channels
    end
    
   rawEEG{1,sessionNum} = eeg_and_time; %% [eeg time]; units? --> match with timestamps *
   
  %% speed + accel + hd
  numLeds = 2; 
  % parse and medfilt pos vector to eliminate major outliers
  t = smoothPos(:,1); % in seconds
  x = medfilt1(smoothPos(:,2)); 
  x2 = medfilt1(smoothPos(:,4));
  y = medfilt1(smoothPos(:,3)); 
  y2 = medfilt1(smoothPos(:,5));
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
  
  %% For stuff JS already generated 
  
  % for fmEvents.mat
  if animal_fmEvents{1,sessionNum}{1,1} == "OFS"
        fmEvents{1,sessionNum} = 'OFS';
  else
      eventStruct = load(animal_fmEvents{1,sessionNum}{1,1});
      fmEvents{1,sessionNum} = eventStruct;
  end
  
  % for behaviour.mat
  if isempty(animal_behaviourlist{1,sessionNum}{1,1})
      behaviour{1,sessionNum} = 'Empty';
  else
      behaviourMat = load(animal_behaviourlist{1,sessionNum}{1,1});
      behaviour{1,sessionNum} = behaviourMat;
  end
  
  
end % end of sessions loop

  [unitInfo] = pullLabNotes(LabNotes_flnm, 4);
  [uniqueID, neuronID, sessType] = getUniqueID(unitInfo, animal_tfilelist);

    %% Bin stuff
    spikeTrain = cell(1, totalSessions);
    for sessionNum = 1:totalSessions
        if ~isempty(pos{1,sessionNum}) && ~isempty(spikeTimes{1,sessionNum}) % do we have everything we need?
            spikeTrain{1,sessionNum} = binSpikes(pos{1,sessionNum}(:,1), spikeTimes{1,sessionNum}); 
        else
            spikeTrain{1,sessionNum} = []; % no tracker file/spikeTimes for that session...
        end
    end
    
    %% Remove cells that don't have a unique cell ID 
    SpikeTrain = cell(1,totalSessions); % capital S
    SpikeTimes = cell(1,totalSessions);
    UniqueID = cell(1,totalSessions);
    UnitID = cell(1,totalSessions);
    
    for sessionNum = 1:totalSessions
        STRN = spikeTrain{1,sessionNum};
        STM = spikeTimes{1,sessionNum};
        U = uniqueID{1,sessionNum};
        ID = unitID{1,sessionNum};
        keep = find(~cellfun(@isempty,U)); % find cells w/ unique IDs
        SpikeTrain{1,sessionNum}  = STRN(keep);
        SpikeTimes{1,sessionNum} = STM(keep);
        UniqueID{1,sessionNum} = U(keep);
        UnitID{1,sessionNum} = ID(keep);
    end

    %% Save all variables to a single .matfile

     save D:\Data\Dataset\24116.mat SpikeTimes SpikeTrain UnitID rawEEG pos speed accel hd sessInfo behaviour fmEvents unitInfo UniqueID neuronID sessType

