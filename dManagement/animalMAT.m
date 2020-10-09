% make matfile for an animal (AXONA ONLY)

% get length of sessions
totalSessions = length(animal_tfilelist);
totalEEGChannels = 4;

% create cell array, where each cell is a session.
cellSpikes = cell(1,totalSessions);
unitID = cell(1,totalSessions);
eeg_and_time = cell(1,totalEEGChannels);
eeg = cell(1,totalSessions);
posLED = cell(1,totalSessions);
pixLED = cell(1,totalSessions);

% loop through the filelists for each session
for sessionNum = 1:totalSessions
    [S, ID] = LoadSpikes(animal_tfilelist{1,sessionNum});
    eegSession = animal_EEGfilelist{1,sessionNum};
    [led_pos, led_pix, header] = read_pos_file(animal_posfilelist{1,sessionNum});
        
    cellSpikes{1, sessionNum} = S;
    unitID{1,sessionNum} = ID;
    posLED{1,sessionNum} = led_pos;
    pixLED{1,sessionNum} = led_pix;
    
    for eegChan = 1:totalEEGChannels
        eeg_and_time{1,eegChan} = read_eeg_file(eegSession{1,eegChan});
    end
    
    EEG{1,sessionNum} = eeg_and_time;
end



% save here-- D:\Data\24116\matFiles