function [info] = getSessInfo(posFilePath)
%SESSINFO Reads session information from .pos file into a struct, info.
%   INPUT:
%   posFilePath:    File path for .pos file for a single session, 
%                   for example: posFilePath = animal_posfilelist{1,1}{1,1}; 
%                   animal_posfilelist can be generated with the sessionDirectory.m
%                   function.

%   OUTPUT: 
%   info:           Struct with a bunch of fields. (Haven't decided whether to
%                   add/remove any of the current fields.)

%   DEPENDENCIES:
%   mTintCore:      A set of functions for loading *Axona* data and cut files
%                   into MATLAB (by Daniel Manson, UCL). 
%                   Can be downloaded at: https://wiki.ucl.ac.uk/display/Hippo/mTintCore

%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath("C:\Users\17145\Documents\github_local\egoHC\mTintCore")); % Add mTintCore to path

[header, ~] = read_binary_data(posFilePath,'pos'); % grab position header

    info.EEGsamplesPerPos =     header.KeyValue("EEG_samples_per_position");
    info.sampleRate =           header.KeyValue("sample_rate");
    info.bytesPerTS =           header.KeyValue("bytes_per_timestamp");
    info.bytes_per_coord =      header.KeyValue("bytes_per_coord");
    info.pixels_per_metre =     header.KeyValue("pixels_per_metre");
    info.num_pos_samples =      header.KeyValue("num_pos_samples");
    info.duration =             header.KeyValue("duration");
    info.trial_date  =          header.KeyValue("trial_date");
    info.trial_time =           header.KeyValue("trial_time");
    info.window_min_x =         header.KeyValue("window_min_x");
    info.window_max_x =         header.KeyValue("window_max_x");
    info.window_min_y =         header.KeyValue("window_min_y");
    info.window_max_y =         header.KeyValue("window_max_y");
    info.min_x =                header.KeyValue("min_x");
    info.max_x =                header.KeyValue("max_x");
    info.min_y =                header.KeyValue("min_y");
    info.max_y =                header.KeyValue("max_y");     
end

