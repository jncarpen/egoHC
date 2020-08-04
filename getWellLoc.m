function [hwLoc, rdLoc] = getWellLoc(labNotes, trialType)
%   DESCRIPTION:
%   Find well locations for FM trials (for a single animal).
%
%   INPUTS
%   labNotes:       Rx15 table (R = # rows)
%   trialType:      1xS cell array (S = # sessions)
%
%   OUTPUTS
%   hwLoc:          1xS cell array. For FM sessions, well location for home
%                   well in FM coordinates. For other session types, either
%                   "OF" (open field) or "TS" (training session) will be
%                   specified.
%   rdLoc:          1xS cell array. For FM sessions, a *list* of well
%                   locations will be specified (this will be a 1xT double)
%                   where T is the number of trials in the session.
%
%   Jordan Carpenter, 2020.
%
%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%

hwLoc = cell(1, length(trialType));
rdLoc = cell(1,length(trialType));

for sess = 1:length(trialType)
    if trialType{1,sess} == "FM"
        hwLoc{1,sess} = labNotes.hw(find(labNotes.Sess_ == sess));
        rdWellStr = labNotes.GoalLocations(find(labNotes.Sess_ == sess));
        rdWellVec = str2num(rdWellStr{1,1})'; % convert to vector
        rdLoc{1,sess} = rdWellVec;
    elseif trialType{1,sess} == "OF"
        hwLoc{1,sess} = "OF";
        rdLoc{1,sess} = "OF";
    elseif trialType{1,sess} == "TS"
        hwLoc{1,sess} = "TS";
        rdLoc{1,sess} = "TS";
    else
        hwLoc{1,sess} = "?";
        rdLoc{1,sess} = "?";
    end
end

end

