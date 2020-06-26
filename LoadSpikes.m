% LoadSpikes: from AD Redish's mClust program **

function S = LoadSpikes(tfilelist)
%
% S = LoadSpikes(tfilelist)
%
%   inp: tfilelist is a cellarray of strings, each of which is a
% 	tfile to open.  Note: this is incompatible with version unix3.1.
    
% 	out: Returns a cell array such that each cell contains a ts 
% 	object (timestamps which correspond to times at which the cell fired)

% ADR 1998
%  version L4.0
%  status: PROMOTED
%  ADR 2012/12 removed flags.  Now only returns in seconds.

% Modified in the following ways:
%  Added file names in adjacent cells (i.e. S{cell,2})
%  Jordan Carpenter, June 24, 2020

%-------------------
% Check input type
%-------------------
if isa(tfilelist, 'char')
    tfilelist = {tfilelist};
end
if ~isa(tfilelist, 'cell')
   error('LoadSpikes: tfilelist should be a cell array.');
end

nFiles = length(tfilelist);

%--------------------
% Read files
%--------------------

fprintf(2, 'Reading %d files.', nFiles);

% for each tfile
% first read the header, then read a tfile 
% note: uses the bigendian modifier to ensure correct read format.

S = cell(nFiles, 2);
for iF = 1:nFiles
	tfn = tfilelist{iF};
    [fd, fn, xt] = fileparts(tfn);
    
	if ~isempty(tfn)
		tfp = fopen(tfn, 'rb','b');
		if (tfp == -1)
			warning('MClust:LoadSpikes', 'Could not open tfile %s.', tfn);
		end
		
		MClust.ReadHeader(tfp);    
        switch xt
            case {'.raw64'}
                S{iF,1} = fread(tfp,inf,'uint64')*10000;	%read as 64 bit ints                             
            case {'.raw32'}
                S{iF,1} = fread(tfp,inf,'uint32')*10000;	%read as 32 bit ints
            case {'.t64', '.k64'}
                S{iF,1} = fread(tfp,inf,'uint64');	%read as 64 bit ints
            case {'.t32', '.t', '.k32','.k'}
                S{iF,1} = fread(tfp,inf,'uint32');	%read as 32 bit ints             
            otherwise
                error('MClust::LoadSpikes', 'Unknown t-file extension.');
        end
		fclose(tfp);		

		% Set appropriate time units.
        S{iF,1} = S{iF,1}/10000;
		
        % convert to ts object
		S{iF,1} = ts(S{iF,1});
        
        
        % Add file name next in the adjacent cell 
        S{iF,2} = fn;
      
	end 		% if tfn valid
end		% for all files
fprintf(2,'\n');
