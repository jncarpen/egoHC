function TESTparseFn(input1,varargin)
p = inputParser;
p.addRequired('mainVector',@(x) length(x)>1);
p.addOptional('ntimes',1,@isscalar);
p.addParamValue('title','Default title',@isstr);

p.parse(input1,varargin{:})