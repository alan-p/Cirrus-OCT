function [ data ] = loadPatientSurfaces( varargin )
%loadLayer loads all layer data from an OCT segmentation file
%   Detailed explanation goes here
% INPUTS:
% [varargin] optional - path to folder containing segmentation results
%% Browse for a file
if nargin < 1
    [fname,path]=uigetfile('*.txt');
else
    path=varargin{1};
    fname = dir(fullfile(path,'*surf0.txt'));
    fname=fname.name;
end

basename=fname(1:(length(fname)-5));

surfaces=nan(128,512,11);

for iSurface = 0:10
    newname=[basename,num2str(iSurface),'.txt'];
    surface=readtable(fullfile(path,newname),'Delimiter',' ',...
        'ReadRowNames',false,'ReadVariableNames',0);
    % Reshaping to remove column 512 (NaN) and change to array
    stmp=table2array(surface(:,1:512));
    %coordinates where the OCT is off the screen are very large,
    %find these and replace with a median value from the same column
    %used column not row so don't pass through fovea
    largeVals = find(stmp>65500); 
    [~,y]=ind2sub(size(stmp),largeVals);
    repVals=round(arrayfun(@(x) median(stmp(:,x)),y));
    stmp(largeVals)=repVals;
    
    surfaces(:,:,(iSurface+1))=stmp;
end

data=surfaces;
end

