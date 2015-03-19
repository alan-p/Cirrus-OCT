%% Cirrus_OCT_thicknessIRL
% This function generates an unflattened inner retinal layer 
% (IRL=GCL+IPL+INL), called "thicknessIRL". The dimensions are interpolated
% to achieve 512x512 square shape versus rectangular 512x128.

function [thicknessIRL]=Cirrus_OCT_thicknessIRL(varargin)
%% Browse for a file
if nargin < 1
    [fname,path]=uigetfile('*.txt');
else
    path=varargin{1};
    fname = dir(fullfile(path,'*surf0.txt'));
    fname=fname.name;
end

%% Generating boundaries (surface) of IRL
% Extract the first part of file name
basename=fname(1:(length(fname)-5));

% Create interpolated 512x512 INL-OPL surface (surface5)
newname=[basename,num2str(4),'.txt'];
% Reads file and converts to table
surface5=readtable(fullfile(path,newname),'Delimiter',' ',...
    'ReadRowNames',false,'ReadVariableNames',0);
% Reshaping to remove column 512 (NaN) and change to array
surface5=table2array(surface5(:,1:512));
% Interpolating to create 512x512 square matrix
%surface5=imresize(surface5,[512 512],'bilinear');

% Create interpolated 512x512 NFL-GCL surface (surface2)
newname=[basename,num2str(1),'.txt'];
% Reads file and converts to table
surface2=readtable(fullfile(path,newname),'Delimiter',' ',...
    'ReadRowNames',false,'ReadVariableNames',0);
% Reshaping to remove column 512 (NaN) and change to array
surface2=table2array(surface2(:,1:512));
% Interpolating to create 512x512 square matrix
%surface2=imresize(surface2,[512 512],'bilinear');

%% Creating thicknessIRL
% Subtract surface2 from surface5 to generate thicknessIRL 
thicknessIRL=surface5-surface2;

end