%% Cirrus_OCT_pathlist
% This script creates a variable, called "pathlist", that allows multiple
% participants to be run subsequently and automatically. To add a
% participant, list their 'SubID/subfolder'. Note: The subfolder must 
% contain the corresponding surf0.txt file.

% Address
basepath = 'x:/Diabetes OCT Segmented data/Preliminary Data';
% Directory of participants
pathlist = {'204/E542',...
           '304/E472',...
           '310/E507',...
           '346/E572',...
           '451/E973',...
           '497/E683',...
           '498/E892',...
           '499/E524',...
           '504/E327',...
           '507/E849',...
           '508/E001',...
           '521/E525',...
           '522/E462'};