clc;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(projectRoot);
setup();
cd(fullfile(projectRoot, 'GUI'));
global AllFrameParameters
global AllFramesTable
AllFrameParameters = [];    AllFramesTable = [];

FrameAnalyzer;
