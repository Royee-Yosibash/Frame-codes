function projectRoot = setup()
%SETUP Add the active Matlab project folders to the MATLAB path.
%   The path is resolved from this file so setup does not depend on the
%   current working directory. Existing scripts retain their original
%   working-directory behavior.

projectRoot = fileparts(mfilename('fullpath'));

addpath(projectRoot);
addpath(fullfile(projectRoot, 'FramesTOOLBOX'));
addpath(fullfile(projectRoot, 'GUI'));
addpath(fullfile(projectRoot, 'Coding Scheme'));
addpath(fullfile(projectRoot, 'Results'));

end
