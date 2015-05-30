%addpath ('/home/1016/ma/mmedvin/Documents/MATLAB/SVN/DPM/');

mPath = pwd;
Seps = strfind(pwd,filesep)
mPath(Seps(end):end)='';
addpath(mPath)
feature('getpid')


%RunLaplacian346;
%fprintf('\n\n');
%RunLaplacian351;


%RunExterior;

RunInteriorHomo;


%RunExtSubmarine;
%RunIntHomoSubmarine;
%RunTransReflAboutSubmarine;


%RunTransReflAboutStarShapedBody;