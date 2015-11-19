mPath = pwd;
Seps = strfind(pwd,filesep)
mPath(Seps(end):end)='';
addpath(mPath)
feature('getpid')


RunLapBL53;
RunLapBL54;

RunLapExterior;
RunLapInterior;

RunLapInterior02;
RunLapInterior03;
RunLaplacian01;

RunLaplacian346;
RunLaplacian351;

