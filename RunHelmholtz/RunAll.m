mPath = pwd;
Seps = strfind(pwd,filesep)
mPath(Seps(end):end)='';
addpath(mPath)
feature('getpid')


RunExterior;
RunInterior;
RunInteriorHomo;
RunSimpleTransReflAboutCircle;
RunSimpleTransReflAboutElps;
RunTransReflAboutElps;
RunTransReflAboutStarShapedBody;
RunNested;
RunNestedHomo;
RunSimpleMultScat;
RunMultScat;