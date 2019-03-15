%addpath ('/home/1016/ma/mmedvin/Documents/MATLAB/SVN/DPM/');
addpath ('../');

feature('getpid')


%RunLaplacian346;
%fprintf('\n\n');
%RunLaplacian351;


%RunExterior;

%RunInteriorHomo;


%RunExtSubmarine;
%RunIntHomoSubmarine;
%RunTransReflAboutSubmarine;


%RunTransReflAboutStarShapedBody;

%%%%%%%%%%%%%%%%%%%

%dPhi=0.004;
dK = 0.5;
%Phi = pi*(-0.2:dPhi:0.2);%0;%-10:5:10;% 
K0 = 50:dK:55;%5;%5:7;%50;%
Phi =linspace(-0.2*pi,0.2*pi,400);        
ScattererType = 'ellipse'; %'circle';%'star' , 'ellipse'
%AR=2;
%filename = [pwd filesep 'SAR.mat'];


GridParam=8;
rchoice=[10,40,100];
%SARsimulation([pwd filesep 'SAR_Ellipse_a1_b2_gp8.mat'],2,GridParam,ScattererType, K0, Phi,rchoice);
%SARsimulation([pwd filesep 'SAR_Circle_r1_gp8.mat'],1,GridParam,ScattererType, K0, Phi,rchoice);

%GridParam=10;
%rchoice=[40,160,400];
%SARsimulation([pwd filesep 'SAR_Ellipse_a1_b2_gp10.mat'],2,GridParam,ScattererType, K0, Phi,rchoice);
%SARsimulation([pwd filesep 'SAR_Circle_r1_gp10.mat'],1,GridParam,ScattererType, K0, Phi,rchoice);




% dK = 5;
%K0 = 200:dK:220;        
%GridParam=11;
%rchoice=[80,320,800];
%SARsimulation([pwd filesep 'SAR_k200_Ellipse_a1_b2_gp11.mat'],2,GridParam,ScattererType, K0, Phi,rchoice);
%SARsimulation([pwd filesep 'SAR_k200_Circle_r1_gp11.mat'],1,GridParam,ScattererType, K0, Phi,rchoice);

      
%GridParam=11;
%rchoice=[80,320,800];
%ParSAR([pwd filesep 'SAR_k200_k204_Circle_r1_gp11.mat'],1,GridParam,ScattererType, 200:204, Phi,rchoice);



%GridParam=11;
%rchoice=[80,320,800];
%ParSAR([pwd filesep 'SAR_k400_k409_Circle_r1_gp11.mat'],1,GridParam,ScattererType, 400:409, Phi,rchoice);


%GridParam=11;
%rchoice=[10,20,50];
%ParSAR([pwd filesep 'SAR_k400_k409_Circle_r09_shft03_gp11.mat'],1,GridParam,'shifted', 400:0.5:409, Phi,rchoice);


GridParam=11;
rchoice=[10,20];
K0=400:0.5:440;
FolderName = [pwd filesep 'SAR_k' num2str(K0(1)) '_k' num2str(K0(end)) filesep 'gp' num2str(GridParam) ]
if 7~=exist(FolderName,'dir'), mkdir(FolderName); end

%ParSAR([FolderName filesep 'Circle_r09_shft03_gp11_kf1.1.mat'],1,GridParam,'shifted', 400:0.5:440, Phi,rchoice,1.1);


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
   parpool%('local',9);
end

%ParSAR([FolderName filesep 'Circle_r09_Shft1Slash3_kf_half'],1,GridParam,'shifted', 400:0.5:440, Phi,rchoice,1/2);

GridParam=8;

if 7~=exist([pwd filesep 'tmp'],'dir'), mkdir([pwd filesep 'tmp']); end


P=['tmp' filesep 'Circle_r09_shft03_kf1_1']
%profile on
ParSAR(['tmp' ],1,GridParam,'shifted', 400:410, Phi,rchoice,1.1);

%profsave
%profile off











