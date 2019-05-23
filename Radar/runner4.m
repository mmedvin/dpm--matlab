%addpath ('/home/1016/ma/mmedvin/Documents/MATLAB/SVN/DPM/');
addpath ('/home/mmedvin/SVN/DpmRadar/');

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




%GridParam=11;
%rchoice=[80,320,800];
%ParSAR([pwd filesep 'SAR_k430_k440_Circle_r1_gp11.mat'],1,GridParam,ScattererType, 430:440, Phi,rchoice);


GridParam=11;
rchoice=[10,20,50];
ParSAR([pwd filesep 'SAR_k429.5_k434_Circle_r09_shft03_gp11.mat'],1,GridParam,'shifted', 429.5:0.5:434, Phi,rchoice);
%ParSAR([pwd filesep 'SAR_k434.5_k440_Circle_r09_shft03_gp11.mat'],1,GridParam,'shifted', 434.5:0.5:440, Phi,rchoice);







