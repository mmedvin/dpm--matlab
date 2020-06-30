addpath ('../');

feature('getpid')


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
   parpool%('local',9);
end


dK = 0.5;
K0 = 50:dK:55;%5;%5:7;%50;%
Phi =linspace(-0.2*pi,0.2*pi,400);        
ScattererType = 'shifted';%'ellipse'; %'circle';%'star' , 'ellipse'
rchoice = [10,20,30];


GridParam=8;
if 7~=exist([pwd filesep 'tmp'],'dir'), mkdir([pwd filesep 'tmp']); end
%P=['tmp' filesep 'Circle_r09_shft03_kf1_1']


if 0

	shift=[0,0];
	ParSAR(['tmp' ],1,GridParam,ScattererType,shift, K0, Phi,rchoice,1.1);

	shift = [0,1/3];
	ParSAR(['tmp' ],1,GridParam,ScattererType,shift, K0, Phi,rchoice,1.1);


	shift=[1/3,1/3];
	rchoice = [4,5,6];
	ParSAR(['tmp' ],1,GridParam,ScattererType,shift, K0, Phi,rchoice,1.1);
end



path = [pwd filesep 'ResK50K55' filesep 'ScattData'];

newAdaptor([path filesep 'ParSar_Shifted0.00_0.00gp8kf1.1_k50_k55.mat'],[path filesep  'Adapted' filesep 'Shifted0.00_0.00gp8kf1.1_k50_k55']);
newAdaptor([path filesep 'ParSar_Shifted0.33_0.33gp8kf1.1_k50_k55.mat'],[path filesep  'Adapted' filesep 'Shifted0.33_0.33gp8kf1.1_k50_k55']);
newAdaptor([path filesep 'ParSar_Shifted0.00_0.33gp8kf1.1_k50_k55.mat'],[path filesep  'Adapted' filesep 'Shifted0.00_0.33gp8kf1.1_k50_k55']);

