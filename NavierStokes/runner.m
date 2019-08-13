addpath('../');

feature('getpid')


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool %('local',9);
end

if 1
    
    RN=[10,50,100,400]
    
    parfor i=1:numel(RN)
        
        FileName = ['RN' num2str(RN(i)) '.txt'];
        fID = fopen(FileName,'w');
        
        
        UsingConvectionTerm=Tools.Enums.Bool.No;
        NavierStokesTime(8,RN(i),UsingConvectionTerm,fID);
        
        
        UsingConvectionTerm=Tools.Enums.Bool.Yes;
        NavierStokesTime(8,RN(i),UsingConvectionTerm,fID);
        
        fclose(fID);
    end
    
else
    NavierStokesTime;
end
