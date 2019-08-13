addpath('../');

feature('getpid')


p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    parpool %('local',9);
end

if 1
    ExactChoice = [NavierStokesExact.Exact1Time,NavierStokesExact.Exact2Time];
    RN=[10,50,100,200]
    
    for ECn=1:numel(ExactChoice)
        
        EC = ExactChoice(ECn);
        
        parfor i=1:numel(RN)
            
            FileName = ['RN' num2str(RN(i)) 'EC' num2str(ECn) '.txt'];
            fID = fopen(FileName,'w');
            
            
            UsingConvectionTerm=Tools.Enums.Bool.No;
            NavierStokesTime(8,RN(i),UsingConvectionTerm,fID,EC);
            
            
            UsingConvectionTerm=Tools.Enums.Bool.Yes;
            NavierStokesTime(8,RN(i),UsingConvectionTerm,fID,EC);
            
            fclose(fID);
        end
    end
    
else
    NavierStokesTime;
end
