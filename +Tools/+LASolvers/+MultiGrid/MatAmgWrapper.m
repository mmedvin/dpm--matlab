classdef MatAmgWrapper < handle
	%UNTITLED Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		MyPath='';
		OldPath;
		Options;
		A;
		fid=-1;
	end
	
	methods
		
		function obj = MatAmgWrapper(Params)
			if isunix
				obj.MyPath = [pwd filesep '..' filesep 'ThirdParties' filesep 'bin' filesep 'matamg' filesep 'Linux' ];
			else
				obj.MyPath = [pwd filesep '..' filesep 'ThirdParties' filesep 'bin' filesep 'matamg' ];
			end
			if ~exist([obj.MyPath filesep 'amg.m'],'file')
				error('The path to the MATAMG doesn''t exists');
			end
			

			
			%
			%
			% 			obj.OldPath = 
			addpath(obj.MyPath);
			% 			newPath = path;
			% 			if ~strcmp(oldPath ,newPath)
			% 				obj.MyPath = Path; %save the
			% 			end
			
			obj.RegisterInstance();
			
			obj.Options = amgset;
            obj.Options = amgset(obj.Options,'PrintOnScreen','off', 'MaxCycle',Params.MaxCycle,'PreCond','pcg','TolAMG', Params.tollerance);
			%'IntpType','lramg');%, 'SaveCsn','off','SaveIntp','off','Log','off');

			obj.A = Params.A;
		end
		
		function Instance = RegisterInstance(obj)
			obj.fid = fopen([obj.MyPath filesep 'matamg.semafore'],'w+');
			instances = fread(obj.fid,[1,1]);
			if isempty(instances)
				instances = 0;
			end
			
			Instance = instances + 1;
			fwrite(obj.fid,Instance);
			
		end
		
		function delete(obj)
			%destructor
			if obj.fid > 0
				instances = fread(obj.fid,[1,1]);
				if instances == 1
					fclose(obj.fid);
					delete([obj.MyPath filesep 'matamg.semafore']);
					rmpath(obj.MyPath);
				else
					fwrite(obj.fid,instances - 1);
					fclose(obj.fid);
				end				
			end
		end
		
		function Result = Solve(obj,Rhs)
			Result = amg(obj.A,rand(size(Rhs)),Rhs,obj.Options);
		end
	end
	
end

