classdef AbstractLA_Solver < handle
	%UNTITLED2 Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
	end
	
	methods(Abstract)
		Solve(obj,rhs);		
	end
	
end

