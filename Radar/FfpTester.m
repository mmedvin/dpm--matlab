function FfpTester

    k             = 10;
    
    BC          = Tools.Enums.BoundaryConditions.Dirichlet;%;Neumann;%
    
    for IncidentAngle = [30:20:120]*pi/180
        uc = ScatteringProblem(k,IncidentAngle,BC,Tools.Enums.Scatterer.Circle);
        ue = ScatteringProblem(k,IncidentAngle,BC,Tools.Enums.Scatterer.Ellipse);
        us = ScatteringProblem(k,IncidentAngle,BC,Tools.Enums.Scatterer.StarShaped);
    end
end