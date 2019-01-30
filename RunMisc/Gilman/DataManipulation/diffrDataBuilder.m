function diffrDataBuilder(irs)

    if ~exist('irs', 'var'), irs = 10; end

    addpath('./dpm--matlab-Radar'); 
    %addpath('../matlab')

    inMatFile = './SAR_Ellipse_a1_b2_gp8.mat'; 
    outMatFilePrefix = 'medvinScatt'; 

%    inMatFile = './SAR_Circle_r1_gp8.mat'; 
%    outMatFilePrefix = 'circle_medvinScatt'; 
    
    
    tic
    load(inMatFile,'M');
    toc

    for iir = 1:numel(irs)
        ir = irs(iir); 
        
        extractDataFor_ir(ir, M, outMatFilePrefix); 
    end
    
end
 
function extractDataFor_ir(ir, M, outMatFilePrefix)

    tol = 1e-9; 
     
    spec.k_range = 50 : 0.5 : 55; 
    spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
    %spec.theta = pi * (0 : 0.0025 : 1.999); 
    
    outMatFile = sprintf('%s%s%d', outMatFilePrefix, '_ir_', ir);  
   
    Grid = M.PlrGrid;
    n = numel(Grid.r) - ir;  
    assert(abs(Grid.r(n)   - Grid.r(n-1) - Grid.dx) < tol); 
    assert(abs(Grid.r(n+1) - Grid.r(n)   - Grid.dx) < tol);     

    r = Grid.r(n);     
    theta = Grid.theta;     
    curve.z =  r * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 
  %{
    % pseudo video
    for il = 1:numel(spec.k_range)    
        for im = 1:numel(spec.phi_range)    

            
            % According to the spec, it should be: 
            % scattData = M.scattData{im, il}
            medvinScattData = M.scattData{il, im}; 
            % --- DAMN!
           

            U = medvinScattData.Extu;
          
            
            figure(3); customPlot(Grid, U); 
            pause
        end
    end    
   %} 
    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
            % According to the spec, it should be: 
            % scattData = M.scattData{im, il}
            medvinScattData = M.scattData{il, im}; 
            % --- DAMN!
            
            k = medvinScattData.k;
            phi = medvinScattData.phi;    
            assert(abs(k   - spec.k_range(il))   < tol); 
            assert(abs(phi - spec.phi_range(im)) < tol);             

            U = medvinScattData.Extu;
            u = U(n,:);
            % SECOND ORDER
            %un = (U(n+1,:) - U(n-1,:)) / (2 * Grid.dx); 
            
            % FOURTH ORDER
            un = (8 * (U(n+1,:) - U(n-1,:)) - (U(n+2,:) - U(n-2,:)) ) / (12 * Grid.dx);             
            
            assert(numel(u)  == numel(theta)); 
            assert(numel(un) == numel(theta)); 
            
            scattFields(im,il).value = u; 
            scattFields(im,il).normal_deriv = un;            
            
        end
    end
    toc 

    save(outMatFile, 'curve', 'scattFields', 'ir', 'r')
end
%{
scattField = 

  struct with fields:

           value: [1×113 double]
    normal_deriv: [1×113 double]

curve

curve = 

  struct with fields:

     z: [2×113 double]
    nz: [2×113 double]

%}