function trySARfromMedvin(irs)
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    if ~exist('irs', 'var'), irs = 10; end
    
    %inMatFilePattern = 'circle_medvinScatt_ir_%d.mat';         
    inMatFilePattern = 'medvinScatt_ir_%d.mat';          

    addpath('../matlab')  
    
%{
    figure; 
    ir = irs(1); 
    inMatFile = sprintf(inMatFilePattern, ir);    
    plotSAR_for_ir(ir, inMatFile); 
%}    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.6, 0.8], 'paperpositionmode', 'auto'); 
    
    assert(numel(irs) <= 4); 
   
    
    for iir = 1:numel(irs)
        ir = irs(iir); 
        subplot(2, 2, iir); 
        
        inMatFile = sprintf(inMatFilePattern, ir);    
        plotSAR_for_ir(ir, inMatFile); 
    end
end

function plotSAR_for_ir(ir, inMatFile)    
    
 
    
    [y1, y2] = meshgrid(-1.5 : 0.02 : 1.8, -1.3 : 0.02 : 1.5);      
    
    spec.k_range = 50 : 0.5 : 55; 
    spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
    spec.theta = pi * (0 : 0.0025 : 1.999); 
    
    phi_range = spec.phi_range; 
    k_band = spec.k_range; 
    
    %k_band = k_band(1 : round(numel(k_band) / 2) ); 
    phi_range = phi_range(46 : 56 ); 
    
    load(inMatFile); % created by: save(outMatFile, 'curve', 'scattFields')
        
    % create field for each k and each direction 
    uinf = nan(numel(phi_range), numel(k_band)); 
    for im = 1:numel(phi_range)
        phi_val = phi_range(im); 
        for il = 1:numel(k_band)
            k = k_band(il); 
            scattField = scattFields(im,il); 
            traj = SARUtils.doFFP_curve(k, phi_val, curve, scattField);  
            assert(numel(traj.ffp) == 1); 
            uinf(im, il) = traj.ffp; 
        end
    end
    
    I = SARUtils.buildSARimage(y1, y2, k_band, phi_range, uinf); 

    pcolor(y1, y2, log(abs(I))); 
    
    hold on; plot(r * cos(spec.theta), r * sin(spec.theta), 'w--', 'LineWidth', 3)    
    hold on; plot(r * cos(spec.phi_range), r * sin(spec.phi_range), 'w-', 'LineWidth', 5)
    
    
    hc = colorbar; shading flat;
    ylabel(hc, 'log(magnitude)')
    daspect([1 1 1])
    xlabel('x'); ylabel('y')    
    titleStr = sprintf('SAR image; ir = %d, r = %7.4f \n min = %g, max = %g \n fname = %s', ... 
                   ir, r, min(abs(I(:))), max(abs(I(:))), inMatFile); 
    title(titleStr, 'interpreter', 'none'); 
    
    %print('-djpeg', 'test1_SAR_from_numerics.jpg')
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

