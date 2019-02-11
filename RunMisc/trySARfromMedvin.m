function trySARfromMedvin(irs,inMatFilePrefix)
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 

    if ~exist('inMatFilePrefix', 'var') 
        %inMatFilePattern = 'ellipse_medvinScatt_';
        
        %inMatFilePrefix = 'circle_medvinScatt_';  spec = getSpec();
        
        inMatFilePrefix = 'circle_gp10_'; spec = getSpec();
        
        % inMatFilePrefix = 'ellipse_gp11_'; spec = getSpec_2K();
        
        %inMatFilePrefix = 'circle_gp11_'; spec = getSpec_2K();
    else
        load([inMatFilePrefix 'ind_a.mat'],'spec');
        clear irs
    end
    
    if exist('irs', 'var') 
        processIndexedNum(inMatFilePrefix, irs, spec)
    else
        processIndexedAlph(inMatFilePrefix, spec)        
    end
    
end
    
function processIndexedAlph(inMatFilePrefix, spec)
 
    
    ind_letters = {'a', 'b', 'c'}; 

    ind_pattern = 'ind_%s.mat'; 
    scatterer_pattern = 'iface.mat'; 
    
   % addpath('../matlab')  
        
    figure('units', 'normalized', 'position', [0.02 0.1, 0.6, 0.8], 'paperpositionmode', 'auto'); 
    
    assert(numel(ind_letters) <= 3);   
    
    for ind = 1:numel(ind_letters)

        subplot(2, 2, ind); 
        
        inMatFile = sprintf([inMatFilePrefix, ind_pattern], ind_letters{ind});    
        plotSAR_from_matfile(inMatFile, spec); 
    end
    
    inMatFile = [inMatFilePrefix, scatterer_pattern]; 
    
    if exist(inMatFile, 'file')
        subplot(2, 2, 4);     
        plotSAR_from_matfile(inMatFile, spec);
    end
end

function processIndexedNum(inMatFilePrefix, irs, spec)
       

    ir_pattern = 'ir_%d.mat'; 
    scatterer_pattern = 'iface.mat'; 
    
    addpath('../matlab')  
        
    figure('units', 'normalized', 'position', [0.02 0.1, 0.6, 0.8], 'paperpositionmode', 'auto'); 
    
    assert(numel(irs) <= 3);   
    for iir = 1:numel(irs)
        ir = irs(iir); 
        subplot(2, 2, iir); 
        
        inMatFile = sprintf([inMatFilePrefix, ir_pattern], ir);    
        plotSAR_from_matfile(inMatFile, spec); 
    end
    
    subplot(2, 2, 4);     
    inMatFile = [inMatFilePrefix, scatterer_pattern]; 
    plotSAR_from_matfile(inMatFile, spec);
end

function plotSAR_from_matfile(inMatFile, spec)    

    [y1, y2] = meshgrid(-2.5 : 0.02 : 2.8, -2.3 : 0.02 : 2.5);      
    
    phi_range = spec.phi_range; 
    k_band = spec.k_range; 
    
    %k_band = k_band(1 : round(numel(k_band) / 2) ); 
    
    phi_range = phi_range(46 : 56 ); fprintf('Warning: modified range of phi!\n'); 
    
    
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

    pcolor(y1, y2, log10(abs(I))); 
    
    hold on; plot(r * cos(spec.theta), r * sin(spec.theta), 'w--', 'LineWidth', 3)    
    hold on; plot(r * cos(phi_range), r * sin(phi_range), 'w-', 'LineWidth', 10)
    
    
    hc = colorbar; shading flat;
    ylabel(hc, 'log10(magnitude)')
    daspect([1 1 1])
    xlabel('x'); ylabel('y')    
    titleStr = sprintf('SAR image; r = %7.4f \n min = %g, max = %g \n fname = %s', ... 
                       r, min(abs(I(:))), max(abs(I(:))), inMatFile); 
    title(titleStr, 'interpreter', 'none'); 
    
    %print('-djpeg', 'test1_SAR_from_numerics.jpg')
end 

% function spec = getSpec_2K()
% 
%     spec.k_range = 200:5:220; 
%     spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
%     spec.theta = pi * (0 : 0.0025 : 1.999); 
% 
% end
% 
% function spec = getSpec()
% 
%     spec.k_range = 50:0.5:55; 
%     spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
%     spec.theta = pi * (0 : 0.0025 : 1.999); 
% 
% end