% calls: 
%trySARfromMedvin(nan, 'circle_gp10_')
%trySARfromMedvin(nan, 'circle_k400n_')
%trySARfromMedvin(nan, 'Circler09Shifted03_gp11_')
%trySARfromMedvin(nan, 'Circler09Shifted03_gp11_kf_half_')

function trySARfromMedvin(irs,inMatFilePrefix)

    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 

    restoredefaultpath;
    addpath('../minus_i')   
    %addpath('../plus_i')   
    addpath('./dpm--matlab-Radar'); 
    
    if ~exist('inMatFilePrefix', 'var') 
        %inMatFilePattern = 'ellipse_medvinScatt_';
        %inMatFilePrefix = 'circle_medvinScatt_';  spec = getSpec();
        
        inMatFilePrefix = 'circle_gp10_'; spec = getSpec();
        
        %inMatFilePrefix = 'ellipse_gp11_'; spec = getSpec_2K(); % too many ambiguities     
        % inMatFilePrefix = 'circle_gp11_'; spec = getSpec_2K();
    else
        fprintf('Reading spec from a large file\n'); 
        tic
        load([inMatFilePrefix 'ind_a.mat'],'spec');
        toc
        clear irs
    end
    
    if exist('irs', 'var') 
        processIndexedNum(inMatFilePrefix, irs, spec)
    else
        processIndexedAlph(inMatFilePrefix, spec)        
    end
%{
    inMatFile = '../Analytic/AnalyticEllipse.mat'; 
    % spec should be read from the matfile; used to be: spec = getSpec(); 
    load(inMatFile); 
    plotSAR_from_matfile(inMatFile, spec);
%} 
end
    
function processIndexedAlph(inMatFilePrefix, spec)
 
    
    ind_letters = {'a'}; %{'a', 'b', 'c'}; 

    ind_pattern = 'ind_%s.mat'; 
    scatterer_pattern = 'iface.mat'; 
    
   % addpath('../matlab')  
        
    figure('units', 'normalized', 'position', [0.02 0.1, 0.6, 0.8], 'paperpositionmode', 'auto'); 
    
    assert(numel(ind_letters) <= 3);   
  
    for ind = 1:numel(ind_letters)

        %subplot(2, 2, ind); 
        
        inMatFile = sprintf([inMatFilePrefix, ind_pattern], ind_letters{ind});    
        plotSAR_from_matfile(inMatFile, spec); 
    end
       
    
    inMatFile = [inMatFilePrefix, scatterer_pattern]; 
    
    if exist(inMatFile, 'file')
    %    subplot(2, 2, 4);     
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

    [y2D.y1, y2D.y2] = meshgrid(-2.5 : 0.02 : 2.8, -2.3 : 0.02 : 2.5);      
    

    k_band = spec.k_range; 
    
    %k_band = k_band(1 : round(numel(k_band) / 2) ); 
   
    %phi_inc = spec.phi_range(abs(spec.phi_range) < 0.05 * pi); fprintf('Warning: modified range of phi!\n'); 
    %phi_inc = spec.phi_range(spec.phi_range > 0.1 * pi); fprintf('Warning: modified range of phi!\n');     
    phi_inc = spec.phi_range(spec.phi_range > 0.16 * pi); fprintf('Warning: modified range of phi!\n');    
    %phi_inc = spec.phi_range(spec.phi_range < -0.16 * pi); fprintf('Warning: modified range of phi!\n');     
    %phi_inc = spec.phi_range(spec.phi_range < 0); fprintf('Warning: modified range of phi!\n'); 
    
    %phi_inc = spec.phi_range(46 : 56 ); fprintf('Warning: modified range of phi!\n'); 
    %phi_inc = spec.phi_range(1 : 26 ); fprintf('Warning: modified range of phi!\n'); 
    
    %phi_inc = spec.phi_range; % full aperture

    
    tic
    load(inMatFile); % created by: save(outMatFile, 'curve', 'scattFields')
    toc   
    
    isShiftedCircle = strcmp(inMatFile, 'Circler09Shifted03_gp11_iface.mat') ... 
                   || strcmp(inMatFile, 'Circler09Shifted03_gp11_kf_half_iface.mat'); 
    
    if isShiftedCircle 
        warning('HACK: correcting interface for a shifted case '); 
        curve.z(1,:) = curve.z(1,:) + 1/3; 
        curve.z(2,:) = curve.z(2,:) + 1/3; 
    end
    % BEFORE EXPERIMENT  
    phi_refl = phi_inc + pi; 
 
    % AFTER EXPERIMENT
    % phi_refl = phi_inc; 
 
    uinf = nan(numel(phi_refl), numel(k_band));     
    for im = 1:numel(phi_refl)
        
        if rem(im, 40) == 1
            fprintf('im = %d out of %d\n', im, numel(phi_refl)); 
        end
        
        phi_val = phi_refl(im); 
        for il = 1:numel(k_band)
            k = k_band(il); 
            scattField = scattFields(im,il); 
            %traj = SARUtils.doFFP_curve(k, phi_val, curve, scattField);  
            traj = mi__SARUtils.doFFP_curve(k, phi_val, curve, scattField); 
            assert(numel(traj.ffp) == 1); 
            uinf(im, il) = traj.ffp; 
        end
    end

    %I = SARUtils.buildSARimage(y2D.y1, y2D.y2, k_band, phi_refl, uinf); 
    I = mi__SARUtils.buildSARimage(y2D, k_band, phi_refl, uinf); 

    pcolor(y2D.y1, y2D.y2, log10(abs(I))); 
 %   warning('restore log!')
 %   pcolor(y1, y2, (abs(I))); 
        
    hc = colorbar; shading flat;
    ylabel(hc, 'log10(magnitude)')
    daspect([1 1 1])
    xlabel('x'); ylabel('y')    
    if exist('r', 'var')
        hold on; plot(r * cos(spec.theta), r * sin(spec.theta), 'w--', 'LineWidth', 3)    
        hold on; plot(1.5 * r * cos(phi_refl), 1.5 * r * sin(phi_refl), 'w-', 'LineWidth', 10)
    
        titleStr = sprintf('SAR image; r = %7.4f \n min = %g, max = %g \n fname = %s', ... 
                           r, min(abs(I(:))), max(abs(I(:))), inMatFile); 
        title(titleStr, 'interpreter', 'none'); 
    end
    
    if isShiftedCircle  
        r2=0.9;
        th = 0:0.02:2*pi;
        hold on; plot(r2*cos(th)+1/3,r2*sin(th)+1/3, 'm--');
    end
    
    %print('-djpeg', 'test1_SAR_from_numerics.jpg')
end 

 function spec = getSpec_2K()

     spec.k_range = 200:5:220; 
     spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
     spec.theta = pi * (0 : 0.0025 : 1.999); 
 
 end
% 
% function spec = getSpec()
% 
%     spec.k_range = 50:0.5:55; 
%     spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
%     spec.theta = pi * (0 : 0.0025 : 1.999); 
% 
% end