function tryFFPfromMedvin

    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 

    addpath('../minus_i')  
    
    
    inMatFilePrefix = 'circle_gp10_'; spec = getSpec();               % wide PSF
    %inMatFilePrefix = 'circle_gp11_'; spec = getSpec_2K();            % many stripes
    %inMatFilePrefix = 'Circler09Shifted03_gp11_'; spec = nan; 
    %inMatFilePrefix = 'Circler09Shifted03_gp11_kf_half_'; spec = nan; 

    %inMatFilePrefix = 'ellipse_gp11_ind_a.mat'; spec = getSpec_2K(); 
         
    %phiIncToPlot = [-0.2, 0, 0.2] * pi; 
    phiIncToPlot = [-0.1, 0, 0.1] * pi;     
    
    ilToPlot = [1, 3, 5]; 
    
    plotFFPs(inMatFilePrefix, spec, phiIncToPlot, ilToPlot);
end

function plotFFPs(inMatFilePrefix, spec, phiIncToPlot, ilToPlot)

    phi_refl = (0: 0.02: 2) * pi; 

    if contains(inMatFilePrefix, 'ellipse_gp11')
        inMatFile = inMatFilePrefix;  
    else
        if ~isstruct(spec)
            fprintf('Reading spec from a large file\n'); 
            tic
            load([inMatFilePrefix 'ind_a.mat'],'spec');
            toc
            clear irs
        end

        scatterer_pattern = 'iface.mat'; 
        inMatFile = [inMatFilePrefix, scatterer_pattern];     

    end
    
    tic
    load(inMatFile); % created by: save(outMatFile, 'curve', 'scattFields')
    toc       
   
    if    strcmp(inMatFile, 'Circler09Shifted03_gp11_iface.mat') ... 
       || strcmp(inMatFile, 'Circler09Shifted03_gp11_kf_half_iface.mat')
        warning('HACK: correcting interface for a shifted case ');
        curve.z(1,:) = curve.z(1,:) + 1/3; 
        curve.z(2,:) = curve.z(2,:) + 1/3; 
    end
    
    phi_inc = spec.phi_range; 
    k_band = spec.k_range;     

    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.6], 'name', inMatFile);     
    
    
    cmap_lines = colormap(lines); 
    numPhiPlots = numel(phiIncToPlot); 
    for iphiToPlot = 1:numPhiPlots
        
        phiInc_prescribed = phiIncToPlot(iphiToPlot); 
         
        [~, im] = min(abs(phi_inc - phiInc_prescribed));
        if numel(im) > 1 
            im = im(1); 
        end      
        phi_inc_val = phi_inc(im); 
        
        subplot(1, numPhiPlots, iphiToPlot); 
        hil = nan(1, numel(ilToPlot)); 
        for ili = 1:numel(ilToPlot)
            il = ilToPlot(ili);
            
            k = k_band(il); 
            scattField = scattFields(im,il); 
            traj = mi__SARUtils.doFFP_curve(k, phi_refl, curve, scattField);  

            hold on; hil(ili) = plot(phi_refl, abs(traj.ffp), '-.'); 
            set(hil(ili), 'Color', cmap_lines(ili, :));
            set(hil(ili), 'DisplayName', sprintf('k = %6.2f', k));  
            
        end       
        

        hold on; hi = plot(mod(phi_inc_val     , 2*pi), 0, 'bs', 'DisplayName', 'incident');         
        hold on; hs = plot(mod(phi_inc_val + pi, 2*pi), 0, 'ks', 'DisplayName', 'specular');         

        legend([hi, hs, hil], 'location', 'southoutside'); 
        
        ylabel('abs(FFP)')
        title(sprintf('phi inc = %3.1f * pi', phi_inc_val / pi))
    end

end

function spec = getSpec_2K()

     spec.k_range = 200:5:220; 
     spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
     %spec.theta = pi * (0 : 0.0025 : 1.999); 
 
end