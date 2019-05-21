function tryFFPfromMedvin
    %inMatFilePrefix = 'circle_gp10_'; spec = getSpec();               % wide PSF
    %inMatFilePrefix = 'circle_gp11_'; spec = getSpec_2K();            % many stripes
    %inMatFilePrefix = 'Circler09Shifted03_gp11_'; spec = nan; 
    inMatFilePrefix = 'Circler09Shifted03_gp11_kf_half_'; spec = nan; 
    
    phiIncToPlot = [-0.2, 0, 0.2] * pi; 
       
    ilToPlot = [1, 3, 5]; 
    
    plotFFPs(inMatFilePrefix, spec, phiIncToPlot, ilToPlot);
end

function plotFFPs(inMatFilePrefix, spec, phiIncToPlot, ilToPlot)

    phi_refl = (0: 0.02: 2) * pi; 

    if ~isstruct(spec)
        fprintf('Reading spec from a large file\n'); 
        tic
        load([inMatFilePrefix 'ind_a.mat'],'spec');
        toc
        clear irs
    end

    scatterer_pattern = 'iface.mat'; 
    inMatFile = [inMatFilePrefix, scatterer_pattern];     
    
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

    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.4], 'name', inMatFile);     
    
    
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
        for ili = 1:numel(ilToPlot)
            il = ilToPlot(ili);
            
            k = k_band(il); 
            scattField = scattFields(im,il); 
            traj = SARUtils.doFFP_curve(k, phi_refl, curve, scattField);  

            hold on; h = plot(phi_refl, abs(traj.ffp), '-.'); 
            set(h, 'Color', cmap_lines(ili, :));
            set(h, 'DisplayName', sprintf('phi inc = %3.1f * pi', phi_inc_val / pi));  
            
        end       
        

        hold on; hi = plot(mod(phi_inc_val     , 2*pi), 0, 'bs', 'DisplayName', 'incident');         
        hold on; hs = plot(mod(phi_inc_val + pi, 2*pi), 0, 'ks', 'DisplayName', 'specular');         

        legend([hi, hs]); 
        
        title(sprintf('phi inc = %3.1f * pi', phi_inc_val / pi))
    end

end

 function spec = getSpec_2K()

     spec.k_range = 200:5:220; 
     spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
     spec.theta = pi * (0 : 0.0025 : 1.999); 
 
 end