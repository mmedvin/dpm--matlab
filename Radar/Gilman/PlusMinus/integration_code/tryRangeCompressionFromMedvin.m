function tryRangeCompressionFromMedvin

    %inMatFilePrefix = 'circle_gp10_'; spec = getSpec();               % wide PSF
    %inMatFilePrefix = 'circle_gp11_'; spec = getSpec_2K();            % many stripes
    %inMatFilePrefix = 'Circler09Shifted03_gp11_'; spec = nan; 
    inMatFilePrefix = 'Circler09Shifted03_gp11_kf_half_'; spec = nan; 
    
    phiIncToPlot = [-0.2, 0, 0.2] * pi; 
    
    %phiToPlot = 0 * pi; 
    
    plotRangeProfiles(inMatFilePrefix, spec, phiIncToPlot);
end

function plotRangeProfiles(inMatFilePrefix, spec, phiIncToPlot)

    [y1, y2] = meshgrid(-2:0.05:2.1, -2.3:0.05:2.1); % for 2D 
    yRange = -4:0.02:4;                              % for 1D 
    
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
    
    % create field for each k and each _specified_ direction (unlike trySARfromMedvin.m)

    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.4], 'name', inMatFile); 
    numPhiPlots = numel(phiIncToPlot); 
    cmap_pcolor = colormap; 
    cmap_lines = colormap(lines); 
    for iphiToPlot = 1:numPhiPlots
        
        phiInc_prescribed = phiIncToPlot(iphiToPlot); 
         
        [~, im] = min(abs(phi_inc - phiInc_prescribed));
        if numel(im) > 1 
            im = im(1); 
        end    
        phi_refl = phi_inc(im) + pi; 
        
        uinf = nan(1,numel(k_band));         
        for il = 1:numel(k_band)
            k = k_band(il); 
            scattField = scattFields(im,il); 
            traj = SARUtils.doFFP_curve(k, phi_refl, curve, scattField);  
            assert(numel(traj.ffp) == 1); 
            uinf(il) = traj.ffp; 
        end
 
        % what is xhat? is it the argument of FFP, i.e., the direction of _scattered_ field
        x_hat = [cos(phi_refl); sin(phi_refl)];           
        
        I2D = SARUtils.buildRangeCompressionImg2D(y1, y2, k_band, x_hat, uinf);         
 
        subplot(1, numPhiPlots+1, iphiToPlot); 
        colormap(cmap_pcolor); 
        pcolor(y1, y2, log10(abs(I2D))); 

        hc = colorbar; shading flat;
        ylabel(hc, 'log10(magnitude)')
        title(sprintf('im = %d, phi (refl) = %g', im, phi_refl)); 
        
        %%%%%%%%%%%%
        
        subplot(1, numPhiPlots+1, numPhiPlots+1); 
        I1D = SARUtils.buildRangeCompressionImg1D(yRange, k_band, x_hat, uinf); 
        hold on; h = plot(yRange, abs(I1D), '-.');   
        set(h, 'Color', cmap_lines(iphiToPlot, :));
        set(h, 'DisplayName', sprintf('phi (refl) = %3.1f * pi', phi_refl / pi));         
        
    end    
    
    subplot(1, numPhiPlots+1, numPhiPlots+1); 
    legend('location', 'NorthEast')
        
end


function spec = getSpec_2K()
 
    spec.k_range = 200:5:220; 
    spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
    spec.theta = pi * (0 : 0.0025 : 1.999); 

end

