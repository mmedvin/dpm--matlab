function mi__trySARfromVaryingSector
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    phi_range = pi * (1 + (-0.12 : 0.004 : 0.08)); 
    k_band = 50 : 1 : 55; 
  
    R = 5.4; 
    rCenter = [1.2; 1.7]; 
    
    % ??? sharp transition about endx = 0.5 (btw. 0.4 and 0.6) - put assert to catch it 
    theta = pi * (1 + (-0.08 : 0.0025 : 0.2)); 
    reflectivity_floor = 0.07; % (epsilon-1)/4
    reflectivity_peak_height = 0.2; % 0.2; 
    peak_theta_pos = pi * (1 + [-0.02, 0.04, 0.12]); 
    peak_theta_width = pi * 0.05; 

    [y2D.y1, y2D.y2] = meshgrid(-6 : 0.1 : -2, -2 : 0.1 : 4);     

    [rowsTheta, colsTheta] = size(theta); assert(rowsTheta == 1); 
    curve.z = repmat(rCenter, 1, colsTheta) + R * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 
    
    reflector.curve = curve; 
    
    tmp_refl = reflectivity_floor * ones(1, colsTheta); 
    for ipeak = 1:numel(peak_theta_pos) 
        pos_theta = peak_theta_pos(ipeak); 
        tmp_refl = tmp_refl + reflectivity_peak_height * ... 
                   triangularPulse(pos_theta - peak_theta_width/2, pos_theta + peak_theta_width/2, theta); 
    end
    reflector.coeff = tmp_refl;      
       
    % create field for each k and each direction 
    uinf = nan(numel(phi_range), numel(k_band)); 
    for im = 1:numel(phi_range)
        phi_xhat_val = phi_range(im); 
        for il = 1:numel(k_band)
            k = k_band(il);             
            isCheckShallow = true; % for varying incident angle, we don't want incidence get shallow
            scattField = mi__SARUtils.getScattField(k, phi_xhat_val, reflector, isCheckShallow); % use phi instead of phi_xhat here to always get specular reflection 
            traj = mi__SARUtils.doFFP_curve(k, phi_xhat_val, curve, scattField);  
            assert(numel(traj.ffp) == 1); 
            uinf(im, il) = traj.ffp; 
        end
    end
    
    I = mi__SARUtils.buildSARimage(y2D, k_band, phi_range, uinf); 

    figure('units', 'normalized', 'position', [0.02 0.1, 0.5, 0.4], 'paperpositionmode', 'auto'); 
    subplot(121)
    mi__SARUtils.displayReflector(reflector); 
    hc = colorbar; 
    ylabel(hc, 'reflectivity')    
    plotForPhi(rCenter, R, min(phi_range)); 
    plotForPhi(rCenter, R, max(phi_range)); 
    xlim([-6, -2]); ylim([-2, 4]); grid on
    xlabel('x'); ylabel('y')
    title('scatterer and specular limits')
    
    subplot(122) 
    pcolor(y2D.y1, y2D.y2, abs(I)); 
    hc = colorbar; shading flat;
    ylabel(hc, 'magnitude')
    daspect([1 1 1])
    xlabel('x'); ylabel('y')    
    title('SAR image')

    print('-djpeg', 'SAR_extendedScatt.jpg')
end 

function plotForPhi(rCenter, R, phi)
    assert(numel(phi) == 1); 
    coord = rCenter + R * [cos(phi); sin(phi)]; 
    hold on; plot(coord(1), coord(2), 'ko', 'MarkerFaceColor', 'k')
end

