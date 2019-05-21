function trySpecularFromVaryingSector

    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 

    k = 50; 
    R = 5.4; 
    rCenter = [1.2; 1.7]; 
    
    % ??? sharp transition about endx = 0.5 (btw. 0.4 and 0.6) - put assert to catch it 
    endx = 0.2;
    theta = pi * (1 + (-0.1 : 0.0005 : endx)); 
    reflectivity_floor = 0.07; % (epsilon-1)/4
    reflectivity_peak_height = 0.2; 
    peak_theta_pos = pi * (1 + [0.05, 0.17]); 
    peak_theta_width = pi * 0.05; 

    % cannot view from all angles because shallow incidence is prohibited in the Born-Kirchhoff model 
    %pi * (0 : 0.02 : 2); % where we calculate FFP
    min_phi = min(theta) - 0.3 * (max(theta) - min(theta)); 
    max_phi = max(theta) + 0.3 * (max(theta) - min(theta)); 
    phi = min_phi : 0.01 : max_phi;  

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
    


    specular_ffp = nan(size(phi)); 
    for iphi = 1:numel(phi)
        phi_val = phi(iphi); 
        scattField = SARUtils.getScattField(k, phi_val, reflector); % use phi instead of phi_xhat here to always get specular reflection 
        traj = SARUtils.doFFP_curve(k, phi_val, curve, scattField);  
        assert(numel(traj.ffp) == 1); 
        specular_ffp(iphi) = traj.ffp; 
    end
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.4]); 
    subplot(121)
    SARUtils.displayReflector(reflector); 
    title('reflector')
    
    subplot(122)
    plot(phi, abs(specular_ffp), 'bx')
    hold on; plot(peak_theta_pos, zeros(size(peak_theta_pos)), 'rs'); 
    hold on; plot([min(theta), max(theta)], [0, 0], 'rs', 'MarkerFaceColor', 'r');     
    legend('abs(ffp)', sprintf('reflectivity \nmaxima'), sprintf('reflector \nlimits'))
    title('specular reflection')
end


