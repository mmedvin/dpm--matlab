function tryFFPfromShiningLineSegment
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 

    k = 50; 
    
    zStart = [5; 0]; zEnd   = [0; 2]; 
    %zStart = [0; 0]; zEnd   = [3; 0]; 
    %zStart = [0; 0]; zEnd   = [0; 3]; 
    zN = 700; 
  
    phi = pi * (0 : 0.001 : 2); % where we calculate FFP

    z1 = linspace(zStart(1), zEnd(1), zN); 
    z2 = linspace(zStart(2), zEnd(2), zN);     
    curve.z = [z1; z2]; 
       
    n_z = SARUtils.getUnitVectorPerpTo(zEnd - zStart); 
    curve.nz =  repmat(n_z, 1, zN); 
    
    scattField.value = ones(1, zN);  
    scattField.normal_deriv = 1i * k .* scattField.value; 
    
    traj = SARUtils.doFFP_curve(k, phi, curve, scattField); 
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.4]); 
    subplot(131)
    plot(phi, abs(traj.ffp), 'b.-'); title('abs') 
    
    subplot(132)
    plot(phi, angle(traj.ffp), 'bx'); title('angle')     
    
    [~, ind] = max(traj.ffp); 
    
    rngToPlot = max(1, ind - 50) : min (ind + 50, numel(phi)); 
 
    subplot(133)
    plot(phi(rngToPlot), abs(traj.ffp(rngToPlot)), '.-'); 
    title(sprintf('abs peak, nz = (%5.3f, %5.3f)', n_z(1), n_z(2)));  
    
    n_z_angle = mod(atan2(n_z(2), n_z(1)), 2*pi); % want the result in the interval (0, 2*pi)
    hold on; plot(n_z_angle, 0, 'rs')
    legend('diagram', 'normal to line')
    
    suptitle('Shining line segment')
end


