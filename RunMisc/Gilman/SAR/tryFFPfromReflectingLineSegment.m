function tryFFPfromReflectingLineSegment
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 


    k = 50; 
    phi_xhat = 1.45 * pi; % incident wave
    
    %zStart = [5; 0]; zEnd   = [0; 2]; 
    %zStart = [0; 0]; zEnd   = [3; 0]; 
    zStart = [0; 0]; zEnd   = [0; 3]; % grazing incidence
    zN = 700; 
    
    A = 0.13; % (epsilon-1)/4
    
    phi = pi * (0 : 0.001 : 2); % where we calculate FFP

    z1 = linspace(zStart(1), zEnd(1), zN); 
    z2 = linspace(zStart(2), zEnd(2), zN);     
    curve.z = [z1; z2]; 
    
    n_z = SARUtils.getUnitVectorPerpTo(zEnd - zStart);  
    curve.nz =  repmat(n_z, 1, zN); 
    
    xhat = [cos(phi_xhat); sin(phi_xhat)];
    k_inc = - k * xhat; 
    % sum(..., 1) works as columnwise dot(,), i.e., dot(,) at each z 
    % BTW dot here OK because args are real
    % same as: q_inc_z = k * sum(repmat(xhat, 1, colsTheta) .* nz, 1);    
    q_inc_z = - sum(repmat(k_inc, 1, zN) .* curve.nz, 1); 
                     

    
    incWavePhase = sum(repmat(k_inc, 1, zN) .* curve.z, 1); 
    scattField.value = - A * k^2 * q_inc_z.^(-2) .* exp(1i * incWavePhase); 
    scattField.normal_deriv = 1i * q_inc_z .* scattField.value; 
    
    %{
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.6]); 
    subplot(221); plot(1:zN, abs(scattField.value), 'bx'); title('abs value')
    subplot(222); plot(1:zN, angle(scattField.value), 'bx'); title('angle value')    
    subplot(223); plot(1:zN, abs(scattField.normal_deriv), 'bx'); title('abs normal deriv')
    subplot(224); plot(1:zN, angle(scattField.normal_deriv), 'bx'); title('angle normal deriv')    
    %}
   
    
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
    title(sprintf('abs peak, phi hat = %6.3f * pi', phi_xhat/pi));       

    n_z_angle = atan2(n_z(2), n_z(1)); 
    geometric_refl_angle = mod(2 * n_z_angle - phi_xhat, 2*pi); % want the result in the interval (0, 2*pi)
    hold on; plot(geometric_refl_angle, 0, 'rs')
    legend('diagram', sprintf('specular \n refl angle'))    
    
end


