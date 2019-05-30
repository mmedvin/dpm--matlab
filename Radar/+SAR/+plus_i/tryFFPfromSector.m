function tryFFPfromSector
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 

    k = 50; 
    % TODO: what is this? -phi_inc?
    phi_xhat = 1.15 * pi; 
    R = 5.4; 
    rCenter = [1.2; 1.7]; 
    
    % ??? sharp transition about endx = 0.5 (btw. 0.4 and 0.6) - put assert to catch it 
    endx = 0.4;
    theta = pi * (1 + (-0.2 : 0.0005 : endx)); 
    reflectivity_const = 0.07; % (epsilon-1)/4
    
    phi = pi * (0 : 0.02 : 2); % where we calculate FFP
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.4]); 
    subplot(121)
    traj1layer = tryOneLayer(k, phi_xhat, R, rCenter, theta, reflectivity_const, phi); 
    plot(phi, abs(traj1layer.ffp), 'bx'); title('one layer') 
    hold on; plot( 2*[min(theta), max(theta)] - phi_xhat, [0,0], 'rs'); % geometric boundaries of specular reflection 
    legend('ffp', 'geom. range')
    
    subplot(122)
    traj2layers = tryTwoLayers(k, phi_xhat, R, rCenter, theta, reflectivity_const, phi); 
    plot(phi, abs(traj2layers.ffp), 'bx'); title('two layers')     
end

function traj = tryOneLayer(k, phi_xhat, R, rCenter, theta, reflectivity_const, phi)
    [rowsTheta, colsTheta] = size(theta); assert(rowsTheta == 1); 
    curve.z = repmat(rCenter, 1, colsTheta) + R * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 
    
    reflector.curve = curve; 
    reflector.coeff = reflectivity_const * ones(1, colsTheta); 
    
    TODO: FIX! 
    scattField = SARUtils.getScattField(k, phi_xhat, reflector); 
    
    %{
    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.4]); 
    subplot(131); plot(theta, abs(scattField.value), 'bx'); title('value')
    subplot(132); plot(theta, abs(scattField.normal_deriv), 'bx'); title('normal deriv')
    subplot(133); plot(theta, q_inc_z, 'bx'); title('q inc')
    %}
    
    traj = SARUtils.doFFP_curve(k, phi, curve, scattField); 
end

function traj = tryTwoLayers(k, phi_xhat, R, rCenter, theta, reflectivity_const, phi)
    
    [rowsTheta, colsTheta] = size(theta); assert(rowsTheta == 1); 
    curve.z = repmat(rCenter, 1, colsTheta) + R * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 
   
    %adding second layer
    curve.z (:,(colsTheta + 1) : (2*colsTheta)) =  curve.z (:,1:colsTheta); 
    curve.nz(:,(colsTheta + 1) : (2*colsTheta)) = -curve.nz(:,1:colsTheta);   
    
    reflector.curve = curve; 
    reflector.coeff = reflectivity_const * ones(1, 2*colsTheta); % note "2" here!     
 
    scattField = SARUtils.getScattField(k, phi_xhat, reflector);     
    
    %{
    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.4]); 
    subplot(131); plot([theta, theta], abs(scattField.value), 'bx'); title('value')
    subplot(132); plot([theta, theta], abs(scattField.normal_deriv), 'bx'); title('normal deriv')
    subplot(133); plot([theta, theta], q_inc_z, 'bx'); title('q inc')
    %}
    
    traj = SARUtils.doFFP_curve(k, phi, curve, scattField); 
    

end