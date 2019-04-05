function sigmoid_based

    R = 2.45;
    
    limits.theta_m = pi/4; 
    limits.theta_c = 2 * pi / 5;    
    
    Pthree = [-20, 70, -84, 35, 0, 0, 0, 0]; % such that sigmoid(x) is polyval(Pthree,x); 
    
    theta = (-1 : 0.005 : 1) * pi; 
    
    assert(limits.theta_m < limits.theta_c); 
    assert(limits.theta_c < pi/2);  
    limits.theta_s = 2 * limits.theta_c - limits.theta_m; 
    
    d = R * cos(limits.theta_c);  
    
    vt = vartheta(theta, Pthree, limits); 
    r = d ./ cos(vt); 
    x = - r .* cos(theta); % minus because I want the flat region face waves with k_x > 0
    y =   r .* sin(theta); 
    
    figure('units', 'normalized', 'position', [0.1 0.1 0.8 0.5]);  
    
    subplot(121)
    plot(theta, vt, 'b.'); grid on
    title('vartheta')
    
    subplot(122); 
             hs = plot(x, y, 'b-', 'DisplayName', 'target shape', 'LineWidth', 4); grid on  
    hold on; hcirc = plot(R * cos(theta), R * sin(theta), 'm-', 'DisplayName', 'prescribed circle');  
    hold on; hchrd = plot([-d, -d], [R * sin(limits.theta_c), - R * sin(limits.theta_c)], 'r--', 'DisplayName', 'full chord');     
    
    % mark all angles
    % use "minus" in the first argument because of the minus in "x = - r .* cos(theta);"
    hold on; plot([0, -d], [0,  d/cos(limits.theta_m)], 'm--')    
    hold on; plot([0, -d], [0, -d/cos(limits.theta_m)], 'm--')
    hold on; plot([0, -R * cos(limits.theta_s)], [0,  R * sin(limits.theta_s)], 'k--'); 
    hold on; plot([0, -R * cos(limits.theta_s)], [0, -R * sin(limits.theta_s)], 'k--');
    hold on; plot([0, -R * cos(limits.theta_c)], [0,  R * sin(limits.theta_c)], 'c--'); 
    hold on; plot([0, -R * cos(limits.theta_c)], [0, -R * sin(limits.theta_c)], 'c--');
    
    hl = legend([hs, hcirc, hchrd]);
    set(hl, 'location', 'EastOutside')
    title('target shape') 
    
    axis equal
    
    print -djpeg sigmoid_chord.jpg
end

function ret = vartheta(theta, P, limits)

    mask_slide  = (abs(theta) <= limits.theta_m); 

    mask_neg_foot  = (theta <= - limits.theta_s);  
    mask_neg_curve = ((theta < 0) & (~mask_slide) & (~mask_neg_foot)); 
    
    mask_pos_foot  = (theta >=   limits.theta_s);  
    mask_pos_curve = ((theta > 0) & (~mask_slide) & (~mask_pos_foot)); 
    
    % sanity check: make sure all masks are not intersecting and don't leave anything out 
    assert(sum(mask_slide) + sum(mask_neg_foot) + sum(mask_neg_curve) + sum(mask_pos_foot) + sum(mask_pos_curve) == numel(theta)); 
    assert(all(mask_slide | mask_neg_foot | mask_neg_curve | mask_pos_foot | mask_pos_curve)); 
    
    D = 2 * (limits.theta_c - limits.theta_m); 
    Q = polyint(P);

    T_D = @(x) D * polyval(Q, x/D); 
    
    ret = nan(size(theta)); 
    ret(mask_slide) = theta(mask_slide); 
    ret(mask_neg_foot) = - limits.theta_c; 
    ret(mask_pos_foot) =   limits.theta_c; 
    
    ret(mask_neg_curve) = - limits.theta_c + T_D(  theta(mask_neg_curve) + limits.theta_s); 
    ret(mask_pos_curve) =   limits.theta_c - T_D(- theta(mask_pos_curve) + limits.theta_s); 
    
    % a couple more (paranoic) sanity checks
    assert(~any(isnan(ret)));         % did not leave anything out
    assert(all( abs(ret) < pi/2 ));   % cos(ret) > 0 
end