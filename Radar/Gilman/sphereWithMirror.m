function sphereWithMirror
    set(0, 'defaultLineLineWidth', 2);    
    set(0, 'defaultLineMarkerSize', 30);     
    set(0, 'defaultAxesFontSize', 30);

    figure('units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.7]);     

    %%%%%%%%

    subplot(221)
    arg = -1.6:0.001: 1.5; 
    step_halfw = 1.3; 
    
    plot(arg, hump(arg, step_halfw), 'b');    
    hold on; plot(arg, smooth_step(arg, step_halfw), 'r'); 
    title('hump and step')
    
    hold on; plot(  [step_halfw, step_halfw], [0, 1], 'm--');     
    hold on; plot(- [step_halfw, step_halfw], [0, 1], 'm--');     
    grid on

    %%%%%%%%%
    subplot(222)
    
    arg = - 3.5 : 0.001 : 3.7; 
    plateau_halfw = 1.6; 
    step_halfw = 0.74; 
    
    [val, curveMax] = slide(arg, plateau_halfw, step_halfw); 

    plot(arg, val, 'k.-'); 
    %title(sprintf('slide, curveMax = %g', curveMax))
    title('slide')
    
    foot_lim = plateau_halfw + 2 * step_halfw; % TODO: refactor    
    levelHeight = plateau_halfw + curveMax;    
    % turns out that curveMax == step_halfw - let's check it!
    hold on; plot(  [-foot_lim, foot_lim], [ levelHeight - step_halfw,  levelHeight - step_halfw], 'k--');        
    hold on; plot(  [-foot_lim, foot_lim], [-levelHeight + step_halfw, -levelHeight + step_halfw], 'k--'); 
    
    hold on; plot(  [plateau_halfw, plateau_halfw], [-levelHeight, levelHeight], 'b--'); 
    hold on; plot(- [plateau_halfw, plateau_halfw], [-levelHeight, levelHeight], 'b--');
    
    % verticals
    hold on; plot(  [foot_lim, foot_lim], [-levelHeight,  levelHeight], 'm--');     
    hold on; plot(- [foot_lim, foot_lim], [-levelHeight,  levelHeight], 'm--');     

    % horizontals
    hold on; plot(  [-foot_lim, foot_lim], [ levelHeight,  levelHeight], 'm--');        
    hold on; plot(  [-foot_lim, foot_lim], [-levelHeight, -levelHeight], 'm--');        
    
    axis equal; % this is important - the slope should be 45 deg
    grid on; 
    
    yl = ylim; 
    hold on; plot(yl * 0.9, yl * 0.9, 'r-')     
    
    %%%%%%%%%
    subplot(223)
   
    theta = (0:0.02:2) * pi; 
    
    theta_plateau_center = pi; 
    theta_plateau_halfw  = pi/4; 
    theta_step_halfw     = pi/5; 
    
    assert(theta_plateau_halfw + theta_step_halfw < pi/2); 
          
    theta_truncated = slide(theta - theta_plateau_center, theta_plateau_halfw, theta_step_halfw);
    
    plot(theta, theta_truncated, 'b')
    title(sprintf('truncated angle\n halfwidths: plateau = %g * pi, step = %g * pi', ... 
                   theta_plateau_halfw / pi, theta_step_halfw / pi)); 
               
    theta_eff_max = theta_plateau_halfw + theta_step_halfw; % property of the "slide" 
    
    hold on; plot([0, 2*pi], [ theta_eff_max,  theta_eff_max], 'k--'); 
    hold on; plot([0, 2*pi], [-theta_eff_max, -theta_eff_max], 'k--'); 
    
    axis equal; % this is important - the slope should be 45 deg
    grid on
    
    %%%%%%%%%
    subplot(224)       
    
    distToMirror = cos(theta_plateau_halfw); 
    r_const = distToMirror / cos(theta_eff_max);     
    r_flat =  distToMirror ./ cos(theta_truncated); 
    
    shape_flat  = [r_flat .* cos(theta);  r_flat .* sin(theta)];         
    ht = plot(shape_flat(1,:), shape_flat(2,:),  'm--', 'DisplayName', 'target shape');  
    title(sprintf('circle with mirror\n dist = %5.3f, radius = %6.3f', distToMirror, r_const));  
    
    yl = [min(shape_flat(2,:)), max(shape_flat(2,:))]; 
    hold on; 
    hs = plot(- [distToMirror, distToMirror], yl, 'k--', 'DisplayName', 'straight line'); % minus because the angles are centered around pi 
    
    hold on; 
    hc = plot(r_const * cos(theta), r_const * sin(theta),  'b--', 'DisplayName', 'circle');  
    legend([ht, hc, hs])
    
    grid on
    axis equal    
    
    print('-djpeg', 'smooth_shape.jpg')
end

function [ret, curveMax]  = slide(x, plateau_halfw, step_halfw)

    ngrid = 131; 
    arg = linspace(- 1.1 * step_halfw, 1.1 * step_halfw, ngrid); 
    smooth_step_interp = griddedInterpolant(arg, smooth_step(arg, step_halfw), 'cubic', 'none'); 
       
    smooth_curve = @(upperLim) arrayfun(@(uL) integral(@(x) smooth_step_interp(x), -step_halfw, uL), upperLim);

    curveMax = smooth_curve(step_halfw); 
    
    
    foot_lim = plateau_halfw + 2 * step_halfw; % TODO: refactor        
    
    mask_plateau  = (abs(x) <= plateau_halfw); 

    mask_neg_foot  = (x <= - foot_lim);  
    mask_neg_curve = ((x < 0) & (~mask_plateau) & (~mask_neg_foot)); 
    
    mask_pos_foot  = (x >=   foot_lim);  
    mask_pos_curve = ((x > 0) & (~mask_plateau) & (~mask_pos_foot)); 
    
    assert(sum(mask_plateau) + sum(mask_neg_foot) + sum(mask_neg_curve) + sum(mask_pos_foot) + sum(mask_pos_curve) == numel(x)); 
    assert(all(mask_plateau | mask_neg_foot | mask_neg_curve | mask_pos_foot | mask_pos_curve)); 
            
    ret = nan(size(x)); 
    ret(mask_neg_foot)  = - plateau_halfw - curveMax; 
    ret(mask_neg_curve) = - plateau_halfw - curveMax + smooth_curve(  x(mask_neg_curve) + plateau_halfw + step_halfw);    
    
    ret(mask_pos_foot)  =   plateau_halfw + curveMax;         
    ret(mask_pos_curve) =   plateau_halfw + curveMax - smooth_curve(- x(mask_pos_curve) + plateau_halfw + step_halfw);    
    
    ret(mask_plateau) = x(mask_plateau); 
   
    assert(all(~isnan(ret))); 
end

% https://math.stackexchange.com/questions/328868/how-to-build-a-smooth-transition-function-explicitly
% answer by user1337
function ret = smooth_step(x, step_halfw)

    intOfHump = @(upperLim) integral(@(x) hump(x, step_halfw), -step_halfw, upperLim);

    A = intOfHump(step_halfw); 
    
    mask_need_integral = (abs(x) < step_halfw); 
    
    ret(~mask_need_integral & (x < 0)) = 0; 
    ret(~mask_need_integral & (x > 0)) = 1;     
    ret( mask_need_integral) = (1/A) * arrayfun(@(x_) intOfHump(x_), x(mask_need_integral)); 
    
end

function ret = hump(x, halfw)
    
    assert(halfw > 0); 
    mask_gt_w = (abs(x) >= halfw); 
    x(mask_gt_w) = 0; 

    in_exp = 1 ./ (x.^2 - halfw^2); 
    assert(all(in_exp < 0)); 
    ret = exp(in_exp);     
    ret(mask_gt_w) = 0; 
end
