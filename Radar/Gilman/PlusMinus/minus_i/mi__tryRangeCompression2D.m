function mi__tryRangeCompression2D
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    testSixDirections(); 
    testHalfLambda()
end


function testSixDirections()
    
    k_band = 50:0.1:55; 
    pointScattererSep   = 2.3; 
    pointScattererAngle = deg2rad(40); 
    
    [y2D.y1, y2D.y2] = meshgrid(-5:0.2:5, -3:0.2:3); 
    
    pointScatterers(1) = mi__SARUtils.createPointWithAmplitude(2.3, ... 
              pointScattererSep * [cos(pointScattererAngle); sin(pointScattererAngle)]); 
    pointScatterers(2) = mi__SARUtils.createPointWithAmplitude(1.7, [0; 0]); 
    pointScatterers(3) = mi__SARUtils.createPointWithAmplitude(1.1, [-2.5; 1.6]);  
    
    % the first two form a pair of scatterers, the third is an extra scatterer
    pointColors = {'m', 'm', 'w'};  
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.6]); 
    
    subplot(231); 
    plotRangeImg(y2D, pi, k_band, pointScatterers, pointColors); 
    title('horiz')   
    
    subplot(232); 
    plotRangeImg(y2D, pi + pointScattererAngle + pi/3, k_band, pointScatterers, pointColors); 
    title('parallel + pi/3')  
    
    subplot(233); 
    plotRangeImg(y2D, pi + pointScattererAngle + pi/7, k_band, pointScatterers, pointColors); 
    title('parallel + pi/7')    
    
    subplot(234); 
    plotRangeImg(y2D, pi + pointScattererAngle, k_band, pointScatterers, pointColors); 
    title('parallel')      

    subplot(235); 
    plotRangeImg(y2D, pi + pointScattererAngle - pi/2, k_band, pointScatterers, pointColors); 
    title('perp -')      

    subplot(236); 
    plotRangeImg(y2D, pi + pointScattererAngle + pi/2, k_band, pointScatterers, pointColors); 
    title('perp +')  

    suptitle('The direction is relative to the pair of ''m'' scatterers')
end


function testHalfLambda()

    k_band = 50:0.1:52; 
    [y2D.y1, y2D.y2] = meshgrid(-5:0.2:5, -3:0.2:3);     
    
    k0 = mean(k_band); 
    
    pointSources(1) = mi__SARUtils.createPointWithAmplitude(1.1, [-2.5; 1.6]); 
    anotherPos = pointSources(1).pos + (1/4) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = mi__SARUtils.createPointWithAmplitude(pointSources(1).ampl, anotherPos);  
    pointColors = {'m', 'm'}; % there will always be just two scatterers
                                           
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.2]); 

    subplot(131)
    anotherPos = pointSources(1).pos + (1/4) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = mi__SARUtils.createPointWithAmplitude(pointSources(1).ampl, anotherPos);  
    plotRangeImg(y2D, pi, k_band, pointSources, pointColors); 
    title(sprintf('separation: 1/4 times lambda \n (!!!NOTE: yields low return!)'))
 
    subplot(132)
    anotherPos = pointSources(1).pos + (1/2) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = mi__SARUtils.createPointWithAmplitude(pointSources(1).ampl, anotherPos);  
    plotRangeImg(y2D, pi, k_band, pointSources, pointColors);     
    title(sprintf('separation: 1/2 times lambda \n (two returns in phase)'))

    subplot(133)
    anotherPos = pointSources(1).pos + (1) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = mi__SARUtils.createPointWithAmplitude(pointSources(1).ampl, anotherPos);  
    plotRangeImg(y2D, pi, k_band, pointSources, pointColors);     
    title(sprintf('separation: 1  times lambda \n (two returns in phase)'))    
end

function plotRangeImg(y2D, phi_refl, k_band, pointScatterers, pointColors) 
    % what is xhat? is it the argument of FFP, i.e., the direction of _scattered_ field
    x_hat = [cos(phi_refl); sin(phi_refl)];     
    uinf = mi__SARUtils.create_uinf_fromPoints(k_band, pointScatterers, phi_refl);

    I = mi__SARUtils.buildRangeCompressionImg2D(y2D, k_band, x_hat, uinf); 

    pcolor(y2D.y1, y2D.y2, abs(I)); colorbar
    
    for ij = 1:numel(pointScatterers)
        ps = pointScatterers(ij); 
        
        hold on; hs = plot(ps.pos(1), ps.pos(2), 's', 'MarkerSize', 12); 
        set(hs, 'MarkerEdgeColor', pointColors{ij}, 'MarkerFaceColor', pointColors{ij}); 
    end
end 


