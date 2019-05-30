function tryRangeCompression2D
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    testSixDirections(); 
    %testHalfLambda()
end


function testSixDirections()
    
    k_band = 50:0.1:55; 
    pointScattererSep   = 2.3; 
    pointScattererAngle = deg2rad(40); 
    
    pointScatterers(1) = SARUtils.createPointSrc(2.3, ... 
              pointScattererSep * [cos(pointScattererAngle); sin(pointScattererAngle)]); 
    pointScatterers(2) = SARUtils.createPointSrc(1.7, [0; 0]); 

    pointScatterers(3) = SARUtils.createPointSrc(1.1, [-2.5; 1.6]);     
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.6]); 
    
    subplot(231); 
    plotRangeImg(pi, k_band, pointScatterers); 
    title('horiz')   
    
    subplot(232); 
    plotRangeImg(pi + pointScattererAngle + pi/3, k_band, pointScatterers); 
    title('parallel + pi/3')  
    
    subplot(233); 
    plotRangeImg(pi + pointScattererAngle + pi/7, k_band, pointScatterers); 
    title('parallel + pi/7')    
    
    subplot(234); 
    plotRangeImg(pi + pointScattererAngle, k_band, pointScatterers); 
    title('parallel')      

    subplot(235); 
    plotRangeImg(pi + pointScattererAngle - pi/2, k_band, pointScatterers); 
    title('perp -')      

    subplot(236); 
    plotRangeImg(pi + pointScattererAngle + pi/2, k_band, pointScatterers); 
    title('perp +')  

end


function testHalfLambda()

    k_band = 50:0.1:52; 
    
    k0 = mean(k_band); 
    
    pointSources(1) = SARUtils.createPointSrc(1.1, [-2.5; 1.6]); 
    anotherPos = pointSources(1).pos + (1/4) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = SARUtils.createPointSrc(pointSources(1).ampl, anotherPos);  
                                           
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.2]); 

    subplot(131)
    anotherPos = pointSources(1).pos + (1/4) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = SARUtils.createPointSrc(pointSources(1).ampl, anotherPos);  
    plotRangeImg(pi, k_band, pointSources); 
    title(sprintf('separation: 1/4 times lambda \n (yields low return!)'))
 
    subplot(132)
    anotherPos = pointSources(1).pos + (1/2) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = SARUtils.createPointSrc(pointSources(1).ampl, anotherPos);  
    plotRangeImg(pi, k_band, pointSources);     
    title(sprintf('separation: 1/2 times lambda \n (two returns in phase)'))

    subplot(133)
    anotherPos = pointSources(1).pos + (1) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = SARUtils.createPointSrc(pointSources(1).ampl, anotherPos);  
    plotRangeImg(pi, k_band, pointSources);     
    title(sprintf('separation: 1  times lambda \n (two returns in phase)'))    
end

function plotRangeImg(phi_refl, k_band, pointSources) 

    [y1, y2] = meshgrid(-5:0.2:5, -3:0.2:3); 
    
    % what is xhat? is it the argument of FFP, i.e., the direction of _scattered_ field
    x_hat = [cos(phi_refl); sin(phi_refl)];     
    
    % create field for each k 
    uinf = zeros(size(k_band)); 
    for il = 1:numel(uinf)
        for ij = 1:numel(pointSources)
            ps = pointSources(ij); 
            % dot(,) below is OK because both args are real
            uinf(il) = uinf(il) + ps.ampl * exp(-2i * k_band(il) * dot(x_hat, ps.pos)); 
        end
    end
    
    I = SARUtils.buildRangeCompressionImg2D(y1, y2, k_band, x_hat, uinf); 

    pcolor(y1, y2, abs(I)); colorbar
    
    for ij = 1:numel(pointSources)
        ps = pointSources(ij); 
        
        hold on; plot(ps.pos(1), ps.pos(2), 'ms', 'MarkerSize', 20); 
    end
end 

