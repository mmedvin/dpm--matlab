function tryRangeCompression2D
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    testSixDirections(); 
    testHalfLambda()
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
    title('1/4')
 
    subplot(132)
    anotherPos = pointSources(1).pos + (1/2) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = SARUtils.createPointSrc(pointSources(1).ampl, anotherPos);  
    plotRangeImg(pi, k_band, pointSources);     
    title('1/2')

    subplot(133)
    anotherPos = pointSources(1).pos + (1) * (2*pi/k0) * [1; 0]; 
    pointSources(2) = SARUtils.createPointSrc(pointSources(1).ampl, anotherPos);  
    plotRangeImg(pi, k_band, pointSources);     
    title('11')    
end



function testSixDirections()
    
    k_band = 50:0.1:55; 
    pointSourceSep   = 2.3; 
    pointSourceAngle = deg2rad(40); 
    
    pointSources(1) = SARUtils.createPointSrc(2.3, ... 
              pointSourceSep * [cos(pointSourceAngle); sin(pointSourceAngle)]); 
    pointSources(2) = SARUtils.createPointSrc(1.7, [0; 0]); 

    pointSources(3) = SARUtils.createPointSrc(1.1, [-2.5; 1.6]);     
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.6]); 
    
    subplot(231); 
    plotRangeImg(pi, k_band, pointSources); 
    title('horiz')   
    
    subplot(232); 
    plotRangeImg(pi + pointSourceAngle + pi/3, k_band, pointSources); 
    title('parallel + pi/3')  
    
    subplot(233); 
    plotRangeImg(pi + pointSourceAngle + pi/7, k_band, pointSources); 
    title('parallel + pi/7')    
    
    subplot(234); 
    plotRangeImg(pi + pointSourceAngle, k_band, pointSources); 
    title('parallel')      

    subplot(235); 
    plotRangeImg(pi + pointSourceAngle - pi/2, k_band, pointSources); 
    title('perp -')      

    subplot(236); 
    plotRangeImg(pi + pointSourceAngle + pi/2, k_band, pointSources); 
    title('perp +')  

end

function plotRangeImg(phi, k_band, pointSources) 

    %y2D = - [yRange * cos(phi); yRange * sin(phi)]; % when the source is "left and below", we want both coords positive

    
    [y1, y2] = meshgrid(-5:0.2:5, -3:0.2:3); 
    
    x_hat = [cos(phi); sin(phi)];     
    
    % create field for each k 
    uinf = zeros(size(k_band)); 
    for il = 1:numel(uinf)
        for ij = 1:numel(pointSources)
            ps = pointSources(ij); 
            % dot(,) below is OK because both args are real
            uinf(il) = uinf(il) + ps.ampl * exp(-2i * k_band(il) * dot(x_hat, ps.pos)); 
        end
    end
    
    % create image
    [rowsI, colsI] = size(y1); 
    I = zeros(rowsI, colsI); 

    for il = 1:numel(uinf)
        dotprod = x_hat(1) * y1 + x_hat(2) * y2; 
        I = I + uinf(il) * exp(2i * k_band(il) * dotprod); 
    end    

    pcolor(y1, y2, abs(I)); colorbar
    
    for ij = 1:numel(pointSources)
        ps = pointSources(ij); 
        
        hold on; plot(ps.pos(1), ps.pos(2), 'ks'); 
    end
end 

