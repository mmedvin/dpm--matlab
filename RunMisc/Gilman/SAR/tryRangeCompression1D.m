function tryRangeCompression1D
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    yRange = -4:0.02:4; 
    k_band = 50:0.1:55; 
    pointSourceSep   = 2.3; 
    pointSourceAngle = deg2rad(40); 
    
    pointSources(1) = SARUtils.createPointSrc(2.3, ... 
              pointSourceSep * [cos(pointSourceAngle); sin(pointSourceAngle)]); 
    pointSources(2) = SARUtils.createPointSrc(1.7, [0; 0]); 

    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.6]); 
    
    subplot(231); 
    plotRangeImg(yRange, pi, k_band, pointSources); 
    title('horiz')   
    
    subplot(232); 
    plotRangeImg(yRange, pi + pointSourceAngle + pi/3, k_band, pointSources); 
    title('parallel + pi/3')  
    
    subplot(233); 
    plotRangeImg(yRange, pi + pointSourceAngle + pi/7, k_band, pointSources); 
    title('parallel + pi/7')    
    
    subplot(234); 
    plotRangeImg(yRange, pi + pointSourceAngle, k_band, pointSources); 
    title('parallel')      

    subplot(235); 
    plotRangeImg(yRange, pi + pointSourceAngle - pi/2, k_band, pointSources); 
    title('perp -')      

    subplot(236); 
    plotRangeImg(yRange, pi + pointSourceAngle + pi/2, k_band, pointSources); 
    title('perp +')   
end

function plotRangeImg(yRange, phi, k_band, pointSources) 

    x_hat = [cos(phi); sin(phi)]; 
    
    y2D = - [yRange * cos(phi); yRange * sin(phi)]; % when the source is "left and below", we want both coords positive

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
    I = zeros(size(yRange)); 

    for iy = 1:numel(I)
        % dot(,) below is OK because both args are real        
        dotprod = dot(x_hat, y2D(:,iy)); 
        for il = 1:numel(uinf)
            I(iy) = I(iy) + uinf(il) * exp(2i * k_band(il) * dotprod); 
        end
    end
    

    plot(yRange, abs(I), 'b.'); 
    
    for ij = 1:numel(pointSources)
        ps = pointSources(ij); 
        
        %dist = ps.pos(1); 
        
        dist = ps.pos(1) * cos(phi+pi) ... 
             + ps.pos(2) * sin(phi+pi); 
        
        hold on; plot(dist, 0, 'rx'); 
    end
end 

