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

    figTitle = sprintf('First point at (%4.1f,%4.1f)', pointSources(1).pos(1), pointSources(1).pos(2)); 
    figure('units', 'normalized', 'position', [0.02 0.1, 0.9, 0.6], 'name', figTitle); 
    
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

function plotRangeImg(yRange, phi_refl, k_band, pointSources) 

    % phi is the direction towards the source

    % what is xhat? is it the argument of FFP, i.e., the direction of _reflected_ field, 
    % which is the same as direction TO the source
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

    I = SARUtils.buildRangeCompressionImg1D(yRange, k_band, x_hat, uinf); 

    plot(yRange, abs(I), 'b.'); 
    xlabel('range distance')
    
    % we measure distance _from the source_, i.e., in the direction of the _incident wave_
    phi_inc = phi_refl + pi; 
    for ij = 1:numel(pointSources)
        ps = pointSources(ij); 
        
        % we can also do: dot(ps.pos, -k_hat) --- note minus here! 
        dist = ps.pos(1) * cos(phi_inc) ... 
             + ps.pos(2) * sin(phi_inc); 
        
        hold on; plot(dist, 0, 'rx'); 
    end
end 
