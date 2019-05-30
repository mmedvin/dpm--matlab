function mi__tryRangeCompression1D
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    y_r = -4:0.02:4; 
    k_band = 50:0.1:55; 
    pointSourceSep   = 2.3; 
    pointSourceAngle = deg2rad(40); 
    
    pointSources(1) = mi__SARUtils.createPointWithAmplitude(2.3, ... 
              pointSourceSep * [cos(pointSourceAngle); sin(pointSourceAngle)]); 
    pointSources(2) = mi__SARUtils.createPointWithAmplitude(1.7, [0; 0]); 

    figTitle = sprintf('First point at (%4.1f,%4.1f)', pointSources(1).pos(1), pointSources(1).pos(2)); 
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.6], 'name', figTitle); 
    
    xhat_angles_with_descrs = ... 
        {{pi, 'horizontal'}, ... 
         {pi + pointSourceAngle + pi/3,   'parallel + pi/3'}, ... 
         {pi + pointSourceAngle - 4*pi/5, 'parallel - 4*pi/5'}, ... 
         {pi + pointSourceAngle,          'parallel'}, ... 
         {pi + pointSourceAngle + pi/2,   'perp-1: parallel + pi/2'}, ... 
         {pi + pointSourceAngle - pi/2,   'perp-2: parallel - pi/2'}}; 
     
    for istr = 1:numel(xhat_angles_with_descrs)
        subplot(2, 3, istr); 
        xhat_angle_with_descr = xhat_angles_with_descrs{istr}; 
        xhat_angle = xhat_angle_with_descr{1}; 
        plotRangeImg(y_r, xhat_angle, k_band, pointSources); 
        title(sprintf('%d: %6.2f deg. -- %s', istr, rad2deg(xhat_angle), xhat_angle_with_descr{2})); 
    end

    suptitle('NOTE: y_r = (xhat,y) INCREASES towards the source') 
    
    %%%%%%%%%%%%
    figure('units', 'normalized', 'position', [0.02 0.1, 0.4, 0.6], 'name', 'scatterers and angles'); 
    
    for ipnt = 1:numel(pointSources)
        pos = pointSources(ipnt).pos; 
        hold on; plot(pos(1), pos(2), 'r*'); 
        text(pos(1), pos(2), sprintf(' point %d: ampl: %g', ipnt, pointSources(ipnt).ampl), 'FontSize', 24);
        grid on; axis equal
    end
     
    
    for istr = 1:numel(xhat_angles_with_descrs)
        xhat_angle_with_descr = xhat_angles_with_descrs{istr};    
        xhat_angle = xhat_angle_with_descr{1};
        hold on; quiver(0, 0, cos(xhat_angle), sin(xhat_angle), 0); % The final zero in the quiver call turns off the automatic scaling.
        text(cos(xhat_angle), sin(xhat_angle), sprintf('%d', istr), 'FontSize', 18, 'HorizontalAlignment', 'right'); 
    end
    
    xl = xlim; yl = ylim; 
    xlim([xl(1) - 0.5, xl(2) + 0.5]); 
    ylim([yl(1) - 0.5, yl(2) + 0.5]); 
    

end

function plotRangeImg(y_r, phi_refl, k_band, pointSources)     
    uinf = mi__SARUtils.create_uinf_fromPoints(k_band, pointSources, phi_refl); 
    I = mi__SARUtils.buildRangeCompressionImg1D(y_r, k_band, uinf); 

    plot(y_r, abs(I), 'b.'); 
    xlabel('y_r')
    grid on
    
    % see "-i10": positive y_r are closer to the source 
    for ij = 1:numel(pointSources)
        ps = pointSources(ij); 
        
        dist = ps.pos(1) * cos(phi_refl) ... 
             + ps.pos(2) * sin(phi_refl); 
        
        hold on; plot(dist, 0, 'rx'); 
    end
end 

