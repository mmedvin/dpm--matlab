function trySARfromPointScatterers
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    % settings as of 2018
    %phi_range = pi + (-0.05 : 0.005 : 0.05); 
    %k_band = 50 : 0.25 : 52;     
    %plotLims.rng = -5:0.1:5; plotLims.az  = -3:0.1:3; 
    
    phi_range = pi + (-0.05 : 0.005 : 0.05); 
    k_band = 50 : 0.25 : 52; 
    plotLims.rng = -15:0.1:14; plotLims.az  = -13:0.1:14.3;
    
    
    pointSources(1) = SARUtils.createPointSrc(2.3, [2; 0.5]); 
    pointSources(2) = SARUtils.createPointSrc(1.7, [-1; -1]); 
    pointSources(3) = SARUtils.createPointSrc(1.1, [-2.5; 1.6]);     
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.8,0.6], 'paperpositionmode', 'auto');    
    plotFullImg(phi_range, k_band, pointSources, plotLims); 


end

function plotFullImg(phi_range, k_band, pointSources, plotLims) 

    assert(isvector(phi_range)); 
    assert(isvector(k_band));    
   
    [y1, y2] = meshgrid(plotLims.rng, plotLims.az); 
    
    % create field for each k 
    uinf = zeros(numel(phi_range), numel(k_band)); 
    for im = 1:numel(phi_range)
        phi = phi_range(im); 
        x_hat = [cos(phi); sin(phi)]; 
        for il = 1:numel(k_band)
            for ij = 1:numel(pointSources)
                ps = pointSources(ij); 
                % dot(,) below is OK because both args are real
                uinf(im, il) = uinf(im, il) + ps.ampl * exp(-2i * k_band(il) * dot(x_hat, ps.pos)); 
            end
        end
    end
    
    I = SARUtils.buildSARimage(y1, y2, k_band, phi_range, uinf); 


    %pcolor(y1, y2, abs(I)); colorbar
    contourf(y1, y2, abs(I)); 
    title('SAR image; crosses indicate point scatterers') 
    hc = colorbar; 
    ylabel(hc, 'magnitude')  
    daspect([1 1 1])
    xlabel('x'); ylabel('y')
    
    for ij = 1:numel(pointSources)
        ps = pointSources(ij); 
        indicatePoint(ps.pos(1), ps.pos(2)); 
    end
    print('-djpeg', 'SAR_pointScatt.jpg')
end 

function indicatePoint(x, y) 

    barLength = 0.3; 
    hold on; plot([x-barLength, x+barLength], [y, y], 'r'); 
    hold on; plot([x, x], [y-barLength, y+barLength], 'r'); 
    
end


