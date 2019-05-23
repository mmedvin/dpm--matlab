function AnalyticData(AR,ispec)
% AR is an aspect ratio of an ellipse
    if nargin==0, AR=2; end
    if nargin < 2
        spec = getSpec();
    else
        spec=ispec;
    end
    
    
    ellipse.a=1;%major axis
    ellipse.b=ellipse.a/AR;%minor axis
 
    theta = spec.theta; % elliptical angle
    
    x = ellipse.a * cos(theta);
    y = ellipse.b * sin(theta);

    % xi is elliptical radius
    % focal distance is f=sqrt(a^2-b^2)
    % hence a=f cosh(xi) and b=f sinh(xi), 
    dxdxi = ellipse.b * cos(theta); %dx/dxi
    dydxi = ellipse.a * sin(theta); %dy/dxi
    
    curve.z =  [x; y];
    %hxi = sqrt(	(ellipse.a).^2 * sin(theta).^2 + (ellipse.b).^2 * cos(theta).^2	); the same as sqrt(dxdxi.^2 + dydxi.^2)
    curve.nz = [ dxdxi; dydxi]./sqrt(dxdxi.^2 + dydxi.^2); %unit normal vector

    
    for im = 1:numel(spec.phi_range)
        for il = 1:numel(spec.k_range)

            k = spec.k_range(il);
            phi = spec.phi_range(im); %incident angle
            
            PlaneWave = -exp( 1i.* k .* (x.*cos(phi) + y.*sin(phi)) ); 
            
            if 1 % dirichlet
                scattFields(im,il).value = PlaneWave;
                scattFields(im,il).normal_deriv = ell_du_dn_exact_SM(ellipse.a,ellipse.b,k,theta,phi)';
            else %neumann             
                scattFields(im,il).value = ell_hard_exact_Sol(ellipse.a,ellipse.b,k,phi,theta)';%phi and theta not as in previouse one! 
                scattFields(im,il).normal_deriv = 1i .* k * (dxdxi.*cos(phi) + dydxi.*sin(phi)).*  PlaneWave;
            end
            
        end
    end
    
    
    save('AnalyticEllipse', 'curve', 'scattFields', 'ellipse')
    
end

function spec = getSpec()
    
    spec.k_range = 50:0.5:55;
    spec.phi_range = pi * (-0.2 : 0.004 : 0.2);
    spec.theta = pi * (0 : 0.0025 : 1.999);
    
end