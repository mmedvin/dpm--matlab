% TODO: arrange it all as actual tests? 
function mi__tryPointInCircle
%     set(0, 'defaultLineLineWidth', 2); 
%     set(0, 'defaultLineMarkerSize', 15);     
%     set(0, 'defaultAxesFontSize', 20); 
    
    %circularMesh = struct('rCenter', [1.2; 1.7], 'R', 1.3);  
    circularMesh = struct('rCenter', [1.2; 1.7], 'R', 3);  

    Ntheta = 1500;        % # of DIFFERENT angles in the mesh     
    Nphi = 1300;          % # of directions of outgoing waves

    k = 4.3; % common freq
    twoPointSetup.ampls = [2.6, 2.3]; 
    twoPointSetup.firstPointShift = [0.3; 0.7]; 
    twoPointSetup.rectShiftCoeffs = [0.2, 0.3];
    twoPointSetup.rectSizeCoeffs = [1.7, 1];
    twoPointSetup.Nvert = 400; 
    twoPointSetup.Nhoriz = 700; 
    twoPointSetup.radiusMultFactor = 1.8; 
    
    phi = 2 * pi * (0:(Nphi-1)) / Nphi;
    
    
    % SOME TESTS
    %warning('testFFPsInsideAndOutsideCircle and testFFPfromPlaneWave disabled!'); 
    testFFPsInsideAndOutsideCircle(circularMesh, Ntheta, k, phi); 
    testFFPfromPlaneWave(circularMesh, Ntheta, k, phi);      
    testTwoPoints(circularMesh, Ntheta, k, phi, twoPointSetup); 
    
    warning('Lengthy graphics disabled')
    %warning('Some obscure graphics follows, takes several minutes to build')    
    %plotNormDerivativeAtMesh(mesh, field);    %!!! LENGTHY
    %displayFieldAndNormDeriv(src1, mesh, k);  %!!! LENGTHY 
    
    %%%%%%%%%%%
    

end
    
function corners = createCorners(circularMesh, twoPointSetup)
    % right, left, top, bottom
    corners.RT = circularMesh.rCenter(1) + (twoPointSetup.rectShiftCoeffs(1) + twoPointSetup.rectSizeCoeffs(1)) * circularMesh.R; 
    corners.LT = circularMesh.rCenter(1) + (twoPointSetup.rectShiftCoeffs(1) - twoPointSetup.rectSizeCoeffs(1)) * circularMesh.R;  
    
    corners.TP = circularMesh.rCenter(2) + (twoPointSetup.rectShiftCoeffs(2) + twoPointSetup.rectSizeCoeffs(2)) * circularMesh.R;     
    corners.BT = circularMesh.rCenter(2) + (twoPointSetup.rectShiftCoeffs(2) - twoPointSetup.rectSizeCoeffs(2)) * circularMesh.R; 
end

function rectCurve = createRectangularCurve(corners, twoPointSetup)
    % avoid repeating coords in the corners by adding one extra point and the taking centers of the intervals 
    horizOneExtra = linspace(corners.LT, corners.RT, (twoPointSetup.Nhoriz + 1));
    horizIncr = (horizOneExtra(1:end-1) + horizOneExtra(2:end))/2; 
    horizDecr = fliplr(horizIncr); 
    
    vertOneExtra = linspace(corners.BT, corners.TP, (twoPointSetup.Nvert + 1)); 
    vertIncr  = (vertOneExtra(1:end-1) + vertOneExtra(2:end))/2;
    vertDecr  = fliplr(vertIncr); 
        
    % bottom, right, top, left: coordinates in CCW direction
    bottom = [horizIncr;  ones(1, twoPointSetup.Nhoriz) * corners.BT]; 
    right  = [ones(1, twoPointSetup.Nvert) * corners.RT; vertIncr]; 
    top    = [horizDecr; ones(1, twoPointSetup.Nhoriz) * corners.TP]; 
    left   = [ones(1, twoPointSetup.Nvert) * corners.LT; vertDecr];
    
    % normals 
    bottom_nz = repmat([0; -1], [1, twoPointSetup.Nhoriz]); 
    top_nz    = repmat([0;  1], [1, twoPointSetup.Nhoriz]);     
    
    left_nz   = repmat([-1; 0], [1, twoPointSetup.Nvert]);     
    right_nz  = repmat([ 1; 0], [1, twoPointSetup.Nvert]); 
    
    % assemble in CCW order
    rectCurve.z  = [bottom, right, top, left]; 
    rectCurve.nz = [bottom_nz, right_nz, top_nz, left_nz];
    
end

function displayRectCurve(rectCurve, twoPointSetup)
    
    [rows_z, cols_z] = size(rectCurve.z);    assert(rows_z == 2); 
    [rows_nz, cols_nz] = size(rectCurve.nz); assert(rows_nz == 2);
    
    assert(cols_z == cols_nz); 
    assert(cols_z == 2 * (twoPointSetup.Nvert + twoPointSetup.Nhoriz)); 
        
    figure('units', 'normalized', 'position', [0.2 0.2 0.7 0.3]); 
    
    subplot(121)
    hold on; plot(1:cols_z, rectCurve.z(1,:), 'b.'); 
    hold on; plot(1:cols_z, rectCurve.z(2,:), 'm.'); 
    hold on; plot(1:cols_z, rectCurve.nz(1,:), 'k-'); 
    hold on; plot(1:cols_z, rectCurve.nz(2,:), 'k--'); 
    grid on
    legend('x', 'y', 'nx', 'ny', 'location', 'northwest')
    
    subplot(122)
    plot(rectCurve.z(1,:), rectCurve.z(2,:), 'b.')
    
    mask = 1 : 20 : cols_z; 
    hold on; 
    quiver(rectCurve.z(1,mask), rectCurve.z(2,mask), rectCurve.nz(1,mask), rectCurve.nz(2,mask)); 
    grid on; axis equal
    
    sgtitle('Coordinates and shape of rectangular contour')
end

function circularCurve = createCircularCurve(circularMesh, Ntheta)
    theta = 2 * pi * (0:(Ntheta-1)) / Ntheta; 
    circularCurve.z  = circularMesh.rCenter + circularMesh.R * [cos(theta); sin(theta)]; 
    circularCurve.nz = [cos(theta); sin(theta)]; 
end

function traj = doFFP_curve(k, phi, curve, field)
  
    traj = SAR.minus_i.mi__SARUtils.doFFP_curve(k, phi, curve, field); 
end

function assertPointWithinRect(src, corners)
    assert( (corners.LT < src.pos(1)) && (src.pos(1) < corners.RT) && ... 
            (corners.BT < src.pos(2)) && (src.pos(2) < corners.TP)); 
end

function testTwoPoints(circularMesh, Ntheta, k, phi, twoPointSetup)

    assert(sum(twoPointSetup.firstPointShift.^2) < circularMesh.R^2); % we want two points inside the circle
    src1 = SAR.minus_i.mi__SARUtils.createPointWithAmplitude(twoPointSetup.ampls(1), circularMesh.rCenter + twoPointSetup.firstPointShift); 
    src2 = SAR.minus_i.mi__SARUtils.createPointWithAmplitude(twoPointSetup.ampls(2), circularMesh.rCenter);     
    
    % "Regular curve" 
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.85])
    
    circularCurve = createCircularCurve(circularMesh, Ntheta);   
    field = getFieldFromTwoPoints(src1, src2, circularCurve, k);     
    traj = doFFP_curve(k, phi, circularCurve, field); 
    plotTraj(traj, 1, 'circle'); 

    %%% "Big curve" %%%%%%%%%%%%%%%%%%
    
    circularMeshXX = circularMesh; 
    circularMeshXX.R = twoPointSetup.radiusMultFactor * circularMeshXX.R; 
    
    circularCurveXX = createCircularCurve(circularMeshXX, Ntheta);   
    fieldXX = getFieldFromTwoPoints(src1, src2, circularCurveXX, k);     
    trajXX = doFFP_curve(k, phi, circularCurveXX, fieldXX); 
    plotTraj(trajXX, 2, 'circleXX'); 
                 
                  
    %%%%%% Rectangular contour %%%%%%%%%%%%%%%
    corners = createCorners(circularMesh, twoPointSetup); 
    assertPointWithinRect(src1, corners); 
    assertPointWithinRect(src2, corners);     
    rectCurve = createRectangularCurve(corners, twoPointSetup); 

    field_rect = getFieldFromTwoPoints(src1, src2, rectCurve, k);     
    traj_rect = doFFP_curve(k, phi, rectCurve, field_rect); 
    plotTraj(traj_rect, 3, 'rectangle');     
                                
    sgtitle(sprintf('circular: Ntheta = %d, rCenter = (%g, %g), R = %g (circleXX: R * %g) \n rectangular, horiz = (%g--%g), vert = (%g--%g)', ... 
                      Ntheta, circularMesh.rCenter(1), circularMesh.rCenter(2), circularMesh.R, ... 
                      twoPointSetup.radiusMultFactor, ... 
                      corners.LT, corners.RT, corners.BT, corners.TP)); 
    %%%%%%%                            
    diffCircRel = max(abs(traj.ffp - trajXX.ffp)) / max(abs(traj.ffp)); 
    diffRectRel = max(abs(traj.ffp - traj_rect.ffp)) / max(abs(traj.ffp)); 
    assert( (diffCircRel < 0.02) && (diffRectRel < 0.02))
    fprintf('Relative difference in FFP: for circles: %g, circle vs. rectangle: %6.4f\n', diffCircRel, diffRectRel); 
    %%%%%%%%
    
    
    % Extra test...
    displayRectCurve(rectCurve, twoPointSetup); 
    
    figure,plot(traj.phi_refl,traj.ffp,'bo',trajXX.phi_refl,trajXX.ffp,'r+',traj_rect.phi_refl,traj_rect.ffp,'g^')
    legend('circle','circleXX','rectangle')
    
    assert(all(traj.phi_refl - trajXX.phi_refl   <1e-10))
    assert(all(traj.phi_refl - traj_rect.phi_refl<1e-10))
    norm(traj.ffp - trajXX.ffp   ,inf)
    norm(traj.ffp - traj_rect.ffp,inf)
    
end


function testFFPsInsideAndOutsideCircle(circularMesh, Ntheta, k, phi)
    fprintf('TEST: testFFPsInsideAndOutsideCircle ... '); 
    
    circularCurve = createCircularCurve(circularMesh, Ntheta); 

    % when the source is just 0.03 inside of the sphere, the 
    % field is still pretty big 
    testFFPfromPointInside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R - 0.03) * [0; 1]); 
    testFFPfromPointInside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R - 0.03) * [1; 0]); 
    testFFPfromPointInside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R - 0.03) * [cos(1); sin(1)]); 
    testFFPfromPointInside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R - 0.03) * [cos(2); sin(2)]);  
    testFFPfromPointInside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R - 0.03) * [cos(3); sin(3)]);     
    
    % when the source is just 0.03 outside of the sphere, the 
    % field drops to near zero
    testFFPfromPointOutside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R + 0.03) * [0; 1]); 
    testFFPfromPointOutside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R + 0.03) * [1; 0]); 
    testFFPfromPointOutside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R + 0.03) * [cos(1); sin(1)]); 
    testFFPfromPointOutside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R + 0.03) * [cos(3); sin(3)]);  
    testFFPfromPointOutside(circularCurve, k, phi, circularMesh.rCenter + (circularMesh.R + 0.03) * [cos(5); sin(5)]); 
       
    fprintf(' SUCCESS!\n'); 
end

function testFFPfromPlaneWave(circularMesh, Ntheta, k, phi)
    fprintf('TEST: testFFPfromPlaneWave ... '); 
    
    circularCurve = createCircularCurve(circularMesh, Ntheta);     
    
    % when we use the trace of a plane wave as a source, 
    % the resulting FFP is zero (can be shown analytically)
    planeWave.ampl = 3.3; 
    planeWave.k = k;  
    planeWave.dirAngle = deg2rad(32); 
    planeWave.phase = deg2rad(58); 

    field = getPlaneWaveField(planeWave, circularCurve);   
    traj = doFFP_curve(planeWave.k, phi, circularCurve, field); 
    assert(all( abs(traj.ffp) < 1e-6 )); 
    
    fprintf(' SUCCESS!\n');     
end

function testFFPfromPointInside(circularCurve, k, phi, pointPos)
    src_inside = SAR.minus_i.mi__SARUtils.createPointWithAmplitude(2.3, pointPos); 
    field = getFieldFromPoint(src_inside, circularCurve, k);  
    traj = doFFP_curve(k, phi, circularCurve, field); 
    assert(all( abs(traj.ffp) > 0.1 )); % replaced 1e2 with 0.1 due to implementation of curve_dl in mi__SARUtils.doFFP_curve
end

function testFFPfromPointOutside(circularCurve, k, phi, pointPos)
    src_outside = SAR.minus_i.mi__SARUtils.createPointWithAmplitude(2.3, pointPos); 
    field = getFieldFromPoint(src_outside, circularCurve, k);   
    traj = doFFP_curve(k, phi, circularCurve, field); 
    assert(all( abs(traj.ffp) < 1e-6 )); % no field from a point outside
end

function plotNormDerivativeAtMesh(mesh, field)

    figure; 
    plot(mesh.theta, real(field.normal_deriv), 'r-.', ... 
         mesh.theta, imag(field.normal_deriv), 'b-.'); 
    legend('real', 'imag'); title('normal deriv')
end

function displayFieldAndNormDeriv(pointSrc, mesh, k)
    
    Ntheta = numel(mesh.theta); 
    R_range = 0.1 : 0.1 : 4; 
    NR = numel(R_range); 
        
    x = nan(NR * Ntheta, 1);     
    y = nan(NR * Ntheta, 1); 
    v = nan(NR * Ntheta, 1); 
    vn = nan(NR * Ntheta, 1);     
    
    for ii = 1:NR
        
        tmp_mesh = createMesh(mesh.rCenter, R_range(ii), Ntheta); % add pointSrc
        field = getFieldFromPoint(pointSrc, tmp_mesh, k);  
       
        point_range = ((ii-1)*Ntheta + 1) : (ii*Ntheta);

        x(point_range) = mesh.rCenter(1) + tmp_mesh.R * cos(tmp_mesh.theta); 
        y(point_range) = mesh.rCenter(2) + tmp_mesh.R * sin(tmp_mesh.theta);      
        v(point_range) = field.value; 
        vn(point_range) = field.normal_deriv;         
    end
    

    [xq,yq] = meshgrid(-max(R_range) : 0.1 : max(R_range)); 
    
    %%%%%%%%
    figure('name', 'fun', 'units', 'normalized', 'position', [0.02 0.1, 0.7, 0.4])
    subplot(131)
    realF = scatteredInterpolant(x,y,real(v), 'linear', 'none');      
    pcolor(xq, yq, realF(xq, yq)); shading flat; colorbar; daspect([1 1 1]); title('Real')
  
    subplot(132)
    imagF = scatteredInterpolant(x,y,real(v), 'linear', 'none');      
    pcolor(xq, yq, imagF(xq, yq)); shading flat; colorbar; daspect([1 1 1]); title('Imag')    
    
    subplot(133)
    absF = scatteredInterpolant(x,y,log(abs(v)), 'linear', 'none');      
    pcolor(xq, yq, absF(xq, yq)); shading flat; colorbar; daspect([1 1 1]); title('logAbs')       
    
    %%%%%%%
    figure('name', 'normal deriv', 'units', 'normalized', 'position', [0.02 0.1, 0.7, 0.4])
    subplot(131)
    realF = scatteredInterpolant(x,y,real(vn), 'linear', 'none');      
    pcolor(xq, yq, realF(xq, yq)); shading flat; colorbar; daspect([1 1 1]); title('Real')
  
    subplot(132)
    imagF = scatteredInterpolant(x,y,real(vn), 'linear', 'none');      
    pcolor(xq, yq, imagF(xq, yq)); shading flat; colorbar; daspect([1 1 1]); title('Imag')    
    
    subplot(133)
    absF = scatteredInterpolant(x,y,log(abs(vn)), 'linear', 'none');      
    pcolor(xq, yq, absF(xq, yq)); shading flat; colorbar; daspect([1 1 1]); title('logAbs')        
    
    testMeshAndNormals(mesh); 
    testMeshAndField(mesh, field); 
end   


function plotTraj(traj, plotRow, name)
   
    subplot(3, 2, (plotRow - 1) * 2 + 1); 
    plot(traj.phi_refl, abs(traj.ffp), 'b-.'); grid on; title(['abs(FFP): ' name]); 
    subplot(3, 2, (plotRow - 1) * 2 + 2);
    plot(traj.phi_refl, angle(traj.ffp), 'b-.'); grid on; title(['phase(FFP): ' name]);     

end


function field = getFieldFromTwoPoints(src1, src2, curve, k)

    fld1 = getFieldFromPoint(src1, curve, k); 
    fld2 = getFieldFromPoint(src2, curve, k); 

    field.value = fld1.value + fld2.value;
    field.gradient = fld1.gradient + fld2.gradient;
    field.normal_deriv = fld1.normal_deriv + fld2.normal_deriv;    
end

function field = getFieldFromPoint(pointSrc, curve, k)

    [rows, Npoints] = size(curve.z); assert(rows == 2); 
    srcToCurveVects = curve.z - repmat(pointSrc.pos, 1, Npoints); % ??? IT USED TO WORK WITHOUT REPMAT 
    distToSrc  = sqrt(sum(srcToCurveVects.^2,1));
    srcToCurveUnitVects = srcToCurveVects ./ [distToSrc; distToSrc]; % I could do just ./ distToSrc

    besselKind = SAR.minus_i.mi__SARUtils.besselKind; 
   
    field.value    = (1i/4) * pointSrc.ampl * besselh(0, besselKind, k * distToSrc);    
    deriv          = (1i/4) * pointSrc.ampl * besselh(1, besselKind, k * distToSrc) * (-k);  % see DLMF 10.6.3: derivative of (J,Y,H)_0 is -(J,Y,H)_1 
    field.gradient = [deriv .* srcToCurveUnitVects(1,:); ... 
                      deriv .* srcToCurveUnitVects(2,:)]; 
    field.normal_deriv = sum(field.gradient .* curve.nz, 1); 
end

function field = getPlaneWaveField(planeWave, curve)
    [rows, Npoints] = size(curve.z); assert(rows == 2); 

    kIncVector = planeWave.k * [cos(planeWave.dirAngle); sin(planeWave.dirAngle)]; 
    kz = sum(repmat(kIncVector, [1, Npoints]) .* curve.z, 1);   
       
    field.value    = planeWave.ampl * exp(-1i * (kz + planeWave.phase));  % PMI see -i1,(a)
    field.gradient = -1i * kIncVector * field.value;                      % PMI         % this relation only holds for plane wave 
    field.normal_deriv = sum(field.gradient .* curve.nz, 1); 
end
    
function testMeshAndField(mesh, field)

    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.7]); 
    subplot(231); plotMeshAndMult(mesh, real(field.value)); title('Re(field)'); 
    subplot(234); plotMeshAndMult(mesh, imag(field.value)); title('Im(field)'); 

    subplot(232); plotMeshAndMult(mesh, real(field.gradient(1,:))); title('Re(grad x)'); 
    subplot(235); plotMeshAndMult(mesh, imag(field.gradient(1,:))); title('Im(grad x)');     
 
    subplot(233); plotMeshAndMult(mesh, real(field.gradient(2,:))); title('Re(grad y)'); 
    subplot(236); plotMeshAndMult(mesh, imag(field.gradient(2,:))); title('Im(grad y)');       
    
end
    
function plotMeshAndMult(mesh, mult) 
    for it = 1:numel(mesh.theta)
        multval = mult(it); 
        

        meshx = mesh.rCenter(1) + mesh.R * cos(mesh.theta(it)); 
        meshy = mesh.rCenter(2) + mesh.R * sin(mesh.theta(it));         
        hold on; plot(meshx, meshy, 'b*'); 
        
        hold on; quiver(meshx, meshy, multval * mesh.normals(1, it), multval * mesh.normals(2, it), ... 
                         'color', 'b', 'ShowArrowHead', 'off'); 
                      
    end 
    daspect([1 1 1])

end 

function testMeshAndNormals(mesh)

%    figure; 
%    hold on; plot(mesh.theta, real(field.value), 'r-.', ... 
%                  mesh.theta, imag(field.value), 'b-.'); 

    figure; 
    for it = 1:numel(mesh.theta)
    
        meshx = mesh.rCenter(1) + mesh.R * cos(mesh.theta(it)); 
        meshy = mesh.rCenter(2) + mesh.R * sin(mesh.theta(it));         
        hold on; plot(meshx, meshy, 'b*'); 
        
        hold on; quiver(meshx, meshy, mesh.normals(1, it), mesh.normals(2, it)); 
                      
    end 
    daspect([1 1 1])
end