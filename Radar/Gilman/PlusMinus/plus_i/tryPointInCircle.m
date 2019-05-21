% TODO: arrange it all as actual tests? 
function tryPointInCircle
    set(0, 'defaultLineLineWidth', 2); 
    set(0, 'defaultLineMarkerSize', 15);     
    set(0, 'defaultAxesFontSize', 20); 
    
    rCenter = [1.2; 1.7]; % center of the circle where we set up a mesh 
    R = 1.3;              % radius ot the circle 
    Ntheta = 1500;        % # of DIFFERENT angles in the mesh 
    
    Nphi = 1300;    % # of directions of outgoing waves

    k = 4.3; % common freq
    
    phi = 2 * pi * (0:(Nphi-1)) / Nphi;
    mesh = createMesh(rCenter, R, Ntheta); 

    % SOME TESTS
    testAllFFPs(mesh, k, phi); 
    %fprintf('Tests disabled!\n'); 
    
    % I could have made tests from these as well:
    %
    %plotNormDerivativeAtMesh(mesh, field);  
    %displayFieldAndNormDeriv(pointSrc, mesh);
    
    
    src1 = SARUtils.createPointSrc(2.6, rCenter + [0.3; 0.7]); 
    src2 = SARUtils.createPointSrc(2.3, rCenter + [0; 0]); 
    field = getFieldFromTwoPoints(src1, src2, mesh, k);     
    
    traj = doFFP_circle(k, phi, mesh, field); 
    plotTraj(traj); 
    suptitle(sprintf('Ntheta = %d, rCenter = (%g, %g)', Ntheta, rCenter(1), rCenter(2))); 
 
    %%%%%%%%%%%
end

function testAllFFPs(mesh, k, phi)
    % when the source is just 0.03 inside of the sphere, the 
    % field is still pretty big 
    testFFPfromPointInside(mesh, k, phi, mesh.rCenter + (mesh.R - 0.03) * [0; 1]); 
    testFFPfromPointInside(mesh, k, phi, mesh.rCenter + (mesh.R - 0.03) * [1; 0]); 
    testFFPfromPointInside(mesh, k, phi, mesh.rCenter + (mesh.R - 0.03) * [cos(1); sin(1)]); 
    testFFPfromPointInside(mesh, k, phi, mesh.rCenter + (mesh.R - 0.03) * [cos(2); sin(2)]);  
    testFFPfromPointInside(mesh, k, phi, mesh.rCenter + (mesh.R - 0.03) * [cos(3); sin(3)]);     
    
    % when the source is just 0.03 outside of the sphere, the 
    % field drops to near zero
    testFFPfromPointOutside(mesh, k, phi, mesh.rCenter + (mesh.R + 0.03) * [0; 1]); 
    testFFPfromPointOutside(mesh, k, phi, mesh.rCenter + (mesh.R + 0.03) * [1; 0]); 
    testFFPfromPointOutside(mesh, k, phi, mesh.rCenter + (mesh.R + 0.03) * [cos(1); sin(1)]); 
    testFFPfromPointOutside(mesh, k, phi, mesh.rCenter + (mesh.R + 0.03) * [cos(3); sin(3)]);  
    testFFPfromPointOutside(mesh, k, phi, mesh.rCenter + (mesh.R + 0.03) * [cos(5); sin(5)]); 
    
    testFFPfromPlaneWave(mesh, k, phi);     
end

function testFFPfromPlaneWave(mesh, k, phi)

    % when we use the trace of a plane wave as a source, 
    % the resulting FFP is zero (can be shown analytically)
    planeWave.ampl = 3.3; 
    planeWave.k = k;  
    planeWave.dirAngle = deg2rad(32); 
    planeWave.phase = deg2rad(58); 

    field = getPlaneWaveField(planeWave, mesh);   
    traj = doFFP_circle(planeWave.k, phi, mesh, field); 
    assert(all( abs(traj.ffp) < 1e-6 )); 
end

function testFFPfromPointInside(mesh, k, phi, pointPos)

    src_inside = SARUtils.createPointSrc(2.3, pointPos); 
    field = getFieldFromPoint(src_inside, mesh, k);  
    traj = doFFP_circle(k, phi, mesh, field); 
    assert(all( abs(traj.ffp) > 0.1 )); % replaced 1e2 with 0.1 due to implementation of curve_dl in SARUtils.doFFP_curve
end

function testFFPfromPointOutside(mesh, k, phi, pointPos)

    src_outside = SARUtils.createPointSrc(2.3, pointPos); 
    field = getFieldFromPoint(src_outside, mesh, k);   
    traj = doFFP_circle(k, phi, mesh, field); 
    assert(all( abs(traj.ffp) < 1e-6 )); 
end

function plotNormDerivativeAtMesh(mesh, field)

    figure; 
    plot(mesh.theta, real(field.normal_deriv), 'r-.', ... 
         mesh.theta, imag(field.normal_deriv), 'b-.'); 
    legend('real', 'imag'); title('normal deriv')
end

function displayFieldAndNormDeriv(pointSrc, mesh)
    
    Ntheta = numel(mesh.theta); 
    R_range = 0.1 : 0.1 : 4; 
    NR = numel(R_range); 
        
    x = nan(NR * Ntheta, 1);     
    y = nan(NR * Ntheta, 1); 
    v = nan(NR * Ntheta, 1); 
    vn = nan(NR * Ntheta, 1);     
    
    for ii = 1:NR
        tmp_mesh = createMesh(R_range(ii), Ntheta); % add pointSrc
        field = getFieldFromPoint(pointSrc, tmp_mesh, k);  
       
        point_range = ((ii-1)*Ntheta + 1) : (ii*Ntheta);
        ADD rCenter!!
        x(point_range) = tmp_mesh.R * cos(tmp_mesh.theta); 
        y(point_range) = tmp_mesh.R * sin(tmp_mesh.theta);      
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
    
    %testMeshAndNormals(mesh); 
    %testMeshAndField(mesh, field); 
end

function mesh = createMesh(rCenter, R, Ntheta)

    mesh.rCenter = rCenter; 
    mesh.R = R;
    
    mesh.theta = 2 * pi * (0:(Ntheta-1)) / Ntheta; 
    mesh.normals = [cos(mesh.theta); sin(mesh.theta)]; 
end
    
function traj = doFFP_circle(k, phi, mesh, field)
    %fprintf('why displacement does not just add a phase?\n'); 
    %curve.z  = mesh.rCenter + mesh.R * [cos(mesh.theta); sin(mesh.theta)]; 
    curve.z  = [200;-350] + mesh.R * [cos(mesh.theta); sin(mesh.theta)]; 
    %curve.z  = mesh.R * [cos(mesh.theta); sin(mesh.theta)]; 
    curve.nz = mesh.normals;     
    traj = SARUtils.doFFP_curve(k, phi, curve, field); 
end


function plotTraj(traj)
    
    figure('units', 'normalized', 'position', [0.02 0.1, 0.7, 0.4])
    subplot(121); plot(traj.phi, abs(traj.ffp), 'b-.'); grid on; title('abs')
    subplot(122); plot(traj.phi, angle(traj.ffp), 'b-.'); grid on; title('phase')    

end


function field = getFieldFromTwoPoints(src1, src2, mesh, k)

    fld1 = getFieldFromPoint(src1, mesh, k); 
    fld2 = getFieldFromPoint(src2, mesh, k); 

    field.value = fld1.value + fld2.value;
    field.gradient = fld1.gradient + fld2.gradient;
    field.normal_deriv = fld1.normal_deriv + fld2.normal_deriv;    
end

function field = getFieldFromPoint(pointSrc, mesh, k)

    meshCoords = mesh.rCenter + mesh.R * [cos(mesh.theta); sin(mesh.theta)]; 
    srcToMeshVects = meshCoords - repmat(pointSrc.pos, 1, numel(mesh.theta)); % ??? IT USED TO WORK WITHOUT REPMAT 
    distToSrc  = sqrt(sum(srcToMeshVects.^2,1));
    srcToMeshUnitVects = srcToMeshVects ./ [distToSrc; distToSrc]; % I could do just ./ distToSrc

    field.value    = (1i/4) * pointSrc.ampl * besselh(0, k * distToSrc); 
    deriv          = (1i/4) * pointSrc.ampl * besselh(1, k * distToSrc) * (-k);  
    field.gradient = [deriv .* srcToMeshUnitVects(1,:); ... 
                      deriv .* srcToMeshUnitVects(2,:)]; 
    field.normal_deriv = sum(field.gradient .* mesh.normals, 1); 
end

function field = getPlaneWaveField(planeWave, mesh)

    kIncVector = planeWave.k * [cos(planeWave.dirAngle); sin(planeWave.dirAngle)]; 
    % dot(,) below is OK because both args are real    
    kz = dot(kIncVector, mesh.rCenter) ... 
       + planeWave.k * mesh.R * cos(planeWave.dirAngle - mesh.theta);    

    
    field.value    = planeWave.ampl * exp(1i * (kz + planeWave.phase));  
    field.gradient = 1i * kIncVector * field.value; % only for plane wave 
    field.normal_deriv = sum(field.gradient .* mesh.normals, 1); 
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
        
        ADD rCenter
        meshx = mesh.R * cos(mesh.theta(it)); 
        meshy = mesh.R * sin(mesh.theta(it));         
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
        ADD rCenter        
        meshx = mesh.R * cos(mesh.theta(it)); 
        meshy = mesh.R * sin(mesh.theta(it));         
        hold on; plot(meshx, meshy, 'b*'); 
        
        hold on; quiver(meshx, meshy, mesh.normals(1, it), mesh.normals(2, it)); 
                      
    end 
    daspect([1 1 1])
end