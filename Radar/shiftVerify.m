function shiftVerify
    tol = 1e-5; 

    setup.phi_refl = linspace(0, 2*pi, 987); 
    
    % from my spec (?)
    setup.shiftVec = [0.3; 0]; % this is how (x,y) are arranged
    
    % from MM's files
    setup.k           = 5;
    setup.IncAng_windConventions      = pi/5; % Medvinsky-style: "northerly" blows toward south
    
    load('./CenteredCircle.mat');     
    load('./ShiftedCircle.mat'); 

    checkNormalLengthAndAngle(CenteredCircle, tol); 
    checkNormalLengthAndAngle(ShiftedCircle, tol);     
    
    checkShiftOfCoords(CenteredCircle, ShiftedCircle, setup.shiftVec, tol); 
    % THIS SHOULD FAIL: checkShiftOfCoords(ShiftedCircle, CenteredCircle, setup.shiftVec, tol); 

    checkShiftOfField(CenteredCircle, ShiftedCircle, setup, tol); 
    
    checkShiftOfFFP(CenteredCircle, ShiftedCircle, setup, tol);
    
    fprintf('SUCCESS!\n'); 
    %{
    mask = 1:20:300; % 999; 

    figure; 
    plot(CenteredCircle.curve.z(1,mask), CenteredCircle.curve.z(2,mask), 'rx', ... 
         ShiftedCircle.curve.z(1,mask), ShiftedCircle.curve.z(2,mask), 'bo') 
    title('z')
    axis equal

    figure; 
    plot(CenteredCircle.curve.nz(1,mask), CenteredCircle.curve.nz(2,mask), 'rx', ... 
         ShiftedCircle.curve.nz(1,mask), ShiftedCircle.curve.nz(2,mask), 'bo') 
    title('nz')
    axis equal
    %}
    
end

function checkShiftOfFFP(CenteredCircle, ShiftedCircle, setup, tol)
    
    centeredFFP = SAR.minus_i.mi__SARUtils.doFFP_curve(setup.k, setup.phi_refl, CenteredCircle.curve, CenteredCircle.field); 
    shiftedFFP  = SAR.minus_i.mi__SARUtils.doFFP_curve(setup.k, setup.phi_refl, ShiftedCircle.curve, ShiftedCircle.field); 
    
    xhat = [cos(setup.phi_refl); sin(setup.phi_refl)];    
    k_inc = get_k_inc(setup); 
    
    % the second term is just a scalar - see checkShiftOfField below
    n_refl = numel(setup.phi_refl); 
    shiftMat = repmat(setup.shiftVec, 1, n_refl); 
    mult = exp(1i*( setup.k * dot(xhat, shiftMat, 1) - dot(k_inc, setup.shiftVec, 1) )); 
    
    approxShiftedFFP.ffp = centeredFFP.ffp .* mult; 
    diff = approxShiftedFFP.ffp - shiftedFFP.ffp; 
    
    assert(all( abs(diff) < tol )); 
    
    figure('units', 'normalized', 'position', [0.1 0.1 0.6 0.8]); 
    plot(setup.phi_refl, real(approxShiftedFFP.ffp), 'rx', setup.phi_refl, real(shiftedFFP.ffp), 'bo')
    
end

function checkShiftOfField(CenteredCircle, ShiftedCircle, setup, tol)

    k_inc = get_k_inc(setup); 
   
    mult = exp(-1i * dot(k_inc, setup.shiftVec, 1)); % see "-i5", C = e^{-i(k^i,R)}
    assert(numel(mult) == 1); 
    
    approxShiftedCircle.field.value        = CenteredCircle.field.value * mult; 
    approxShiftedCircle.field.normal_deriv = CenteredCircle.field.normal_deriv * mult;     
    
    field_diff = approxShiftedCircle.field.value - ShiftedCircle.field.value; 
    deriv_diff = approxShiftedCircle.field.normal_deriv - ShiftedCircle.field.normal_deriv; 
    
    assert(all( abs(field_diff) < tol )); 
    assert(all( abs(deriv_diff) < tol ));     

end

function checkNormalLengthAndAngle(S, tol)

    normal_length_lazy = sqrt(S.curve.nz(1,:).^2 + S.curve.nz(2,:).^2);     
    normal_length = sqrt(sum(S.curve.nz .^2, 1));      

    assert(all( abs(normal_length - normal_length_lazy) < tol )); 
    assert(all( abs(normal_length - 1) < tol )); 

    % see SARUtils.calculate_dl
    curve_z_prev = circshift(S.curve.z,  1, 2); 
    curve_z_next = circshift(S.curve.z, -1, 2); 

    dl = curve_z_next - curve_z_prev; 
    dl_length = sqrt(sum(dl.^2, 1));
    
    tangent = dl ./ repmat(dl_length, 2, 1); 
    
    assert(all( abs(sum(S.curve.nz .* tangent, 1)) < tol )); 
end

function checkShiftOfCoords(CenteredCircle, ShiftedCircle, shiftVec, tol)

    observedShift = ShiftedCircle.curve.z - CenteredCircle.curve.z; 
    
    [rows, cols] = size(observedShift); assert(rows == 2); 
    shiftMat = repmat(shiftVec, 1, cols); % prescribed 
    diff = observedShift - shiftMat; 
    assert(all( abs(diff(:)) < tol )); 

end

function k_inc = get_k_inc(setup)
    incAng = setup.IncAng_windConventions - pi; 
    k_inc = setup.k * [cos(incAng); sin(incAng)]; 
end