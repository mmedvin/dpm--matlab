classdef mi__SARUtils
    
    properties (Constant)          
        besselKind = 2; % PMI
    end
    
    methods (Static)
        function pointSrc = createPointWithAmplitude(ampl, pos)
            pointSrc.ampl = ampl; 
            pointSrc.pos = pos; 
        end        
        

        function uinf = create_uinf_fromPoints(k_band, pointSources, phi_refl)
            % what is xhat? is it the argument of FFP, i.e., the direction of _reflected_ field, 
            % which is the same as direction TO the source
            x_hat = [cos(phi_refl); sin(phi_refl)]; 
            uinf = zeros(size(k_band)); 
            for il = 1:numel(uinf)
                for ij = 1:numel(pointSources)
                    ps = pointSources(ij); 
                    % dot(,) below is OK because both args are real            
                    uinf(il) = uinf(il) + ps.ampl * exp(2i * k_band(il) * dot(x_hat, ps.pos)); % PMI see -i10,(a)
                end
            end
        end
        
        function I = buildRangeCompressionImg1D(y_r, k_band, uinf)
            % x_hat is the argument of FFP, so its direction is that of the _reflected_ wave, 

            % create image for these 2D
            I = zeros(size(y_r)); 
            for il = 1:numel(uinf)
                I = I + uinf(il) * exp(- 2i * k_band(il) * y_r); % PMI  see -i10,(b)
            end
        end
        
        function I = buildRangeCompressionImg2D(y2D, k_band, x_hat, uinf)
            % x_hat is the argument of FFP, so its direction is that of the _reflected_ wave

            I = zeros(size(y2D.y1));  

            for il = 1:numel(uinf)
                dotprod = x_hat(1) * y2D.y1 + x_hat(2) * y2D.y2; 
                I = I + uinf(il) * exp(-2i * k_band(il) * dotprod); %PMI see -i10,(c)
            end  
        end
        
        % phi_range corresponds to scattered field - see expression for x_hat
        function I = buildSARimage(y2D, k_band, phi_refl_range, uinf)

            I = zeros(size(y2D.y1)); 

            for im = 1:numel(phi_refl_range)
                phi_refl = phi_refl_range(im); 
                x_hat = [cos(phi_refl); sin(phi_refl)]; % scattered wave - don't need a minus

                for il = 1:numel(k_band)
                    dotprod = x_hat(1) * y2D.y1 + x_hat(2) * y2D.y2; 
                    I = I + uinf(im, il) * exp(-2i * k_band(il) * dotprod); % PMI see -i10,(c)
                end
            end  
        end
        
        function traj = doFFP_curve(k, phi_refl, curve, field)
            
            traj.phi_refl = phi_refl; 
            traj.ffp = nan(size(traj.phi_refl));    

            [rows_z, cols_z] = size(curve.z); assert(rows_z == 2); 

            curve_dl = mi__SARUtils.calculate_dl(curve.z); % the length element along the curve
            
            for iphi = 1:numel(traj.phi_refl)
                phi = traj.phi_refl(iphi); 

                % incident wave propagates from "-hat{x}"
                xhat = [cos(phi); sin(phi)]; 
                k_xhat_z = k * sum(repmat(xhat, 1, cols_z) .* curve.z, 1); 

                expInKirchhoffKernel.value = exp(1i * k_xhat_z);                             % PMI see -i4,(a)
                expInKirchhoffKernel.gradient = 1i * k * xhat * expInKirchhoffKernel.value;  % PMI see -i4,(b)            
                
                expInKirchhoffKernel.normal_deriv = sum(expInKirchhoffKernel.gradient .* curve.nz, 1);

                % don't use dot(,) for complex numbers! it makes conj() of the first argument! 
                sum_for_ffp = sum(curve_dl .* (  field.value .* expInKirchhoffKernel.normal_deriv ... 
                                               - field.normal_deriv .* expInKirchhoffKernel.value)); 

                % there is a factor of exp(1i * pi/4) / sqrt(8 * pi * k) in front of the integral 
                % I retain 1/sqrt(k) because I use different k's in SAR 
                traj.ffp(iphi) = sum_for_ffp / sqrt(k); % see -i4,(c)
            end
        end    
 
        function [scattField, q_inc_z] = getScattField(k, phi_xhat, reflector, isCheckShallow)
            [rows_z, cols_z] = size(reflector.curve.z); assert(rows_z == 2);
            xhat = [cos(phi_xhat); sin(phi_xhat)];
            
            k_inc = - k * xhat; 

            % BTW dot here OK because args are real    
            % same as: q_inc_z = k * sum(repmat(xhat, 1, colsTheta) .* nz, 1);   
            curve = reflector.curve; 
            q_inc_z = - sum(repmat(k_inc, 1, cols_z) .* curve.nz, 1); 
            if isCheckShallow
                assert(all( abs(q_inc_z / k) > (1/3) ), 'incident wave too shallow'); 
            end
                
            k_inc_times_z = sum(repmat(k_inc, 1, cols_z) .* curve.z, 1); 
             
            scattField.value = - k^2 * q_inc_z.^(-2) .* reflector.coeff .* exp(-1i * k_inc_times_z); % PMI
            scattField.normal_deriv = mi__SARUtils.scattFieldValueToNormalDeriv(scattField.value, q_inc_z); 
    
        end        
        
        function ret = scattFieldValueToNormalDeriv(value, qz)
            assert(all(imag(qz) == 0)); 
            ret = -1i * qz .* value; % PMI, see -i1,(a)
        end
        
        function n_z = getUnitVectorPerpTo(vect_in)
            % use 3D cross product to calculate normal vector 
            % - it comes rotated counterclockwise from (zEnd - zStart), 
            % then drop the third coord

            vectCross = cross([0;0;1], [vect_in; 0]); 
            n_z = [vectCross(1); vectCross(2)]/norm(vectCross);  
        end   
        
        function displayReflector(reflector)
            % see  https://www.mathworks.com/matlabcentral/answers/5042-how-do-i-vary-color-along-a-2d-line
            x_curve = reflector.curve.z(1,:);
            y_curve = reflector.curve.z(2,:);
            z_curve = zeros(size(x_curve)); % this is a trick! 
            col_curve = reflector.coeff;  % This is the color, vary with x in this case.
            surface([x_curve;x_curve],[y_curve;y_curve],[z_curve;z_curve],[col_curve;col_curve],...
                    'facecol','no',...
                    'edgecol','interp',...
                    'linew',4);
            daspect([1 1 1]); colorbar; 
            cl = caxis; caxis([0, cl(2)]); 
        end        
    end
       
    methods (Static, Access = private)
        function ret = isLoop(curve_z)
            rel_threshold = 0.1;  
            rel_diff_threshold = 0.3;
            minPointsInCurve = 10; 
            
            [rows, cols] = size(curve_z); assert(rows == 2); 
            assert(cols > minPointsInCurve, sprintf('Too few points in curve: %d', cols)); 
            dist_first = norm(curve_z(:,2)   - curve_z(:,1)); 
            dist_last  = norm(curve_z(:,end) - curve_z(:,(end-1))); 
            dist_chord = norm(curve_z(:,end) - curve_z(:,1)); 
           
            rel_first = dist_first / dist_chord; 
            rel_last  = dist_last  / dist_chord;             
            rel_diff_first = abs(dist_first - dist_chord) / dist_chord; 
            rel_diff_last  = abs(dist_last  - dist_chord) / dist_chord;      
           
            if     ((rel_diff_first < rel_diff_threshold) && (rel_diff_last < rel_diff_threshold)) 
                % the chord is about the same as the first and the last distances
                ret = 1; 
            elseif ((rel_first < rel_threshold) && (rel_last < rel_threshold))
                % the chord is much bigger than the first and the last distances
                ret = 0;
            else
                error('Uncertain curve type'); 
            end
            % TODO: need a test for it!!! 
        end
        
        function dl = calculate_dl(curve_z)
                       
            curve_z_prev = circshift(curve_z,  1, 2); 
            curve_z_next = circshift(curve_z, -1, 2); 
        
            vect_to_prev = curve_z_prev - curve_z; 
            vect_to_next = curve_z_next - curve_z;
            
            dist_to_prev = sqrt(sum(vect_to_prev.^2, 1)); 
            dist_to_next = sqrt(sum(vect_to_next.^2, 1)); 
            
            dl = (dist_to_prev + dist_to_next)/2; 
            
            if (~ mi__SARUtils.isLoop(curve_z)) 
                % chord is too big - remove it from dl 
                dl(1)   = dist_to_next(1)/2; 
                dl(end) = dist_to_prev(end)/2; 
            end
        end
        
    end
end
        