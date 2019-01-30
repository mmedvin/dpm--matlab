classdef SARUtils
    
    methods (Static) 
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
        
        function pointSrc = createPointSrc(ampl, pos)
            pointSrc.ampl = ampl; 
            pointSrc.pos = pos; 
        end
        
        function n_z = getUnitVectorPerpTo(vect_in)
            % use 3D cross product to calculate normal vector 
            % - it comes rotated counterclockwise from (zEnd - zStart), 
            % then drop the third coord

            vectCross = cross([0;0;1], [vect_in; 0]); 
            n_z = [vectCross(1); vectCross(2)]/norm(vectCross);  
        end        
        
        function scattField = getScattField(k, phi_xhat, reflector)
            [rows_z, cols_z] = size(reflector.curve.z); assert(rows_z == 2);
            xhat = [cos(phi_xhat); sin(phi_xhat)];
            k_inc = - k * xhat; 
            % sum(..., 1) works as columnwise dot(,), i.e., dot(,) at each z 
            % BTW dot here OK because args are real    
            % same as: q_inc_z = k * sum(repmat(xhat, 1, colsTheta) .* nz, 1);   
            curve = reflector.curve; 
            q_inc_z = - sum(repmat(k_inc, 1, cols_z) .* curve.nz, 1); 
            assert(all( abs(q_inc_z / k) > (1/3) ), 'incident wave too shallow'); 


            incWavePhase = sum(repmat(k_inc, 1, cols_z) .* curve.z, 1); 
            scattField.value = - k^2 * q_inc_z.^(-2) .* reflector.coeff .* exp(1i * incWavePhase); 
            scattField.normal_deriv = 1i * q_inc_z .* scattField.value; 
        end
        
        function traj = doFFP_curve(k, phi, curve, field)

            traj.phi = phi; 
            traj.ffp = nan(size(traj.phi));    

            [rows_z, cols_z] = size(curve.z); assert(rows_z == 2); 

            for iphi = 1:numel(traj.phi)
                phi = traj.phi(iphi); 

                % incident wave propagates from "-hat{x}"
                xhat_vector = [cos(phi); sin(phi)]; 
                k_xhat_z = k * sum(repmat(xhat_vector, 1, cols_z) .* curve.z, 1); 

                % outgoing wave is also a plane wave 
                outgoingExp.value = exp(-1i * k_xhat_z); 
                outgoingExp.gradient = -1i * k * [cos(phi); sin(phi)] * outgoingExp.value;
                outgoingExp.normal_deriv = sum(outgoingExp.gradient .* curve.nz, 1);

                % don't use dot(,) for complex numbers! it makes conj() of the first argument! 
                ffp = sum(field.value .* outgoingExp.normal_deriv ... 
                        - field.normal_deriv .* outgoingExp.value); 

                traj.ffp(iphi) = ffp; 
                %{ 
                DEBUG
                if iphi < 5 
                    fprintf('Start k_xhat_z: %g, diffs = ', k_xhat_z(1)); 
                    fprintf('%8.5f ', k_xhat_z(1:5) - k_xhat_z(1)); 
                    fprintf('\n Mean curve = (%8.5f, %8.5f), result: %8.5f + i %8.5f\n', ... 
                                mean(curve.z(1)), mean(curve.z(2)), real(ffp), imag(ffp)); 
                end
                %}
            end

        end
        
        function I = buildSARimage(y1, y2, k_band, phi_range, uinf)

            [rowsI, colsI] = size(y1); 
            I = zeros(rowsI, colsI); 

            for im = 1:numel(phi_range)
                phi = phi_range(im); 
                x_hat = [cos(phi); sin(phi)]; % scattered wave - don't need a minus
                for il = 1:numel(k_band)
                    dotprod = x_hat(1) * y1 + x_hat(2) * y2; 
                    I = I + uinf(im, il) * exp(2i * k_band(il) * dotprod); 
                end
            end  

        end

    end

end