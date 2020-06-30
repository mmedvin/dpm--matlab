function newAdaptor(inMatFile,outMatFilePrefix)
    
    tol = 1e-9;
    
    %addpath('./dpm--matlab-Radar');
    %addpath('../matlab')

    %M.scattData = [M1.scattData M2.scattData M3.scattData M4.scattData];
    %also need to combine spec.k_range and meta
    
    if 0 %nargin==0
        inMatFile = '../SAR_Circle_r1_gp10.mat';
        outMatFilePrefix = 'circle_gp10';
        spec = getSpec();
        
        %inMatFile = './SAR_k200_Ellipse_a1_b2_gp11.mat';
        %outMatFilePrefix = 'ellipse_gp11';
        %spec = getSpec_2K();
        
        %     inMatFile = './SAR_k200_Circle_r1_gp11.mat';
        %     outMatFilePrefix = 'circle_gp11';
        %     spec = getSpec_2K();
        
    end
    
    tic
    %load(inMatFile, 'M', 'a', 'b');
    E = Extractor(inMatFile);
    fprintf('Finished readling %s.', inMatFile);
    toc
    
    ind_letters = {'a', 'b', 'c'};
    
    for ind = 1:numel(ind_letters)
        extractDataFor_ind(ind, ind_letters, E, outMatFilePrefix, tol);
    end
    
    
    if (E.a == E.b) 
        r = E.a;
        theta = linspace(0,2*pi,200);
        extractDataFor_scatterer(E, outMatFilePrefix, r, theta, tol);
    end
    
end
 
function extractDataFor_ind(ind, ind_letters, E, outMatFilePrefix, tol)
    
    outMatFile = sprintf('%s%s%s.mat', outMatFilePrefix, '_ind_', ind_letters{ind});  
       
    spec = E.spec;
    

    tmp  = zeros(1,numel(E.theta));
    scattFields(numel(spec.phi_range),numel(spec.k_range)) = struct('value',tmp,'normal_deriv',tmp);
    
    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
         [scattFields(im,il),k,phi] =  E.FieldOnCircles(im,il,ind);
            
         assert(abs(k   - spec.k_range(il))   < tol);
         assert(abs(phi - spec.phi_range(im)) < tol);
                       
        end
    end
    toc 
    
    curve = E.Circle(ind);
    r = E.r(ind);

    save(outMatFile, 'curve', 'scattFields', 'r','spec','-v7.3');
end

function extractDataFor_scatterer(E, outMatFilePrefix, r,theta, tol)
   
    outMatFile = sprintf('%s%s.mat', outMatFilePrefix, '_iface');  
   
    spec = E.spec;
    
    tmp  = zeros(1,numel(theta));
    scattFields(numel(spec.phi_range),numel(spec.k_range)) = struct('value',tmp,'normal_deriv',tmp);  
    
    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
           [scattFields(im,il),k,phi] = E.FieldOnScatterer(im,il,theta);
            
           assert(abs(k   - spec.k_range(il))   < tol);
           assert(abs(phi - spec.phi_range(im)) < tol);
           
        end
    end
    toc 

    curve = E.Curve(theta);
    
    save(outMatFile, 'curve', 'scattFields', 'r','-v7.3')
end

function spec = getSpec_2K()

    spec.k_range = 200:5:220; 
    spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
    spec.theta = pi * (0 : 0.0025 : 1.999); 

end

function spec = getSpec()

    spec.k_range = 50:0.5:55; 
    spec.phi_range = pi * (-0.2 : 0.004 : 0.2); 
    spec.theta = pi * (0 : 0.0025 : 1.999); 

end
