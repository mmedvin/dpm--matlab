function newAdaptor(inMatFile,outMatFilePrefix)
    
    tol = 1e-9;
    
    %addpath('./dpm--matlab-Radar');
    %addpath('../matlab')

    %M.scattData = [M1.scattData M2.scattData M3.scattData M4.scattData];
    %also need to combine spec.k_range and meta
    
    if nargin==0
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
    load(inMatFile, 'M', 'a', 'b');
    fprintf('Finished readling %s.', inMatFile);
    toc
    
    ind_letters = {'a', 'b', 'c'};
    
    for ind = 1:numel(ind_letters)
        extractDataFor_ind(ind, ind_letters, M, outMatFilePrefix, tol);
    end
    
    
    if (a == b);
        r = a;
        extractDataFor_scatterer(M, outMatFilePrefix, r, tol);
    end
    
end
 
function extractDataFor_ind(ind, ind_letters, M, outMatFilePrefix, tol)
    
    outMatFile = sprintf('%s%s%s', outMatFilePrefix, '_ind_', ind_letters{ind});  
       
    r = M.scattData{1,1}.r(ind);
    spec = M.spec;
    theta = spec.theta;     
    curve.z =  r * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 

    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
            medvinScattData = M.scattData{im, il}; 
            
            k = medvinScattData.k;
            phi = medvinScattData.phi;    
%             assert(abs(k   - spec.k_range(il))   < tol); 
%             assert(abs(phi - spec.phi_range(im)) < tol);     
            
            u  = medvinScattData.u {ind}; 
            un = medvinScattData.un{ind}; 
            
            scattFields(im,il).value = u; 
            scattFields(im,il).normal_deriv = un;            
            
        end
    end
    toc 

    save(outMatFile, 'curve', 'scattFields', 'r','spec','-v7.3');
end

function extractDataFor_scatterer(M, outMatFilePrefix, r, tol)
   
    outMatFile = sprintf('%s%s', outMatFilePrefix, '_iface');  
   
    spec = M.spec;
    theta = spec.theta;%(0 : 0.0025 : 1.999);%spec.theta;     
    curve.z =  r * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 

    Basis = M.Basis;    
    if ~isfield(Basis,'AddParams')
        Basis.AddParams=[];
    end
    
    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
            medvinScattData = M.scattData{im, il};
            
            k = medvinScattData.k;
            phi = medvinScattData.phi;    
            assert(abs(k   - spec.k_range(il))   < tol); 
            assert(abs(phi - spec.phi_range(im)) < tol);     
            
            cn0 = medvinScattData.cn0;
            cn1 = medvinScattData.cn1;

            u = 0;
            un = 0;

            for j=1:numel(Basis.Indices0)
                uj = Basis.Handle(theta,Basis.Indices0(j),Basis.AddParams);
                u = u + cn0(j).*uj.xi0;
            end

            for j=1:numel(Basis.Indices1)
                uj = Basis.Handle(theta,Basis.Indices1(j),Basis.AddParams);
                un = un + cn1(j).*uj.xi0;
            end         
            
            assert(numel(u)  == numel(theta)); 
            assert(numel(un) == numel(theta)); 
            
            scattFields(im,il).value = u; 
            scattFields(im,il).normal_deriv = un;            
            
        end
    end
    toc 

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
