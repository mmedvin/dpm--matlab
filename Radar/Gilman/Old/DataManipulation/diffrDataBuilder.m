function diffrDataBuilder(irs)

    if ~exist('irs', 'var'), irs = 10; end
    tol = 1e-9; 

    addpath('./dpm--matlab-Radar'); 
    %addpath('../matlab')

%    inMatFile = './SAR_Ellipse_a1_b2_gp8.mat'; 
%    outMatFilePrefix = 'ellipse_medvinScatt'; 

    inMatFile = './SAR_Circle_r1_gp8.mat'; 
    outMatFilePrefix = 'circle_medvinScatt'; 

    
    tic
    load(inMatFile, 'M', 'a', 'b');
    fprintf('Finished readling long file. ');
    toc

    for iir = 1:numel(irs)
        ir = irs(iir); 
        
        extractDataFor_ir(ir, M, outMatFilePrefix, tol); 
    end

    assert(a == b); 
    r = a; 
    extractDataFor_scatterer(M, outMatFilePrefix, r, tol); 
end
 
function extractDataFor_ir(ir, M, outMatFilePrefix, tol)
   
    spec = getSpec(); 
    
    outMatFile = sprintf('%s%s%d', outMatFilePrefix, '_ir_', ir);  
   
    Grid = M.PlrGrid;
    n = numel(Grid.r) - ir;  
    assert(abs(Grid.r(n)   - Grid.r(n-1) - Grid.dx) < tol); 
    assert(abs(Grid.r(n+1) - Grid.r(n)   - Grid.dx) < tol);     

    r = Grid.r(n);     
    theta = Grid.theta;     
    curve.z =  r * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 

    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
            % According to the spec, it should be: 
            % scattData = M.scattData{im, il}
            medvinScattData = M.scattData{il, im}; 
            % --- DAMN!
            
            k = medvinScattData.k;
            phi = medvinScattData.phi;    
            assert(abs(k   - spec.k_range(il))   < tol); 
            assert(abs(phi - spec.phi_range(im)) < tol);     
            
            U = medvinScattData.Extu;

            u = U(n,:);
            % SECOND ORDER
            %un = (U(n+1,:) - U(n-1,:)) / (2 * Grid.dx); 
            
            % FOURTH ORDER
            un = (8 * (U(n+1,:) - U(n-1,:)) - (U(n+2,:) - U(n-2,:)) ) / (12 * Grid.dx);             
            
            assert(numel(u)  == numel(theta)); 
            assert(numel(un) == numel(theta)); 
            
            scattFields(im,il).value = u; 
            scattFields(im,il).normal_deriv = un;            
            
        end
    end
    toc 

    save(outMatFile, 'curve', 'scattFields', 'ir', 'r')
end

function extractDataFor_scatterer(M, outMatFilePrefix, r, tol)
   
    spec = getSpec(); 
    
    outMatFile = sprintf('%s%s', outMatFilePrefix, '_iface');  
   
    theta = spec.theta;     
    curve.z =  r * [cos(theta); sin(theta)]; 
    curve.nz = [cos(theta); sin(theta)]; 

    Basis = M.Basis;    
    if ~isfield(Basis,'AddParams')
        Basis.AddParams=[];
    end
    
    for im = 1:numel(spec.phi_range)    
        for il = 1:numel(spec.k_range)
            
            % According to the spec, it should be: 
            % scattData = M.scattData{im, il}
            medvinScattData = M.scattData{il, im}; 
            % --- DAMN!
            
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

    save(outMatFile, 'curve', 'scattFields', 'r')
end



