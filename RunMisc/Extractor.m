function [u,un] = Extractor(filename,im,il,theta)
    load(filename,'M');
    
    Basis = M.Basis;
    cn0 = M.scattData{im,il}.cn0;
    cn1 = M.scattData{im,il}.cn1;
    
    if ~isfield(Basis,'AddParams')
        Basis.AddParams=[];
    end
    
    u=0;
    un=0;
    
    for j=1:numel(Basis.Indices0)
        uj = Basis.Handle(theta,Basis.Indices0(j),Basis.AddParams);
        u = u + cn0(j).*uj.xi0;
    end

    for j=1:numel(Basis.Indices1)
        uj = Basis.Handle(theta,Basis.Indices1(j),Basis.AddParams);
        un = un + cn1(j).*uj.xi0;
    end

end