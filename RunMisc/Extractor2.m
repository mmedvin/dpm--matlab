function [u,un,k,phi,r,theta] = Extractor2(filename,im,il,rind)

    tic
    load(filename,'M');
    toc
    
    k=M.scattData{im,il}.k;
    phi=M.scattData{im,il}.phi;
    
    Grid = M.PlrGrid;
    n=numel(Grid.r)-rind;
    r=Grid.r(n);
    theta=Grid.theta;

    U = M.scattData{im,il}.Extu;
    u = U(n,:);
    
    un = (U(n+1,:) - U(n-1,:)) / (2 * Grid.dx);
    toc
end

%{

scattField = 

  struct with fields:

           value: [1×113 double]
    normal_deriv: [1×113 double]

curve

curve = 

  struct with fields:

     z: [2×113 double]
    nz: [2×113 double]

%}