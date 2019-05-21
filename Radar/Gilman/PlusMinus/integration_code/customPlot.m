function customPlot(PlrGrid, Extu)
            R  = ones(size(PlrGrid.R)).*NaN;
            Th = ones(size(PlrGrid.Theta)).*NaN;
%            Nm = ExtPrb.Scatterer.Nm;
%            R(Nm) = PlrGrid.R(Nm);
%            Th(Nm)= PlrGrid.Theta(Nm);
            
            Nth = PlrGrid.Ny; 
            R(:,Nth)=R(:,Nth-1);
            Th(:,Nth)=2*pi;
            
            XExt = (PlrGrid.R .* cos(PlrGrid.Theta));
            YExt = (PlrGrid.R .* sin(PlrGrid.Theta));
            
%            tExtu= ones(size(PlrGrid.R)).*NaN;
%            tExtu(Nm)=Extu(Nm);
%            tExtu(:,Nth)=tExtu(:,1);
 
%            pcolor([XExt],[YExt],abs([tExtu]));

            pcolor(XExt,YExt,abs(full(Extu))); 
            shading flat
end