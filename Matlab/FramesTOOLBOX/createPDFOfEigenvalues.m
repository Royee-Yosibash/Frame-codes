function [xHist,PDF] = createPDFOfEigenvalues(eigenVect,nBins, Manova_Gamma)
%CREATEHISTOFEIGENVALUES Summary of this function goes here
%   Detailed explanation goes here
maxHist = max(2*abs(eigenVect));
xHist = linspace(0, maxHist, nBins);
dxHist = abs(xHist(2)-xHist(1));


if ~isempty(Manova_Gamma)
   deltLocation = 1/Manova_Gamma;
   eigenVectNew = eigenVect;
   ind = eigenVectNew >= (deltLocation*(1-dxHist));
   probDelta = sum(ind)/numel(eigenVect);
   eigenVectNew(ind) = [];
   PDF = hist(eigenVectNew, xHist);
   [PDF] = normPDF(xHist,PDF); 
   PDF = PDF * (1-probDelta);
   
   dif = abs( xHist - deltLocation );
    [~,IndMin] = min(dif);
    
    PDF(IndMin) = PDF(IndMin) + probDelta; 
   
else
    PDF = hist(eigenVect, xHist);
    [PDF] = normPDF(xHist,PDF);
end
end

