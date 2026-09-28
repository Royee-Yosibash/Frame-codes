function [hFig] = plotFramePDF(eigenVect, nBins, beta, gamma, titleStr)
% beta = ((nWorkers - iDeleted)/mComputations)^-1
% gamma = mComputations/nWorkers;

[lambda,PDF] = createPDFOfEigenvalues(eigenVect,nBins);

hFig = figure('units','normalized','outerposition',[0 0 1 1]);   hold on;
plot(lambda, PDF, 'b*', 'DisplayName','Gram matrix \lambda hist'); grid on;
ylabel('Probabilty'); xlabel('\lambda (eigenvalues of A = transpose(H)*H)');

meanEig = sum(lambda.*PDF);
plot(meanEig, 0, 'bx', 'markersize', 20, 'DisplayName','Mean of Gram matrix \lambda');
plot(beta, 0, 'g+', 'markersize', 20, 'DisplayName','\beta');

[f] = marchenkoPasturPDF(lambda, gamma); % This is originaly beta, but due to gamma being defined as m/n this is the one
% meanF = sum(lambda.*f*beta);
% plot(lambda*beta, f, 'r', 'DisplayName','marchenko pastur dist');
meanF = sum(lambda.*f);
plot(lambda, f, 'r', 'DisplayName','marchenko pastur dist');
plot(meanF, 0, 'rx', 'markersize', 20, 'DisplayName','Mean of marchenko pastur dist');
ymaxPlot = find(PDF>1e-4);  ymaxPlot = ymaxPlot(end);
xlim([0, lambda(ymaxPlot)]);
title(titleStr);    legend();

%% Display the results
disp(['beta = ', num2str(beta), '. gamma = ', num2str(gamma)]);
disp(['The mean of the eigenvalues is: ', num2str(meanEig)]);
disp(['The mean of the MP distribution is: ', num2str(meanF)]);

end

