close all; clear all; clc
nWorkers = 1000;
WorkerRedundancyNeeded  = 500;
mComputations = nWorkers - WorkerRedundancyNeeded;
gamma = mComputations/nWorkers;

codeType = 'RS';
% codeType = 'BPF';
Code = getCode(nWorkers, mComputations, codeType, 'None');

DFT = dftmtx(nWorkers);
IDFT = inv(DFT);
H = Code * IDFT;

normH = zeros(1,nWorkers);
for i=1:nWorkers
    normH(i) = norm(H(:,i));
end

idx = abs(H) < 1e-10;
H(idx) = 0;


H_tag = DFT * ctranspose(Code);

normH_tag = zeros(1,nWorkers);
for i=1:nWorkers
    normH_tag(i) = norm(H_tag(i,:));
end


idx = abs(H_tag) < 1e-10;
H_tag(idx) = 0;

figure();
bar(normH); 
xlabel('Norm element index');
ylabel('Power');
title(['Code type: ', codeType, '. n= ', num2str(nWorkers), ' m= ', num2str(mComputations)]);
