clear clc
addpath(genpath(cd));
%% Indian_pines ’LANDSAT' RANK=100 add noise 25dB 
SRI = cell2mat(struct2cell(load('Indian_pines.mat')));
SRI(:,:,[104:108 150:163 220]) = []; %Regions of water absorption
SRI(1,:,:) = []; SRI(:,1,:) = [];
Pm = spectral_deg(SRI,"LANDSAT");
MSI = tmprod(SRI,Pm,3);
d1 = 4; d2 = 4; q = 9;
[P1,P2] = spatial_deg(SRI, q, d1, d2);
HSI = tmprod(tmprod(SRI,P1,1),P2,2);
[I, J, K] = size(SRI);
S1=reshape(SRI,[I,K*J])'; % mode 1 matricization 
for k=1:size(HSI,3) %加噪声
    HSI(:,:,k) = awgn(HSI(:,:,k),25,'measured'); 
end
for k=1:size(MSI,3)
    MSI(:,:,k) = awgn(MSI(:,:,k),25,'measured');
end
%% The proposed Super-resolution method
HSI(HSI<0) = 0;
MSI(MSI<0) = 0;
%%  参数
par.beta=1e-4;%-7
par.lambda1 = 1; par.lambda2 = 1e1;%1
par.lambda3 = 1e-3; par.lambda4 = 1e4;  
par.mu1=1e-2; par.mu2=1e-2; par.mu3=1e-2;
t=clock;
[A,B,C]= UTV_SR1(HSI, MSI, P1, P2, Pm, par);%主函数
t1=etime(clock,t);
%%
SRI_hat1 = cpdgen({A,B,C});
d1=4; d2=4;
err1 = cell2mat(compute_metrics(SRI,SRI_hat1,d1,d2));
tab = ["Algorithm" "R-SNR" "CC" "SAM" "ERGAS" "ssim" "uiqi" "times" ;...
     "Indian pines" err1 t1];
tab

figure
subplot(2,2,1); imagesc(SRI(:,:,31)); title('Groundtruth SRI'); colorbar; lim = caxis; axis off
subplot(2,2,2); imagesc(SRI_hat1(:,:,31)); title('Ours'); colorbar; caxis(lim); axis off
subplot(2,2,3); imagesc(SRI_hat1(:,:,31)); title('Ours'); colorbar; caxis(lim); axis off
subplot(2,2,4); imagesc(SRI_hat1(:,:,31)); title('Ours'); colorbar; caxis(lim); axis off
saveas(gcf,'fig1','png');