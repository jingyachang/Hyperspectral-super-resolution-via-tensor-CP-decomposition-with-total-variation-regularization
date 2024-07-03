function uiqi = uiqi(SRI,SRI_hat1)
q_band   = zeros(1, size(SRI_hat1,3));
for idx1 = 1:size(SRI_hat1,3)
    q_band(idx1)=img_qi(SRI(:,:,idx1), SRI_hat1(:,:,idx1), 10);
end
uiqi = mean(q_band);
end