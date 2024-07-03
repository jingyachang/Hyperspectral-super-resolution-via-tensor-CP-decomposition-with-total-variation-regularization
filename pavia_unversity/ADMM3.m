function AA=ADMM3(Y,E0, B,Z, P,A,beta,mu3,lambda1,lambda4,tol)
AA=zeros(size(E0)); %AA初始矩阵，mu1,mu2,lambda1分别是beta，罚参数，正则化参数。
Dd1=zeros(size(E0,1)-1,size(E0,2));%拉格朗日乘子L3
Dy=diffy(E0);% 差分矩阵,如文献中定义的，差分矩阵
relerr=[];% 相对误差矩阵
fprintf('#ADMM_3iter    relerr\n')
fprintf('===================================\n');
for i=1:400
    A_old=AA; %初始化
%     Gaussian = fspecial('gaussian',[3 3],0.5); %创建二位滤波器
%     AA = conv2(AA,Gaussian,'same');  %高斯卷积并返回一个与Bn相同大小的中心矩阵
%     [~,Gy]=gradient(AA);
%     AA = 1./(1+0.1*Gy).*AA;
    Vv1=soft(Dy*AA-Dd1,lambda4/2*mu3); %V3的更新
    AA=GE3(Y,E0,AA, B,Z, P,A,beta,mu3,lambda1,Vv1,Dd1,Dy);%调用函数
    %Dd1=Dd1+mu3*(Dy*AA-Vv1);
    Dd1=Dd1-(Dy*AA-Vv1);
    relerr(i)=abs(norm(AA(:)-A_old(:))/norm(A_old)+eps); % 相对误差。
    relerr(1)=1;
    fprintf('%3d \t %25.24f\n',i,relerr(i));
    if relerr(i)<tol
        %fprintf('i=%d, relerr=%1.4f \n',i, relerr);
        break
    end
end
end
