function AA=ADMM1(Y,E0, B,Z, P,A,beta,mu1,lambda1,lambda2,tol)
AA=zeros(size(E0)); %AA初始化矩阵，mu1,mu2,lambda1分别是beta，罚参数，正则化参数。
Dd1=zeros(size(E0,1)-1,size(E0,2));%拉格朗日乘子L3
Dy=diffy(E0);%差分矩阵,如文献中定义的，差分矩阵
relerr=[];%相对误差矩阵
fprintf('#ADMM1_iter    relerr\n')
fprintf('===================================\n');
for i=1:400%40
    A_old=AA; 
%     Gaussian = fspecial('gaussian',[3 3],0.5); %创建二位滤波器
%     AA = conv2(AA,Gaussian,'same');  %高斯卷积并返回一个与AA相同大小的中心矩阵
%     [~,Gy]=gradient(AA);
%     AA = 1./(1+0.1*Gy).*AA;
    Vv1=soft(Dy*AA-Dd1,lambda2/2*mu1);%V1的更新
    AA=GE(Y,E0,AA,B,Z, P,A,beta,mu1,lambda1,Vv1,Dd1,Dy);%调用CG函数
    %Dd1=Dd1+mu1*(Dy*AA-Vv1);
    Dd1=Dd1-(Dy*AA-Vv1);
    relerr(i)=abs(norm(AA(:)-A_old(:))/norm(A_old)+eps);%相对误差。
    relerr(1)=1;
    fprintf('%3d \t %25.24f\n',i,relerr(i));
    if relerr(i)<tol
        %fprintf('i=%d, relerr=%1.4f \n',i, relerr);
        break
    end
end
end