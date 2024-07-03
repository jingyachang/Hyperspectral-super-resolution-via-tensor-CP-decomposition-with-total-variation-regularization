function AA=ADMM1(Y,E0, B,Z, P,A,beta,mu1,lambda1,lambda2,tol)
AA=zeros(size(E0)); %AA��ʼ������mu1,mu2,lambda1�ֱ���beta�������������򻯲�����
Dd1=zeros(size(E0,1)-1,size(E0,2));%�������ճ���L3
Dy=diffy(E0);%��־���,�������ж���ģ���־���
relerr=[];%���������
fprintf('#ADMM1_iter    relerr\n')
fprintf('===================================\n');
for i=1:400%40
    A_old=AA; 
%     Gaussian = fspecial('gaussian',[3 3],0.5); %������λ�˲���
%     AA = conv2(AA,Gaussian,'same');  %��˹���������һ����AA��ͬ��С�����ľ���
%     [~,Gy]=gradient(AA);
%     AA = 1./(1+0.1*Gy).*AA;
    Vv1=soft(Dy*AA-Dd1,lambda2/2*mu1);%V1�ĸ���
    AA=GE(Y,E0,AA,B,Z, P,A,beta,mu1,lambda1,Vv1,Dd1,Dy);%����CG����
    %Dd1=Dd1+mu1*(Dy*AA-Vv1);
    Dd1=Dd1-(Dy*AA-Vv1);
    relerr(i)=abs(norm(AA(:)-A_old(:))/norm(A_old)+eps);%�����
    relerr(1)=1;
    fprintf('%3d \t %25.24f\n',i,relerr(i));
    if relerr(i)<tol
        %fprintf('i=%d, relerr=%1.4f \n',i, relerr);
        break
    end
end
end