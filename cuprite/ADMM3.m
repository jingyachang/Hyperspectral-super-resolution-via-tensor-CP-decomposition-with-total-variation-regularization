function AA=ADMM3(Y,E0, B,Z, P,A,beta,mu3,lambda1,lambda4,tol)
AA=zeros(size(E0)); %AA��ʼ����mu1,mu2,lambda1�ֱ���beta�������������򻯲�����
Dd1=zeros(size(E0,1)-1,size(E0,2));%�������ճ���L3
Dy=diffy(E0);% ��־���,�������ж���ģ���־���
relerr=[];% ���������
fprintf('#ADMM_3iter    relerr\n')
fprintf('===================================\n');
for i=1:400
    A_old=AA; %��ʼ��
%     Gaussian = fspecial('gaussian',[3 3],0.5); %������λ�˲���
%     AA = conv2(AA,Gaussian,'same');  %��˹���������һ����Bn��ͬ��С�����ľ���
%     [~,Gy]=gradient(AA);
%     AA = 1./(1+0.1*Gy).*AA;
    Vv1=soft(Dy*AA-Dd1,lambda4/2*mu3); %V3�ĸ���
    AA=GE3(Y,E0,AA, B,Z, P,A,beta,mu3,lambda1,Vv1,Dd1,Dy);%���ú���
    %Dd1=Dd1+mu3*(Dy*AA-Vv1);
    Dd1=Dd1-(Dy*AA-Vv1);
    relerr(i)=abs(norm(AA(:)-A_old(:))/norm(A_old)+eps); % �����
    relerr(1)=1;
    fprintf('%3d \t %25.24f\n',i,relerr(i));
    if relerr(i)<tol
        %fprintf('i=%d, relerr=%1.4f \n',i, relerr);
        break
    end
end
end
