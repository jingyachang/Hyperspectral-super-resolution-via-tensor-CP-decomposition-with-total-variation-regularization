function E =GE3(Y,E0,E,X1,Z, P,X2,beta,mu,lambda1,V111,D111,Dy )
maxIter =40;
PPP=P'*P;
A1=X2'*X2;
B1=X1'*X1+beta*eye(size(X1'*X1,1));%B*B'=X3*X3'
H=Y*X1+beta*E0+lambda1*P'*Z*X2+mu*Dy'*(V111+D111);%P'*Z*A'=D3'*C(3)*Y3'��A=������Y3��B=������X3��Y=������B(3).
%C(3)��ʾ�����ͼ��,B(3)��ʾ�߹���ͼ��
r0=H-E*B1-lambda1*PPP*E*A1-mu*(Dy)'*Dy*E;%�������ұ߼�ȥ��ߡ�
p0=r0;
for i=1:maxIter
    pp= (p0*B1+lambda1*PPP*p0*A1+mu*(Dy)'*Dy*p0);%
    pp1=p0(:)'*pp(:);
   
    a=(r0(:)')*r0(:)/pp1;
    E=E+a*p0;
    r1=r0-a*pp;
    
    b1=(r1(:)'*r1(:))/(r0(:)'*r0(:));
    p1=r1+b1*p0;
    p0=p1;
    r0=r1;
end


