function [A,B,C] = UTV_SR1(HSI,MSI,P1,P2,Pm,par)
beta=par.beta;
mu1=par.mu1;
mu2=par.mu2;
mu3=par.mu3;
lambda1=par.lambda1;
lambda2=par.lambda2;
lambda3=par.lambda3;
lambda4=par.lambda4;
%%
[Ih,Jh,K]=size(HSI);
[I,J,Km]=size(MSI);
HSI1=reshape(HSI,[Ih,Jh*K]);
HSI2=reshape(permute(HSI,[2 1 3]),[Jh,K*Ih]);%permute置换数组维度
H3=reshape(HSI,[Ih*Jh,K]);
HSI3=H3';
MSI1=reshape(MSI,[I,J*Km]);
MSI2=reshape(permute(MSI,[2 1 3]),[J,Km*I]);
MSI3=reshape(MSI,[I*J,Km])';
% %%  simulate LR-HSI
% %HSI;
% HSI1=Unfold(HSI,size(HSI),1);
% HSI2=Unfold(HSI,size(HSI),2);
% HSI3=Unfold(HSI,size(HSI),3);%与H3之间差别是一个转置的关系。
% [Ih,Jh,K]=size(HSI);
% H3=reshape(HSI,[Ih*Jh,K]);%420*204
% %%  simulate HR-MSI
% %MSI;
% MSI1=Unfold(MSI,size(MSI),1);
% MSI2=Unfold(MSI,size(MSI),2);
% MSI3=Unfold(MSI,size(MSI),3);
%% TenRec (initialization algorithm)
t_rank=400;
%rank of the tensor -- change for different noise levels
%使用第二个数据集的时候t_rank取60,70，R-SNR的值保持不变，且取60的时候是最大的。
%取50会变化，取80会报错，取90以及更大会变化
maxit=25;
%[A,B,C,~,~,~]=TenRec(MSI,H3,maxit,t_rank,P1,P2);%调用TenRec用于CG初始值。
load('A');
load('B');
load('C');
%% 设置第一个数据集的初始值
% load('A1')
% A_hat=A1;
% load('B1')
% B_hat=B1;
% load('C1')
% C_hat=C1;
%% 设置第二个数据集的初始值
% load('A11')
% A=A11;
% load('B11')
% B=B11;
% load('C11')
% C=C11;
%% 设置第三个数据集的初始值
% load('A111')
% A_hat=A111;
% load('B111')
% B_hat=B111;
% load('C111')
% C_hat=C111;
%% 设置第四个数据集的初始值
% load('A1111')
% A_hat=A1111;
% load('B1111')
% B_hat=B1111;
% load('C1111')
% C_hat=C1111;
%%
tol=1e-2;
%MAXIT=10;
%cost(1)=inf;
%for i=1:5
%%  update A
X1=khatri_rao(Pm*C,B);%X1的大小为504*50
X2=khatri_rao(C,P2*B);%X2的大小为4284*50
A=ADMM1(MSI1, A, X1, HSI1, P1, X2, beta, mu1, lambda1, lambda2, tol);
%%  update B1
X3=khatri_rao(Pm*C,A);
X4=khatri_rao(C,P1*A);
B=ADMM2(MSI2, B, X3, HSI2, P2, X4, beta, mu2, lambda1, lambda3, tol);
%%  update C1
X5=khatri_rao(P2*B,P1*A);
X6=khatri_rao(B,A);
C=ADMM3(HSI3, C, X5, MSI3, Pm, X6, beta, mu3, lambda1, lambda4, tol);
%% 终止条件
%     H_ferror=norm(HSI3'-X5*C','fro');
%     M_ferror=norm(MSI3'-X6*C'*Pm','fro');
%     cost(iter)=H_ferror+M_ferror;
%     if  abs((cost(iter)-cost(iter-1))/cost(iter-1))<tol
%         %fprintf('iter=%d, cost(iter)=%1.4f \n',iter,cost(iter));
%         break
%     end 
end