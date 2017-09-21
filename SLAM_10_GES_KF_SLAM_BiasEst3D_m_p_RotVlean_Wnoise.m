%% SLAM demonstrator 
% clc
clear all
% close all
addpath Functions

u_d=1;

x_0=0;
y_0=0;
z_0=0;

p_n_0=[x_0;
       y_0;
       z_0;];
psi_0=0;
R_0=Rot(psi_0);


u_0=u_d;
v_0=0;
w_0=u_d*0.1;

nu_0=[u_0;
      v_0;
      w_0];
  
r_0=u_d/75;

omega_0=[0;
         0;
         r_0];
     
 leg=['Comp';'SLAM';'True'];

% PP=[[20,10,0]';
%    [30,20,0]';
%    [-20,0,0]'];
% m_p=length(PP)/3;
p_n=p_n_0;
m_p=4;

PP=zeros(m_p*3,1);

for i=1:m_p
PP(3*i-2:3*i,1)=[rand*50 rand*50 rand*50]'*10-[25 25 25]'*10;
end


time_end=1000;
h=0.1; 

Delta=PP-kron(ones(m_p,1),p_n);

sigma_Li=zeros(3*m_p,1);
sigma_nudot=0.01*0;
sigma_omega=0.001;
sigma_mag=0.065;
sigma_acc=0.1;
sigma_rho=0.005;
sigma_L=0.01;



L_hat=zeros(3*m_p,1);
Q_hat=ones(m_p,1);
b_hat=zeros(3,1);
L_hat_dot_old=zeros(3*m_p,1);




for l=1:m_p
Q_hat(l)=sqrt(Delta(3*l-2:3*l,1)'*Delta(3*l-2:3*l,1));
L_hat(3*l-2:3*l,1)=Delta(3*l-2:3*l,1)/Q_hat(l);
%  L_hat(3*l-2:3*l,1)=[1 0 0];

end

 Q_hat=100*ones(m_p,1); %%%%%%%%%%%%%%%%%
Li_mes_old=L_hat;

Isigma_rho=zeros(3,1);

R_nb=R_0;
p_n=p_n_0;
p_n_old=p_n_0;
nu=nu_0;
omega_bib=omega_0;

% bias=[1.5 2 1]';
bias=[0.8 0.1 -0.5]';


 

%%
% creating system matrices
A_s=[0 zeros(1,m_p)];

for l=1:m_p
   A_s=[A_s;1 zeros(1,m_p)];
end
A_c=kron(A_s,eye(3));



B_s=[1 zeros(1,m_p)]';
B_c=kron(B_s,eye(3));

A_d=expm(A_c*h);

MM=expm([A_c, B_c;
        zeros(3,3*(m_p+2))]*h);

A_d_N=MM(1:3*(m_p+1),1:3*(m_p+1));
B_d=MM(1:3*(m_p+1),3*(m_p+1)+1:3*(m_p+2));
    
C_s=[zeros(m_p,1), eye(m_p)];
C_c=kron(C_s,R_nb');
C_d=C_c;


%% States
X=[p_n;PP];
X_hat=zeros(length(X),1);
X_hat_p=zeros(length(X),1);
% Uncertaintees
% P_hat=diag([0.5 0.5 0.5 100 100 100 100 100 100 100 100 100]);
% P_p=diag([0.5 0.5 0.5 100 100 100 100 100 100 100 100 100]);
% Q=diag([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]);
% R=diag([1 1 0.10 1 1 0.10 1 1 0.10]);
% R=eye(9)*0.1;

P_p=eye(3*(m_p+1))*100;
P_p(1:3,1:3)=diag([0.5 0.5 0.5]);
P_hat=P_p;
Q=eye(3*(m_p+1))*0.1;
R=kron(eye(m_p),diag([0.10 0.10 0.10]));


Delta_Hat=zeros(3*m_p,1);
time=0:h:time_end;
  P_N=zeros(3,length(time));
  ATT=zeros(3,length(time));
  DELTA=zeros(3*m_p,length(time));
  DELTA_HAT=zeros(3*m_p,length(time));
  VV=zeros(3,length(time));
  FF=zeros(3,length(time));
 
  LL_M=zeros(3*m_p,length(time));
  QQ_M=zeros(m_p,length(time));
  VV_M=zeros(3,length(time));
  
  L_bM=zeros(3*m_p,length(time));
  L_HAT=zeros(3*m_p,length(time));
  L_DOT_HAT=zeros(3*m_p,1);
  L_DOT_M=zeros(3*m_p,1);
  
  B_HAT=zeros(3,length(time));
  RHO_HAT=zeros(m_p,length(time));
  
  EULER_M=zeros(3,length(time));
  EULER_HAT=zeros(3,length(time));
  EULER_GRIP_HAT=zeros(3,length(time));

 
  
  B_GRIP_HAT=zeros(3,length(time));

  OMEGA_BIB=zeros(3,length(time));
  OMEGA_M=zeros(3,length(time));
  
  XXHAT=zeros((m_p+1)*3,length(time));
  k=1;
  
  Rot_hat=eye(3);
  Li_hat_dot=zeros(3,1);
  L_dot_hat=zeros(3*m_p,1);
  L_dot_m=zeros(3*m_p,1);

    compas_n=[1 0 0]';
  gravity_m=[0 0 -1]';
  
  %%
  
  rho=5*pi/6;
  theta=0;
  psi=pi*3/4;
  R_x=[1 0 0;0 cos(rho) -sin(rho); 0 sin(rho) cos(rho)];
R_y=[cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta) ];
R_z=[cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0; 0 0 1];

RR=R_z*R_y*R_x;
  %%
  Rot_hat=RR;
  
  Rot_grip_hat=eye(3);
   Rot_grip_hat=RR;
  b_grip_hat=zeros(3,1);
  R2Euler=@(Rot)([atan2(Rot(3,2),Rot(3,3));atan2(-Rot(3,1),sqrt(Rot(3,2)*Rot(3,2)+Rot(3,3)*Rot(3,3)));atan2(Rot(2,1),Rot(1,1))]);
  
for i=0:h:time_end
    
       omega_bib=[r_0*0.1*cos(i*0.005);
         r_0*0.3*sin(i*0.01);
         r_0*log(1+i*0.01)];
    
    R_nb_dot=R_nb*S(omega_bib);
    
    R_nb_new=intEuler(R_nb_dot,R_nb,h);
    [U,LAMBDA,V]=svd(R_nb_new);
    R_nb_new=U*V'; %% Keep the rotation matrix orthogonal
    R_nb=R_nb_new; %%%%%%%%%

    Euler_m=R2Euler(R_nb);
    p_n_dot=R_nb*nu+[0 0.002 0]'*0;
    p_n_new=intEuler(p_n_dot,p_n,h);
    v=(p_n_new-p_n)/h;
    f=(p_n_new-2*p_n+p_n_old)/h^2;
    Somega_m=R_nb'*R_nb_dot;
    
    omega_m=[Somega_m(3,2) Somega_m(1,3) Somega_m(2,1)]'+bias+sigma_omega*randn(3,1);
    
    p_n_old=p_n;
    p_n=p_n_new;
        nu_dot=[0;0;0];
    nu_dot=sigma_nudot*randn(3,1);
    nu_new=intEuler(nu_dot,nu,h);
    nu=nu_new;
    Delta=PP-kron(ones(m_p,1),p_n);
    
    for r=1:m_p
    Q_m(r)=sqrt(Delta((r-1)*3+1:r*3)'*Delta((r-1)*3+1:r*3))+sigma_rho*randn; % varrho
    L_m((r-1)*3+1:r*3,1)=Delta((r-1)*3+1:r*3,1)/Q_m(r)+S(sigma_L*randn(3,1))*Delta((r-1)*3+1:r*3,1)/Q_m(r);      % LOS
    v_b_m=nu;
    end
    LL_M(:,k)=L_m;  %% her kan støy introduseres!!! 
    QQ_M(:,k)=Q_m;  %% her kan støy introduseres!!! 
    VV_M(:,k)=v_b_m;  %% her kan støy introduseres!!! 
    EULER_M(:,k)=Euler_m;
    
    

        
     %vectors
     v_n_1=compas_n;
%      v_n_2=v/norm(v);
     v_n_2=gravity_m;
     v_n_3=S(gravity_m)*compas_n;
%      v_n_3=S(v_n_2)*compas_n;
     
          v_b_1=R_nb'*compas_n+sigma_mag*randn(3,1);
%      v_b_2=nu/norm(nu);
     v_b_2=R_nb'*(gravity_m+f/9.81)+sigma_acc*randn(3,1);
     v_b_3=S(R_nb'*(gravity_m+f/9.81))*R_nb'*compas_n;
%      v_b_3=S(R_nb'*v_n_2)*R_nb'*compas_n;
    
    
    %% Bias Estimator 
    sigma_b=0;
    k_l=1;
    k_b=1;
    k_lib=1/m_p;
%      b_hat=[0 0 1]'; % just for testing
    
    for n=1:m_p
        Li_hat_old=L_hat(3*n-2:3*n,1);
        Li_mes=R_nb'*L_m(3*n-2:3*n,1);
        rho_mes=Q_m(n);
        
        
        k_li=k_l*1;
        sigma_Li=S(Li_mes)*Li_hat_old*k_li;
        sigma_b=sigma_b+S(Li_mes)*Li_hat_old*k_lib;
         sigma_b=sigma_b-S(Li_hat_old)*S(Li_mes)*S(Li_mes)*Li_hat_old*k_lib;

                
        Li_hat_dot=-S(omega_m-b_hat+sigma_Li)*Li_hat_old+1/rho_mes*(Li_mes*Li_mes'-eye(3))*v_b_m;
        Li_hat_dot=-S(omega_m-b_hat+sigma_Li)*Li_hat_old+1/rho_mes*S(Li_hat_old)*S(Li_mes)*v_b_m;
%         Li_hat_dot=-S(omega_m)*Li_hat_old-S(Li_mes)*(b_hat-sigma_Li)+  1/rho_mes*S(Li_mes)*S(Li_mes)*v_b_m; %% LTV system

        L_hat_new=intEuler(Li_hat_dot,Li_hat_old,h);
        L_hat_new=L_hat_new/norm(L_hat_new);
        
        L_hat(3*n-2:3*n,1)=L_hat_new;
        Delta_Hat(3*n-2:3*n,1)=R_nb*L_hat_new*rho_mes;
        Li_mes_dot=(Li_mes-Li_mes_old(3*n-2:3*n,1))/h;
        L_dot_hat(n*3-2:n*3,1)=Li_hat_dot;
        L_dot_m(n*3-2:n*3,1)=Li_mes_dot;
        Li_mes_old(3*n-2:3*n,1)=Li_mes;
    end
    
    %% Rotation Estimate
    
    
    
    b_hat_dot=-k_b*sigma_b;
    b_hat=intEuler(b_hat_dot,b_hat,h);
    

    k_1=0.20;    
    k_2=0.5;
    k_3=0.3;
     

     
    
    %Rot estimation
    sigma_rot=k_1*S(v_b_1)*Rot_hat'*v_n_1+k_2*S(v_b_2)*Rot_hat'*v_n_2+k_3*S(v_b_3)*Rot_hat'*v_n_3;
    Rot_hat_dot=Rot_hat*S(omega_m-b_hat+sigma_rot);
    Rot_hat=intEuler(Rot_hat_dot,Rot_hat,h);
    [U,LAMBDA,V]=svd(Rot_hat);
    Rot_hat=U*V'; %%%%%%%%%% keeping orthogonal
    
   % storing bias callculation
     L_bM(:,k)=Li_mes_old;
     L_HAT(:,k)=L_hat;
     B_HAT(:,k)=b_hat;
     DELTA_HAT(:,k)=Delta_Hat;
%      RHO_HAT(:,k)=Q_hat;
     OMEGA_BIB(:,k)=omega_bib;
     OMEGA_M(:,k)=omega_m;
     L_DOT_HAT(:,k)=L_dot_hat;
     L_DOT_M(:,k)=L_dot_m;
     EULER_HAT(:,k)=R2Euler(Rot_hat);

   %% Mahoney/grip estimator
     k_g_i=0.15;

    
     
     sigma_grip=k_1*S(v_b_1)*Rot_grip_hat'*v_n_1+k_2*S(v_b_2)*Rot_grip_hat'*v_n_2+k_3*S(v_b_3)*Rot_grip_hat'*v_n_3;
     
      Rot_hat_grip_dot=Rot_grip_hat*S(omega_m-b_grip_hat+sigma_grip);
      b_grip_hat_dot=-k_g_i*sigma_grip;
      Rot_grip_hat=intEuler(Rot_hat_grip_dot,Rot_grip_hat,h);
      [U,LAMBDA,V]=svd(Rot_grip_hat);
      Rot_grip_hat=U*V'; %%%%%%%%%% keeping orthogonal
      
      b_grip_hat=intEuler(b_grip_hat_dot,b_grip_hat,h);
      
      Rot_grip_tilde=R_nb*Rot_grip_hat';
      
      EULER_GRIP_HAT(:,k)=R2Euler(Rot_grip_hat);
      B_GRIP_HAT(:,k)=b_grip_hat;


     

    %% Bearing and range measurment
      C_d=kron(C_s,R_nb');

        X=[v;Delta];
        y=C_d*X;
        
              
    %% Kalman Updte
    %State estimation update
    K_k=P_p*C_d'*((C_d*P_p*C_d'+R)\eye(length(R)));
    X_hat_k=X_hat_p+K_k*(y-C_d*X_hat_p); 
    P_hat=(eye(length(K_k))-K_k*C_d)*P_p*(eye(length(K_k))-K_k*C_d)'+K_k*R*K_k';
    %State propagation
   
    X_hat_p=A_d*X_hat_k+B_d*f;
    P_p=A_d*P_hat*A_d'+Q;
    
     
    %%
    
    
    
    
    
    
    
    
%     X_hat_new=A_d*X_hat+B;
    

    
  P_N(:,k)=p_n;  
  ATT(:,k)=R_nb*[1;0;0];
  DELTA(:,k)=Delta;
  FF(:,k)=f;
  VV(:,k)=v;
  XXHAT(:,k)=X_hat_k;
  k=k+1;
end
%%
figure
% plot(P_N(2,:),P_N(1,:))
plot3(P_N(3,:),P_N(2,:),P_N(1,:))
hold on
% quiver(P_N(2,1:10:end),P_N(1,1:10:end),ATT(2,1:10:end)*100,ATT(1,1:10:end)*100)
quiver3(P_N(3,1:10:end),P_N(2,1:10:end),P_N(1,1:10:end),ATT(3,1:10:end)*100,ATT(2,1:10:end)*100,ATT(1,1:10:end)*100)

hold on 
% scatter(PP(2:3:end,1),PP(1:3:end,1))
scatter3(PP(3:3:end,1),PP(2:3:end,1),PP(1:3:end,1))
axis equal


% DELTA_1=DELTA(1:3,:);
% DELTA_2=DELTA(4:6,:);
% DELTA_3=DELTA(7:9,:);

DELTA_1_hat=DELTA_HAT(1:3,:);
% DELTA_2_hat=DELTA_HAT(4:6,:);
% DELTA_3_hat=DELTA_HAT(7:9,:);


for k=1:m_p
eval(['DELTA_' num2str(k) '_hat=DELTA_HAT(' num2str(k*3-2) ':' num2str(k*3) ',:);'])
end

% for k=1:m_p
% eval(['DELTA_' num2str(k) '_hat=DELTA(' num2str(k*3-2) ':' num2str(k*3) ',:);'])
% end

%quiver(P_N(2,1:5000:end),P_N(1,1:5000:end),DELTA_3(2,1:5000:end),DELTA_3(1,1:5000:end),0) %% Paints arrows to target point
res_DH=800;
% res_DH=10;

ploted=1:res_DH:length(P_N);
% ploted=floor(length(P_N)*4/5):res_DH:length(P_N);
% ploted=1:res_DH:floor(length(P_N)*1/5);
% ploted=floor(length(P_N)*1/5):res_DH:floor(length(P_N)*2/5);

% quiver(P_N(2,1:res_DH:end),P_N(1,1:res_DH:end),DELTA_1_hat(2,1:res_DH:end),DELTA_1_hat(1,1:res_DH:end),0) %% Paints arrows to target point
% quiver(P_N(2,1:res_DH:end),P_N(1,1:res_DH:end),DELTA_2_hat(2,1:res_DH:end),DELTA_2_hat(1,1:res_DH:end),0) %% Paints arrows to target point
% quiver(P_N(2,1:res_DH:end),P_N(1,1:res_DH:end),DELTA_3_hat(2,1:res_DH:end),DELTA_3_hat(1,1:res_DH:end),0) %% Paints arrows to target point

% quiver3(P_N(3,1:res_DH:end),P_N(2,1:res_DH:end),P_N(1,1:res_DH:end),DELTA_1_hat(3,1:res_DH:end),DELTA_1_hat(2,1:res_DH:end),DELTA_1_hat(1,1:res_DH:end),0) %% Paints arrows to target point
% quiver3(P_N(3,1:res_DH:end),P_N(2,1:res_DH:end),P_N(1,1:res_DH:end),DELTA_2_hat(3,1:res_DH:end),DELTA_2_hat(2,1:res_DH:end),DELTA_2_hat(1,1:res_DH:end),0) %% Paints arrows to target point
% quiver3(P_N(3,1:res_DH:end),P_N(2,1:res_DH:end),P_N(1,1:res_DH:end),DELTA_3_hat(3,1:res_DH:end),DELTA_3_hat(2,1:res_DH:end),DELTA_3_hat(1,1:res_DH:end),0) %% Paints arrows to target point
%%
  
for k=1:m_p
eval(['quiver3(P_N(3,ploted),P_N(2,ploted),P_N(1,ploted),DELTA_' num2str(k) '_hat(3,ploted),DELTA_' num2str(k) '_hat(2,ploted),DELTA_' num2str(k) '_hat(1,ploted),0) %% Paints arrows to target point'])
end
%%

% for j=1:m_p
% figure(j+1)
% eval(['plot(time(:,1:1000:end),DELTA_' int2str(j) '(:,1:1000:end))'])
% end

% figure(6)
% plot(time(:,1:1000:end),VV(:,1:1000:end))
% 
% figure(7)
% plot(time(:,1:1000:end),FF(:,1:1000:end))

%%
sigmaBias=sqrt(var(B_HAT(1:3,1000*0.1/h:end)'))'
sigmaGripBias=sqrt(var(B_GRIP_HAT(1:3,1000*0.1/h:end)'))'

figure 
subplot(3,1,1)
plot(time(1,1000*0.1/h:end),B_HAT(1,1000*0.1/h:end)','b')
hold on
plot(time(1,1000*0.1/h:end),B_GRIP_HAT(1,1000*0.1/h:end)','g')
title('B_hat_Elias')
xlabel('[s]')

subplot(3,1,2)
plot(time(1,1000*0.1/h:end),B_HAT(2,1000*0.1/h:end)','b')
hold on
plot(time(1,1000*0.1/h:end),B_GRIP_HAT(2,1000*0.1/h:end)','g')
xlabel('[s]')

subplot(3,1,3)
plot(time(1,1000*0.1/h:end),B_HAT(3,1000*0.1/h:end)','b')
hold on
plot(time(1,1000*0.1/h:end),B_GRIP_HAT(3,1000*0.1/h:end)','g')
xlabel('[s]')

%%
figure 
subplot(3,1,1)
plot(time(1,1:end+0*1000*0.1/h),B_GRIP_HAT(1,1:end+0*1000*0.1/h)','g')
hold on
plot(time(1,1:end+0*1000*0.1/h),B_HAT(1,1:end+0*1000*0.1/h)','b')
title('Bias Estimation')
plot(time(1,1:end+0*1000*0.1/h),bias(1)*ones(1,length(B_GRIP_HAT(3,1:end+0*1000*0.1/h)))','r')
legend(leg)
xlabel('[s]')

subplot(3,1,2)
plot(time(1,1:end+0*1000*0.1/h),B_GRIP_HAT(2,1:end+0*1000*0.1/h)','g')
hold on
plot(time(1,1:end+0*1000*0.1/h),B_HAT(2,1:end+0*1000*0.1/h)','b')
plot(time(1,1:end+0*1000*0.1/h),bias(2)*ones(1,length(B_GRIP_HAT(3,1:end+0*1000*0.1/h)))','r')
legend(leg)
xlabel('[s]')

subplot(3,1,3)
plot(time(1,1:end+0*1000*0.1/h),B_GRIP_HAT(3,1:end+0*1000*0.1/h)','g')
hold on
plot(time(1,1:end+0*1000*0.1/h),B_HAT(3,1:end+0*1000*0.1/h)','b')
plot(time(1,1:end+0*1000*0.1/h),bias(3)*ones(1,length(B_GRIP_HAT(3,1:end+0*1000*0.1/h)))','r')
legend(leg)
xlabel('[s]')

Rot_grip_tilde
%%


figure 
subplot(3,1,1)
plot(time(1,2000*0.1/h:10:end),OMEGA_M(1,2000*0.1/h:10:end)'-B_GRIP_HAT(1,2000*0.1/h:10:end)','g')
hold on
plot(time(1,2000*0.1/h:10:end),OMEGA_M(1,2000*0.1/h:10:end)'-B_HAT(1,2000*0.1/h:10:end)','b')
plot(time(1,2000*0.1/h:10:end),OMEGA_BIB(1,2000*0.1/h:10:end)','r-.')
title('\omega estimates end')
legend(leg)
ylabel('[rad/s]')
xlabel('[s]')

subplot(3,1,2)
plot(time(1,2000*0.1/h:10:end),OMEGA_M(2,2000*0.1/h:10:end)'-B_GRIP_HAT(2,2000*0.1/h:10:end)','g')
hold on
plot(time(1,2000*0.1/h:10:end),OMEGA_M(2,2000*0.1/h:10:end)'-B_HAT(2,2000*0.1/h:10:end)','b')
plot(time(1,2000*0.1/h:10:end),OMEGA_BIB(2,2000*0.1/h:10:end)','r-.')
legend(leg)
xlabel('[s]')
ylabel('[rad/s]')

subplot(3,1,3)
plot(time(1,2000*0.1/h:10:end),OMEGA_M(3,2000*0.1/h:10:end)'-B_GRIP_HAT(3,2000*0.1/h:10:end)','g')
hold on
plot(time(1,2000*0.1/h:10:end),OMEGA_M(3,2000*0.1/h:10:end)'-B_HAT(3,2000*0.1/h:10:end)','b')
plot(time(1,2000*0.1/h:10:end),OMEGA_BIB(3,2000*0.1/h:10:end)','r-.')
legend(leg)
xlabel('[s]')
ylabel('[rad/s]')
%%
figure 
subplot(3,1,1)
plot(time(1,1:10:2000),OMEGA_M(1,1:10:2000)'-B_GRIP_HAT(1,1:10:2000)','g')
hold on
plot(time(1,1:10:2000),OMEGA_M(1,1:10:2000)'-B_HAT(1,1:10:2000)','b')
plot(time(1,1:10:2000),OMEGA_BIB(1,1:10:2000)','r-.')
title('\omega estimates start')
legend(leg,'Location','EastOutside')
xlabel('[s]')
ylabel('[rad/s]')

subplot(3,1,2)
plot(time(1,1:10:2000),OMEGA_M(2,1:10:2000)'-B_GRIP_HAT(2,1:10:2000)','g')
hold on
plot(time(1,1:10:2000),OMEGA_M(2,1:10:2000)'-B_HAT(2,1:10:2000)','b')
plot(time(1,1:10:2000),OMEGA_BIB(2,1:10:2000)','r-.')
legend(leg,'Location','EastOutside')
xlabel('[s]')
ylabel('[rad/s]')

subplot(3,1,3)
plot(time(1,1:10:2000),OMEGA_M(3,1:10:2000)'-B_GRIP_HAT(3,1:10:2000)','g')
hold on
plot(time(1,1:10:2000),OMEGA_M(3,1:10:2000)'-B_HAT(3,1:10:2000)','b')
plot(time(1,1:10:2000),OMEGA_BIB(3,1:10:2000)','r-.')
legend(leg,'Location','EastOutside')
xlabel('[s]')
ylabel('[rad/s]')
%   EULER_M=zeros(3,length(time));
%   EULER_HAT=zeros(3,length(time));
%   EULER_GRIP_HAT=zeros(3,length(time));

subplot(3,1,1)
handaxes2 = axes('Position', [0.55 0.805 0.2 0.17]);

    hold on
    %plot(time(1000:n), E(i,1000:n), colors(i,:))
    plot(time(1:200), OMEGA_M(1,1:200)-B_GRIP_HAT(1,1:200), 'g')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), OMEGA_M(1,1:200)-B_HAT(1,1:200), 'b')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), OMEGA_BIB(1,1:200), 'r-.')
    set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')



subplot(3,1,2)
handaxes2 = axes('Position', [0.55 0.505 0.2 0.17]);

    hold on
    %plot(time(1000:n), E(i,1000:n), colors(i,:))
    plot(time(1:200), OMEGA_M(2,1:200)-B_GRIP_HAT(2,1:200), 'g')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), OMEGA_M(2,1:200)-B_HAT(2,1:200), 'b')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), OMEGA_BIB(2,1:200), 'r-.')
    set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')


subplot(3,1,3)
handaxes2 = axes('Position', [0.55 0.205 0.2 0.17]);

    hold on
    %plot(time(1000:n), E(i,1000:n), colors(i,:))
    plot(time(1:200), OMEGA_M(3,1:200)-B_GRIP_HAT(3,1:200), 'g')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), OMEGA_M(3,1:200)-B_HAT(3,1:200), 'b')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), OMEGA_BIB(3,1:200), 'r-.')
    set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')

%%
figure
subplot(3,1,1)
plot(time(1,1:10:end),EULER_GRIP_HAT(1,1:10:end)'*180/pi,'g')
xlabel('[s]')
ylabel('[deg]')

hold on
plot(time(1,1:10:end),EULER_HAT(1,1:10:end)'*180/pi,'b')

plot(time(1,1:10:end),EULER_M(1,1:10:end)'*180/pi,'r-.')
title('Euler angles')
legend(leg,'Location','EastOutside')

% handaxes2 = axes('Position', [0.7 0.655 0.25 0.1]);
% 
%     hold on
%     %plot(time(1000:n), E(i,1000:n), colors(i,:))
%     plot(time(1:250), EULER_GRIP_HAT(1,1:250), 'g')
%     set(handaxes2, 'Box', 'off')
%     plot(time(1:250), EULER_HAT(1,1:250), 'b')
%     set(handaxes2, 'Box', 'off')
%     plot(time(1:250), EULER_M(1,1:250), 'r-.')
%     set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')



subplot(3,1,2)
plot(time(1,1:10:end),EULER_GRIP_HAT(2,1:10:end)'*180/pi,'g')
ylabel('[deg]')

hold on
plot(time(1,1:10:end),EULER_HAT(2,1:10:end)'*180/pi,'b')
plot(time(1,1:10:end),EULER_M(2,1:10:end)'*180/pi,'r-.')
legend(leg,'Location','EastOutside')


subplot(3,1,3)
plot(time(1,1:10:end),EULER_GRIP_HAT(3,1:10:end)'*180/pi,'g')
hold on
plot(time(1,1:10:end),EULER_HAT(3,1:10:end)'*180/pi,'b')
plot(time(1,1:10:end),EULER_M(3,1:10:end)'*180/pi,'r-.')
ylabel('[deg]')

legend(leg,'Location','EastOutside')
%%
subplot(3,1,1)
handaxes2 = axes('Position', [0.55 0.805 0.2 0.17]);

    hold on
    %plot(time(1000:n), E(i,1000:n), colors(i,:))
    plot(time(1:200), EULER_GRIP_HAT(1,1:200)*180/pi, 'g')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), EULER_HAT(1,1:200)*180/pi, 'b')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), EULER_M(1,1:200)*180/pi, 'r-.')
    set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')



subplot(3,1,2)
handaxes2 = axes('Position', [0.55 0.505 0.2 0.17]);

    hold on
    %plot(time(1000:n), E(i,1000:n), colors(i,:))
    plot(time(1:200), EULER_GRIP_HAT(2,1:200)*180/pi, 'g')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), EULER_HAT(2,1:200)*180/pi, 'b')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), EULER_M(2,1:200)*180/pi, 'r-.')
    set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')


subplot(3,1,3)
handaxes2 = axes('Position', [0.55 0.205 0.2 0.17]);

    hold on
    %plot(time(1000:n), E(i,1000:n), colors(i,:))
    plot(time(1:200), EULER_GRIP_HAT(3,1:200)*180/pi, 'g')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), EULER_HAT(3,1:200)*180/pi, 'b')
    set(handaxes2, 'Box', 'off')
    plot(time(1:200), EULER_M(3,1:200)*180/pi, 'r-.')
    set(handaxes2, 'Box', 'off')
%     xlabel('[s]')
%     ylabel('[m]')
%%

