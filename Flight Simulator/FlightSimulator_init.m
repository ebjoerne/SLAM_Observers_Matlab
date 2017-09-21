%% Flight Simulator for an AUV with Camera pointing downwards

clc
clear all
addpath Functions


time_end=10;
h=0.05;
k_max=time_end/h;

leg=['Comp';'SLAM';'True'];
Plot=0;
noise=0;

sigma_Li=0.01*0;
sigma_nudot=0.01*0;
sigma_omega=0.001;




max_estimates=10;% (can be 100)




STORE_P_N=zeros(1,3,k_max);
STORE_P_HAT=zeros(1,3,k_max);
STORE_V_B=zeros(1,3,k_max);
STORE_PSI=zeros(1,3,k_max);
STORE_PSI_HAT=zeros(1,3,k_max);



STORE_L_HAT=zeros(max_estimates,3,k_max);
STORE_D_HAT=zeros(max_estimates,3,k_max);
STORE_B_HAT=zeros(1,3,k_max);
STORE_OMEGA_BIB=zeros(1,3,k_max);
STORE_OMEGA_M=zeros(1,3,k_max);




%starting states point in NED
p_n_0=[0, 0 , -2]';
psi_0=[0, 0, 0]';
v_b_0=[0.1, 0, 0]';
omega_0=[0 , 0, 0.1/7.5]';
R_nb_0=eye(3);
R_bc=[0 0 -1; 0 1 0;1 0 0];
compas_n=[1 0 0]';
gravity_n=[0 0 -1]';
b_true=[0.8 -0.3 -0.5]';

% Estimates Initialization
rho_hat=pi/2*0.1;
theta_hat=pi*3/4*0.1;
psi_hat=pi*2/3*0.1;
R_x=[1 0 0;0 cos(rho_hat) -sin(rho_hat); 0 sin(rho_hat) cos(rho_hat)];
R_y=[cos(theta_hat) 0 sin(theta_hat); 0 1 0; -sin(theta_hat) 0 cos(theta_hat) ];
R_z=[cos(psi_hat) -sin(psi_hat) 0;sin(psi_hat) cos(psi_hat) 0; 0 0 1];
R_nb_hat=R_z*R_y*R_x;
 R2Euler=@(Rot)([atan2(Rot(3,2),Rot(3,3));atan2(-Rot(3,1),sqrt(Rot(3,2)*Rot(3,2)+Rot(3,3)*Rot(3,3)));atan2(Rot(2,1),Rot(1,1))]);


%Creating map
N=1500;
PP=[-25+50*rand(N,1), -25+50*rand(N,1) ,zeros(N,1)];

L_HAT=zeros(N,3);
D_HAT=ones(N,1)*0.1;
b_hat=zeros(1,3)';

k=0;
p_n=p_n_0; p_n_old=p_n_0;
v_b=v_b_0;
R_nb=R_nb_0;

Index_M=[];
Index_C=[];


for i=0:h:time_end
k=k+1;
    %% Kinematics of flight
    omega_bib=omega_0+[0;0.0005*sin(i*0.1);0];
    R_nb_dot=R_nb*S(omega_bib);
    R_nb_new=intEuler(R_nb_dot,R_nb,h);
    [U,LAMBDA,V]=svd(R_nb_new);
    R_nb=U*V'; %% Keep the rotation matrix orthogonal

    Euler_m=R2Euler(R_nb);
    p_n_dot=R_nb*v_b+[0 0.002 0]'*0;
    p_n_new=intEuler(p_n_dot,p_n,h);
    v=(p_n_new-p_n)/h;
    f=(p_n_new-2*p_n+p_n_old)/h^2;
    p_n_old=p_n;
    p_n=p_n_new;
    Somega_m=R_nb'*R_nb_dot;
    
    omega_m=[Somega_m(3,2) Somega_m(1,3) Somega_m(2,1)]'+b_true+sigma_omega*randn(3,1)*noise;
    
    v_b_dot=[0;0;0];
    v_b_dot=sigma_nudot*randn(3,1);
    v_b=intEuler(v_b_dot,v_b,h);
    
    STORE_P_N(:,:,k)=p_n';
    STORE_V_B(:,:,k)=v_b';
    STORE_PSI(:,:,k)=Euler_m;
    STORE_OMEGA_BIB(:,:,k)=omega_m'-b_true';
    STORE_OMEGA_M(:,:,k)=omega_m;
    
    %% Estimation
    % Recording Camera Sensor
    [Index_M, LOS_mesurments, rho_measurments]=Cameramodel(p_n,R_nb,PP);
    Index_M
    %Estimating LOS and inverce range
    k_l=10;
    k_d=5;
    k_b=5;
    k_lib=1/length(Index_M);
    l=0;
    sigma_b=0;
    for j=Index_M'
        l=l+1;
        if (norm(L_HAT(j,:))==0 || L_HAT(j,:)*LOS_mesurments(l,:)'<0.7)
            L_HAT(j,:)=LOS_mesurments(l,:);
            D_HAT(j)=mean(D_HAT(Index_M))*0+0.1;
%             D_HAT(j)=1/rho_measurments(l);
        else
            l_m=LOS_mesurments(l,:)';
            l_hat=L_HAT(j,:)';
            d_hat=D_HAT(j);
%             d_hat=1/rho_measurments(l);
            v_b_m=v_b;
           
            sigma_li=S(l_m')*l_hat*k_l;
%             sigma_li
            if(l_m'*l_hat>0.9)
            %  sigma_b=sigma_b-S(l_hat)*S(l_m)*S(l_m)*l_hat*k_lib*0+S(l_m')*l_hat;
              sigma_b=sigma_b+S(l_m')*l_hat;
            end
             
            Li_hat_dot=-S(omega_m-b_hat+sigma_li)*l_hat+d_hat*S(l_hat)*S(l_m)*v_b_m;
            L_hat_new=intEuler(Li_hat_dot,l_hat,h);
            L_hat_new=L_hat_new/norm(L_hat_new);
            if (isnan(L_hat_new))
               disp('LOS_error')
               break
            end
            % InvRange estimate
            d_hat_dot=d_hat^2*l_m'*v_b_m + k_d/(norm(v_b_m)+0.001)*v_b_m'*S(l_m)*S(l_m)*S(l_hat)*S(l_m)*l_hat;
            if (d_hat>2)
            if(d_hat_dot>=0)
                d_hat_dot=0;
            end
            elseif(d_hat<=0)
                if(d_hat_dot<0)
                d_hat_dot=0;
                end
            end
            d_hat_new=intEuler(d_hat_dot,d_hat,h);
            L_HAT(j,:)=L_hat_new';
            D_HAT(j)=d_hat_new;        
        end    
    end
    b_hat_dot=-sigma_b*k_b;
    b_hat=intEuler(b_hat_dot,b_hat,h);
    STORE_B_HAT(:,:,k)=b_hat';
end



%%
% close all
[Index_M, LOS_mesurments, rho_measurments]=Cameramodel(p_n,R_nb,PP);
STORE_P_N_t=squeeze(STORE_P_N(:,:,1:k));
figure(1)
plot3(STORE_P_N_t(1,:)',STORE_P_N_t(2,:)',-STORE_P_N_t(3,:)')
hold on 
scatter3(PP(:,1),PP(:,2),PP(:,3),'b')
scatter3(PP(Index_M,1),PP(Index_M,2),PP(Index_M,3),'r')


scatter3(p_n(1),p_n(2),-p_n(3),'rx')
axis equal
nn=length(Index_M);

LOS_mesurments_n=R_nb*LOS_mesurments';
LOS_mesurments_n=LOS_mesurments_n';
quiver3(ones(nn,1)*p_n(1),ones(nn,1)*p_n(2),-ones(nn,1)*p_n(3),LOS_mesurments_n(:,1)*5,LOS_mesurments_n(:,2)*5,-LOS_mesurments_n(:,3)*5)

%%
%CIndex=Cameramodel(p_test,PP);
%PPC=PP(CIndex,:);
figure(2) 
plot(squeeze(STORE_B_HAT)')















