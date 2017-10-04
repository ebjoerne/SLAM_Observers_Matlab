%% Flight Simulator for an AUV with Camera pointing downwards

clc
clear all
addpath Functions

time_end=2.5;
h=0.0125;
k_max=time_end/h;
floorAndSealing=1;
comp_filter=0;

leg=['Comp';'SLAM';'True'];
Plot=0;
video_rec=0;
noise=0;
Just_yaw=0;
Just_yaw_cel=0;

sigma_Li=0.01*0;
sigma_nudot=0.01*0;
sigma_omega=0.001;
sigma_compas=0.1;
sigma_Icompass=0.01;

compInoise=zeros(3,1);



max_estimates=10;% (can be 100)




STORE_P_N=zeros(1,3,k_max);
STORE_P_HAT=zeros(1,3,k_max);
STORE_V_B=zeros(1,3,k_max);
STORE_V_N=zeros(1,3,k_max);
STORE_PSI=zeros(1,3,k_max);
STORE_PSI_HAT=zeros(1,3,k_max);



STORE_L_HAT=zeros(max_estimates,3,k_max);
STORE_D_HAT=zeros(max_estimates,3,k_max);
STORE_B_HAT=zeros(1,3,k_max);
STORE_OMEGA_BIB=zeros(1,3,k_max);
STORE_OMEGA_M=zeros(1,3,k_max);

STORE_A_N=zeros(1,3,k_max);
STORE_VEL=zeros(1,3,k_max);
STORE_LAMBDA=zeros(1,1,k_max);



%starting states point in NED
p_n_0=[0, 0 , -2]';
psi_0=[0, 0, 0]';
v_b_0=[0.5, 0, -10]';
omega_0=[0 , 0, 0.1/7.5]'*0;
R_nb_0=eye(3);
R_bc=[0 0 -1; 0 1 0;1 0 0];
compas_n=[1 0 0]';
g_n=[0 0 9.81]';
b_true=[0.8 0.1 -0.5]';

% Estimates Initialization
rho_hat=pi/2*0.1;
theta_hat=pi*3/4*0.1;
psi_hat=pi*2/3*0.1;
R_x=[1 0 0;0 cos(rho_hat) -sin(rho_hat); 0 sin(rho_hat) cos(rho_hat)];
R_y=[cos(theta_hat) 0 sin(theta_hat); 0 1 0; -sin(theta_hat) 0 cos(theta_hat) ];
R_z=[cos(psi_hat) -sin(psi_hat) 0;sin(psi_hat) cos(psi_hat) 0; 0 0 1];
R_nb_hat=R_z*R_y*R_x;
 R2Euler=@(Rot)([atan2(Rot(3,2),Rot(3,3));atan2(-Rot(3,1),sqrt(Rot(3,2)*Rot(3,2)+Rot(3,3)*Rot(3,3)));atan2(Rot(2,1),Rot(1,1))]);
 S=@(nu)([0 -nu(3) nu(2);   nu(3) 0 -nu(1);   -nu(2) nu(1) 0]);

%Creating map
N=1500;
PP=[-25+50*rand(N,1), -25+50*rand(N,1) ,-4*(rand(N,1) < 0.5)*floorAndSealing];

Convergance_iterations=50; % How many iterations before we are confident on convergance

L_HAT=zeros(N,3);
D_HAT=ones(N,1)*0.1;
b_hat=zeros(1,3)';
Converge_Counter=Convergance_iterations*ones(N,1);
sigma_conv=0;
sigma_b_comp=0;
if (video_rec==1)
    framerat_hz=5;
    K_mov=1/h*framerat_hz;
    Movie(k_max/K_mov) = struct('cdata',[],'colormap',[]);
    mov_counter=0;
end

x_vel_hat=zeros(7,1); % X_vel_hat=[a^n, v^n, lambda_vel]
P_vel_0=eye(7,7)*2;
P_vel=P_vel_0;
a_n_hat=[0 0 0]';

R_vel=eye(6,6)*0.5;
Q_vel=eye(7,7)*0.5;

k=0;
p_n=p_n_0; p_n_old=p_n_0;
v_b=v_b_0;
R_nb=R_nb_0;



Index_M=[];
Index_C=[];


for i=0:h:time_end
k=k+1;
    %% Kinematics of flight
    omega_bib=omega_0*0+[0;0.0005*sin(i*0.1);0]*0;
%     v_b=v_b_0+[-cos(i*0.001)*0.05, -cos(i*0.01)*0.02, 0]';

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
    g_m=R_nb'*(g_n+f*0-a_n_hat*0)/norm(g_n+f*0-a_n_hat*0);
    f_m=R_nb'*(g_n+f);
    compInoise=compInoise*0.99+sigma_Icompass*randn(3,1);
    compas_m=R_nb'*compas_n+compInoise*0+sigma_compas*0;
    compas_m=compas_m/norm(compas_m);
    
    v_b_dot=[0;0;9.81];
%     v_b_dot=sigma_nudot*randn(3,1);
    v_b=intEuler(v_b_dot,v_b,h);
    
    STORE_P_N(:,:,k)=p_n';
    STORE_V_B(:,:,k)=v_b';
    STORE_V_N(:,:,k)=v';
    STORE_PSI(:,:,k)=Euler_m;
    STORE_OMEGA_BIB(:,:,k)=omega_m'-b_true';
    STORE_OMEGA_M(:,:,k)=omega_m;
    
    %% Estimation
    % Recording Camera Sensor
    Index_M_old=Index_M;
    [Index_M, LOS_mesurments, rho_measurments]=Cameramodel(p_n,R_nb,PP);
    
    flag = ~ismember(Index_M_old,Index_M);
    index = find(flag);
    if (~isempty(index)&& k>10)
        Index_out=Index_M_old(index); % Indecies of features not seen anymore
        for m=1:length(Index_out)
            r=Index_out(m);
            Converge_Counter(r)=Convergance_iterations;
            D_HAT(r)=0.1;
            L_HAT(r,:)=[0 0 0];
        end
    end
    
    %Estimating LOS and inverce range
    k_l=15;
    k_d=10;
    k_b=10;
    k_lib=1/length(Index_M);
    l=0;
    sigma_b=0;
    
    if(k>30000)
%         g_m
%         g_hat
        sigma_g=S(g_m)*g_hat*2;
        b_hat
        sigma_b=sigma_b+S(g_m)*g_hat*1;
        g_hat_dot=-S(omega_m-b_hat+sigma_g)*g_hat;
        g_hat_new=intEuler(g_hat_dot,g_hat,h);
        g_hat=g_hat_new/norm(g_hat_new);
        if (isnan(g_hat))
           disp('LOS_error')
           break
        end        
        
    else
        g_hat=g_m;
    end
    
    
    for j=Index_M'
        l=l+1;
        if (norm(L_HAT(j,:))==0 || L_HAT(j,:)*LOS_mesurments(l,:)'<0.97*0)
            L_HAT(j,:)=LOS_mesurments(l,:);
            D_HAT(j)=mean(D_HAT(Index_M))*0+0.1;
             D_HAT(j)=1/rho_measurments(l);
        else
            l_m=LOS_mesurments(l,:)';
            l_hat=L_HAT(j,:)';
            d_hat=D_HAT(j);
%             d_hat=1/rho_measurments(l);
            v_b_m=v_b;
           
            sigma_li=S(l_m')*l_hat*k_l;
%             sigma_li
            % Deciding if LOS/Range estimates are to be used in the bias
            % estimation
            if((sigma_conv==0 || Converge_Counter(j)<=0))
            %  sigma_b=sigma_b-S(l_hat)*S(l_m)*S(l_m)*l_hat*k_lib*0+S(l_m')*l_hat;
              sigma_b=sigma_b+S(S(g_m)*l_m)*S(g_m)*l_hat*Just_yaw_cel+(1-Just_yaw_cel)*S(l_m)*l_hat;
              
              if (sigma_conv==0)
                    Converge_Counter(j)=0;
              end
            else 
                Converge_Counter(j)=Converge_Counter(j)-1;
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
    
    if (norm(sigma_b)<10^-2&& k>10)
        sigma_conv=1;
        if (Just_yaw==1)
            Just_yaw_cel=1;
        end
    end
    b_hat_dot=-sigma_b*k_b*(1-comp_filter)-sigma_b_comp*k_b*comp_filter;
    b_hat=intEuler(b_hat_dot,b_hat,h);
    STORE_B_HAT(:,:,k)=b_hat';
    
    
    %% Attitude observer
   k_1=0.5;  v_n_1=g_n/norm(g_n);               v_b_1=g_m;
   k_2=0.5*0.3;   v_n_2=compas_n;          v_b_2=compas_m;
   k_3=0.3*0.1*0;   v_n_3=S(v_n_1)*v_n_2;    v_b_3=S(v_b_1)*v_b_2;

    
    sigma_r=k_1*S(v_b_1)*R_nb_hat'*v_n_1+k_2*S(v_b_2)*R_nb_hat'*v_n_2+k_3*S(v_b_3)*R_nb_hat'*v_n_3;   %%% set to zero for c
    sigma_r;
    R_nb_dot=R_nb_hat*S(omega_m-b_hat+sigma_r);
    R_nb_hat_new=intEuler(R_nb_dot,R_nb_hat,h);
    if (k==500)
        R_nb_hat_new=R_nb;
    end
    
   sigma_b_comp=k_1*S(v_b_1)*R_nb_hat'*v_n_1+k_2*S(v_b_2)*R_nb_hat'*v_n_2+k_3*S(v_b_3)*R_nb_hat'*v_n_3;
    %% Algorithm for orthogonalizing the rotation matrix
    r_1_bar=R_nb_hat_new(:,1)/norm(R_nb_hat_new(:,1));
    r_2_bar=(eye(3)-r_1_bar*r_1_bar')*R_nb_hat_new(:,2);
    r_2_bar=r_2_bar/norm(r_2_bar);
    r_3_bar=S(r_1_bar)*r_2_bar;
    
    R_nb_hat=[r_1_bar r_2_bar r_3_bar];
    [U,LAMBDA,V]=svd(R_nb_hat_new);
    R_nb_hat=U*V'; %%%%%%%%%% keeping orthogonal
    STORE_PSI_HAT(:,:,k)=R2Euler(R_nb_hat);
    STORE_PSI(:,:,k)=R2Euler(R_nb);
    
    
    %% Velocity scaling estimator
 
   if(k>3500) 
    if (norm(v)<10^-5)
        v_n_bar=[0 0 0]';
    else
            v_n_bar=v/norm(v);
    end
%         a_n_m=R_nb*g_m-
    if(k>3) 
        a_n_m=f;
    else
        a_n_m=zeros(3,1);
    end
    
    % Estimation update
    A_vel=zeros(7,7); A_vel(4:6,1:3)=eye(3); A_vel(7,1:3)=v_n_bar';
    A_vel_d=eye(7,7)+h*A_vel;
    C_vel=zeros(6,7); C_vel(1:3,1:3)=eye(3); C_vel(4:6,4:6)=eye(3); C_vel(4:6,7)=-v_n_bar;
    C_vel_d=C_vel;
    % B_vel=0;
    y_m=[R_nb_hat*f_m-g_n;zeros(3,1)];
    y_hat=C_vel*x_vel_hat;
    y_tilde=y_m-y_hat;

    
    S0 = C_vel*P_vel*C_vel' + R_vel;
    K_vel = P_vel*C_vel'*(S0\eye(6,6));
    dx=h*K_vel*y_tilde;
    
    x_vel_hat=x_vel_hat+dx;
    P_vel=(eye(7,7)-K_vel*C_vel)*P_vel;    
    
    STORE_A_N(:,:,k)=x_vel_hat(1:3)';
    STORE_VEL(:,:,k)=x_vel_hat(4:6)';
    STORE_LAMBDA(:,:,k)=x_vel_hat(7)';
    
    a_n_hat=x_vel_hat(1:3);
    
    % Propagation
    x_vel_hat=A_vel_d*x_vel_hat;
    P_vel=A_vel_d*P_vel*A_vel_d'+Q_vel;
    
   end  
    
    
    %% Video Storage
    if (video_rec==1)
    if (mod(k,K_mov)==0)
        mov_counter=mov_counter+1;
        [Index_M, LOS_mesurments, rho_measurments]=Cameramodel(p_n,R_nb,PP);
        STORE_P_N_t=squeeze(STORE_P_N(:,:,1:k));
        figure(1)
        plot3(STORE_P_N_t(1,:)',STORE_P_N_t(2,:)',-STORE_P_N_t(3,:)')
        hold on 
        scatter3(PP(:,1),PP(:,2),-PP(:,3),'y.')
        scatter3(PP(Index_M,1),PP(Index_M,2),-PP(Index_M,3),'r')
        scatter3(p_n(1),p_n(2),-p_n(3),'rx')
        axis equal
        nn=length(Index_M);

        LOS_mesurments_n=R_nb*LOS_mesurments';
        LOS_mesurments_n=LOS_mesurments_n';
        
%         D_HAT(find(D_HAT<0.05))=10;
         D_hat_temp=D_HAT(Index_M);
          D_hat_temp(find(D_hat_temp<0.1))=10;
        
         quiver3(ones(nn,1)*p_n(1),ones(nn,1)*p_n(2),-ones(nn,1)*p_n(3),LOS_mesurments_n(:,1).*1./D_hat_temp,LOS_mesurments_n(:,2).*1./D_hat_temp,-LOS_mesurments_n(:,3).*1./D_hat_temp)
%         quiver3(ones(nn,1)*p_n(1),ones(nn,1)*p_n(2),-ones(nn,1)*p_n(3),LOS_mesurments_n(:,1).*rho_measurments,LOS_mesurments_n(:,2).*rho_measurments,-LOS_mesurments_n(:,3).*rho_measurments)

        drawnow
        Movie(mov_counter)=getframe;
        hold off 
    end
    end
end



%%
 close all
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

figure(3)
plot(squeeze(STORE_PSI)')
hold on
plot(squeeze(STORE_PSI_HAT)','--')

figure(4)

%%
plot(squeeze(STORE_VEL)','--')
hold on
plot(squeeze(STORE_V_N)','-r')
% movie(Movie,2)











