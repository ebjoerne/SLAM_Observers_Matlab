function [Index1, LOS_mesurments rho_measurments]=Cameramodel(Position,R_nb,PP)
    
    Xmin=Position(1,1)-(-Position(3,1));
    Xmax=Position(1,1)+(-Position(3,1));
    Ymin=Position(2,1)-(-Position(3,1));
    Ymax=Position(2,1)+(-Position(3,1));
    n=length(PP);

    IX=find(PP(:,1)>Xmin*ones(n,1) & PP(:,1)<Xmax*ones(n,1));
    IY=find(PP(:,2)>Ymin*ones(n,1) & PP(:,2)<Ymax*ones(n,1));
    
%     Index1 = intersect(IX,IY);
    flag = ismember(IX,IY);
    Index1 = IX(find(flag));
  
    n=length(Index1);
    LOS_mesurments=zeros(n,3);
    rho_measurments=zeros(n,1);
    for i=1:n;
         LOS_n=(PP(Index1(i),:)'-Position)/norm(PP(Index1(i),:)-Position');
         LOS_b=R_nb'*LOS_n;
         LOS_mesurments(i,:)=LOS_b';
         rho_measurments(i,1)=norm(PP(Index1(i),:)-Position');
    end

end