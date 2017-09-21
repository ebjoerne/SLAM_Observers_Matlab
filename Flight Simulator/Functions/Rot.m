function matrix=Rot(N)
% R_nb=Rot(psi)
    yaw=N;
    matrix=[cos(yaw) -sin(yaw) 0;
            sin(yaw)  cos(yaw) 0;
            0         0        1]; 
end
