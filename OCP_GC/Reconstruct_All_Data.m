function alliterdata=Reconstruct_All_Data(w_opt,Options,NMuscle,NParameters,N,d,nq,scaling,jointi,tau_root)


    NParameters = 1;    
    alliterdata.tf_opt = w_opt(1:NParameters);
    starti = NParameters+1;
    alliterdata.a_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    alliterdata.a_col_opt = reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    alliterdata.FTtilde_opt = reshape(w_opt(starti:starti+NMuscle*(N+1)-1),NMuscle,N+1)';
    starti = starti + NMuscle*(N+1);
    alliterdata.FTtilde_col_opt =reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    alliterdata.Qs_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    alliterdata.Qs_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);
    alliterdata.Qdots_opt = reshape(w_opt(starti:starti+nq.all*(N+1)-1),nq.all,N+1)';
    starti = starti + nq.all*(N+1);
    alliterdata.Qdots_col_opt = reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,d*N)';
    starti = starti + nq.all*(d*N);    
    alliterdata.a_a_opt = reshape(w_opt(starti:starti+nq.arms*(N+1)-1),nq.arms,N+1)';
    starti = starti + nq.arms*(N+1);
    alliterdata.a_a_col_opt = reshape(w_opt(starti:starti+nq.arms*(d*N)-1),nq.arms,d*N)';
    starti = starti + nq.arms*(d*N);
    alliterdata.vA_opt = reshape(w_opt(starti:starti+NMuscle*N-1),NMuscle,N)';
    starti = starti + NMuscle*N;
    alliterdata.e_a_opt = reshape(w_opt(starti:starti+nq.arms*N-1),nq.arms,N)';
    starti = starti + nq.arms*N;   
    alliterdata.dFTtilde_col_opt=reshape(w_opt(starti:starti+NMuscle*(d*N)-1),NMuscle,d*N)';
    starti = starti + NMuscle*(d*N);
    alliterdata.qdotdot_col_opt =reshape(w_opt(starti:starti+nq.all*(d*N)-1),nq.all,(d*N))';
    starti = starti + nq.all*(d*N);
    if starti - 1 ~= length(w_opt)
        error('error when extracting results')
    end
    % Combine results at mesh and collocation points
    alliterdata.a_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    alliterdata.a_mesh_col_opt(1:(d+1):end,:)= alliterdata.a_opt;
    alliterdata.FTtilde_mesh_col_opt=zeros(N*(d+1)+1,NMuscle);
    alliterdata.FTtilde_mesh_col_opt(1:(d+1):end,:)= alliterdata.FTtilde_opt;
    alliterdata.Qs_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    alliterdata.Qs_mesh_col_opt(1:(d+1):end,:)= alliterdata.Qs_opt;
    alliterdata.Qdots_mesh_col_opt=zeros(N*(d+1)+1,nq.all);
    alliterdata.Qdots_mesh_col_opt(1:(d+1):end,:)= alliterdata.Qdots_opt;
    alliterdata.a_a_mesh_col_opt=zeros(N*(d+1)+1,nq.arms);
    alliterdata.a_a_mesh_col_opt(1:(d+1):end,:)= alliterdata.a_a_opt;
    for k=1:N
        rangei = k*(d+1)-(d-1):k*(d+1);
        rangebi = (k-1)*d+1:k*d;
        alliterdata.a_mesh_col_opt(rangei,:) = alliterdata.a_col_opt(rangebi,:);
        alliterdata.FTtilde_mesh_col_opt(rangei,:) = alliterdata.FTtilde_col_opt(rangebi,:);
        alliterdata.Qs_mesh_col_opt(rangei,:) = alliterdata.Qs_col_opt(rangebi,:);
        alliterdata.Qdots_mesh_col_opt(rangei,:) = alliterdata.Qdots_col_opt(rangebi,:);
        alliterdata.a_a_mesh_col_opt(rangei,:) = alliterdata.a_a_col_opt(rangebi,:);
    end

    %% Unscale results
    % States at mesh points
    % Qs (1:N-1)
    alliterdata.q_opt_unsc.rad = alliterdata.Qs_opt(1:end-1,:).*repmat(...
        scaling.Qs,size(alliterdata.Qs_opt(1:end-1,:),1),1); 
    alliterdata.q_opt_unsc.rad(:,jointi.knee_ty.r) = alliterdata.Qs_opt(1:end-1,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
    % Convert in degrees
    alliterdata.q_opt_unsc.deg = alliterdata.q_opt_unsc.rad;
    alliterdata.q_opt_unsc.deg(:,[1:3,7:16 20:end]) = alliterdata.q_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qs (1:N)
    alliterdata.q_opt_unsc_all.rad = alliterdata.Qs_opt(1:end,:).*repmat(...
        scaling.Qs,size(alliterdata.Qs_opt(1:end,:),1),1); 
    alliterdata.q_opt_unsc_all.rad(:,jointi.knee_ty.r)=alliterdata.Qs_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
    % Convert in degrees
    alliterdata.q_opt_unsc_all.deg = alliterdata.q_opt_unsc_all.rad;
    alliterdata.q_opt_unsc_all.deg(:,[1:3,7:16 20:end]) = ...
        alliterdata.q_opt_unsc_all.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qdots (1:N-1)
    alliterdata.qdot_opt_unsc.rad = alliterdata.Qdots_opt(1:end,:).*repmat(...
        scaling.Qdots,size(alliterdata.Qdots_opt(1:end,:),1),1);
    % Convert in degrees
    alliterdata.qdot_opt_unsc.deg = alliterdata.qdot_opt_unsc.rad;
    alliterdata.qdot_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
    alliterdata.qdot_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qdots (1:N)
    alliterdata.qdot_opt_unsc_all.rad =alliterdata.Qdots_opt.*repmat(scaling.Qdots,size(alliterdata.Qdots_opt,1),1); 
    % Muscle activations (1:N-1)
    alliterdata.a_opt_unsc = alliterdata.a_opt(1:end-1,:).*repmat(...
        scaling.a,size(alliterdata.a_opt(1:end-1,:),1),size(alliterdata.a_opt,2));
    alliterdata.a_opt_unsc_all = alliterdata.a_opt(1:end,:).*repmat(...
        scaling.a,size(alliterdata.a_opt(1:end,:),1),size(alliterdata.a_opt,2));
    % Muscle-tendon forces (1:N-1)
    alliterdata.FTtilde_opt_unsc = alliterdata.FTtilde_opt(1:end-1,:).*repmat(...
        scaling.FTtilde,size(alliterdata.FTtilde_opt(1:end-1,:),1),1);
    % Muscle-tendon forces
    alliterdata.FTtilde_opt_unsc_all = alliterdata.FTtilde_opt(1:end,:).*repmat(...
        scaling.FTtilde,size(alliterdata.FTtilde_opt(1:end,:),1),1);
    % Arm activations (1:N-1)
    alliterdata.a_a_opt_unsc = alliterdata.a_a_opt(1:end-1,:).*repmat(...
        scaling.a_a,size(alliterdata.a_a_opt(1:end-1,:),1),size(alliterdata.a_a_opt,2));
    % Arm activations (1:N)
    alliterdata.a_a_opt_unsc_all = alliterdata.a_a_opt(1:end,:).*repmat(...
        scaling.a_a,size(alliterdata.a_a_opt(1:end,:),1),size(alliterdata.a_a_opt,2));
    
    % Controls at mesh points
    % Time derivative of muscle activations (states)
    alliterdata.vA_opt_unsc = alliterdata.vA_opt.*repmat(scaling.vA,size(alliterdata.vA_opt,1),size(alliterdata.vA_opt,2));
    tact = 0.015;
    tdeact = 0.06;
    % Get muscle excitations from time derivative of muscle activations
    alliterdata.e_opt_unsc = computeExcitationRaasch(alliterdata.a_opt_unsc,alliterdata.vA_opt_unsc,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    % Arm excitations
    alliterdata.e_a_opt_unsc = alliterdata.e_a_opt.*repmat(scaling.e_a,size(alliterdata.e_a_opt,1),...
        size(alliterdata.e_a_opt,2));
    % State and Controls at collocation points
    % Qs a
    alliterdata.q_col_opt=alliterdata.Qs_col_opt;
    alliterdata.qdot_col_opt=alliterdata.Qdots_col_opt;
    alliterdata.q_col_opt_unsc.rad = alliterdata.q_col_opt(1:end,:).*repmat(...
    scaling.Qs,size(alliterdata.q_col_opt(1:end,:),1),1); 
    alliterdata.q_col_opt_unsc.rad(:,jointi.knee_ty.r) = alliterdata.q_col_opt(1:end,jointi.knee_ty.r)*scaling.knee_ty.b+scaling.knee_ty.a;
    % Convert in degrees
    alliterdata.q_col_opt_unsc.deg = alliterdata.q_col_opt_unsc.rad;
    alliterdata.q_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = alliterdata.q_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Qdots
    alliterdata.qdot_col_opt_unsc.rad = alliterdata.qdot_col_opt(1:end,:).*repmat(...
        scaling.Qdots,size(alliterdata.qdot_col_opt(1:end,:),1),1);
    % Convert in degrees
    alliterdata.qdot_col_opt_unsc.deg = alliterdata.qdot_col_opt_unsc.rad;
    alliterdata.qdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
        alliterdata.qdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Muscle activations
    alliterdata.a_col_opt_unsc = alliterdata.a_col_opt.*repmat(...
        scaling.a,size(alliterdata.a_col_opt,1),size(alliterdata.a_col_opt,2));
    %Muscle-tendon forces
    alliterdata.FTtilde_col_opt_unsc= alliterdata.FTtilde_col_opt.*repmat(...
        scaling.FTtilde,size(alliterdata.FTtilde_col_opt,1),1);
    % Arm activations
    alliterdata.a_a_col_opt_unsc=alliterdata.a_a_col_opt.*repmat(...
        scaling.a_a,size(alliterdata.a_a_col_opt,1),size(alliterdata.a_a_col_opt,2));
    % "Slack" controls at collocation points   
    % Time derivative of Qdots
    alliterdata.qdotdot_col_opt_unsc.rad = ...
        alliterdata.qdotdot_col_opt.*repmat(scaling.Qdotdots,size(alliterdata.qdotdot_col_opt,1),1);
    % Convert in degrees
    alliterdata.qdotdot_col_opt_unsc.deg = alliterdata.qdotdot_col_opt_unsc.rad;
    alliterdata.qdotdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]) = ...
        alliterdata.qdotdot_col_opt_unsc.deg(:,[1:3,7:16 20:end]).*180/pi;
    % Time derivative of muscle-tendon forces
    alliterdata.dFTtilde_col_opt_unsc = alliterdata.dFTtilde_col_opt.*repmat(...
        scaling.dFTtilde,size(alliterdata.dFTtilde_col_opt,1),size(alliterdata.dFTtilde_col_opt,2));

    %% Time grid    
    % Mesh points
    tgrid = linspace(0,alliterdata.tf_opt,N+1);
    dtime = zeros(1,d+1);
    for i=1:4
        dtime(i)=tau_root(i)*((alliterdata.tf_opt-0)/N);
    end
    % Mesh points and collocation points
    alliterdata.tgrid_ext = zeros(1,(d+1)*N+1);
    for i=1:N
        alliterdata.tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
    end
    alliterdata.tgrid_ext(end)=alliterdata.tf_opt(end); 


end