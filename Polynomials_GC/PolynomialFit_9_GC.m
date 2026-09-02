% This function computes the polynomials to approximate muscle-tendon 
% lengths, velocities and moment arms. 
%
% Authors: Original code from Wouter Aerts, adapted by Antoine Falisse
% Date: 12/19/2018
%
function [muscle_spanning_joint_INFO,MuscleInfo] =  PolynomialFit_9_GC(MuscleData)
  
%% Construct the polynomials for the moment arms and muscle length

    %muscle_sel=[];    
   %muscle_sel=34:50; %all right muscles
    muscle_sel=1:50; %all right muscles
%     for m_nr = 1:length(MuscleData.muscle_names)
%         if strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_r') || ...
%                 strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_l') % was _l before
%             muscle_sel = [muscle_sel m_nr];
%         end
%     end

    muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO<=0.0001 & ...
        muscle_spanning_joint_INFO>=-0.0001) = 0;
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO~=0) = 1;

    q_all = MuscleData.q;

    max_order = 9; % maximal polynomial order
    threshold_r = 0.003; % maximal deviation from reference OpenSim data (in m) for rotational degrees of freedom
    threshold_t = 0.03; % maximal deviation from reference OpenSim data (in []) for translational degrees of freedom
    nr_samples = length(q_all(:,1));

    lMT_all_error = zeros(length(muscle_sel), 1);
    DM_all_error = zeros(length(muscle_sel), length(q_all(1,:)));
    order_all = zeros(length(muscle_sel), 1);
    mat_all = cell(length(muscle_sel), max_order);
    diff_mat_q_all = cell(length(muscle_sel), max_order, nr_samples);
    for m_nr=1:length(muscle_sel)
        muscle_index = muscle_sel(m_nr);

        index_dof_crossing_withknee = ...
            find(muscle_spanning_joint_INFO(muscle_index,:)==1);


        nr_dof_crossing = length(index_dof_crossing_withknee);


        lMT = MuscleData.lMT(:,muscle_index);
        dM = zeros(nr_samples, nr_dof_crossing);
        for dof_nr = 1:nr_dof_crossing
            dM(:,dof_nr) = ...
                MuscleData.dM(:,muscle_index,index_dof_crossing_withknee(dof_nr));
        end

        criterion_full_filled = 0;
        order = 3;
        while criterion_full_filled==0
            fprintf('Order: %d\n', order); % Print the current order
            try
            [mat,diff_mat_q] = n_art_mat_9_GC(q_all(:,index_dof_crossing_withknee), order);
            
            mat_all(m_nr, order) = mat;
            
            %not needed?
            for d= 1:size(diff_mat_q, 3)
                diff_mat_q_all(m_nr, order, d)= diff_mat_q(:,:,d);
            end
            catch
                %keyboard;
            end
            
            
            nr_coeffs = length(mat(1,:));

            diff_mat_q_all = zeros(nr_samples*nr_dof_crossing, nr_coeffs);
            for dof_nr = 1:nr_dof_crossing
                diff_mat_q_all(nr_samples*(dof_nr-1)+1:nr_samples*dof_nr,:)= ...
                    -squeeze(diff_mat_q(:,:,dof_nr));
            end

            coeff{order}=[mat ; diff_mat_q_all]\[lMT; dM(:)];
            dM_recon = zeros(nr_samples, nr_dof_crossing);         

            for dof_nr = 1:nr_dof_crossing
                dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff{order};
            end
            lMT_recon=mat*coeff{order};

            lMT_error_rms{order} = sqrt(mean((lMT - lMT_recon).^2));
            dm_error_rms{order} = sqrt(mean((dM - dM_recon).^2));
            % Print the RMS errors
            fprintf('lMT RMS error: %.4f\n', lMT_error_rms{order});
            fprintf('dM RMS error: %.4f\n', max(dm_error_rms{order}));

            % RMS errors should be lower than threshold for muscle-tendon
            % lengths and moment arms. Otherwise, we increase the
            % polynomial order.
            if any(contains(MuscleData.dof_names(index_dof_crossing_withknee),{'knee_tx','knee_ty','knee_tz'}))
                t_pos=find(contains(MuscleData.dof_names(index_dof_crossing_withknee),{'knee_tx','knee_ty','knee_tz'}));
                r_pos=1:length(index_dof_crossing_withknee);
                r_pos(t_pos)=[];
                criterion_full_filled = ...
                    lMT_error_rms{order}<=threshold_r &...
                    max(dm_error_rms{order}(r_pos))<=threshold_r &...
                    max(dm_error_rms{order}(t_pos))<=threshold_t;
            else
                criterion_full_filled = ...
                    lMT_error_rms{order}<=threshold_r & max(dm_error_rms{order})<=threshold_r;
            end
            if order==max_order
                criterion_full_filled = 1;
            end
            if criterion_full_filled==0
                order = order+1;
            end
        end

        MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing_withknee);
        MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names_forMA{muscle_index};
        MuscleInfo.muscle(m_nr).coeff = coeff;
        MuscleInfo.muscle(m_nr).order = order;
        MuscleInfo.muscle(m_nr).lMT_error_rms = lMT_error_rms;
        MuscleInfo.muscle(m_nr).dm_error_rms = dm_error_rms;

        lMT_all_error(m_nr) = lMT_error_rms{order};
        DM_all_error(m_nr, index_dof_crossing_withknee) = dm_error_rms{order};
        order_all(m_nr) = order;      

        clear coeff lMT_error_rms dm_error_rms;

    end
    save('mat_all.mat', 'mat_all');
    save('diff_mat_q_all.mat', 'diff_mat_q_all');

    figure();
    hold on;
    plot(lMT_all_error)
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold_r, threshold_r], 'r', 'linewidth', 2)
    title('RMS error on the approximated muscle-tendon length')
    ylabel('RMS error (m)')

    figure();
    hold on;
    plot(max(DM_all_error, [], 2))
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold_t, threshold_t], 'r', 'linewidth', 2); %though this is the threshold for translational dofs
    title('maximal RMS error on the approximated muscle moment arm')
    ylabel('RMS error (m)'); 

    figure();
    hold on;
    plot(order_all)
    ylim([0 max_order+1])
    xlimits = get(gca, 'XLim');
    plot(xlimits, [max_order, max_order], 'r', 'linewidth', 2)
    title('Order of the polynomial approximation')
    ylabel('Order')

end
