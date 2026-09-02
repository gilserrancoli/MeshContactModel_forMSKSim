function MuscleInfo=PolyFitKnee(MuscleData,MuscleInfo)

for m=1:length(MuscleInfo.muscle)

    if any(contains(MuscleInfo.muscle(m).DOF,'hip'))&&...
            any(contains(MuscleInfo.muscle(m).DOF,'knee'))
        span.hipknee=1;
        span.ankleknee=0;
        span.kneeMonoArt=0;
    elseif any(contains(MuscleInfo.muscle(m).DOF,'ankle'))&&...
            any(contains(MuscleInfo.muscle(m).DOF,'knee'))
        span.hipknee=0;
        span.ankleknee=1;
        span.kneeMonoArt=0;
    elseif any(contains(MuscleInfo.muscle(m).DOF,'knee'))&&...
            (length(MuscleInfo.muscle(m).DOF)==1)
        span.kneeMonoArt=1;
        span.hipknee=0;
        span.ankleknee=0;
    else
        span.kneeMonoArt=0;
        span.hipknee=0;
        span.ankleknee=0;
    end
    
    dM_knee=squeeze(MuscleData.dM(:,m,[4:5 7:9]));
    if span.hipknee
        q_hip=MuscleData.q(:,1:4);
        dep_side=[ones(length(q_hip),1) q_hip q_hip.^2 q_hip.^3];
        for i=1:5
            coeffs(:,i)=dep_side\dM_knee(:,i);
        end
    elseif span.ankleknee
        q_ankle=MuscleData.q(:,4:6);
        dep_side=[ones(length(q_ankle),1) q_ankle q_ankle.^2 q_ankle.^3];
        for i=1:5
            coeffs(:,i)=dep_side\dM_knee(:,i);
        end
    elseif span.kneeMonoArt
        q_knee=MuscleData.q(:,4);
        dep_side=[ones(length(q_knee),1) q_knee q_knee.^2 q_knee.^3];
        for i=1:5
            coeffs(:,i)=dep_side\dM_knee(:,i);
        end
    else
        coeffs=[];
    end
    MuscleInfo.muscle(m).coeffs_knee=coeffs;
    if ~isempty(coeffs)
        dM_rec(:,:,m)=dep_side*coeffs;
        err_dM(:,:,m)=sqrt((dM_rec(:,:,m)-dM_knee).^2);
    end
    
    clear coeffs;
%     clear err_dM;
end
        
        





end