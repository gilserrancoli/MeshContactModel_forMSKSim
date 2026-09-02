data=importdata('dummy_motion_initial.mot');
data_all=data.data(:,1:10);
colheaders=data.colheaders(1:10);
data_all(:,11)=zeros(size(data.data,1),1); %knee add
colheaders{11}='knee_adduction';
data_all(:,12)=zeros(size(data.data,1),1); %knee rot
colheaders{12}='knee_rotation';
data_all(:,13)=data.data(:,11); %knee flexion
colheaders{13}=data.colheaders{11};
data_all(:,14)=zeros(size(data.data,1),1); %knee tx
colheaders{14}='knee_tx';
data_all(:,15)=zeros(size(data.data,1),1)+0.042; %knee ty
colheaders{15}='knee_ty';
data_all(:,16)=zeros(size(data.data,1),1); %knee tz
colheaders{16}='knee_tz';
data_all(:,17:21)=data.data(:,12:16);
colheaders(17:21)=data.colheaders(12:16);
data_all(:,22)=zeros(size(data.data,1),1); %knee add l
colheaders{22}='knee_adduction_l';
data_all(:,23)=zeros(size(data.data,1),1); %knee rot l
colheaders{23}='knee_rotation_l';
data_all(:,24)=data.data(:,17); %knee flexion l
colheaders{24}=data.colheaders{17};
data_all(:,25)=zeros(size(data.data,1),1); %knee tx l
colheaders{25}='knee_tx_l';
data_all(:,26)=zeros(size(data.data,1),1)+0.042; %knee ty l
colheaders{26}='knee_ty_l';
data_all(:,27)=zeros(size(data.data,1),1); %knee tz l
colheaders{27}='knee_tz_l';
data_all=[data_all data.data(:,18:end)];
colheaders=[colheaders data.colheaders(18:end)];
for i=1:length(colheaders)
    colheaders{i}=strrep(colheaders{i},' ','');
end
q.data=data_all;
q.labels=colheaders;

write_motionFile(q,'dummy_motion.mot');