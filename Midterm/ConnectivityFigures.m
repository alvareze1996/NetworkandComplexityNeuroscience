%% Figures for resting state connectivity 

root_folder=['C:\Users\ea472\Documents'];
%%% Load correaltion matrix
load([root_folder filesep 'conn_parcellation_15_ROI.mat'])

ROI_coords = conn_parcellation_15ROI.parcell_coords;
parcell_coords_name = conn_parcellation_15ROI.parcell_coords_name;

%%%load demographics
load([root_folder filesep 'RiskAmb_Luce_subject_level_characteristics.mat']);


%%% Only keep scanned subjects + with resting state file smoothed
all_subjects=table2array(dems_dx_model_drug_data(:,1));
dfolders=dir(['C:\Users\ea472\Box\RiskAmbAnalysis' filesep 'sub*']);
SubID = [];
for s = 1:size(dfolders,1)
    SubSess = extractAfter(dfolders(s).name,'sub-');
    s_name = SubSess(1:8); %8 digit SubID number
    SubID = [SubID; str2num(s_name)]; %ID stored as strings
end
id_scanned_subjects=ismember(all_subjects, SubID);
diagnosis=table2array(dems_dx_model_drug_data(id_scanned_subjects,4)); %'OUD'/'CONTROL'


%%% Delete missing subjects 
    SubID([36 38 46])=[];
    diagnosis([36 38 46])=[];


%% CONNECTIVITY BETWEEN PATIENTS AND CONTROLS
dx_nb=double(diagnosis); dx_nb(dx_nb==1)=0;  dx_nb(dx_nb==2)=1; %CTR=0 / OUD=1

        %CONTROL
        A_ctr_all=mean(conn_parcellation_15ROI.conn_matrix_strength_Fishers_Z(:,:,dx_nb==0),3); %CTR
        median=nanmean(nanmean(A_ctr_all));
        A_ctr=tril(A_ctr_all); %only keep lower half of the square matrix
        Amed_ctr=A_ctr; Amed_id=find(Amed_ctr<median | Amed_ctr>=1); Amed_ctr(Amed_id)=NaN; %connections with Z>=median (median = 0.64)
        Ahigh_ctr=A_ctr; Ahigh_id=find(Ahigh_ctr<1); Ahigh_ctr(Ahigh_id)=NaN; %connections with Z>=1

        [row_med_ctr, col_med_ctr]=find(~isnan(Amed_ctr));
        x_med_ctr=[ROI_coords(row_med_ctr,1) ROI_coords(col_med_ctr,1)];
        y_med_ctr=[ROI_coords(row_med_ctr,2) ROI_coords(col_med_ctr,2)];
        z_med_ctr=[ROI_coords(row_med_ctr,3) ROI_coords(col_med_ctr,3)];

        [row_high_ctr, col_high_ctr]=find(~isnan(Ahigh_ctr));
        x_high_ctr=[ROI_coords(row_high_ctr,1) ROI_coords(col_high_ctr,1)];
        y_high_ctr=[ROI_coords(row_high_ctr,2) ROI_coords(col_high_ctr,2)];
        z_high_ctr=[ROI_coords(row_high_ctr,3) ROI_coords(col_high_ctr,3)];

        %OUD PATIENTS
        A_oud_all=mean(conn_parcellation_15ROI.conn_matrix_strength_Fishers_Z(:,:,dx_nb==1),3); %OUD
        median=nanmean(nanmean(A_oud_all));
        A_oud=tril(A_oud_all); %only keep lower half of the square matrix
        Amed_oud=A_oud; Amed_id=find(Amed_oud<median | Amed_oud>=1); Amed_oud(Amed_id)=NaN; %connections with Z>=median (median = 0.64)
        Ahigh_oud=A_oud; Ahigh_id=find(Ahigh_oud<1); Ahigh_oud(Ahigh_id)=NaN; %connections with Z>=1

        [row_med_oud, col_med_oud]=find(~isnan(Amed_oud));
        x_med_oud=[ROI_coords(row_med_oud,1) ROI_coords(col_med_oud,1)];
        y_med_oud=[ROI_coords(row_med_oud,2) ROI_coords(col_med_oud,2)];
        z_med_oud=[ROI_coords(row_med_oud,3) ROI_coords(col_med_oud,3)];

        [row_high_oud, col_high_oud]=find(~isnan(Ahigh_oud));
        x_high_oud=[ROI_coords(row_high_oud,1) ROI_coords(col_high_oud,1)];
        y_high_oud=[ROI_coords(row_high_oud,2) ROI_coords(col_high_oud,2)];
        z_high_oud=[ROI_coords(row_high_oud,3) ROI_coords(col_high_oud,3)];

                  
        figure %heatmap
            subplot(1,3,1) %CONTROL
            A_ctr_all(isnan(A_ctr_all))=0;
            heatmap(parcell_coords_name,parcell_coords_name,A_ctr_all)
            title([num2str(length(dx_nb)-sum(dx_nb)) ' CONTROL'])
            colormap jet
            caxis([0 1.6])

            subplot(1,3,2) %OUD
            A_oud_all(isnan(A_oud_all))=0;
            heatmap(parcell_coords_name,parcell_coords_name,A_oud_all)
            title([num2str(sum(dx_nb)) 'OUD PATIENTS'])
            colormap jet
            caxis([0 1.6])

            subplot(1,3,3) %CTR-OUD
            A_diffOUD_CTR=A_oud_all-A_ctr_all;
            heatmap(parcell_coords_name,parcell_coords_name,A_diffOUD_CTR,'ColorLimits',[-0.5 0.5])
            title('difference OUD - CTR')
            colormap jet

            sgtitle(['Resting state connectivity - ' num2str(size(ROI_coords,1)) 'ROI - OUD vs CTR'])
