clear all; clc

SUB_G = "all"; % dep, norm, all

cd('D:\')
psyc_path = 'D:\exp_data\over_generalization\psyc\';

path = ".\exp_data\over_generalization\behavior_all\";

save_path = 'D:\exp_data\over_generalization\Analysis\23_06_final\';

cd(path)
sub_list = dir(fullfile('*.csv'));

noPSSERQ_ind = find(strcmp({sub_list.name}, "TOD01.csv"));

%% data reorganize : Anxiety / Depression score

psyc_table_fin = []; psyc_battery_ori = []; psyc_battery_plt = [];  temp_table_psyc = []; ultimate_sub_data = []; gender_table =[];age_table = [];
sub_i = 0;
for subi = 1 : length(sub_list)
    sub_name = sub_list(subi).name;
    sub_path = sub_list(subi).folder;

    %LOAD TABLE
    temp_table = readtable([sub_path '/' sub_name]);

    if subi < noPSSERQ_ind
        select_var_names = ["age", "res_keys", "res_2_keys","Questionnaire","Question","gender"];
    else
        select_var_names = ["age", "res_keys", "res_2_keys","res_3_keys","Questionnaire","Question","gender"];
    end

    questionnaire_ind = find(strcmp(select_var_names,"Questionnaire"));
    gender_ind = find(strcmp(select_var_names,"gender"));
    age_ind = find(strcmp(select_var_names,"age"));
    question_ind = find(strcmp(select_var_names,"Question"));

    temp_ind = [];
    for vari = 1 : length(select_var_names)
        temp_var_name = select_var_names(vari);
        temp_var_ind = find(string(temp_table.Properties.VariableNames) == temp_var_name);
        temp_ind(vari) = temp_var_ind;
    end

    age_table(subi) = table2array(temp_table(1,temp_ind(age_ind)));
    if subi == 132
        age_table(subi) = 56; % 
    end
    gender_table_1 = table2array(temp_table(1,temp_ind(gender_ind)));

    if  strfind(cell2mat(gender_table_1),'Female') == 1
        gender_table(subi) = 2;
    else
        gender_table(subi) = 1;
    end

    psyc_start_ind = find(ismissing(temp_table(:,temp_ind(questionnaire_ind))),10);

    psyc_start_ind_1 = psyc_start_ind(max(find(psyc_start_ind<100)))+1;
    psyc_start_ind_2 = psyc_start_ind(min(find(psyc_start_ind>psyc_start_ind_1)))-1;


    % ATTENTION QUESTION INDEXING
    temp_table_response = temp_table([psyc_start_ind_1:psyc_start_ind_2],temp_ind([questionnaire_ind question_ind]));

    at_ind = [];
    for at_id = 1 : size(temp_table_response,1)
        temp_at_split = split(temp_table_response.Question(at_id),' ');
        if strcmp(temp_at_split{1},'Attention')
            at_ind = [at_ind, at_id];
        end
    end


    % PSYC RESPONSE INDEXING
    if subi < noPSSERQ_ind
        psyc_answer_ind = [find(strcmp(select_var_names,"res_keys")), find(strcmp(select_var_names,"res_2_keys"))];
    else
        psyc_answer_ind = [find(strcmp(select_var_names,"res_keys")), find(strcmp(select_var_names,"res_2_keys")), find(strcmp(select_var_names,"res_3_keys"))];
    end
    
    temp_table_response = horzcat(temp_table_response,temp_table([psyc_start_ind_1:psyc_start_ind_2],temp_ind([psyc_answer_ind])));


    temp_table_response(at_ind,:) = [];


    
    psyc_response_only = [];
    if subi < noPSSERQ_ind
        temp_table_response_final1 = temp_table_response;
        psyc_1 = [min(find(strcmp(temp_table_response.Questionnaire,'PHQ-9')==1)) : max(find(strcmp(temp_table_response.Questionnaire,'BDI-II')==1))];
        psyc_2 = [min(find(strcmp(temp_table_response.Questionnaire,'PANAS')==1)) : max(find(strcmp(temp_table_response.Questionnaire,'PANAS')==1))];
        temp_table_response_final1(psyc_2,3) = temp_table_response_final1(psyc_2,4);
        temp_table_response_final1(:,[4]) = [];
        psyc_response_only = table2array(temp_table_response_final1(:,3));
    else
        temp_table_response_final2 = temp_table_response;
        psyc_1 = [min(find(strcmp(temp_table_response.Questionnaire,'PHQ-9')==1)) : max(find(strcmp(temp_table_response.Questionnaire,'BDI-II')==1))];
        psyc_2 = [min(find(strcmp(temp_table_response.Questionnaire,'PANAS')==1)) : max(find(strcmp(temp_table_response.Questionnaire,'PSS')==1))];
        psyc_3 = [min(find(strcmp(temp_table_response.Questionnaire,'ERQ')==1)) : max(find(strcmp(temp_table_response.Questionnaire,'ERQ')==1))];
        temp_table_response_final2(psyc_2,3) = temp_table_response_final2(psyc_2,4);
        temp_table_response_final2(psyc_3,3) = temp_table_response_final2(psyc_3,5);

        temp_table_response_final2(:,[4 5]) = [];
        psyc_response_only = table2array(temp_table_response_final2(:,3));
    end

    if subi< noPSSERQ_ind
        psyc_battery_plt(subi,:) = psyc_response_only;
    else
        sub_i = sub_i+1;
        psyc_battery_ori(sub_i,:) = psyc_response_only;
    end
end

sub_number =  1 : length(sub_list);
sub_number= sub_number(:);
sub_number = table(sub_number,'VariableNames',"Sub ID");
% CALCULATION
  CESD= []; CESD_dep = []; CESD_pos  = []; CESD_soma  =[]; CESD_interper = []; PHQ9 = []; PHQ_soma = []; PHQ_cog = []; BDI = []; BDI_soma = []; BDI_cog = []; BDI_aff = [];
    GAD_7 =[]; STAI_2 = [];
for group_i = 1 : 2 % 1= pilot 2 = ori
  
    if group_i == 1
        psyc_battery = psyc_battery_plt;
        scale_i = 5;
        PHQ_ind = find(strcmp(temp_table_response_final1.Questionnaire,'PHQ-9'));
        CESD_ind = find(strcmp(temp_table_response_final1.Questionnaire,'CES-D'));
        BDI_ind = find(strcmp(temp_table_response_final1.Questionnaire,'BDI-II'));
        STAI2_ind = find(strcmp(temp_table_response_final1.Questionnaire,'STAI-2'));
        PANAS_ind = find(strcmp(temp_table_response_final1.Questionnaire,'PANAS'));
        
    else
        psyc_battery = psyc_battery_ori;
        PHQ_ind = find(strcmp(temp_table_response_final2.Questionnaire,'PHQ-9'));
        CESD_ind = find(strcmp(temp_table_response_final2.Questionnaire,'CES-D'));
        BDI_ind = find(strcmp(temp_table_response_final2.Questionnaire,'BDI-II'));
        STAI2_ind = find(strcmp(temp_table_response_final2.Questionnaire,'STAI-2'));
        PANAS_ind = find(strcmp(temp_table_response_final2.Questionnaire,'PANAS'));
        GAD_ind = find(strcmp(temp_table_response_final2.Questionnaire,'GAD-7'));
        PSS_ind = find(strcmp(temp_table_response_final2.Questionnaire,'PSS'));
        ERQ_ind = find(strcmp(temp_table_response_final2.Questionnaire,'ERQ'));
        scale_i = 8;
    end

    for scale = 1 : scale_i
        if scale == 1 % PHQ-9
            PHQ9 = psyc_battery(:,PHQ_ind);
            PHQ9 = PHQ9-1;
            PHQ_soma = sum(PHQ9(:,[1 2 3 4 5]),2);
            PHQ_cog = sum(PHQ9(:,[6 7 8 9]),2);
            PHQ9 = sum(PHQ9,2);

        elseif scale == 2 % CES-D
            CESD = psyc_battery(:,CESD_ind);
            CESD = CESD-1;
            for i = [ 4 8 12 16]
                cesd_ind_0 = find(CESD(:,i) == 0);
                cesd_ind_1 = find(CESD(:,i) == 1);
                cesd_ind_2 = find(CESD(:,i) == 2);
                cesd_ind_3 = find(CESD(:,i) == 3);
                CESD(cesd_ind_0,i) = 3;
                CESD(cesd_ind_1,i) = 2;
                CESD(cesd_ind_2,i) = 1;
                CESD(cesd_ind_3,i) = 0;
            end
                CESD_dep = sum(CESD(:,[3 6 14 17 18]),2);
                CESD_pos = sum(CESD(:,[4 12 16]),2);
                CESD_soma = sum(CESD(:,[2 5 7 11 13 20]),2);
                CESD_interper = sum(CESD(:,[15 19]),2);
                CESD = sum(CESD,2);

        elseif scale == 3 % BDI
            BDI = psyc_battery(:, [BDI_ind]);
            BDI = BDI-1;
            BDI_aff = sum(BDI(:,[4 10 12]),2);
            BDI_cog = sum(BDI(:,[2 3 6 8 9 14]),2);
            BDI_soma = sum(BDI(:,[16 17 18 19 20 21]),2);
            BDI = sum(BDI,2);

        elseif scale == 6 % GAD-7
            GAD_7 = psyc_battery(:,GAD_ind);
            GAD_7 = GAD_7-1;
            GAD_7 = sum(GAD_7,2);

        elseif scale == 4 % STAI-2
            STAI_2 = psyc_battery(:,[STAI2_ind]);
            for i = [1 3 6 7 10 13 14 16 19]
                stai_ind_1 = find(STAI_2(:,i) == 1);
                stai_ind_2 = find(STAI_2(:,i) == 2);
                stai_ind_3 = find(STAI_2(:,i) == 3);
                stai_ind_4 = find(STAI_2(:,i) == 4);
                
                STAI_2(stai_ind_4,i) = 1;
                STAI_2(stai_ind_3,i) = 2;
                STAI_2(stai_ind_2,i) = 3;
                STAI_2(stai_ind_1,i) = 4;

            end
            STAI_2 = sum(STAI_2,2);
       
        elseif scale == 5 % PANAS
            PANAS = psyc_battery(:,[PANAS_ind]);
            PANAS_P = PANAS(:,[1 3 5 7 9 11 13 15 17 19]);
            PANAS_P = sum(PANAS_P,2);
            PANAS_N = PANAS(:,[2 4 6 8 10 12 14 16 18 20]);
            PANAS_N = sum(PANAS_N,2);

        elseif scale == 8
            PSS = psyc_battery(:,PSS_ind);
            PSS = PSS-1;
            for i = [4 5 7 8]
                pss_ind_0 = find(PSS(:,i) == 0);
                pss_ind_1 = find(PSS(:,i) == 1);
                pss_ind_2 = find(PSS(:,i) == 2);
                pss_ind_3 = find(PSS(:,i) == 3);
                pss_ind_4 = find(PSS(:,i) == 4);
                PSS(pss_ind_4,i) = 0;
                PSS(pss_ind_3,i) = 1;
                PSS(pss_ind_2,i) = 2;
                PSS(pss_ind_1,i) = 3;
                PSS(pss_ind_0,i) = 4;
            end
            PSS = sum(PSS,2);
        elseif scale == 7
            ERQ = psyc_battery(:,ERQ_ind);
            ERQ_R = sum(ERQ(:,[1 3 5 7 8 10]),2);
            ERQ_S = sum(ERQ(:,[2 4 6 9]),2);
        end
    end
    if group_i == 1
        psyc_table_plt = table(PHQ9,PHQ_cog,PHQ_soma,CESD,CESD_dep,CESD_pos,CESD_soma,CESD_interper,BDI,BDI_aff,BDI_cog,BDI_soma,STAI_2,PANAS_P, PANAS_N,'VariableNames',{'PHQ-9','PHQ_cognitive', 'PHQ_somatic',...
            'CES-D', 'CESD_depression', 'CESD_Positive', 'CESD_somatic', 'CESD_interpersonal', 'BDI-II','BDI_affect','BDI_Cognitive' ,'BDI_somatic' , 'STAI-Y2', 'PANAS_P', 'PANAS_N'});
    else
        psyc_table_ori = table(PHQ9,PHQ_cog,PHQ_soma,CESD,CESD_dep,CESD_pos,CESD_soma,CESD_interper,BDI,BDI_aff,BDI_cog,BDI_soma,STAI_2,PANAS_P, PANAS_N,GAD_7,ERQ_R,ERQ_S,PSS,'VariableNames', {'PHQ-9','PHQ_cognitive', 'PHQ_somatic', ...
            'CES-D', 'CESD_depression', 'CESD_Positive', 'CESD_somatic', 'CESD_interpersonal', 'BDI-II','BDI_affect','BDI_Cognitive' ,'BDI_somatic' , 'STAI-Y2', 'PANAS_P', 'PANAS_N', 'GAD7', 'ERQ-R', 'ERQ-S', 'PSS'});
    end
end


age_table = table(age_table','VariableNames',"age");
gender_table = table(gender_table','VariableNames',"gender");

psyc_table_fin = vertcat(psyc_table_plt,psyc_table_ori(:,[1:15]));
psyc_table_fin = horzcat(age_table,gender_table,psyc_table_fin);


young_age = find(table2array(age_table)<41);
old_age = find(table2array(age_table)>40);
age_index = table(double(table2array(age_table)<41),'VariableNames', "Age index");

% psychology data
aux_func_get_residual = @(mdl) mdl.Residuals.Standardized;
func_get_residual = @(x,y) aux_func_get_residual(fitlm(x(:),y(:)));

index_psyc = [5,17];
data_psyc = table2array(psyc_table_fin(:,index_psyc));
data_psyc = arrayfun(@(i) data_psyc(:,i), 1:2, 'uni', 0);

data_psyc2 = table2array(psyc_table_fin(:,13));
phq_residual = func_get_residual(data_psyc{1}, data_psyc{2});
stai_residual = func_get_residual(data_psyc{2}, data_psyc{1});
stai_residual2 = func_get_residual(data_psyc{2}, data_psyc2);
bdi_residual = func_get_residual(data_psyc2, data_psyc{2});

psyc_residual_table = table(phq_residual,stai_residual,stai_residual2,bdi_residual,'VariableNames', {'PHQ residual', 'STAI residual','STAI residual2','BDI residual'});


psyc_table_fin = horzcat(sub_number,age_index,psyc_table_fin,psyc_residual_table);

young_psyc_data = psyc_table_fin(young_age,:);

%% Valence Arousal rating & Memory data
complete_data = []; emo_ind_rating = [];  emo_ind_ori = []; aro_stim_check = []; sorted_valence = []; sorted_arousal = [];
aro_ind_rating_equal = []; aro_ind_rating =[];

%%%%%
stim_remove_opt = 0;  % not valence matching stim removal
%%%%

complete_data.psyc = psyc_table_fin;

for subi = 1 : length(sub_list)

    sub_name = sub_list(subi).name;
    sub_path = sub_list(subi).folder;

    %LOAD TABLE
    temp_table = readtable([sub_path '/' sub_name]);


    if subi < noPSSERQ_ind
        select_var_names = ["fore_list_real", "back_list_real", "loc_list_real", "valence_resp_response", "arousal_resp_response","stroop_resp_corr", "recog_pic_index", "pic_answer", "loc_answer", "lure_response", "pic_response","confi_response_recog_response", "loc_response", "confi_response_Loc_response","stroop_resp_rt","answer_loc","stroop_resp_corr","stroop_thisIndex"];
        var_indices = find(ismember(temp_table.Properties.VariableNames, select_var_names));
        temp_table = temp_table(:, var_indices);

        % table reconstruction
        indices = find(ismember(temp_table.loc_list_real, [1,2,3,4]));
        encod_st_ind = indices(1);
        indices = find(ismember(temp_table.pic_answer, [1 2 3 4]));
        recog_st_ind = indices(1);
        indices = find(ismember(temp_table.stroop_resp_corr, [0 1]));
        stroop_st_ind = indices(1);
        table_encod = temp_table([encod_st_ind:encod_st_ind+79], [1:5]);  %% 79
        table_encod(find(isnan(table_encod.loc_list_real)),:) = [];
        table_recog = temp_table([recog_st_ind:recog_st_ind+71],[9:16]); %% 71
        table_recog = sortrows(table_recog,"recog_pic_index");
        table_stroop = temp_table([stroop_st_ind:stroop_st_ind+99],[6:8]); %% 99
        table_new = [table_encod, table_recog];

        % FIXING THE RECOGNITION RESPONSE
        fixed_lure_response = table_new.lure_response;
        ind_0to4 = find(fixed_lure_response == 0);
        fixed_lure_response(ind_0to4) = 4;
        table_new.lure_response = fixed_lure_response;

    else % RECOGNITION RT ADD
        select_var_names = ["fore_list_real", "back_list_real", "loc_list_real", "valence_resp_response", "arousal_resp_response","stroop_resp_corr", "recog_pic_index", "pic_answer", "loc_answer", "lure_response", "pic_response","confi_response_recog_response", "loc_response", "confi_response_Loc_response","stroop_resp_rt","answer_loc","pic_resp_time","loc_resp_time","stroop_resp_corr","stroop_thisIndex"];
        var_indices = find(ismember(temp_table.Properties.VariableNames, select_var_names));
        temp_table = temp_table(:, var_indices);

        % table reconstruction
        indices = find(ismember(temp_table.loc_list_real, [1,2,3,4]));
        encod_st_ind = indices(1);
        indices = find(ismember(temp_table.pic_answer, [1 2 3 4]));
        recog_st_ind = indices(1);
        indices = find(ismember(temp_table.stroop_resp_corr, [0 1]));
        stroop_st_ind = indices(1);


        table_encod = temp_table([encod_st_ind:encod_st_ind+79], [1:5]);  %% 79
        table_encod(find(isnan(table_encod.loc_list_real)),:) = [];
        table_recog = temp_table([recog_st_ind:recog_st_ind+71],[9:19]); %% 71
        table_recog = sortrows(table_recog,"recog_pic_index");
        table_stroop = temp_table([stroop_st_ind:stroop_st_ind+99],[6:8]); %% 99


        table_new = [table_encod, table_recog];
    end

    fore_name = arrayfun(@(i) regexprep(table_new.fore_list_real{i}, '(?<=\D)(\d)(?=\D*\.png)', '0$1'),1:size(table_new,1),'uni',0);
table_new(:,1) = fore_name';

    table_new = sortrows(table_new,1);

    % location accuracy recreation
    cor_loc_ind = find(table_new.loc_answer == table_new.loc_response);
    table_new.loc_cor = zeros(size(table_new.loc_response));
    table_new.loc_cor(cor_loc_ind) = 1;

    

    % NOT matching valence stim removal
    if stim_remove_opt == 1
        table_new([4 7 18 21 22 28 31 38 48 49 50 54 56 59 64 71 72],:) = [];
    end

    complete_data.memory{subi} = table_new;
    complete_data.memory_stroop{subi} = table_stroop;

    %%%%%%%%% NEW PARAMETER
    tri_N = size(table_new,1);


    % Stim 별 valence arousal rating GATHERING

    % original emotion
    temp_emo_ind_ori = table_new.fore_list_real;


    for tri_i = 1 : size(table_new.fore_list_real,1)
        if strfind(cell2mat(temp_emo_ind_ori(tri_i,1)), 'neg') == 1
            temp_ind = -1;
%            
        elseif strfind(cell2mat(temp_emo_ind_ori(tri_i,1)), 'neu') == 1
            temp_ind = 0;
%            
        elseif strfind(cell2mat(temp_emo_ind_ori(tri_i,1)), 'pos') == 1
            temp_ind = 1;
%             
        end
        emo_ind_ori (tri_i,subi) = temp_ind;
    end




    complete_data.valence_all(:,subi) = table_new.valence_resp_response;
    complete_data.arousal_all(:,subi) = table_new.arousal_resp_response;

    complete_data.valence{subi,1} = table_new.valence_resp_response(find(emo_ind_ori(:,subi) == -1));
    complete_data.valence{subi,2} = table_new.valence_resp_response(find(emo_ind_ori(:,subi) == 0));
    complete_data.valence{subi,3} = table_new.valence_resp_response(find(emo_ind_ori(:,subi) == 1));

    complete_data.arousal{subi,1} = table_new.arousal_resp_response(find(emo_ind_ori(:,subi) == -1));
    complete_data.arousal{subi,2} = table_new.arousal_resp_response(find(emo_ind_ori(:,subi) == 0));
    complete_data.arousal{subi,3} = table_new.arousal_resp_response(find(emo_ind_ori(:,subi) == 1));

    complete_data.val_mean(subi,1) = mean(complete_data.valence{subi,1});
    complete_data.val_mean(subi,2) = mean(complete_data.valence{subi,2});
    complete_data.val_mean(subi,3) = mean(complete_data.valence{subi,3});

    complete_data.aro_mean(subi,1) = mean(complete_data.arousal{subi,1});
    complete_data.aro_mean(subi,2) = mean(complete_data.arousal{subi,2});
    complete_data.aro_mean(subi,3) = mean(complete_data.arousal{subi,3});



    % emotion redistribution based on individual rating
    temp_emo_ind = zeros(tri_N,1);
    temp_aro_ind = zeros(tri_N,1);

    val_temp = [[1:tri_N]',table_new.valence_resp_response];
    val_temp = sortrows(val_temp,2);
    neg_ind = val_temp(1:max(find(emo_ind_ori(:,subi) == -1)),1);
    neu_ind = val_temp(min(find(emo_ind_ori(:,subi) == 0)):max(find(emo_ind_ori(:,subi) == 0)),1);
    pos_ind = val_temp(min(find(emo_ind_ori(:,subi) == 1)):max(find(emo_ind_ori(:,subi) == 1)),1);

    temp_emo_ind(neg_ind) = -1;
    temp_emo_ind(neu_ind) = 0;
    temp_emo_ind(pos_ind) = 1;

    emo_ind_rating(:,subi) = temp_emo_ind;


    % low high arousal stimuli indexing
    aro_temp = [[1:tri_N]',table_new.arousal_resp_response];
    aro_temp = sortrows(aro_temp,2);
    low_ind = aro_temp(1:ceil(tri_N/2),1);
    %     med_ind = aro_temp(25:48,1);
    high_ind = aro_temp(ceil(tri_N/2)+1:tri_N,1);

    
    temp_aro_ind(low_ind) = -1;
    temp_aro_ind(high_ind) = 1;
    aro_ind_rating(:,subi) = temp_aro_ind;


    aro_stim_check{subi} = table(table_new.fore_list_real, emo_ind_ori(:,subi), table_new.arousal_resp_response, aro_ind_rating(:,subi), 'VariableNames', {'pic_id','emotion_index','Arousal','Arousal_ind'});


    sub_emo_ind = emo_ind_ori(:,subi);
    aro_temp_equal = [[1:tri_N]', sub_emo_ind ,table_new.arousal_resp_response, zeros(tri_N,1)];
    neg_ind = sortrows(aro_temp_equal(sub_emo_ind == -1,:),3);
    neg_ind(1:ceil(length(neg_ind)/2),4) = -1;
    neg_ind(ceil(length(neg_ind)/2):length(neg_ind),4) = 1;
    neg_median(subi,1) = median(neg_ind(:,3));

    neu_ind = sortrows(aro_temp_equal(sub_emo_ind == 0,:),3);
    neu_ind(1:ceil(length(neu_ind)/2),4) = -1;
    neu_ind(ceil(length(neu_ind)/2):length(neg_ind),4) = 1;
    neu_median(subi,1) = median(neu_ind(:,3));

    pos_ind = sortrows(aro_temp_equal(sub_emo_ind == 1,:),3);
    pos_ind(1:ceil(length(pos_ind)/2),4) = -1;
    pos_ind(ceil(length(neg_ind)/2):length(pos_ind),4) = 1;
    pos_median(subi,1) = median(pos_ind(:,3));


    aro_temp_equal = [neg_ind;neu_ind;pos_ind];
    aro_temp_equal = sortrows(aro_temp_equal,1);

    aro_ind_rating_equal(:,subi) = aro_temp_equal(:,4);
    aro_stim_check_equal{subi} = table(table_new.fore_list_real, aro_temp_equal(:,2), aro_temp_equal(:,3), aro_temp_equal(:,4), 'VariableNames', {'pic_id','emotion_index','Arousal','Arousal_ind'});
    aro_stim_check_equal{subi} = sortrows(aro_stim_check_equal{subi},1);
end


low_aro_mean_neg = []; high_aro_mean_neg = [];
for subi = 1 : size(aro_stim_check_equal,2)
    low_aro_mean_neg(subi,1) = mean(aro_stim_check_equal{subi}.Arousal(aro_stim_check_equal{subi}.emotion_index == -1 & aro_stim_check_equal{subi}.Arousal_ind == -1));
    high_aro_mean_neg(subi,1) = mean(aro_stim_check_equal{subi}.Arousal(aro_stim_check_equal{subi}.emotion_index == -1 & aro_stim_check_equal{subi}.Arousal_ind == 1));
end





save([save_path 'Behavior_complete_set_23_06.mat'],'complete_data','-mat');


sorted_arousal = sorted_arousal;
sorted_valence = sorted_valence;


%% young group extract

excludeIndex = old_age; % Old group exclusion

complete_data_young = struct();

for fn = fieldnames(complete_data)'
    fieldName = fn{1};

    currentFieldData = complete_data.(fieldName);


    if iscell(currentFieldData)
        if fieldName == "valence" || fieldName == "arousal"
            currentFieldData(excludeIndex,:) = [];
        else
            currentFieldData(:, excludeIndex) = [];
        end

        complete_data_young.(fieldName) = currentFieldData;
    elseif istable(currentFieldData)

        currentFieldData(excludeIndex, :) = [];

        complete_data_young.(fieldName) = currentFieldData;
    elseif ismatrix(currentFieldData)

                if fieldName == "valence_all" || fieldName == "arousal_all"
                    currentFieldData(: , excludeIndex) = [];
                else
                    currentFieldData(excludeIndex, :) = [];
                end

        complete_data_young.(fieldName) = currentFieldData;
    end
end

%% Memory performance based subject rejection

% input_data = complete_data_young.memory;
input_data = complete_data.memory;

input_data2 = complete_data;
% input_data2 = complete_data_young;
% input_data = input_data(1,young_age);
complete_data_fixed = struct();

Memory_performance = [];
for subi = 1 : length(input_data)

    temp_memory_data = input_data{subi};
    temp_memory_data = sortrows(temp_memory_data,1);

    temp_recognition_response = temp_memory_data.lure_response;  % 1 Full memory 2 Event Hit 3 Background Hit 4 Lure
    temp_location_response(subi) = sum(temp_memory_data.loc_response./temp_memory_data.loc_answer ==1) / tri_N; 
    Memory_performance(subi) = length(find(temp_recognition_response ~= 4)) / length (temp_recognition_response);
% Memory_performance(subi) = length(find(ismember(find(temp_recognition_response ~= 4) , find(temp_recognition_response ~= 3)))) / length (temp_recognition_response);
end

% MEMORY PORFORMANCE DISTRIBUTION
figure_JK(1)
histogram(Memory_performance,length(input_data))

figure_JK(2)
% performance_outlier = ~isoutlier(Memory_performance);
% performance_outlier = .50 < Memory_performance & .25 < temp_location_response;
performance_outlier = .50 < Memory_performance ;
Memory_performance = Memory_performance(performance_outlier);
histogram(Memory_performance,length(Memory_performance))


% SUBEJECT DATA Exclusion 
excludeIndex = find(performance_outlier == 0); % 

complete_data_fixed = struct();

for fn = fieldnames(input_data2)'
    fieldName = fn{1};

    currentFieldData = input_data2.(fieldName);

    if iscell(currentFieldData)
        if fieldName == "valence" || fieldName == "arousal"
            currentFieldData(excludeIndex,:) = [];
        else
            currentFieldData(:, excludeIndex) = [];
        end

        complete_data_fixed.(fieldName) = currentFieldData;
    elseif istable(currentFieldData)

        currentFieldData(excludeIndex, :) = [];

        complete_data_fixed.(fieldName) = currentFieldData;
    elseif ismatrix(currentFieldData)

                if fieldName == "valence_all" || fieldName == "arousal_all"
                    currentFieldData(: , excludeIndex) = [];
                else
                    currentFieldData(excludeIndex, :) = [];
                end

        complete_data_fixed.(fieldName) = currentFieldData;
    end
end

length(Memory_performance)
%% Arousal control

arousal_controlled_fore_list = [];
arousal_neg = [];  arousal_pos = [];  pos_limit_outlier = []; neg_limit_outlier = [];
for sbj_i = 1:length(complete_data_fixed.memory)
    memory = complete_data_fixed.memory{sbj_i};

    temp = memory.fore_list_real;
    trial_type = [];
    for trial_i = 1:length(temp)
        if contains(temp{trial_i},'neg')
            trial_type(trial_i) = -1;
        elseif contains(temp{trial_i},'neu')
            trial_type(trial_i) = 0;
        elseif contains(temp{trial_i},'pos')
            trial_type(trial_i) = 1;
        end
    end

    temp_arou = memory.arousal_resp_response;

    arousal_neg = temp_arou(trial_type==-1)
    arousal_pos = temp_arou(trial_type==1)

%     pos_limit = find(arousal_pos < prctile(arousal_neg, ))
    pos_limit = find(arousal_pos < min(arousal_neg))
    neg_limit = find(arousal_neg > max(arousal_pos))

    arousal_neg(neg_limit) = [];
    arousal_pos(pos_limit) = [];

    pos_limit_count(sbj_i) = length(pos_limit);
    neg_limit_count(sbj_i) = length(neg_limit);

    if pos_limit_count(sbj_i) > 10
        pos_limit_outlier =[pos_limit_outlier, sbj_i] ;
    end
    if neg_limit_count(sbj_i) > 10
        neg_limit_outlier = [neg_limit_outlier, sbj_i];
    end
        
    control_check(sbj_i) = ttest2(arousal_pos, arousal_neg);

    pos_limit = pos_limit + length(neg_ind) + length(neu_ind);

    temp(neg_limit) = arrayfun(@(x) strrep(temp(x), 'neg', 'none'), neg_limit)
    temp(pos_limit) = arrayfun(@(x) strrep(temp(x), 'pos', 'none'), pos_limit)
     
    arousal_controlled_fore_list{sbj_i} = temp;
    
end
arousal_ouliter = unique([pos_limit_outlier, neg_limit_outlier]);



excludeIndex = arousal_ouliter; % 

complete_data_fixed_arousal = struct();

for fn = fieldnames(complete_data_fixed)'
    fieldName = fn{1};

    currentFieldData = complete_data_fixed.(fieldName);


    if iscell(currentFieldData)
        if fieldName == "valence" || fieldName == "arousal"
            currentFieldData(excludeIndex,:) = [];
        else
            currentFieldData(:, excludeIndex) = [];
        end

        complete_data_fixed_arousal.(fieldName) = currentFieldData;
    elseif istable(currentFieldData)

        currentFieldData(excludeIndex, :) = [];

        complete_data_fixed_arousal.(fieldName) = currentFieldData;
    elseif ismatrix(currentFieldData)

                if fieldName == "valence_all" || fieldName == "arousal_all"
                    currentFieldData(: , excludeIndex) = [];
                else
                    currentFieldData(excludeIndex, :) = [];
                end

        complete_data_fixed_arousal.(fieldName) = currentFieldData;
    end
end

arousal_controlled_fore_list(excludeIndex) =[];



%% STIMULI AROUSAL CONTROL _ SUBJECT WISE(WEIRD TRIAL REJECTION) + Emotion indexing

input_data = complete_data_fixed;
sub_num = length(complete_data_fixed.valence);
sorted_valence = complete_data_fixed.valence_all;
sorted_arousal = complete_data_fixed.arousal_all;
% Weird catch

neg_weird = zeros(1,length(neg_ind)); neu_weird = zeros(1,length(neu_ind)); pos_weird = zeros(1,length(pos_ind));
neg_range = 1:length(neg_ind);
neu_range = length(neg_ind)+1:length(neg_ind)+length(neu_ind);
pos_range = length(neg_ind)+length(neu_ind)+1:length(neg_ind)+length(neu_ind)+length(pos_ind);

i = 0;

% ORIGINAL EMOTION INDEXING
emo_ind_ori = [];
emo_ind_ori(neg_range) = -1;
emo_ind_ori(neu_range) = 0;
emo_ind_ori(pos_range) = 1;
emo_ind_ori2 = ones(tri_N,sub_num) .* emo_ind_ori' ;
%%%%

% SUBJECT LEVEL INDEXING
emo_ind_ori_sub = ones(tri_N,sub_num) .* emo_ind_ori';

emo_ind_val_control = ones(1,tri_N)*3; % SETTING THE EXCLUSION TRIAL AS 3       EMOTION INDEX = -1 0 1

    neu_val = sorted_valence(neu_range,:)';
    
    neg_val = sorted_valence(neg_range,:)';
    
    pos_val = sorted_valence(pos_range,:)';

    neg_aro = sorted_arousal(neg_range,:)';
    neu_aro = sorted_arousal(neu_range,:)';
    pos_aro = sorted_arousal(pos_range,:)';
    
    for emo_i = 1 : 3
        if emo_i == 1 % neg
            tri_n_emo = length(neg_ind);
        elseif emo_i == 2 % neu
            tri_n_emo = length(neu_ind);
        else % pos
            tri_n_emo = length(pos_ind);
        end
        for tri_i = 1 : tri_n_emo

            if emo_i == 1
                valaro_out_neg_sub{tri_i} = find(isoutlier(neg_val(:,tri_i))==1 | isoutlier(neg_aro(:,tri_i))==1);
            elseif emo_i == 2
                valaro_out_neu_sub{tri_i} = find(isoutlier(neu_val(:,tri_i))==1 | isoutlier(neu_aro(:,tri_i))==1);
            else
                valaro_out_pos_sub{tri_i} = find(isoutlier(pos_val(:,tri_i))==1 | isoutlier(pos_aro(:,tri_i))==1);
            end

            emo_ind_ori_sub(tri_i,valaro_out_neg_sub{tri_i}) = 3;
            emo_ind_ori_sub(tri_i+neg_ind,valaro_out_neu_sub{tri_i}) = 3;
            emo_ind_ori_sub(tri_i+neg_ind+neu_ind,valaro_out_pos_sub{tri_i}) = 3;
        end
    end
    trial_count = [];
    for subi = 1 : size(emo_ind_ori_sub,2)
        neg_tri = emo_ind_ori_sub(neg_range, subi);
        neu_tri = emo_ind_ori_sub(neu_range, subi);
        pos_tri = emo_ind_ori_sub(pos_range, subi);



        trial_count{1}(subi) = length(find((neg_tri~= 3)==1));
        trial_count{2}(subi) = length(find((neu_tri~= 3)==1));
        trial_count{3}(subi) = length(find((pos_tri~= 3)==1));
    end

    length(find(trial_count{2}<12))
out_sub_index = find(trial_count{2}<12) % LACK OF NEUTRAL TRIALS

emo_ind_ori_sub(:,out_sub_index) = [];




% SUBEJECT DATA Exclusion 
excludeIndex = out_sub_index; % 

complete_data_fixed2 = struct();

for fn = fieldnames(complete_data_young)'
    fieldName = fn{1};

    currentFieldData = complete_data_fixed.(fieldName);


    if iscell(currentFieldData)
        if fieldName == "valence" || fieldName == "arousal"
            currentFieldData(excludeIndex,:) = [];
        else
            currentFieldData(:, excludeIndex) = [];
        end

        complete_data_fixed2.(fieldName) = currentFieldData;
    elseif istable(currentFieldData)

        currentFieldData(excludeIndex, :) = [];

        complete_data_fixed2.(fieldName) = currentFieldData;
    elseif ismatrix(currentFieldData)

                if fieldName == "valence_all" || fieldName == "arousal_all"
                    currentFieldData(: , excludeIndex) = [];
                else
                    currentFieldData(excludeIndex, :) = [];
                end

        complete_data_fixed2.(fieldName) = currentFieldData;
    end
end


%% AROUSAL VALENCE PLOT
input_opt = 1 ; % 1 original 2 arousal control

if input_opt == 1
    input_data = complete_data_fixed;
elseif input_opt == 2
    input_data = complete_data_fixed_arousal;
end
figure_JK()
data_val = {input_data.val_mean(:,1),input_data.val_mean(:,2),input_data.val_mean(:,3)};
ms_bar(data_val)
stats= anova([data_val{1},data_val{2},data_val{3}]);
multcompare(stats)

xticklabels(["Neg","Neu","Pos"])
title("Valence")

figure_JK()
data_aro = {input_data.aro_mean(:,1),input_data.aro_mean(:,2),input_data.aro_mean(:,3)};
ms_bar(data_aro)
stats= anova([data_aro{1},data_aro{2},data_aro{3}]);
multcompare(stats)

xticklabels(["Neg","Neu","Pos"])
title("Arousal")

val_table = table(data_val{1},data_val{2},data_val{3},'VariableNames',["Neg_val","Neu_val","Pos_val"]);

aro_table = table(data_aro{1},data_aro{2},data_aro{3},'VariableNames',["Neg_aro","Neu_aro","Pos_aro"]);
group_table = ([ones(size(data_aro{1})), ones(size(data_aro{1}))*2 , ones(size(data_aro{1}))*3]);


%% Memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arousal_opt = 0;   % 1 = arousal control      0 = no control
group_opt = 0; % 0 = all     1 = normal     2 = anxiety only     3 = both
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PSYC GROUP INDEXING
if arousal_opt == 0
    complete_psyc_dep = complete_data_fixed.psyc.("BDI-II");
    complete_psyc_anx = complete_data_fixed.psyc.("STAI-Y2");
elseif arousal_opt == 1
    complete_psyc_dep = complete_data_fixed_arousal.psyc.("BDI-II");
    complete_psyc_anx = complete_data_fixed_arousal.psyc.("STAI-Y2");
end

 anxiety_group = (complete_psyc_anx > 37 & complete_psyc_dep < 20); % 1
 normal_group = (complete_psyc_anx < 38 & complete_psyc_dep < 20); % 2
 comorbid_group = (complete_psyc_anx >= 38 & complete_psyc_dep >= 20); % 3

dep_rank_1 = complete_psyc_dep <= prctile(complete_psyc_dep,33.3); % 4
dep_rank_2 = (complete_psyc_dep <= prctile(complete_psyc_dep,66.6) & complete_psyc_dep > prctile(complete_psyc_dep,33.3));% 5
dep_rank_3 = (complete_psyc_dep <= prctile(complete_psyc_dep,100) & complete_psyc_dep > prctile(complete_psyc_dep,66.6));% 6

anx_rank_1 = complete_psyc_anx <= prctile(complete_psyc_anx,33.3);% 7
anx_rank_2 = (complete_psyc_anx <= prctile(complete_psyc_anx,66.6) & complete_psyc_anx > prctile(complete_psyc_anx,33.3));% 8
anx_rank_3 = (complete_psyc_anx <= prctile(complete_psyc_anx,100) & complete_psyc_anx > prctile(complete_psyc_anx,66.6));% 9

%%%%%%%%%%%%%%%%%%%%%%%%%%%
if arousal_opt == 1
    data_ori = complete_data_fixed_arousal;
elseif arousal_opt == 0
    data_ori = complete_data_fixed;
    if group_opt == 1
        group_ind = normal_group;
    elseif group_opt == 2
        group_ind = anxiety_group;
    elseif group_opt == 3
        group_ind = comorbid_group;
    elseif group_opt == 4
        group_ind = dep_rank_1;
    elseif group_opt == 5
        group_ind = dep_rank_2;
    elseif group_opt == 6
        group_ind = dep_rank_3;
    elseif group_opt == 7
        group_ind = anx_rank_1;
    elseif group_opt == 8
        group_ind = anx_rank_2;
    elseif group_opt == 9
        group_ind = anx_rank_3;
    elseif group_opt == 10
        group_ind = bias_rank_1;
    elseif group_opt == 11
        group_ind = bias_rank_2;
    elseif group_opt == 12
        group_ind = bias_rank_3;
    end

end


% Grouping index change
anxiety_group = (complete_psyc_anx > 37 & complete_psyc_dep < 10)*2;
normal_group = (complete_psyc_anx < 38 & complete_psyc_dep < 10);
comorbid_group = (complete_psyc_anx >= 38 & complete_psyc_dep >= 10)*3;

% dep_rank_1 = complete_psyc_dep <= prctile(complete_psyc_dep,33.3);
% dep_rank_2 = (complete_psyc_dep <= prctile(complete_psyc_dep,66.6) & complete_psyc_dep > prctile(complete_psyc_dep,33.3))*2;
% dep_rank_3 = (complete_psyc_dep <= prctile(complete_psyc_dep,100) & complete_psyc_dep > prctile(complete_psyc_dep,66.6))*3;
dep_rank_1 = complete_psyc_dep <= prctile(complete_psyc_dep,50);
dep_rank_3 = (complete_psyc_dep <= prctile(complete_psyc_dep,100) & complete_psyc_dep > prctile(complete_psyc_dep,50))*2;

anx_rank_1 = complete_psyc_anx <= prctile(complete_psyc_anx,33.3);
anx_rank_2 = (complete_psyc_anx <= prctile(complete_psyc_anx,66.6) & complete_psyc_anx > prctile(complete_psyc_anx,33.3))*2;
anx_rank_3 = (complete_psyc_anx <= prctile(complete_psyc_anx,100) & complete_psyc_anx > prctile(complete_psyc_anx,66.6))*3;

anx_rank_4 = complete_psyc_anx <= prctile(complete_psyc_anx,50);
anx_rank_5 = (complete_psyc_anx <= prctile(complete_psyc_anx,100) & complete_psyc_anx > prctile(complete_psyc_anx,50))*2;




group_all_together = array2table(sum([anxiety_group,normal_group,comorbid_group],2),'VariableNames',"3Group_index");
% depression_rank = array2table(sum([dep_rank_1,dep_rank_2,dep_rank_3],2),'VariableNames',"dep_Group_index");
depression_rank = array2table(sum([dep_rank_1,dep_rank_3],2),'VariableNames',"dep_Group_index");

anxiety_rank = array2table(sum([anx_rank_1,anx_rank_2,anx_rank_3],2),'VariableNames',"anx_Group_index");
anxiety_rank2 = array2table(sum([anx_rank_4,anx_rank_5],2),'VariableNames',"anx_Group_index2");

% psychological variable
index_psyc = [13,17]; % 5 13
data_psyc = table2array(data_ori.psyc(:,index_psyc));
data_psyc = arrayfun(@(i) data_psyc(:,i), 1:2, 'uni', 0);
data_psyc_norm = {data_psyc{1}/63, data_psyc{2}/80 };
name_psyc = data_ori.psyc.Properties.VariableNames(index_psyc);
age = data_ori.psyc.age;
sex = data_ori.psyc.gender==1;

% temp = {};
% for i = 1:2
%     temp{i} = data_psyc(:,i);
% end

% memory variables
name_memory = {'full','fore','back','lure','loc'};
name_memory = [name_memory, cellfun(@(x) [x, '_neg'], name_memory, 'uni', 0), ...
               cellfun(@(x) [x, '_neu'], name_memory, 'uni', 0), cellfun(@(x) [x, '_pos'], name_memory, 'uni', 0)];
name_confidence = {'Pic_confi','pic_confi_full','pic_confi_fore','pic_confi_back','pic_confi_lure','loc_confi_correct'};
name_confidence = [name_confidence, cellfun(@(x) [x, '_neg'], name_confidence, 'uni', 0), ...
               cellfun(@(x) [x, '_neu'], name_confidence, 'uni', 0), cellfun(@(x) [x, '_pos'], name_confidence, 'uni', 0)];

name_rt = {'RT_pic','RT_pic_full','RT_pic_fore','RT_pic_back','RT_pic_lure', 'RT_loc_correct'};
name_rt = [name_rt, cellfun(@(x) [x, '_neg'], name_rt, 'uni', 0), ...
               cellfun(@(x) [x, '_neu'], name_rt, 'uni', 0), cellfun(@(x) [x, '_pos'], name_rt, 'uni', 0)];

name_stroop_rt = {'SRT_congruent','SRT_incongruent', 'SRT_difference'};

data_memory = {}; data_confidence = {}; data_rt = {}; loc_memory_check_emo =[];
i = 0; 
for sbj_i = 1 : length(data_ori.memory)
    memory = data_ori.memory{sbj_i};
    memory_stroop = data_ori.memory_stroop{sbj_i};
    response = memory.lure_response;
    location_response = (memory.loc_response == memory.loc_answer);
    confidence_pic_response = memory.confi_response_recog_response;
    confidence_loc_response = memory.confi_response_Loc_response;
    stroop_rt = memory_stroop.stroop_resp_rt;
    
    
    if  sbj_i > noPSSERQ_ind
        response_time_pic = (memory.pic_resp_time);
        response_time_loc = (memory.loc_resp_time);
    end
    
    if arousal_opt == 1
        temp = arousal_controlled_fore_list{sbj_i};
    else
        temp = memory.fore_list_real;
    end

    trial_type = 3 * ones(length(temp),1);
    for trial_i = 1:length(temp)
        if contains(temp{trial_i},'neg')
            trial_type(trial_i) = -1;
        elseif contains(temp{trial_i},'neu')
            trial_type(trial_i) = 0;
        elseif contains(temp{trial_i},'pos')
            trial_type(trial_i) = 1;
        end
    end

    response = response(:);
    trial_type = trial_type(:);

    % Recognition
    temp1 = arrayfun(@(i) sum(response==i)/length(response), 1:4); % All
    temp2 = arrayfun(@(i) sum(response==i & trial_type==-1)/sum(trial_type==-1), 1:4); % Negative
    temp3 = arrayfun(@(i) sum(response==i & trial_type==0)/sum(trial_type==0), 1:4); % Neutral
    temp4 = arrayfun(@(i) sum(response==i & trial_type==1)/sum(trial_type==1), 1:4); % Positive

    % Location
    temp5 = arrayfun(@(i) sum(location_response==i)/length(location_response), 1);
    temp6 = arrayfun(@(i) sum(location_response==i & trial_type==-1)/sum(trial_type==-1), 1);
    temp7 = arrayfun(@(i) sum(location_response==i & trial_type==0)/sum(trial_type==0), 1);
    temp8 = arrayfun(@(i) sum(location_response==i & trial_type==1)/sum(trial_type==1), 1);

%     temp111 = arrayfun(@(i) sum(location_response==1 & response == 1 & trial_type == i)/sum(response==1 & trial_type==i), [-1 0 1]);
%     temp222 = arrayfun(@(i) sum(location_response==1 & response == 2 & trial_type == i)/sum(response==2 & trial_type==i), [-1 0 1]);
%     temp333 = arrayfun(@(i) sum(location_response==1 & response == 3 & trial_type == i)/sum(response==3 & trial_type==i), [-1 0 1]);
%     temp444 = arrayfun(@(i) sum(location_response==1 & response == 4 & trial_type == i)/sum(response==4 & trial_type==i), [-1 0 1]);
    temp111 = arrayfun(@(i) sum(location_response==i & response == 1)/sum(response==1 ), 1);
    temp222 = arrayfun(@(i) sum(location_response==i & response == 2)/sum(response==2), 1);
    temp333 = arrayfun(@(i) sum(location_response==i & response == 3)/sum(response==3), 1);
    temp444 = arrayfun(@(i) sum(location_response==i & response == 4)/sum(response==4), 1);


    loc_1_emo = arrayfun(@(i) sum(location_response==1 & response == 1 & (trial_type==i))/sum(response==1 & (trial_type==i)), [-1 0 1]);
    loc_2_emo = arrayfun(@(i) sum(location_response==1 & response == 2 & (trial_type==i))/sum(response==2 & (trial_type==i)), [-1 0 1]);
    loc_3_emo = arrayfun(@(i) sum(location_response==1 & response == 3 & (trial_type==i))/sum(response==3 & (trial_type==i)), [-1 0 1]);
    loc_4_emo = arrayfun(@(i) sum(location_response==1 & response == 4 & (trial_type==i))/sum(response==4 & (trial_type==i)), [-1 0 1]);

    

    % Confidence : picture
    temp18 = nanmean(confidence_pic_response); % All
    temp9 = arrayfun(@(i) nanmean(confidence_pic_response(response==i)), 1:4); % All by options
    temp17_1 = arrayfun(@(i) nanmean(confidence_pic_response(trial_type == i)), -1 ); % Negative
    temp10 = arrayfun(@(i) nanmean(confidence_pic_response(response==i & trial_type == -1)), 1:4); % Negative by options
    temp17_2 = arrayfun(@(i) nanmean(confidence_pic_response(trial_type == i)), 0 ); % Neutral
    temp11 = arrayfun(@(i) nanmean(confidence_pic_response(response==i & trial_type == 0)), 1:4); % Neutral by options
    temp17_3 = arrayfun(@(i) nanmean(confidence_pic_response(trial_type == i)), 1 ); % Positive
    temp12 = arrayfun(@(i) nanmean(confidence_pic_response(response==i & trial_type == 1)), 1:4); % Positive by options

    % Confidence : Location
    temp13 = arrayfun(@(i) nanmean(confidence_loc_response(location_response==i)), 1);
    temp14 = arrayfun(@(i) nanmean(confidence_loc_response(location_response==i & trial_type == -1)), 1);
    temp15 = arrayfun(@(i) nanmean(confidence_loc_response(location_response == i & trial_type == 0)), 1);
    temp16 = arrayfun(@(i) nanmean(confidence_loc_response(location_response == i & trial_type == 1)), 1);

    % STROOP RT
    stroop_correct = memory_stroop.stroop_resp_corr;
    stroop_congruent = mod(memory_stroop.stroop_thisIndex+1,2);
    stroop_rt = zscore(stroop_rt);
    stroop_rt(isoutlier(stroop_rt)) = nan;

    temp30 = arrayfun(@(i) nanmean(stroop_rt(~isnan(stroop_rt(stroop_correct==i & stroop_congruent==1)))), 1);
    temp31 = arrayfun(@(i) nanmean(stroop_rt(~isnan(stroop_rt(stroop_correct==i & stroop_congruent==0)))), 1);
    temp32 = temp31-temp30;
    


    if  sbj_i > noPSSERQ_ind
        % Response time : picture
        temp20 = mean(response_time_pic);
        temp22 = arrayfun(@(i) nanmean(response_time_pic(response==i)), 1:4 );
        temp21_1 = arrayfun(@(i) nanmean(response_time_pic(trial_type == i)), -1);
        temp23 = arrayfun(@(i) nanmean(response_time_pic(response==i & trial_type == -1)), 1:4);
        temp21_2 = arrayfun(@(i) nanmean(response_time_pic(trial_type == i)), 0);
        temp24 = arrayfun(@(i) nanmean(response_time_pic(response==i & trial_type == 0)), 1:4);
        temp21_3 = arrayfun(@(i) nanmean(response_time_pic(trial_type == i)), 1);
        temp25 = arrayfun(@(i) nanmean(response_time_pic(response==i & trial_type == 1)), 1:4);
        % Response time : location
        temp26 = arrayfun(@(i) nanmean(response_time_loc(location_response == i)), 1);
        temp27 = arrayfun(@(i) nanmean(response_time_loc(location_response == i & trial_type == -1)), 1);
        temp28 = arrayfun(@(i) nanmean(response_time_loc(location_response == i & trial_type == 0)), 1);
        temp29 = arrayfun(@(i) nanmean(response_time_loc(location_response == i & trial_type == 1)), 1);

        data_rt{sbj_i} = [temp20, temp22, temp26, temp21_1,temp23,temp27, temp21_2, temp24, temp28, temp21_3, temp25,temp29];

    end

    data_memory{sbj_i} = [temp1,temp5, temp2,temp6, temp3,temp7,temp4,temp8];
    data_confidence{sbj_i} = [temp18, temp9, temp13, temp17_1, temp10, temp14, temp17_2, temp11, temp15, temp17_3 temp12, temp16];
    data_stroop{sbj_i} = [temp30, temp31, temp32];
    loc_memory_check{sbj_i} = [temp111, temp222];
    loc_memory_check_emo{sbj_i} = [loc_1_emo, loc_2_emo,loc_3_emo,loc_4_emo];
    loc_memory_check_emo_FONLY{sbj_i} = [loc_1_emo, loc_2_emo];

end

% metric variables
name_metric = {'fore-full(all)','fore-full(neg)','fore-full(neu)','fore-full(pos)', ...
               'fore-back(all)','fore-back(neg)','fore-back(neu)','fore-back(pos)', ...
               'fore-lure(all)','fore-lure(neg)','fore-lure(neu)','fore-lure(pos)', ...
               'back-lure(all)','back-lure(neg)','back-lure(neu)','back-lure(pos)', ...
               'fore/fore+back(all)','fore/fore+back(neg)','fore/fore+back(neu)','fore/fore+back(pos)', ...
               'back/fore+back(all)','back/fore+back(neg)','back/fore+back(neu)','back/fore+back(pos)', ...
               ...
               'fore-full(neg-neu)','fore-full(pos-neu)' ...
               'fore-back(neg-neu)','fore-back(pos-neu)' ...
               'fore-lure(neg-neu)','fore-lure(pos-neu)' ...
               'fore-back-lure(neg-neu)','fore-back-lure(pos-neu)', ...
               ...
               'full(neg-neu)', 'full(pos-neu)', ...
               'fore(neg-neu)', 'fore(pos-neu)', ...
               'back(neg-neu)', 'back(pos-neu)', ...
               'lure(neg-neu)', 'lure(pos-neu)',...
               'Location(neg-neu)', 'Location(pos-neu)'};
name_metric_rt = {'rt neg-neu','rt pos-neu',...
    'rt full neg-neu','rt fore neg-neu','rt back neg-neu','rt lure neg-neu',...
    'rt full pos-neu','rt fore pos-neu','rt back pos-neu','rt lure pos-neu'};


data_metric = {};
i = 0;
for sbj_i = 1:length(data_memory)
    
    temp_sbj = []; temp_rt_sbj = [];

    % fore-full
    offset = 5;
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(1);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(1+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(1+offset*2);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*3) - data_memory{sbj_i}(1+offset*3);

    % fore-back
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(3);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(3+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(3+offset*2);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*3) - data_memory{sbj_i}(3+offset*3);

    % fore-lure
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(4);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(4+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(4+offset*2);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*3) - data_memory{sbj_i}(4+offset*3);

    % back-lure
    temp_sbj(end+1) = data_memory{sbj_i}(3) - data_memory{sbj_i}(4);
    temp_sbj(end+1) = data_memory{sbj_i}(3+offset) - data_memory{sbj_i}(4+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(3+offset*2) - data_memory{sbj_i}(4+offset*2);
    temp_sbj(end+1) = data_memory{sbj_i}(3+offset*3) - data_memory{sbj_i}(4+offset*3);

    % fore/fore+back
    temp_sbj(end+1) = data_memory{sbj_i}(2) / (data_memory{sbj_i}(2) + data_memory{sbj_i}(3));
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) / (data_memory{sbj_i}(2+offset) + data_memory{sbj_i}(3+offset));
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) / (data_memory{sbj_i}(2+offset*2) + data_memory{sbj_i}(3+offset*2));
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*3) / (data_memory{sbj_i}(2+offset*3) + data_memory{sbj_i}(3+offset*3));
%     temp_sbj(end+1) = mean([data_memory{sbj_i}(2+offset*2) / (data_memory{sbj_i}(2+offset*2) + data_memory{sbj_i}(3+offset*2)), data_memory{sbj_i}(2+offset*3) / (data_memory{sbj_i}(2+offset*3) + data_memory{sbj_i}(3+offset*3))]);
        

    % back/fore+back
    temp_sbj(end+1) = data_memory{sbj_i}(3) / (data_memory{sbj_i}(2) + data_memory{sbj_i}(3));
    temp_sbj(end+1) = data_memory{sbj_i}(3+offset) / (data_memory{sbj_i}(2+offset) + data_memory{sbj_i}(3+offset));
    temp_sbj(end+1) = data_memory{sbj_i}(3+offset*2) / (data_memory{sbj_i}(2+offset*2) + data_memory{sbj_i}(3+offset*2));
    temp_sbj(end+1) = data_memory{sbj_i}(3+offset*3) / (data_memory{sbj_i}(2+offset*3) + data_memory{sbj_i}(3+offset*3));

    % fore-full (vs. neurtral)
    offset= 6;
    temp_sbj(end+1) = temp_sbj(2) - temp_sbj(3);
    temp_sbj(end+1) = temp_sbj(4) - temp_sbj(3);

    % fore-back (vs. neurtral)
    temp_sbj(end+1) = temp_sbj(2+offset) - temp_sbj(3+offset);
    temp_sbj(end+1) = temp_sbj(4+offset) - temp_sbj(3+offset);

    % fore-lure (vs. neurtral)    
    temp_sbj(end+1) = temp_sbj(2+offset*2) - temp_sbj(3+offset*2);
    temp_sbj(end+1) = temp_sbj(4+offset*2) - temp_sbj(3+offset*2);

    % fore-back-lure (vs. neurtral)
    temp_sbj(end+1) = temp_sbj(2+offset*3) - temp_sbj(3+offset*3);
    temp_sbj(end+1) = temp_sbj(4+offset*3) - temp_sbj(3+offset*3);

    % full (vs. neutral)
    temp_sbj(end+1) = data_memory{sbj_i}(6) - data_memory{sbj_i}(11);
    temp_sbj(end+1) = data_memory{sbj_i}(16) - data_memory{sbj_i}(11);

    % fore (vs. neutral)
    temp_sbj(end+1) = data_memory{sbj_i}(6+1) - data_memory{sbj_i}(11+1);
    temp_sbj(end+1) = data_memory{sbj_i}(16+1) - data_memory{sbj_i}(11+1);

    % back (vs. neutral)
    temp_sbj(end+1) = data_memory{sbj_i}(6+2) - data_memory{sbj_i}(11+2);
    temp_sbj(end+1) = data_memory{sbj_i}(16+2) - data_memory{sbj_i}(11+2);

    % lure (vs. neutral)
    temp_sbj(end+1) = data_memory{sbj_i}(6+3) - data_memory{sbj_i}(11+3);
    temp_sbj(end+1) = data_memory{sbj_i}(16+3) - data_memory{sbj_i}(11+3);

    % location (vs. neutral)
    temp_sbj(end+1) = data_memory{sbj_i}(6+4) - data_memory{sbj_i}(11+4);
    temp_sbj(end+1) = data_memory{sbj_i}(16+4) - data_memory{sbj_i}(11+4);


    % RESPONSE TIME METRIC
%     temp_rt_sbj(end+1) = 
    data_metric{sbj_i} = temp_sbj;
%     data_metric_confidence{sbj_i} = temp_sbj_confidence;
end

% post processing
target = data_memory;
data_memory = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

target = data_confidence;
data_confidence = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

target = loc_memory_check;
loc_memory = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

target = loc_memory_check_emo_FONLY;
loc_memory_emo = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

target = data_rt;
target(1:noPSSERQ_ind) = [];
data_rt = arrayfun(@(i) cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);
data_rt_none = zeros(1,noPSSERQ_ind);
data_rt_none(:) = nan;
data_rt = arrayfun(@(i) [data_rt_none data_rt{i}], 1:length(data_rt), 'UniformOutput', false);

target = data_metric;
data_metric = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

target = data_stroop;
data_stroop = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

data_psyc
name_psyc

data_memory
name_memory

data_metric
name_metric

data_confidence
name_confidence

data_rt
name_rt

data_stroop
name_stroop_rt

fprintf("\n Mean age = %.3f \n Age std = %.3f",mean(age), std(age))
length(age)
sum(sex)
% Make the table for bar plot and anova
psyc = data_ori.psyc;
if group_opt ~= 0
    data_memory = cellfun(@(x) x(group_ind), data_memory,'UniformOutput',0);
    data_confidence = cellfun(@(x) x(group_ind), data_confidence,'UniformOutput',0);
    data_rt = cellfun(@(x) x(group_ind), data_rt,'UniformOutput',0);
    data_metric = cellfun(@(x) x(group_ind), data_metric,'UniformOutput',0);
    data_stroop = cellfun(@(x) x(group_ind), data_stroop,'UniformOutput',0);
    psyc = psyc(group_ind,:);
end
% cell array를 matrix로 변환합니다.
matrix_memory = cell2mat(cellfun(@(x) x', data_memory, 'UniformOutput', false));
matrix_confidence = cell2mat(cellfun(@(x) x', data_confidence, 'UniformOutput', false));
matrix_rt = cell2mat(cellfun(@(x) x', data_rt, 'UniformOutput', false));
matrix_metric = cell2mat(cellfun(@(x) x', data_metric, 'UniformOutput', false));
matrix_stroop = cell2mat(cellfun(@(x) x', data_stroop, 'UniformOutput', false));

% matrix를 table로 변환합니다.
table_memory = array2table(matrix_memory, 'VariableNames', name_memory);
table_rt = array2table(matrix_rt, 'VariableNames', name_rt);
table_confidence = array2table(matrix_confidence, 'VariableNames', name_confidence);
table_metric = array2table(matrix_metric, 'VariableNames', name_metric);
table_stroop = array2table(matrix_stroop, 'VariableNames', name_stroop_rt);


if group_opt == 0 
    table_psyc_memory_metric = horzcat(psyc,table_memory,table_metric, table_confidence, table_rt,table_stroop,anxiety_rank,anxiety_rank2,depression_rank,group_all_together,val_table,aro_table);
else
    table_psyc_memory_metric = horzcat(psyc,table_memory,table_metric, table_confidence, table_rt,table_stroop,val_table,aro_table);
end

size(table_psyc_memory_metric)

writetable(table_psyc_memory_metric,[save_path '2306_behavior_table.csv'])


figure_JK()
ms_bar(loc_memory,[.7,.7,.7;.5,.5,.5;.5,.5,.5;.5,.5,.5;])
xticklabels({"Full hit", "Foreground hit"})
title("Correct location rate per options")
ylabel("Portion (%)")
[h ,p, ~, stats] = ttest(loc_memory{1},loc_memory{2})


figure_JK()
ms_bar(loc_memory_emo,[.7,.7,.7;.5,.5,.5;.5,.5,.5;.5,.5,.5;])
xticklabels({"C(Neg)", "FH(Neg)","C(NEU)", "FH(NEU)","C(POS)", "FH(POS)"})
title("Correct location rate per options")
ylabel("Portion (%)")
[h p] = ttest(loc_memory{1},loc_memory{2})

%% partial correlation

aux_func_get_residual = @(mdl) mdl.Residuals.Standardized;
func_get_residual = @(x,y) aux_func_get_residual(fitlm(x(:),y(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psyc_i = 2;
data_opt = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if data_opt == 1
    data_target = data_memory; name_target = name_memory; target_metric_list = 1:20;   sub_row_n = 4; sub_col_n=5; fig_size = [2000 1300];
    % data_target = data_memory; name_target = name_memory; target_metric_list = 1:20;   sub_row_n = 2; sub_col_n=4; fig_size = [1581,576];
else
    %     data_target = data_metric; name_target = name_metric;  target_metric_list = 1:24;   sub_row_n = 6; sub_col_n=4; fig_size = [1400 1100];
    data_target = data_metric; name_target = name_metric;  target_metric_list = 17:24;   sub_row_n = 2; sub_col_n=4; fig_size = [1581,576];

end
% data_target = data_metric; name_target = name_metric;  target_metric_list = 17:20;  sub_row_n = 1; sub_col_n=4; fig_size = [1740 270];
% data_target = data_metric; name_target = name_metric;  target_metric_list = 25:34;  sub_row_n = 5; sub_col_n=2; fig_size = [700 967];

opt_partial_corr = 3;

corr_type = 'pearson';
% corr_type = 'spearman';
% corr_type = 'kendall';


flag_list = {};
flag = true(1,length(data_target{1})); flag_list{end+1} = flag;
% flag = data_memory{4}<median(data_memory{4});
% flag = data_memory{4}>median(data_memory{4});
% flag = sex==1; flag_list{end+1} = flag;
% flag = sex==0; flag_list{end+1} = flag;
% flag = data_psyc{psyc_i} > median(data_psyc{psyc_i}); flag_list{end+1} = flag;
% flag = data_psyc{psyc_i} < median(data_psyc{psyc_i}); flag_list{end+1} = flag;
% flag_name_list = {'all','male','female','high anxiety','low anxiety'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for psyc_i = 1:2
% for flag_i = 1:length(flag_list)
%     flag = flag_list{flag_i};

figure('position',[50 50 fig_size]); set(gcf,"Color",[1 1 1])
for metric_i = 1:length(target_metric_list)
    
    if opt_partial_corr == 1 % no partial correlation
%         data1 = data_psyc{psyc_i}(flag);
        data1 = data_psyc{psyc_i}(flag);
        data2 = data_target{target_metric_list(metric_i)}(flag);
    elseif opt_partial_corr == 2 % psyc variable only partial correlation
%         data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag), data_psyc{psyc_i}(flag));
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag), data_psyc{psyc_i}(flag));
        data2 = data_target{target_metric_list(metric_i)}(flag);
    elseif opt_partial_corr == 3 % psyc + memory all partial correlation
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag) , data_psyc{psyc_i}(flag));
        data2 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag), data_target{target_metric_list(metric_i)}(flag));
    elseif opt_partial_corr == 4 % psyc variable only partial correlation (from ORIGINAL DATA)
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}, data_psyc{psyc_i});
        data1 = data1(flag);
        data2 = data_target{target_metric_list(metric_i)}(flag);
    elseif opt_partial_corr == 5 % psyc + memory all partial correlation (from ORIGINAL DATA)
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}, data_psyc{psyc_i});
        data2 = func_get_residual(data_psyc{abs(psyc_i-2)+1}, data_target{target_metric_list(metric_i)});
        data1 = data1(flag);
        data2 = data2(flag);
    end


    if opt_partial_corr == 1
        label1 = [name_psyc{psyc_i}] ;
    else
        label1 = [name_psyc{psyc_i} ' residual'] ;
    end
        label2 = name_target{target_metric_list(metric_i)};

    subplot(sub_row_n, sub_col_n, metric_i)
%     data2(data2 == 1) = nan;
%     data1(data2 == 1) = nan;

    [r,p] = jh_regress(data1, data2','on','type',corr_type);
    xlabel(label1,'interpreter','none');
    ylabel(label2,'interpreter','none');

    if psyc_i == 1
        xlim([-2.2 3])
    else
        xlim([-2.2 2.5])
    end
    if p < 0.001
        title(sprintf('r = %.3f, p < 0.001(***)', r),'color','k')
    elseif p < 0.01
        title(sprintf('r = %.3f, p < 0.01(**)', r),'color','k')
    elseif p < 0.05
        title(sprintf('r = %.3f, p < 0.05(*)', r),'color','k')
    else
        title(sprintf('r = %.3f, p = %.3f', r, p))
    end

end

% sgtitle(flag_name_list{flag_i});

% end
% end

%%  BAR GRAPH
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag_list = {};
flag = true(1,length(data_target{1})); flag_list{end+1} = flag;
% flag = sex==1; flag_list{end+1} = flag;
% flag = sex==0; flag_list{end+1} = flag;

psyc_i = 1;  % 1 phq 2 stai
memory_metric_opt = 2;

% flag = data_psyc{psyc_i} > median(data_psyc{psyc_i}); flag_list{end+1} = flag; % High group
% flag = data_psyc{psyc_i} < median(data_psyc{psyc_i}); flag_list{end+1} = flag; % Low group
% flag_name_list = {'all','male','female','high anxiety','low anxiety'};

if memory_metric_opt == 2
    fig_position = [2,42,2555,881];
else
    fig_position = [2,42,2555,769];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idx = [2,3,4,6,7,8,10,11,12,14,15,16];
offset = 1; 
if memory_metric_opt == 1
    idx = {[1 2 3 4],[6 11 16] [6+offset 11+offset 16+offset],[6+(offset*2) 11+(offset*2) 16+(offset*2)],[6+(offset*3) 11+(offset*3) 16+(offset*3)],...
        [6+(offset*4) 11+(offset*4) 16+(offset*4)]};
else
    idx = {[17 21],[18 19 20] [22 23 24],[41 42]};
end
% title_name = {"Memory response rate","Full memory","Foreground hit","Background hit","Lure","Location"};


color_option = [.6 .6 .6;.6 .6 .6;.6 .6 .6;.6 .6 .6;];
for bar_idx = 1 : length(idx)
    if memory_metric_opt == 1
        data = data_memory(idx{bar_idx});
        label = name_memory(idx{bar_idx});
    else
        data = data_metric(idx{bar_idx});
        label = name_metric(idx{bar_idx});
    end

data = cellfun(@(x) x(flag), data, 'uni', 0);

% [err,avg] = ms_mean_err(data);


if bar_idx == 1
    figure_JK(1)
    ms_bar(data,color_option);
elseif bar_idx == 2
    figure_JK(2,fig_position) ;
    subplot(2,4, bar_idx-1)
    ms_bar(data)
    if memory_metric_opt == 1
        ylim([0 .6])
    else
        ylim([.6 .9])
    end
elseif bar_idx == 6
    figure_JK(2,fig_position) ;
    subplot(2,4, bar_idx-1)
    ms_bar(data)
    ylim([.3 .4])
else
    figure_JK(2,fig_position) ;
    subplot(2,4, bar_idx-1)
    ms_bar(data)
    ylim([0 .6])
    figure_JK(2,fig_position)
end
% figure;
hold on
xticks(1:length(idx))
xticklabels(cellfun(@(x) replace(x,'_',' '),label,'uni',0))
% title(sprintf("%s", title_name{bar_idx}))

ylabel("Response rate(%)")
xtickangle(30)

p = [];
for i = 1:length(data)
    for j = 1:length(data)
        [~,p(i,j)] = ttest(data{i},data{j});
    end
end
p(p>0.05) = nan;
figure
heatmap(p);
end

%% Anxiety Depression distribution
data1 = data_psyc{1};
data2 = data_psyc{2};
corr_type = "pearson";
figure_JK;
[r,p] = jh_regress(data2, data1,'on','type', corr_type); hold on
title(sprintf('Anxiety & Depression Distribution \n r = %.3f, p = %.3f', r, p))
% title(sprintf('Anxiety & Depression Distribution'))

xline(45, '--k','linewidth',2)
yline(20, '--k','linewidth',2)
ylabel('BDI-II')
xlabel('STAI-T')



hold on
data1 = data_psyc{1};
data2 = data_psyc{2};
% [r,p] = jh_regress(data1, data2, corr_type);
% title(sprintf('r = %.3f, p = %.3f', r, p))

% % scatter(data1,data2(sex==1),'filled','b')
% scatter(data2(sex==1),data1(sex==1),'filled','b'); hold on
% scatter(data2(sex==0),data1(sex==0),'filled','r')
% scatter(data1,data2,'filled','b')
xline(45, '--k','linewidth',1.5)
yline(20, '--k','linewidth',1.5)
% xline(median(data1), '--k','linewidth',1.5)
% yline(median(data2), '--k','linewidth',1.5)

figure_JK()
histogram(data1,20)
xlabel("BDI-II distribution")
ylabel("Count")
xline(19, '--k','linewidth',2)
set(gca, "Linewidth", 2, 'fontsize', 18)
title("BDI-II distribution")
box off


figure_JK()
histogram(data2,25)
xlabel("STAI-trait distribution")
ylabel("Count")
xline(45, '--k','linewidth',2)
set(gca, "Linewidth", 2, 'fontsize', 18)
title("STAI-trait distribution")
box off





%% clustering based approach

foreground_bias_performance = [table_psyc_memory_metric.("fore/fore+back(pos)"),table_psyc_memory_metric.("fore/fore+back(neu)")];
outlier1 = find(isoutlier(foreground_bias_performance(:,1)));
outlier2 = find(isoutlier(foreground_bias_performance(:,2)));
% outlier_ind = unique([outlier1; outlier2]);

% foreground_bias_performance(outlier_ind,:) = nan;
[idx , centeroids] = kmeans(foreground_bias_performance,3,'Start', 'plus');

numIterations = 10000;  K = 3;
data = foreground_bias_performance;
minWithinClusterVariance = Inf;
bestClusterResult = [];

for i = 1:numIterations
    [idx, centroids] = kmeans(data, K,'Start', 'plus');
    D = pdist2(data, centroids);
    withinClusterVariance = sum((idx - 1) .* diag(D)' .^ 2);
    
    if withinClusterVariance < minWithinClusterVariance
        minWithinClusterVariance = withinClusterVariance;
        bestClusterResult = idx;
        bestcenter = centroids;
    end
end

% 최적의 클러스터링 결과를 bestClusterResult에 저장

cluster1 = foreground_bias_performance(bestClusterResult == 1, :);
cluster2 = foreground_bias_performance(bestClusterResult == 2, :);
cluster3 = foreground_bias_performance(bestClusterResult == 3, :);

% Scatter plot으로 클러스터 표시
figure_JK(1)
scatter(cluster1(:, 1), cluster1(:, 2), 'r', 'filled');
hold on;
scatter(cluster2(:, 1), cluster2(:, 2), 'g', 'filled');
scatter(cluster3(:, 1), cluster3(:, 2), 'b', 'filled');
hold off;

% 클러스터 중심점 표시
hold on;
scatter(bestcenter(:, 1), bestcenter(:, 2), 100, 'k', 'filled', 'd');
hold off;
xlabel("neg bias"); ylabel("neu bias")
title('K-Means 클러스터링 결과');
legend('1', '2', '3', 'center');


bias_rank_1 = bestClusterResult == 2;
bias_rank_2 = bestClusterResult == 3;
bias_rank_3 = bestClusterResult == 1;

bias_group_ind = array2table(bestClusterResult,'VariableNames',"bias_group_ind")

%%
group_opt = 4;

  if group_opt == 1
        group_ind = normal_group;
    elseif group_opt == 2
        group_ind = anxiety_group;
    elseif group_opt == 3
        group_ind = comorbid_group;
    elseif group_opt == 4
        group_ind = dep_rank_1;
    elseif group_opt == 5
        group_ind = dep_rank_2;
    elseif group_opt == 6
        group_ind = dep_rank_3;
    elseif group_opt == 7
        group_ind = anx_rank_1;
    elseif group_opt == 8
        group_ind = anx_rank_2;
    elseif group_opt == 9
        group_ind = anx_rank_3;
  elseif group_opt == 10
      group_ind = bias_rank_1;
  elseif group_opt == 11
      group_ind = bias_rank_2;
  elseif group_opt == 12
      group_ind = bias_rank_3;
  end


%   table_psyc_memory_metric = horzcat(table_psyc_memory_metric,bias_group_ind);

  table_psyc_memory_metric2 = table_psyc_memory_metric(group_ind,:);

  writetable(table_psyc_memory_metric,[save_path '2306_behavior_table.csv'])

% 
% 
%         dep_res = func_get_residual(data_psyc{abs(1-2)+1}(flag), data_psyc{1}(flag));
%         anx_res = func_get_residual(data_psyc{abs(2-2)+1}(flag), data_psyc{2}(flag));
% 
% [h p] = ranksum(dep_res(bias_rank_1), dep_res(bias_rank_3))
% 
























%% AROUSAL BASED ANALYSIS ONCE

%%%%%%%%%%%%%% 

%% Memory - Arousal version

    data_ori = complete_data_fixed;

% psychological variable
index_psyc = [5,17];
data_psyc = table2array(data_ori.psyc(:,index_psyc));
data_psyc = arrayfun(@(i) data_psyc(:,i), 1:2, 'uni', 0);
name_psyc = data_ori.psyc.Properties.VariableNames(index_psyc);
age = data_ori.psyc.age;
sex = data_ori.psyc.gender==1;

% temp = {};
% for i = 1:2
%     temp{i} = data_psyc(:,i);
% end

% memory variables
name_memory = {'full','fore','back','lure','loc'};
name_memory = [name_memory, cellfun(@(x) [x, '_High'], name_memory, 'uni', 0), ...
               cellfun(@(x) [x, '_Low'], name_memory, 'uni', 0)];

data_memory = {};
for sbj_i = 1:length(data_ori.memory)
    memory = data_ori.memory{sbj_i};
    response = memory.lure_response;
    location_response = (memory.loc_response == memory.loc_answer);

    temp = memory.arousal_resp_response;
    temp_sort = sort(temp,"descend");

    trial_type = 2 * ones(length(temp),1);
    for trial_i = 1:length(temp)
        if temp(trial_i) > temp_sort(round(0.5 * numel(temp_sort)))
            trial_type(trial_i) = 1;
        elseif temp(trial_i) < temp_sort(round(0.5 * numel(temp_sort)))
            trial_type(trial_i) = -1;
        end
    end


    response = response(:);
    trial_type = trial_type(:);

    % recognition
    temp1 = arrayfun(@(i) sum(response==i)/length(response), 1:4); % All
    temp2 = arrayfun(@(i) sum(response==i & trial_type==1)/sum(trial_type==1), 1:4); % HIGH AROUSAL
    temp3 = arrayfun(@(i) sum(response==i & trial_type==-1)/sum(trial_type==-1), 1:4); % LOW AROUSAL
        % Location
    temp5 = arrayfun(@(i) sum(location_response==i)/length(location_response), 1);
    temp6 = arrayfun(@(i) sum(location_response==i & trial_type==1)/sum(trial_type==1), 1); % HIGH AROUSAL
    temp7 = arrayfun(@(i) sum(location_response==i & trial_type==-1)/sum(trial_type==-1), 1); % LOW AROUSAL
    

    data_memory{sbj_i} = [temp1,temp5,temp2,temp6,temp3,temp7];
end

% metric variables
name_metric = {'fore-full(all)','fore-full(high)','fore-full(low)', ...
               'fore-back(all)','fore-back(high)','fore-back(low)', ...
               'fore-lure(all)','fore-lure(high)','fore-lure(low)',...
               'fore-back-lure(all)','fore-back-lure(high)','fore-back-lure(low)', ...
               ...
               'fore-full(high-low)', ...
               'fore-back(high-low)',...
               'fore-lure(high-low)',...
               'fore-back-lure(high-low)', ...
               ...
               'full(high-low)',  ...
               'fore(high-low)',  ...
               'back(high-low)',  ...
               'lure(high-low)', ...
               'Location(high-low)'};

data_metric = {};
for sbj_i = 1:length(data_memory)
    temp_sbj = [];

    % fore-full
    offset = 5;
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(1);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(1+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(1+offset*2);
    
    % fore-back
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(3);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(3+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(3+offset*2);
    
    % fore-lure
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(4);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(4+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(4+offset*2);
    
    % fore-back-lure
    temp_sbj(end+1) = data_memory{sbj_i}(2) - data_memory{sbj_i}(3) - data_memory{sbj_i}(4);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset) - data_memory{sbj_i}(3+offset) - data_memory{sbj_i}(4+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(2+offset*2) - data_memory{sbj_i}(3+offset*2) - data_memory{sbj_i}(4+offset*2);

    % fore-back-lure (high - low)
    offset = 2;
    temp_sbj(end+1) = temp_sbj(2) - temp_sbj(3);
    temp_sbj(end+1) = temp_sbj(2+offset) - temp_sbj(3+offset);
    temp_sbj(end+1) = temp_sbj(2+offset*2) - temp_sbj(3+offset*2);
    temp_sbj(end+1) = temp_sbj(2+offset*3) - temp_sbj(3+offset*3);


    % fore-back-lure
    offset = 1;
    temp_sbj(end+1) = data_memory{sbj_i}(6) - data_memory{sbj_i}(11);
    temp_sbj(end+1) = data_memory{sbj_i}(6+offset) - data_memory{sbj_i}(11+offset);
    temp_sbj(end+1) = data_memory{sbj_i}(6+offset*2) - data_memory{sbj_i}(11+offset*2);
    temp_sbj(end+1) = data_memory{sbj_i}(6+offset*3) - data_memory{sbj_i}(11+offset*3);
    
   % location
   offset = 5;
   temp_sbj(end+1) = data_memory{sbj_i}(offset *2) - data_memory{sbj_i}(offset *3);

    data_metric{sbj_i} = temp_sbj;
end

% post processing
target = data_memory;
data_memory = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

target = data_metric;
data_metric = arrayfun(@(i)  cellfun(@(x) x(i), target)  , 1:length(target{1}), 'uni',0);

data_psyc
name_psyc

data_memory
name_memory

data_metric
name_metric


% Make the table for bar plot and anova

% cell array를 matrix로 변환합니다.
matrix_memory = cell2mat(cellfun(@(x) x', data_memory, 'UniformOutput', false));
matrix_metric = cell2mat(cellfun(@(x) x', data_metric, 'UniformOutput', false));

% matrix를 table로 변환합니다.
table_memory = array2table(matrix_memory, 'VariableNames', name_memory);
table_metric = array2table(matrix_metric, 'VariableNames', name_metric);

table_psyc_memory_metric = horzcat(data_ori.psyc,table_memory,table_metric);

writetable(table_psyc_memory_metric,[save_path '2306_behavior_table.csv'])


%% partial correlation : Arousal

aux_func_get_residual = @(mdl) mdl.Residuals.Standardized;
func_get_residual = @(x,y) aux_func_get_residual(fitlm(x(:),y(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psyc_i = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_target = data_memory; name_target = name_memory; target_metric_list = 1:15;   sub_row_n = 3; sub_col_n=5; fig_size = [1400 1200];
% data_target = data_metric; name_target = name_metric;  target_metric_list = 1:12;   sub_row_n = 4; sub_col_n=3; fig_size = [1400 1100];
% data_target = data_metric; name_target = name_metric;  target_metric_list = 17:24;  sub_row_n = 4; sub_col_n=2; fig_size = [700 967];
% data_target = data_metric; name_target = name_metric;  target_metric_list = 25:34;  sub_row_n = 5; sub_col_n=2; fig_size = [700 967];

opt_partial_corr = 2;

corr_type = 'pearson';
% corr_type = 'spearman';
% corr_type = 'kendall';


flag_list = {};
flag = true(1,length(data_target{1})); flag_list{end+1} = flag;
% flag = data_memory{4}<median(data_memory{4});
% flag = data_memory{4}>median(data_memory{4});
% flag = sex==1; flag_list{end+1} = flag;
% flag = sex==0; flag_list{end+1} = flag;
% flag = data_psyc{psyc_i} > median(data_psyc{psyc_i}); flag_list{end+1} = flag;
% flag = data_psyc{psyc_i} < median(data_psyc{psyc_i}); flag_list{end+1} = flag;
% flag_name_list = {'all','male','female','high anxiety','low anxiety'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for psyc_i = 1:2
% for flag_i = 1:length(flag_list)
%     flag = flag_list{flag_i};

figure('position',[50 50 fig_size]); set(gcf,"Color",[1 1 1])
for metric_i = 1:length(target_metric_list)
    
    if opt_partial_corr == 1 % no partial correlation
        data1 = data_psyc{psyc_i}(flag);
        data2 = data_target{target_metric_list(metric_i)}(flag);
    elseif opt_partial_corr == 2 % psyc variable only partial correlation
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag), data_psyc{psyc_i}(flag));
        data2 = data_target{target_metric_list(metric_i)}(flag);
    elseif opt_partial_corr == 3 % psyc + memory all partial correlation
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag) , data_psyc{psyc_i}(flag));
        data2 = func_get_residual(data_psyc{abs(psyc_i-2)+1}(flag), data_target{target_metric_list(metric_i)}(flag));
    elseif opt_partial_corr == 4 % psyc variable only partial correlation (from ORIGINAL DATA)
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}, data_psyc{psyc_i});
        data1 = data1(flag);
        data2 = data_target{target_metric_list(metric_i)}(flag);
    elseif opt_partial_corr == 5 % psyc + memory all partial correlation (from ORIGINAL DATA)
        data1 = func_get_residual(data_psyc{abs(psyc_i-2)+1}, data_psyc{psyc_i});
        data2 = func_get_residual(data_psyc{abs(psyc_i-2)+1}, data_target{target_metric_list(metric_i)});
        data1 = data1(flag);
        data2 = data2(flag);
    end


    label1 = [name_psyc{psyc_i} ' residual'] ;
    label2 = name_target{target_metric_list(metric_i)};

    subplot(sub_row_n, sub_col_n, metric_i)
    [r,p] = jh_regress(data1, data2, corr_type);
    xlabel(label1,'interpreter','none');
    ylabel(label2,'interpreter','none');

    if p < 0.05
        title(sprintf('r = %.3f, p = %.3f', r, p),'color','r')
    elseif p < 0.1
        title(sprintf('r = %.3f, p = %.3f', r, p),'color','b')
    else
        title(sprintf('r = %.3f, p = %.3f', r, p))
    end

end