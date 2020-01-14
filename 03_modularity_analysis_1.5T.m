clear all
close all

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/'));
P='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/';
RDIR='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/02_modularity/'; 
MATDIR='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/01_matfiles/';
aal_data  = SurfStatReadData('aal_both_rsl_final.txt');

save_file = 0;
printfigs = 0; % don't flag on unless you really want to save the figures as .png, otherwise your precious time is gone...
lesion_exclusion = 0;

%% define colormaps
for define_colormaps = 1
    
    blackblue = [zeros(1,3)*0.8;
        zeros(127,1)   (0:126)'/127   ones(127,1)];
    blackblue = flipud(blackblue);
    blue      = [ones(1,3)*0.8;
        zeros(127,1)   (0:126)'/127   ones(127,1)];
    blue      = flipud(blue);
    
    red       = [ones(1,3)*0.8; ...
        ones(64,1) linspace(0,253,64)'/254 zeros(64,1);...
        ones(64,1) ones(64,1) linspace(0,253,64)'/254];
    
end

%% load stereotaxic surfaces
for load_surfaces = 1
    
    S = SurfStatAvSurf({[ P 'surf_reg_model_left.obj' ], [ P 'surf_reg_model_right.obj' ]});
    S_left  = SurfStatReadSurf([ P 'surf_reg_model_left.obj' ]);
    S_right = SurfStatReadSurf([ P 'surf_reg_model_right.obj' ]);
    
    mask = SurfStatMaskCut(S);
    mask_l_r = SurfStatMaskCut(S_left);
    Sinf = SurfStatInflate(S, 0.75);
    
end

%% read controls demographics
for read_controls = 1
    
    fid = fopen('DemographicNGrouping_for_controls_Data.csv');
    C = textscan(fid, '%s%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    Codes_cont                      = C{1}(:, 1);    
    Path_cont                       = C{1}(:, 2);
    Age_cont                        = C{2};                 % mean/sd: 30.9/11.0
    Sex_cont                        = C{3}(:, 1);           % male/female: 16/25    
    fclose(fid);
    
end

%% read FCD Type-I demographics
for read_FCD_Type_I = 1
    
    fid = fopen('DemographicNGrouping_for_FCD_Type-I_Data.csv');
    C = textscan(fid, '%s%f%s%d%d%d%d%s%d%d','Delimiter',',','CollectOutput', 1);
    Codes_FCD_Type_I                = C{1};
    Age_FCD_Type_I                  = C{2};                 % mean/sd: 29.2/8.7
    Sex_FCD_Type_I                  = C{3};                 % male/female: 7/6
    Left_FCD_Type_I                 = C{4}(:, 1);           % 7
    Right_FCD_Type_I                = C{4}(:, 2);           % 6
    Pre_Central_FCD_Type_I          = C{4}(:, 3);           % 10
    Post_Central_FCD_Type_I         = C{4}(:, 4);           % 3
    Path_FCD_Type_I                 = C{5};    
    Duration_FCD_Type_I             = double(C{6}(:, 1));   % mean/sd: 17.1/10.2
    fclose(fid);
    
end

%% read FCD Type-II demographics
for read_FCD_Type_II = 1
    
    fid = fopen('DemographicNGrouping_for_FCD_Type-II_Data.csv');
    C = textscan(fid, '%s%f%s%d%d%d%d%s%d%d%d%d','Delimiter',',','CollectOutput', 1);
    Codes_FCD_Type_II               = C{1};
    Age_FCD_Type_II                 = C{2};                 % mean/sd: 27.8/10.1
    Sex_FCD_Type_II                 = C{3};                 % male/female: 16/14        
    Left_FCD_Type_II                = C{4}(:, 1);           % 15
    Right_FCD_Type_II               = C{4}(:, 2);           % 15
    Pre_Central_FCD_Type_II         = C{4}(:, 3);           % 16
    Post_Central_FCD_Type_II        = C{4}(:, 4);           % 14    
    Path_FCD_Type_II                = C{5};
    Duration_FCD_Type_II            = double(C{6}(:, 1));   % mean/sd: 20.9/11.9
    Operation_FCD_Type_II           = double(C{6}(:, 2));   % 30
    Scanner_FCD_Type_II             = double(C{6}(:, 3));
    Lesion_Volume                   = double(C{6}(:, 4));    
    
    LabelPath_FCD = '/host/gypsy/local_raid//seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label/';
    fclose(fid); 
    
end

%% read HET demographics
for read_HET = 1
    
    fid = fopen('DemographicNGrouping_for_HET_Data.csv');
    C = textscan(fid, '%s%f%s%d%d%s%s%d%s','Delimiter',',','CollectOutput', 1);
    Inclusion                       = C{4}(:, 1);
    Codes_HET                       = C{1}(Inclusion>0);
    Age_HET                         = C{2}(Inclusion>0);    % mean/sd: 28.1/11.5
    Sex_HET                         = C{3}(Inclusion>0);    % male/female: 11/16
    Left_HET                        = C{4}(Inclusion>0, 2) == 1; 
    Right_HET                       = C{4}(Inclusion>0, 2) == 2;
    Bilateral_HET                   = C{4}(Inclusion>0, 2) == 3;  % Left/Right/Bilateral: 6/10/11  
    Path_HET                        = C{5}(Inclusion>0, 1);
    Subtype_HET                     = C{5}(Inclusion>0, 2); % PVNH/SUBC: 25/2
    Lesion_label_HET                = C{6}(Inclusion>0);
    
    LabelPath_HET = '/host/gypsy/local_raid//seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label/';
    fclose(fid);
    
end

%% read PMG demographics
for read_PMG = 1
    
    fid = fopen('DemographicNGrouping_for_PMG_Data.csv');
    C = textscan(fid, '%s%f%s%d%d%s%s%s','Delimiter',',','CollectOutput', 1);
    Inclusion                       = C{4}(:, 1);
    Codes_PMG                       = C{1}(Inclusion>0);
    Age_PMG                         = C{2}(Inclusion>0);     % mean/sd: 28.8/11.2
    Sex_PMG                         = C{3}(Inclusion>0);     % male/female: 12/9      
    Left_PMG                        = C{4}(Inclusion>0, 2) == 1; 
    Right_PMG                       = C{4}(Inclusion>0, 2) == 2;
    Bilateral_PMG                   = C{4}(Inclusion>0, 2) == 3;  % Left/Right/Bilateral: 6/6/9  
    Path_PMG                        = C{5}(Inclusion>0, 1);
    
    fclose(fid); 
    
end

%% load matfiles
for load_files = 1
    
    if(lesion_exclusion)
        load( [ MATDIR '/01_generate_association_matrix_MCD_maskingout_lesion.mat' ] );        
%         load( [ MATDIR '/03_modularity_analysis_result_MCD_maskingout_lesion.mat' ] );
    else
        load( [ MATDIR '/01_generate_association_matrix_MCD1.mat' ] );
        load( [ MATDIR '/03_modularity_analysis_result_MCD.mat' ] );     
    end
    
end

%% compute optimized modularity (i.e., Q, Ci)
for compute_modularity = 1
            
    % since the Louvain algorithm is heuristic, every times it runs, it gives slighltly different values
    % so here we ran 100 times and extracted consistent modular structures among iterations.
    iteration_num = 100;
    interval_thres = 0.01;
    myspar = 0.11:interval_thres:0.4;    
    
    num_of_ROIs = length(Anatomical_Label_structs_idx);
    
    % Control
    for control = 1
        
        Q_set_uw_pos_cont = zeros(1, length(myspar));
        Ci_set_uw_pos_cont = zeros(length(myspar), num_of_ROIs);
                
        for i = 1 : length(myspar)
            i
            wmatrix = threshold_proportional(adjacent_matrix_control, myspar(i));
            highest_Q = 0;
            
            for j = 1 : iteration_num
                
                [Ci Q] = modularity_louvain_und(double(wmatrix>0));
                
                if(highest_Q < Q)
                    highest_Q = Q;
                    Q_set_uw_pos_cont(i) = Q;
                    Ci_set_uw_pos_cont(i, :) = Ci;
                end                
                
            end            
            
        end
        
    end
    
    % FCD Type I
    for FCD_Type_I = 1
        
        Q_set_uw_pos_FCD_Type_I = zeros(1, length(myspar));
        Ci_set_uw_pos_FCD_Type_I = zeros(length(myspar), num_of_ROIs);
                
        for i = 1 : length(myspar)
            i
            wmatrix = threshold_proportional(adjacent_matrix_FCD_Type_I, myspar(i));
            highest_Q = 0;
            
            for j = 1 : iteration_num
                
                [Ci Q] = modularity_louvain_und(double(wmatrix>0));
                
                if(highest_Q < Q)
                    
                    highest_Q = Q;
                    Q_set_uw_pos_FCD_Type_I(i) = Q;
                    Ci_set_uw_pos_FCD_Type_I(i, :) = Ci;
                    
                end                
                
            end            
            
        end
        
    end
    
    % FCD Type II
    for FCD_Type_II = 1
        
        Q_set_uw_pos_FCD_Type_II = zeros(1, length(myspar));
        Ci_set_uw_pos_FCD_Type_II = zeros(length(myspar), num_of_ROIs);
                
        for i = 1 : length(myspar)
            i
            wmatrix = threshold_proportional(adjacent_matrix_FCD_Type_II, myspar(i));
            highest_Q = 0;
            
            for j = 1 : iteration_num
                
                [Ci Q] = modularity_louvain_und(double(wmatrix>0));
                
                if(highest_Q < Q)
                    
                    highest_Q = Q;
                    Q_set_uw_pos_FCD_Type_II(i) = Q;
                    Ci_set_uw_pos_FCD_Type_II(i, :) = Ci;
                    
                end                
                
            end            
            
        end
        
    end
    
    % HET
    for HET = 1
        
        Q_set_uw_pos_HET = zeros(1, length(myspar));
        Ci_set_uw_pos_HET = zeros(length(myspar), num_of_ROIs);
                
        for i = 1 : length(myspar)
            i
            wmatrix = threshold_proportional(adjacent_matrix_HET, myspar(i));
            highest_Q = 0;
            
            for j = 1 : iteration_num
                
                [Ci Q] = modularity_louvain_und(double(wmatrix>0));
                
                if(highest_Q < Q)
                    
                    highest_Q = Q;
                    Q_set_uw_pos_HET(i) = Q;
                    Ci_set_uw_pos_HET(i, :) = Ci;
                    
                end                
                
            end            
            
        end
        
    end
    
    % PMG
    for PMG = 1
        
        Q_set_uw_pos_PMG = zeros(1, length(myspar));
        Ci_set_uw_pos_PMG = zeros(length(myspar), num_of_ROIs);
                
        for i = 1 : length(myspar)
            i
            wmatrix = threshold_proportional(adjacent_matrix_PMG, myspar(i));
            highest_Q = 0;
            
            for j = 1 : iteration_num
                
                [Ci Q] = modularity_louvain_und(double(wmatrix>0));
                
                if(highest_Q < Q)
                    
                    highest_Q = Q;
                    Q_set_uw_pos_PMG(i) = Q;
                    Ci_set_uw_pos_PMG(i, :) = Ci;
                    
                end                
                
            end            
            
        end
        
    end
            
    % The moddiff_sparsity function estimates modularity and its significance with respect to random network
    for compute_significance = 1
        
        parpool(24);
                
        randomization = 1000;
        num_of_iteration_mod = 10;
        interval_thres = 0.01;
        myspar = 0.11:interval_thres:0.4;
        option = [0.04]; % FCD Type-I: [0.04, 0.04]
        
        ms_FCD_Type_I  = moddiff_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_FCD_Type_I_res_val,    myspar(1), myspar(end), interval_thres, randomization, num_of_iteration_mod, option);
        ms_FCD_Type_II = moddiff_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_FCD_Type_II_res_val,   myspar(1), myspar(end), interval_thres, randomization, num_of_iteration_mod, option);
        ms_HET         = moddiff_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_HET_res_val,           myspar(1), myspar(end), interval_thres, randomization, num_of_iteration_mod, option);
        ms_PMG         = moddiff_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_PMG_res_val,           myspar(1), myspar(end), interval_thres, randomization, num_of_iteration_mod, option);
        ms_patient1    = moddiff_sparsity_parfor(asso_mat_FCD_Type_II_res_val,  asso_mat_HET_res_val,           myspar(1), myspar(end), interval_thres, randomization, num_of_iteration_mod, option);
        ms_patient2    = moddiff_sparsity_parfor(asso_mat_FCD_Type_II_res_val,  asso_mat_PMG_res_val,           myspar(1), myspar(end), interval_thres, randomization, num_of_iteration_mod, option);
                
        ms_FCD_Type_I.aQ.aQ1 = (sum(ms_FCD_Type_I.Q_set_uw_pos_group1) -  sum(ms_FCD_Type_I.Q_set_uw_pos_group1([1 end]))/2)*0.01;
        ms_FCD_Type_I.aQ.aQ2 = (sum(ms_FCD_Type_I.Q_set_uw_pos_group2) -  sum(ms_FCD_Type_I.Q_set_uw_pos_group2([1 end]))/2)*0.01;
        ms_FCD_Type_II.aQ.aQ1  = (sum(ms_FCD_Type_II.Q_set_uw_pos_group1) -  sum(ms_FCD_Type_II.Q_set_uw_pos_group1([1 end]))/2)*0.01;
        ms_FCD_Type_II.aQ.aQ2  = (sum(ms_FCD_Type_II.Q_set_uw_pos_group2) -  sum(ms_FCD_Type_II.Q_set_uw_pos_group2([1 end]))/2)*0.01;
        ms_HET.aQ.aQ1 = (sum(ms_HET.Q_set_uw_pos_group1) -  sum(ms_HET.Q_set_uw_pos_group1([1 end]))/2)*0.01;
        ms_HET.aQ.aQ2 = (sum(ms_HET.Q_set_uw_pos_group2) -  sum(ms_HET.Q_set_uw_pos_group2([1 end]))/2)*0.01;
        ms_PMG.aQ.aQ1  = (sum(ms_PMG.Q_set_uw_pos_group1) -  sum(ms_PMG.Q_set_uw_pos_group1([1 end]))/2)*0.01;
        ms_PMG.aQ.aQ2  = (sum(ms_PMG.Q_set_uw_pos_group2) -  sum(ms_PMG.Q_set_uw_pos_group2([1 end]))/2)*0.01;
        
        ms_patient1.aQ.aQ1  = (sum(ms_patient1.Q_set_uw_pos_group1) -  sum(ms_patient1.Q_set_uw_pos_group1([1 end]))/2)*0.01;
        ms_patient1.aQ.aQ2  = (sum(ms_patient1.Q_set_uw_pos_group2) -  sum(ms_patient1.Q_set_uw_pos_group2([1 end]))/2)*0.01;
        ms_patient2.aQ.aQ1  = (sum(ms_patient2.Q_set_uw_pos_group1) -  sum(ms_patient2.Q_set_uw_pos_group1([1 end]))/2)*0.01;
        ms_patient2.aQ.aQ2  = (sum(ms_patient2.Q_set_uw_pos_group2) -  sum(ms_patient2.Q_set_uw_pos_group2([1 end]))/2)*0.01;
        
        ms_FCD_Type_I_rand.aQ.aQ1 = (sum(ms_FCD_Type_I.Q_set_uw_pos_rand_group1, 2) -  sum(ms_FCD_Type_I.Q_set_uw_pos_rand_group1, 2)/2)*0.01;
        ms_FCD_Type_I_rand.aQ.aQ2 = (sum(ms_FCD_Type_I.Q_set_uw_pos_rand_group2, 2) -  sum(ms_FCD_Type_I.Q_set_uw_pos_rand_group2, 2)/2)*0.01;
        ms_FCD_Type_II_rand.aQ.aQ1  = (sum(ms_FCD_Type_II.Q_set_uw_pos_rand_group1, 2) -  sum(ms_FCD_Type_II.Q_set_uw_pos_rand_group1, 2)/2)*0.01;
        ms_FCD_Type_II_rand.aQ.aQ2  = (sum(ms_FCD_Type_II.Q_set_uw_pos_rand_group2, 2) -  sum(ms_FCD_Type_II.Q_set_uw_pos_rand_group2, 2)/2)*0.01;
        ms_HET_rand.aQ.aQ1 = (sum(ms_HET.Q_set_uw_pos_rand_group1, 2) -  sum(ms_HET.Q_set_uw_pos_rand_group1, 2)/2)*0.01;
        ms_HET_rand.aQ.aQ2 = (sum(ms_HET.Q_set_uw_pos_rand_group2, 2) -  sum(ms_HET.Q_set_uw_pos_rand_group2, 2)/2)*0.01;
        ms_PMG_rand.aQ.aQ1  = (sum(ms_PMG.Q_set_uw_pos_rand_group1, 2) -  sum(ms_PMG.Q_set_uw_pos_rand_group1, 2)/2)*0.01;
        ms_PMG_rand.aQ.aQ2  = (sum(ms_PMG.Q_set_uw_pos_rand_group2, 2) -  sum(ms_PMG.Q_set_uw_pos_rand_group2, 2)/2)*0.01;
        
        ms_patient1_rand.aQ.aQ1  = (sum(ms_patient1.Q_set_uw_pos_rand_group1, 2) -  sum(ms_patient1.Q_set_uw_pos_rand_group1, 2)/2)*0.01;
        ms_patient1_rand.aQ.aQ2  = (sum(ms_patient1.Q_set_uw_pos_rand_group2, 2) -  sum(ms_patient1.Q_set_uw_pos_rand_group2, 2)/2)*0.01;
        ms_patient2_rand.aQ.aQ1  = (sum(ms_patient2.Q_set_uw_pos_rand_group1, 2) -  sum(ms_patient2.Q_set_uw_pos_rand_group1, 2)/2)*0.01;
        ms_patient2_rand.aQ.aQ2  = (sum(ms_patient2.Q_set_uw_pos_rand_group2, 2) -  sum(ms_patient2.Q_set_uw_pos_rand_group2, 2)/2)*0.01;
        
        deltaaQ_FCD_Type_I  = ms_FCD_Type_I.aQ.aQ1 - ms_FCD_Type_I.aQ.aQ2;
        deltaaQ_FCD_Type_II = ms_FCD_Type_II.aQ.aQ1 - ms_FCD_Type_II.aQ.aQ2;
        deltaaQ_HET         = ms_HET.aQ.aQ1 - ms_HET.aQ.aQ2;
        deltaaQ_PMG         = ms_PMG.aQ.aQ1 - ms_PMG.aQ.aQ2;        
        
        deltaaQ_patient1 = ms_patient1.aQ.aQ2 - ms_patient1.aQ.aQ1;
        deltaaQ_patient2 = ms_patient2.aQ.aQ2 - ms_patient2.aQ.aQ1;
        
        deltaaQ_FCD_Type_I_rand  = ms_FCD_Type_I_rand.aQ.aQ1 - ms_FCD_Type_I_rand.aQ.aQ2;
        deltaaQ_FCD_Type_II_rand = ms_FCD_Type_II_rand.aQ.aQ1 - ms_FCD_Type_II_rand.aQ.aQ2;
        deltaaQ_HET_rand         = ms_HET_rand.aQ.aQ1 - ms_HET_rand.aQ.aQ2;
        deltaaQ_PMG_rand         = ms_PMG_rand.aQ.aQ1 - ms_PMG_rand.aQ.aQ2;
        deltaaQ_patient1_rand = ms_patient1_rand.aQ.aQ1 - ms_patient1_rand.aQ.aQ2;
        deltaaQ_patient2_rand = ms_patient2_rand.aQ.aQ1 - ms_patient2_rand.aQ.aQ2;
        
        if(deltaaQ_FCD_Type_I > 0) sign_flag = 1; else sign_flag = -1; end
        ms_FCD_Type_I.aQ.aQp = (sum(sign_flag*deltaaQ_FCD_Type_I_rand > sign_flag*deltaaQ_FCD_Type_I)+1)/randomization;
        
        if(deltaaQ_FCD_Type_II > 0) sign_flag = 1; else sign_flag = -1; end
        ms_FCD_Type_II.aQ.aQp = (sum(sign_flag*deltaaQ_FCD_Type_II_rand > sign_flag*deltaaQ_FCD_Type_II)+1)/randomization;
        
        if(deltaaQ_HET > 0) sign_flag = 1; else sign_flag = -1; end
        ms_HET.aQ.aQp = (sum(sign_flag*deltaaQ_HET_rand > sign_flag*deltaaQ_HET)+1)/randomization;
        
        if(deltaaQ_PMG > 0) sign_flag = 1; else sign_flag = -1; end
        ms_PMG.aQ.aQp = (sum(sign_flag*deltaaQ_PMG_rand > sign_flag*deltaaQ_PMG)+1)/randomization;        
        
        if(deltaaQ_patient1 > 0) sign_flag = 1; else sign_flag = -1; end
        ms_patient1.aQ.aQp = (sum(sign_flag*deltaaQ_patient1_rand > sign_flag*deltaaQ_patient1)+1)/randomization;
        
        if(deltaaQ_patient2 > 0) sign_flag = 1; else sign_flag = -1; end
        ms_patient2.aQ.aQp = (sum(sign_flag*deltaaQ_patient2_rand > sign_flag*deltaaQ_patient2)+1)/randomization;
        
        if(save_file)
            if(lesion_exclusion)
                save([ MATDIR '03_modularity_analysis_result_MCD_maskingout_lesion.mat' ], ...
                    'Q_set_uw_pos_cont', 'Q_set_uw_pos_FCD_Type_I', 'Q_set_uw_pos_FCD_Type_II', 'Q_set_uw_pos_HET', 'Q_set_uw_pos_HET', ...
                    'Ci_set_uw_pos_cont', 'Ci_set_uw_pos_FCD_Type_I', 'Ci_set_uw_pos_FCD_Type_II', 'Ci_set_uw_pos_HET', 'Ci_set_uw_pos_PMG', ...
                    'ms_FCD_Type_I', 'ms_FCD_Type_II', 'ms_HET', 'ms_PMG', 'ms_patient1', 'ms_patient2');
            else
                save([ MATDIR '03_modularity_analysis_result_MCD.mat' ], ...
                    'Q_set_uw_pos_cont', 'Q_set_uw_pos_FCD_Type_I', 'Q_set_uw_pos_FCD_Type_II', 'Q_set_uw_pos_HET', 'Q_set_uw_pos_HET', ...
                    'Ci_set_uw_pos_cont', 'Ci_set_uw_pos_FCD_Type_I', 'Ci_set_uw_pos_FCD_Type_II', 'Ci_set_uw_pos_HET', 'Ci_set_uw_pos_PMG', ...
                    'ms_FCD_Type_I', 'ms_FCD_Type_II', 'ms_HET', 'ms_PMG', 'ms_patient1', 'ms_patient2');
            end
        end
        
    end        
    
    % visualization
    for visualization = 1
        
        if(lesion_exclusion)
            postfix = '_lesion_exclued';
            OUTPATH = [ OUTDIR '/module_projection_lesion_exclusion/' ];
        else
            postfix = '';
            OUTPATH = [ OUTDIR '/module_projection/' ];
        end
        
        randomization = 1000;
        num_of_iteration_mod = 10;
        interval_thres = 0.01;
        myspar = 0.11:interval_thres:0.4;
        
        groups    = {'FCD_Type_I', 'FCD_Type_II', 'HET', 'PMG'};
        group_str = {'FCD Type I', 'FCD Type II', 'HET', 'PMG'};
        figure; plot(1:length(myspar), [ Q_set_uw_pos_cont; Q_set_uw_pos_FCD_Type_I; Q_set_uw_pos_FCD_Type_II; Q_set_uw_pos_HET; Q_set_uw_pos_PMG ]);
        legend('control', group_str);
        color_lines = get(gca, 'ColorOrder');
        
        for modularity = 1
            
            sig_pos  = [ -0.15 -0.15 -0.15 -0.25 ];
            ylim_set = [ -0.2 0.2; -0.2 0.2; -0.2 0.2; -0.3 0.3; ];
            
            for g = 1 : 4
                
                colorline_idx = g + 1;
                eval(['group = ms_' groups{g} ';' ]);
                FDRp = FDR(group.Qgroup12_diff_p, 0.05);
                figure; hold on;
                plot(myspar*100, group.Q_set_uw_pos_rand_group1-group.Q_set_uw_pos_rand_group2, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ group.lower95Q12diff; group.upper95Q12diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, group.rerandQ12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', group.rerandQ12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, group.Q_set_uw_pos_group1-group.Q_set_uw_pos_group2, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', group.Q_set_uw_pos_group1-group.Q_set_uw_pos_group2, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(~isempty(FDRp))
                        if(group.Qgroup12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos(g), '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        else
                            if(group.Qgroup12_diff_p(i)<=0.05)
                                text(myspar(i)*100, sig_pos(g), 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                            end
                        end
                    else
                        if(group.Qgroup12_diff_p(i)<=0.05)
                            text(myspar(i)*100, sig_pos(g), 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                end
                
                ylim(ylim_set(g, :));
                
            end            
            
        end
        
        for modular_structure = 1
            
            % visualization of community strucutre at each density and group
            if printfigs == 1
                for i = 1 : length(myspar)
                    project_detection_community(Ci_set_uw_pos_cont(i, :), num2str(Q_set_uw_pos_cont(i)), AAL_surf_data_both, S, 1);
                    exportfigbo(gcf,[ OUTPATH 'uw_pos_sparce_NC_' num2str(i) '_' num2str(myspar(i)) '_modularity' postfix '.png'], 'png',6); close(gcf);
                    
                    project_detection_community(Ci_set_uw_pos_NLES(i, :), num2str(Q_set_uw_pos_NLES(i)), AAL_surf_data_both, S, 1);
                    exportfigbo(gcf,[ OUTPATH 'uw_pos_sparce_NLES_' num2str(i) '_' num2str(myspar(i)) '_modularity' postfix '.png'], 'png',6); close(gcf);
                    
                    project_detection_community(Ci_set_uw_pos_FCD(i, :), num2str(Q_set_uw_pos_FCD(i)), AAL_surf_data_both, S, 1);
                    exportfigbo(gcf,[ OUTPATH 'uw_pos_sparce_FCD_' num2str(i) '_' num2str(myspar(i)) '_modularity' postfix '.png'], 'png',6); close(gcf);
                end
            end            
           
            %% AAL indexing
            starting_point_L_index = 0;
            reordered_ROI_index = [];
            for j = 1 : size(Anatomical_Label_structs, 1)
                for k = 1 : size(Anatomical_Label_structs(j).AAL_L_index, 2)
                    if(all(isnan(Anatomical_Label_structs(j).centroids_L_index{k})))
                        disp([ 'Left ' Anatomical_Label_structs(j).subarea{k} ' : unassigned' ]);
                    else
                        num_of_ROIs = Anatomical_Label_structs(j).AAL_L_index(k);
                        starting_point_L_index = starting_point_L_index + 1;
                        disp([ 'Left ' Anatomical_Label_structs(j).subarea{k} ' : ' num2str(starting_point_L_index) ]);
                        reordered_ROI_index = [ reordered_ROI_index; num_of_ROIs; ];
                    end
                end
            end
            
            starting_point_R_index = starting_point_L_index;
            for j = 1 : size(Anatomical_Label_structs, 1)
                for k = 1 : size(Anatomical_Label_structs(j).AAL_L_index, 2)
                    if(all(isnan(Anatomical_Label_structs(j).centroids_R_index{k})))
                        disp([ 'Right ' Anatomical_Label_structs(j).subarea{k} ' : unassigned' ]);
                    else
                        num_of_ROIs = Anatomical_Label_structs(j).AAL_R_index(k);
                        starting_point_R_index = starting_point_R_index + 1;
                        disp([ 'Right ' Anatomical_Label_structs(j).subarea{k} ' : ' num2str(starting_point_R_index) ]);
                        reordered_ROI_index = [ reordered_ROI_index; num_of_ROIs; ];
                    end
                end
            end
            
            starting_point_L_index = 0;
            count = 0;
            for j = 1 : size(Anatomical_Label_structs, 1)
                for k = 1 : size(Anatomical_Label_structs(j).AAL_L_index, 2)
                    if(all(isnan(Anatomical_Label_structs(j).centroids_L_index{k})))
                        disp([ 'Left ' Anatomical_Label_structs(j).subarea{k} ' : unassigned' ]);
                    else
                        count = count + 1;
                        Anatomical_Label{count} = ['L_' Anatomical_Label_structs(j).subarea{k} ':' num2str(Anatomical_Label_structs(j).AAL_L_index(k)) ];
                    end
                end
            end
            
            starting_point_R_index = starting_point_L_index;
            for j = 1 : size(Anatomical_Label_structs, 1)
                for k = 1 : size(Anatomical_Label_structs(j).AAL_L_index, 2)
                    if(all(isnan(Anatomical_Label_structs(j).centroids_R_index{k})))
                        disp([ 'Right ' Anatomical_Label_structs(j).subarea{k} ' : unassigned' ]);
                    else
                        count = count + 1;
                        Anatomical_Label{count} = ['R_' Anatomical_Label_structs(j).subarea{k} ':' num2str(Anatomical_Label_structs(j).AAL_R_index(k)) ];
                    end
                end
            end
            
            %% Individual density modular organization
            parcel=SurfStatReadData1('aal_both_rsl_final.txt');
            AAL_ROI_index = unique(parcel); AAL_ROI_index = AAL_ROI_index(2:end);
            nodePos=NOELNetVisuNodePosition(S, parcel, 'n');
            noel_colormap=jet(255);
            
            % 1) before making the color of community structures comparable between groups
            i = 9; % 13%, the most conservative one among thresholds where the group shows significant difference of modularity
            
            % control, intra-community
            wmatrix = threshold_proportional(adjacent_matrix_control, myspar(i));
            f1 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_cont(i, :), wmatrix>0, nodePos, reordered_ROI_index, ...
                S, 'black', 'colormap', noel_colormap, 'links_visible', 0.30, 'link_width', 0.3, 'surf_opacity', 0);            
            % control, inter-community
            f2 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_cont(i, :), wmatrix>0, nodePos, reordered_ROI_index, ...
                S, 'black', 'colormap', noel_colormap, 'links_visible', 0.30, 'link_width', 0.1, 'show_modular_link', 'inter', [ones(255, 1) ones(255, 1) ones(255, 1)]);
            
            if printfigs == 1
                exportfigbo(f1,[OUTPATH 'control_intra_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
                exportfigbo(f2,[OUTPATH 'control_inter_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
            end
            
            % NLES, intra-community
            wmatrix = threshold_proportional(adjacent_matrix_NLES, myspar(i));
            f1 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_NLES(i, :), wmatrix>0, nodePos, reordered_ROI_index, ...
                S, 'black', 'colormap', noel_colormap, 'links_visible', 0.13, 'link_width', 0.3, 'surf_opacity', 0);    
            
            % NLES, inter-community
            f2 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_NLES(i, :), wmatrix>0, nodePos, reordered_ROI_index, ...
                S, 'black', 'colormap', noel_colormap, 'links_visible', 0.30, 'link_width', 0.1, 'show_modular_link', 'inter', [ones(255, 1) ones(255, 1) ones(255, 1)]);
            
            if printfigs == 1
                exportfigbo(f1,[OUTPATH 'NLES_intra_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
                exportfigbo(f2,[OUTPATH 'NLES_inter_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
            end
            
            % FCD, intra-community
            wmatrix = threshold_proportional(adjacent_matrix_FCD, myspar(i));
            f1 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_FCD(i, :), wmatrix>0, nodePos, reordered_ROI_index, ...
                S, 'black', 'colormap', noel_colormap, 'links_visible', 0.15, 'link_width', 0.3, 'surf_opacity', 0);               
            
            % FCD, inter-community
            f2 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_FCD(i, :), wmatrix>0, nodePos, reordered_ROI_index, ...
                S, 'black', 'colormap', noel_colormap, 'links_visible', 0.30, 'link_width', 0.1, 'show_modular_link', 'inter', [ones(255, 1) ones(255, 1) ones(255, 1)]);
            
            if printfigs == 1
                exportfigbo(f1,[OUTPATH 'FCD_intra_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
                exportfigbo(f2,[OUTPATH 'FCD_inter_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
            end        
            
            % 2) after making the color of community structures comparable between groups            
            i = 9;
            backgroun_temp = 'black';
            
            % Control
            % yellow: 4
            % dark red: 6
            % blue: 1
            % light green: 3
            % light blue: 2
            % red : 5
            
            % NLES
            % orange 3 --> light blue 2
            % light blue 1 --> blue 1
            % dark red 4 --> dark red 6
            % light green 2 --> light green 3
            Ci_set_uw_pos_NLES_temp = Ci_set_uw_pos_NLES(i, :);
            Ci_set_uw_pos_NLES_temp(Ci_set_uw_pos_NLES(i, :) == 3) = 2;
            Ci_set_uw_pos_NLES_temp(Ci_set_uw_pos_NLES(i, :) == 1) = 6;
            Ci_set_uw_pos_NLES_temp(Ci_set_uw_pos_NLES(i, :) == 4) = 1;
            Ci_set_uw_pos_NLES_temp(Ci_set_uw_pos_NLES(i, :) == 2) = 3;
            wmatrix = threshold_proportional(adjacent_matrix_NLES, myspar(i));
            f1 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_NLES_temp, wmatrix>0, nodePos, reordered_ROI_index, ...
                S, backgroun_temp, 'colormap', noel_colormap, 'links_visible', 0.13, 'link_width', 0.3, 'surf_opacity', 0);
            
            f2 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_NLES_temp, wmatrix>0, nodePos, reordered_ROI_index, ...
                S, backgroun_temp, 'colormap', noel_colormap, 'links_visible', 0.13, 'link_width', 0.1, 'show_modular_link', 'inter', [ones(255, 1) ones(255, 1) ones(255, 1)]);
            
            if printfigs == 1
                exportfigbo(f1,[OUTPATH 'NLES_intra_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
                exportfigbo(f2,[OUTPATH 'NLES_inter_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
            end            
            
            % FCD
            % cyan 2 --> yellow 4
            % yellow 3 --> light blue 2
            % dark red 5 --> blue 1
            % blue 1 --> light green 3
            % orange 4 --> dark red 6
            Ci_set_uw_pos_FCD_temp = Ci_set_uw_pos_FCD(i, :);
            Ci_set_uw_pos_FCD_temp(Ci_set_uw_pos_FCD(i, :) == 2) = 4;
            Ci_set_uw_pos_FCD_temp(Ci_set_uw_pos_FCD(i, :) == 3) = 2;
            Ci_set_uw_pos_FCD_temp(Ci_set_uw_pos_FCD(i, :) == 5) = 1;
            Ci_set_uw_pos_FCD_temp(Ci_set_uw_pos_FCD(i, :) == 1) = 3;
            Ci_set_uw_pos_FCD_temp(Ci_set_uw_pos_FCD(i, :) == 4) = 6;
            wmatrix = threshold_proportional(adjacent_matrix_FCD, myspar(i));
            f1 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_FCD_temp, wmatrix>0, nodePos, reordered_ROI_index, ...
                S, backgroun_temp, 'colormap', noel_colormap, 'links_visible', 0.13, 'link_width', 0.3, 'surf_opacity', 0);
            
            f2 = figure; NOELNetVisuModulesOnSurf(Ci_set_uw_pos_FCD_temp, wmatrix>0, nodePos, reordered_ROI_index, ...
                S, backgroun_temp, 'colormap', noel_colormap, 'links_visible', 0.13, 'link_width', 0.1, 'show_modular_link', 'inter', [ones(255, 1) ones(255, 1) ones(255, 1)]);
            
            if printfigs == 1
                exportfigbo(f1,[OUTPATH 'FCD_intra_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
                exportfigbo(f2,[OUTPATH 'FCD_inter_modular_organization_' num2str(myspar(i)) '.png'], 'png', 6); close(f);
            end
            
        end
        
    end    
    
end

%% compute scaled inclusivity
for compute_SI = 1
    
    % since the Louvain algorithm is heuristic, every times it runs, it gives slighltly different values
    % so here we ran 100 times and extracted consistent modular structures among iterations.
    iteration_num = 100;
    interval_thres = 0.01;
    myspar = 0.11:interval_thres:0.4;    
    
    num_of_ROIs = length(Anatomical_Label_structs_idx);
    
    for community_structure = 1
        
        % Control
        for control = 1
            
            modularity_projected_map_cont = zeros(length(myspar), size(S.coord, 2));
            figure;
            for i = 1 : length(myspar)
                
                modularity_projected_map = zeros(1, size(S.coord, 2));
                Ci_temp = Ci_set_uw_pos_cont(i, :);
                
                for j = 1 : max(Ci_temp)
                    
                    AAL_temp = Anatomical_Label_structs_idx(Ci_temp == j);
                    modularity_projected_map(find(ismember(aal_data, AAL_temp))) = j;
                    
                end
                
                modularity_projected_map_cont(i, :) = modularity_projected_map;
                SurfStatView(modularity_projected_map_cont(i, :), S);
                export_fig([RDIR '/01_control/modularity_' num2str(myspar(i)) ], '-m1', '-png');
                
            end
            
            close(gcf);
            
        end
        
        % FCD Type I
        for FCD_Type_I = 1
            
            modularity_projected_map_FCD_Type_I = zeros(length(myspar), size(S.coord, 2));
            figure;
            for i = 1 : length(myspar)
                
                modularity_projected_map = zeros(1, size(S.coord, 2));
                Ci_temp = Ci_set_uw_pos_FCD_Type_I(i, :);
                
                for j = 1 : max(Ci_temp)
                    
                    AAL_temp = Anatomical_Label_structs_idx(Ci_temp == j);
                    modularity_projected_map(find(ismember(aal_data, AAL_temp))) = j;
                    
                end
                
                modularity_projected_map_FCD_Type_I(i, :) = modularity_projected_map;
                SurfStatView(modularity_projected_map_FCD_Type_I(i, :), S);
                export_fig([RDIR '/02_FCD_Type_I/modularity_' num2str(myspar(i)) ], '-m1', '-png');
                
            end
            
            close(gcf);
            
        end
        
        % FCD Type II
        for FCD_Type_II = 1
            
            modularity_projected_map_FCD_Type_II = zeros(length(myspar), size(S.coord, 2));
            figure;
            for i = 1 : length(myspar)
                
                modularity_projected_map = zeros(1, size(S.coord, 2));
                Ci_temp = Ci_set_uw_pos_FCD_Type_II(i, :);
                
                for j = 1 : max(Ci_temp)
                    
                    AAL_temp = Anatomical_Label_structs_idx(Ci_temp == j);
                    modularity_projected_map(find(ismember(aal_data, AAL_temp))) = j;
                    
                end
                
                modularity_projected_map_FCD_Type_II(i, :) = modularity_projected_map;
                SurfStatView(modularity_projected_map_FCD_Type_II(i, :), S);
                export_fig([RDIR '/03_FCD_Type_II/modularity_' num2str(myspar(i)) ], '-m1', '-png');
                
            end
            
            close(gcf);
            
        end
        
        % HET
        for HET = 1
            
            modularity_projected_map_HET = zeros(length(myspar), size(S.coord, 2));
            figure;
            for i = 1 : length(myspar)
                
                modularity_projected_map = zeros(1, size(S.coord, 2));
                Ci_temp = Ci_set_uw_pos_HET(i, :);
                
                for j = 1 : max(Ci_temp)
                    
                    AAL_temp = Anatomical_Label_structs_idx(Ci_temp == j);
                    modularity_projected_map(find(ismember(aal_data, AAL_temp))) = j;
                    
                end
                
                modularity_projected_map_HET(i, :) = modularity_projected_map;
                SurfStatView(modularity_projected_map_HET(i, :), S);
                export_fig([RDIR '/04_HET/modularity_' num2str(myspar(i)) ], '-m1', '-png');
                
            end
            
            close(gcf);
            
        end
        
        % PMG
        for PMG = 1
            
            modularity_projected_map_PMG = zeros(length(myspar), size(S.coord, 2));
            figure;
            for i = 1 : length(myspar)
                
                modularity_projected_map = zeros(1, size(S.coord, 2));
                Ci_temp = Ci_set_uw_pos_PMG(i, :);
                
                for j = 1 : max(Ci_temp)
                    
                    AAL_temp = Anatomical_Label_structs_idx(Ci_temp == j);
                    modularity_projected_map(find(ismember(aal_data, AAL_temp))) = j;
                    
                end
                
                modularity_projected_map_PMG(i, :) = modularity_projected_map;
                SurfStatView(modularity_projected_map_PMG(i, :), S);
                export_fig([RDIR '/05_PMG/modularity_' num2str(myspar(i)) ], '-m1', '-png');
                
            end
            
            close(gcf);
            
        end
        
    end
    
    for compute_SI = 1;
        
        parpool(24);
                
        randomization = 1000;
        boot_iteration = 100;
        num_of_iteration_mod = 10;
        interval_thres = 0.01;
        myspar = 0.11:interval_thres:0.4;
        
        ss_FCD_Type_I  = scaled_inclusivity_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_FCD_Type_I_res_val,    myspar(1), myspar(end), interval_thres, randomization, boot_iteration, num_of_iteration_mod, true);
        ss_FCD_Type_II = scaled_inclusivity_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_FCD_Type_II_res_val,   myspar(1), myspar(end), interval_thres, randomization, boot_iteration, num_of_iteration_mod, true);
        ss_HET         = scaled_inclusivity_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_HET_res_val,           myspar(1), myspar(end), interval_thres, randomization, boot_iteration, num_of_iteration_mod, true);
        ss_PMG         = scaled_inclusivity_sparsity_parfor(asso_mat_cont_res_val,         asso_mat_PMG_res_val,           myspar(1), myspar(end), interval_thres, randomization, boot_iteration, num_of_iteration_mod, true);
        ss_patient1    = scaled_inclusivity_sparsity_parfor(asso_mat_FCD_Type_II_res_val,  asso_mat_HET_res_val,           myspar(1), myspar(end), interval_thres, randomization, boot_iteration, num_of_iteration_mod, true);
        ss_patient2    = scaled_inclusivity_sparsity_parfor(asso_mat_FCD_Type_II_res_val,  asso_mat_PMG_res_val,           myspar(1), myspar(end), interval_thres, randomization, boot_iteration, num_of_iteration_mod, true);
        
        if(save_file)
            
        
    end
        
    end

end