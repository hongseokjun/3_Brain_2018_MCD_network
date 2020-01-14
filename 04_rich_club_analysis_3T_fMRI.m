clear all
close all
lesion_exclusion = 0;

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/'));
P='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/';
MATDIR='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/01_matfiles/';
OUTPATH='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/00_association_matrix/';
printfigs = 0; % don't flag on unless you really want to save the figures as .png, otherwise your precious time is gone...
CIVET_ver1 = 1; % 1 -> new data processed by 1.2.1 | 2 -> old data processed by 1.2.0
CIVET_ver2 = 2; % 1 -> HET data processed after the removal of nodules | 2 -> HET data processed before the removal of nodules

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
    
    fid = fopen('DemographicNGrouping_for_controls_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    Codes_cont                      = C{1}(:, 1);    
    Path_cont                       = C{1}(:, 2);
    Age_cont                        = C{2};                 % mean/sd:     29.2/7.1
    Sex_cont                        = C{3}(:, 1);           % male/female:   18/15    
    fclose(fid);
    
end

%% read FCD Type-II demographics
for read_FCD_Type_II = 1
    
    fid = fopen('DemographicNGrouping_for_FCD_Type-II_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%f%s%d%d%d%d%s%d%d%d%d','Delimiter',',','CollectOutput', 1);
    Codes_FCD_Type_II               = C{1};
    Age_FCD_Type_II                 = C{2};                 % mean/sd: 27.3/8.6
    Sex_FCD_Type_II                 = C{3};                 % male/female: 14/13        
    Left_FCD_Type_II                = C{4}(:, 1);           % 13
    Right_FCD_Type_II               = C{4}(:, 2);           % 14
    Pre_Central_FCD_Type_II         = C{4}(:, 3);           % 17
    Post_Central_FCD_Type_II        = C{4}(:, 4);           % 10    
    Path_FCD_Type_II                = C{5};
    Duration_FCD_Type_II            = double(C{6}(:, 1));   % mean/sd: 20.9/11.9
    Operation_FCD_Type_II           = double(C{6}(:, 2));   % 27
    Scanner_FCD_Type_II             = double(C{6}(:, 3));
    Lesion_Volume                   = double(C{6}(:, 4));    
    
    LabelPath_FCD = '/host/gypsy/local_raid/seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label_3T/';
    fclose(fid); 
    
end

%% read HET demographics
for read_HET = 1
    
    fid = fopen('DemographicNGrouping_for_HET_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%f%s%d%s%s%d','Delimiter',',','CollectOutput', 1);
    Inclusion                       = C{4}(:, 1);
    Codes_HET                       = C{1}(Inclusion>0);
    Age_HET                         = C{2}(Inclusion>0);    % mean/sd: 31.4/10.6
    Sex_HET                         = C{3}(Inclusion>0);    % male/female: 8/6
    Left_HET                        = C{4}(Inclusion>0, 1) == 1; 
    Right_HET                       = C{4}(Inclusion>0, 1) == 2;
    Bilateral_HET                   = C{4}(Inclusion>0, 1) == 3;  % Left/Right/Bilateral: 4/4/6  
    Path_HET                        = C{5}(Inclusion>0, 1);
    Subtype_HET                     = C{5}(Inclusion>0, 2); % PVNH/SUBC: 13/1
    res_type_HET                    = C{6}(Inclusion>0);    % 1: low res (4x4x4), 2: high res (2x2x4)
    
    LabelPath_HET = '/host/gypsy/local_raid/seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label/';
    fclose(fid);
    
end

%% read PMG demographics
for read_PMG = 1
    
    fid = fopen('DemographicNGrouping_for_PMG_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%f%s%d%s%s%d','Delimiter',',','CollectOutput', 1);
    Inclusion                       = C{4}(:, 1);
    Codes_PMG                       = C{1}(Inclusion>0);
    Age_PMG                         = C{2}(Inclusion>0);     % mean/sd: 29.5/11.4
    Sex_PMG                         = C{3}(Inclusion>0);     % male/female: 6/5      
    Left_PMG                        = C{4}(Inclusion>0, 1) == 1; 
    Right_PMG                       = C{4}(Inclusion>0, 1) == 2;
    Bilateral_PMG                   = C{4}(Inclusion>0, 1) == 3;  % Left/Right/Bilateral: 2/1/8  
    Path_PMG                        = C{5}(Inclusion>0, 1);
    res_type_PMG                    = C{6}(Inclusion>0);    % 1: low res (4x4x4), 2: high res (2x2x4)
    
    fclose(fid);
    
end

%% load matfiles
for load_files = 1
    
    if(lesion_exclusion)
        
%        load( [ MATDIR '/01_generate_association_matrix_MCD_maskingout_lesion.mat' ] );        
%        load( [ MATDIR '/03_modularity_analysis_result_MCD_maskingout_lesion.mat' ] );

    else
        
        load( [ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_GSR_Age_Sex1.mat' ] );
        load( [ MATDIR '/04_rich_club_analysis_result_MCD_rsfMRI.mat' ] );
        load( [ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_GSR_Age_Sex1_micro.mat' ] );       
        load( [ MATDIR '/04_rich_club_analysis_result_MCD_rsfMRI_micro.mat' ] );

    end
    
end

%% compute rich club coefficient
for compute_rich_club = 1
            
    % since the Louvain algorithm is heuristic, every times it runs, it gives slighltly different values
    % so here we ran 100 times and extracted consistent modular structures among iterations.   
    num_of_ROIs = length(Anatomical_Label_structs_idx);
            
    % The moddiff_sparsity function estimates modularity and its significance with respect to random network
    for compute_significance = 1
        
        for original_AAL = 1
            
            num_of_cont         = length(Codes_cont);
            num_of_FCD_Type_II  = length(Codes_FCD_Type_II);
            num_of_HET          = length(Codes_HET);
            num_of_PMG          = length(Codes_PMG);
            type_mat            = 2;
            
            if(type_mat == 1)
                mat_cont            = zadjacent_matrix_control1;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II1;
                mat_HET             = zadjacent_matrix_HET1;
                mat_PMG             = zadjacent_matrix_PMG1;
            elseif(type_mat == 2)
                mat_cont            = zadjacent_matrix_control2;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II2;
                mat_HET             = zadjacent_matrix_HET2;
                mat_PMG             = zadjacent_matrix_PMG2;
            elseif(type_mat == 3)
                mat_cont            = zadjacent_matrix_control3;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II3;
                mat_HET             = zadjacent_matrix_HET3;
                mat_PMG             = zadjacent_matrix_PMG3;
            elseif(type_mat == 4)
                mat_cont            = zadjacent_matrix_control;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II;
                mat_HET             = zadjacent_matrix_HET;
                mat_PMG             = zadjacent_matrix_PMG;
            end
            
            parpool(32);
            iternation_num = 1000;
            interval_thres = 0.01;
            randwiring = 100;
            myspar = 0.05:interval_thres:0.4;
            AAL_index_iter = length(Anatomical_Label_structs_idx);
            mat_thres = 0.2;
            
            % 0) binarized averaged matrix construction
            averaged_mat_cont        = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            averaged_mat_FCD_Type_II = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            averaged_mat_HET         = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            averaged_mat_PMG         = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            
            for i = 1 : length(myspar)
                
                %% control
                mat_cont_temp = mat_cont*0;
                for j = 1 : num_of_cont
                    temp = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    mat_cont_temp(:, :, j) = temp > 0;
%                     mat_cont_temp(:, :, j) = temp;
                end
                averaged_mat_cont(:, :, i) = (sum(mat_cont_temp, 3)/num_of_cont)>mat_thres;
%                 averaged_mat_cont(:, :, i) = mean(mat_cont_temp, 3)>mat_thres;
                
                %% FCD Type II
                mat_FCD_Type_II_temp = mat_FCD_Type_II*0;
                for j = 1 : num_of_FCD_Type_II
                    temp = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    mat_FCD_Type_II_temp(:, :, j) = temp > 0;
%                     mat_FCD_Type_II_temp(:, :, j) = temp;
                end
                averaged_mat_FCD_Type_II(:, :, i) = (sum(mat_FCD_Type_II_temp, 3)/num_of_FCD_Type_II)>mat_thres;
%                 averaged_mat_FCD_Type_II(:, :, i) = mean(mat_FCD_Type_II_temp, 3)>mat_thres;
                
                %% HET
                mat_HET_temp = mat_HET*0;
                for j = 1 : num_of_HET
                    temp = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    mat_HET_temp(:, :, j) = temp > 0;
%                     mat_HET_temp(:, :, j) = temp;
                end
                averaged_mat_HET(:, :, i) = (sum(mat_HET_temp, 3)/num_of_HET)>mat_thres;
%                 averaged_mat_HET(:, :, i) = mean(mat_HET_temp, 3)>mat_thres;
                
                %% PMG
                mat_PMG_temp = mat_PMG*0;
                for j = 1 : num_of_PMG
                    temp = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    mat_PMG_temp(:, :, j) = temp > 0;
%                     mat_PMG_temp(:, :, j) = temp;
                end
                averaged_mat_PMG(:, :, i) = (sum(mat_PMG_temp, 3)/num_of_PMG)>mat_thres;
%                 averaged_mat_PMG(:, :, i) = mean(mat_PMG_temp, 3)>mat_thres;
                
            end
            
            for rc_comp_individual = 1
                
                rc_uw_pos_cont         = cell(length(myspar), num_of_cont);
                rc_uw_pos_FCD_Type_II  = cell(length(myspar), num_of_FCD_Type_II);
                rc_uw_pos_HET          = cell(length(myspar), num_of_HET);
                rc_uw_pos_PMG          = cell(length(myspar), num_of_PMG);
                
                for i = 1 : length(myspar)
                    
                    %% control
                    parfor j = 1 : num_of_cont
                        wmatrix                          = threshold_proportional(mat_cont(:, :, j), myspar(i));
                        rc_uw_pos_cont{i, j}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_cont{i, j}             = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                    %% FCD Type II
                    parfor j = 1 : num_of_FCD_Type_II
                        wmatrix                          = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                        rc_uw_pos_FCD_Type_II{i, j}      = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_FCD_Type_II{i, j}      = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                    %% HET
                    parfor j = 1 : num_of_HET
                        wmatrix                         = threshold_proportional(mat_HET(:, :, j), myspar(i));
                        rc_uw_pos_HET{i, j}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_HET{i, j}             = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                    %% PMG
                    parfor j = 1 : num_of_PMG
                        wmatrix                          = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                        rc_uw_pos_PMG{i, j}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_PMG{i, j}             = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                end
                
                min_cont        = zeros(length(myspar), 1);
                min_FCD_Type_II = zeros(length(myspar), 1);
                min_HET         = zeros(length(myspar), 1);
                min_PMG         = zeros(length(myspar), 1);
                for i = 1 : length(myspar)
                    
                    min_val = 100;
                    for j = 1 : num_of_cont
                        temp = rc_uw_pos_cont{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_cont(i) = min_val;
                    
                    min_val = 100;
                    for j = 1 : num_of_FCD_Type_II
                        temp = rc_uw_pos_FCD_Type_II{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_FCD_Type_II(i) = min_val;
                    
                    min_val = 100;
                    for j = 1 : num_of_HET
                        temp = rc_uw_pos_HET{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_HET(i) = min_val;
                    
                    min_val = 100;
                    for j = 1 : num_of_PMG
                        temp = rc_uw_pos_PMG{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_PMG(i) = min_val;
                    
                end
                
                min_k = min([ min_cont min_FCD_Type_II min_HET min_PMG ], [], 2);
                
                rc_uw_pos_cont_final        = cell(length(myspar), 1);
                rc_uw_pos_FCD_Type_II_final = cell(length(myspar), 1);
                rc_uw_pos_HET_final         = cell(length(myspar), 1);
                rc_uw_pos_PMG_final         = cell(length(myspar), 1);
                
                for i = 1 : length(myspar)
                    
                    for j = 1 : num_of_cont
                        rc_uw_pos_cont_final{i} = [ rc_uw_pos_cont_final{i}; rc_uw_pos_cont{i, j}.normR(1:min_k(i)); ];
                    end
                    for j = 1 : num_of_FCD_Type_II
                        rc_uw_pos_FCD_Type_II_final{i} = [ rc_uw_pos_FCD_Type_II_final{i}; rc_uw_pos_FCD_Type_II{i, j}.normR(1:min_k(i)); ];
                    end
                    for j = 1 : num_of_HET
                        rc_uw_pos_HET_final{i} = [ rc_uw_pos_HET_final{i}; rc_uw_pos_HET{i, j}.normR(1:min_k(i)); ];
                    end
                    for j = 1 : num_of_PMG
                        rc_uw_pos_PMG_final{i} = [ rc_uw_pos_PMG_final{i}; rc_uw_pos_PMG{i, j}.normR(1:min_k(i)); ];
                    end
                    
                end
                
                rc_uw_pos_group = cell(length(myspar), 1);
                
                for i = 1 : length(myspar)
                    
                    rc_uw_pos_group{i} = [ mean(rc_uw_pos_cont_final{i}, 1); mean(rc_uw_pos_FCD_Type_II_final{i}, 1); mean(rc_uw_pos_HET_final{i}, 1); mean(rc_uw_pos_PMG_final{i}, 1) ];
                    
                end
            end
            
            for rc_comp_group_mat = 1
                
                rc_uw_pos_cont         = cell(length(myspar), 1);
                rc_uw_pos_FCD_Type_II  = cell(length(myspar), 1);
                rc_uw_pos_HET          = cell(length(myspar), 1);
                rc_uw_pos_PMG          = cell(length(myspar), 1);
                
                parfor i = 1 : length(myspar)
                    
                    %% control
                    wmatrix                       = threshold_proportional(mean(mat_cont, 3), myspar(i));
                    rc_uw_pos_cont{i}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_cont{i}           = hrich_club_wu(wmatrix, 1, randwiring);

                    %% FCD Type II
                    wmatrix                        = threshold_proportional(mean(mat_FCD_Type_II, 3), myspar(i));
                    rc_uw_pos_FCD_Type_II{i}       = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_FCD_Type_II{i}     = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% HET
                    wmatrix                     = threshold_proportional(mean(mat_HET, 3), myspar(i));
                    rc_uw_pos_HET{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_HET{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% PMG
                    wmatrix                     = threshold_proportional(mean(mat_PMG, 3), myspar(i));
                    rc_uw_pos_PMG{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_PMG{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                end
                
            end
            
            for rc_comp_group_mat_threshold = 1
                
                rc_uw_pos_cont         = cell(length(myspar), 1);
                rc_uw_pos_FCD_Type_II  = cell(length(myspar), 1);
                rc_uw_pos_HET          = cell(length(myspar), 1);
                rc_uw_pos_PMG          = cell(length(myspar), 1);
                
                parfor i = 1 : length(myspar)
                    
                    %% control
                    wmatrix                       = averaged_mat_cont(:, :, i);
                    rc_uw_pos_cont{i}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_cont{i}           = hrich_club_wu(wmatrix, 1, randwiring);

                    %% FCD Type II
                    wmatrix                        = averaged_mat_FCD_Type_II(:, :, i);
                    rc_uw_pos_FCD_Type_II{i}       = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_FCD_Type_II{i}     = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% HET
                    wmatrix                        = averaged_mat_HET(:, :, i);
                    rc_uw_pos_HET{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_HET{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% PMG
                    wmatrix                        = averaged_mat_PMG(:, :, i);
                    rc_uw_pos_PMG{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_PMG{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                end
                
                randomization = 1000;
                brandnetwork = 1;
                interval_thres = 0.01;
                myspar = 0.05:interval_thres:0.36;
                option.bagging = 1;
                
                [ rc_FCD_Type_II rand_rc_FCD_Type_II ] = rcdiff_sparsity_parfor_rsfMRI(mat_cont, mat_FCD_Type_II, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);  figure; pause(5);
                [ rc_HET         rand_rc_HET         ] = rcdiff_sparsity_parfor_rsfMRI(mat_cont, mat_HET, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);          figure; pause(5);
                [ rc_PMG         rand_rc_PMG         ] = rcdiff_sparsity_parfor_rsfMRI(mat_cont, mat_PMG, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);          figure; pause(5);
                
            end
            
            if(visualization == 1)
                
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/03_rich_club_analysis/';
                
                for i = 9: -1 : 9
                    
                    k = length(rc_uw_pos_cont{i}.normR);
                    
                    temp_cont = rc_uw_pos_cont{i}.R;
                    if(length(rc_uw_pos_FCD_Type_II{i}.R) >= k)
                        temp_FCD_Type_II = rc_uw_pos_FCD_Type_II{i}.R(1:k);
                    else
                        temp_FCD_Type_II = [ rc_uw_pos_FCD_Type_II{i}.R zeros(1, k-length(rc_uw_pos_FCD_Type_II{i}.R)) ];
                    end                        
                    if(length(rc_uw_pos_HET{i}.R) >= k)
                        temp_HET = rc_uw_pos_HET{i}.R(1:k);
                    else
                        temp_HET = [ rc_uw_pos_HET{i}.R zeros(1, k-length(rc_uw_pos_HET{i}.R)) ];
                    end
                    if(length(rc_uw_pos_PMG{i}.R) >= k)
                        temp_PMG = rc_uw_pos_PMG{i}.R(1:k);
                    else
                        temp_PMG = [ rc_uw_pos_PMG{i}.R zeros(1, k-length(rc_uw_pos_PMG{i}.R)) ];
                    end                    
                    
                    figure; hold on;
                    plot(1:k, temp_cont, '--', 'LineWidth', 2);
                    plot(1:k, temp_FCD_Type_II, '--', 'LineWidth', 2);
                    plot(1:k, temp_HET, '--', 'LineWidth', 2);
                    plot(1:k, temp_PMG, '--', 'LineWidth', 2);
                    legend('control', 'FCD Type-II', 'HET', 'PMG', 'Location', 'best');
                    
                    color_lines = get(gca, 'ColorOrder');
                    templine = get(gca, 'Children');
                    set(templine(3), 'Color', color_lines(3, :));
                    set(templine(2), 'Color', color_lines(4, :));
                    set(templine(1), 'Color', color_lines(5, :));
                    
                    temp_cont = rc_uw_pos_cont{i}.normR;
                    if(length(rc_uw_pos_FCD_Type_II{i}.normR) >= k)
                        temp_FCD_Type_II = rc_uw_pos_FCD_Type_II{i}.normR(1:k);
                    else
                        temp_FCD_Type_II = [ rc_uw_pos_FCD_Type_II{i}.normR zeros(1, k-length(rc_uw_pos_FCD_Type_II{i}.normR)) ];
                    end
                    if(length(rc_uw_pos_HET{i}.normR) >= k)
                        temp_HET = rc_uw_pos_HET{i}.normR(1:k);
                    else
                        temp_HET = [ rc_uw_pos_HET{i}.normR zeros(1, k-length(rc_uw_pos_HET{i}.normR)) ];
                    end
                    if(length(rc_uw_pos_PMG{i}.normR) >= k)
                        temp_PMG = rc_uw_pos_PMG{i}.normR(1:k);
                    else
                        temp_PMG = [ rc_uw_pos_PMG{i}.normR zeros(1, k-length(rc_uw_pos_PMG{i}.normR)) ];
                    end
                    
                    plot(1:k, temp_cont, 'Color', color_lines(1, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(1, :), 'LineWidth', 2);
                    plot(1:k, temp_FCD_Type_II, 'Color', color_lines(3, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(3, :), 'LineWidth', 2);
                    plot(1:k, temp_HET, 'Color', color_lines(4, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(4, :), 'LineWidth', 2);
                    plot(1:k, temp_PMG, 'Color', color_lines(5, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(5, :), 'LineWidth', 2);
                    plot(0:k+1, ones(1, k+2), ':', 'Color', [0.6 0.6 0.6]);
                    
                end
                xlim([0 21]); ylim([0 2.5]);
                exportfigbo(gcf,[OUTPATH 'rich_club_coefficient2_rsfMRI.png'], 'png', 13); close(gcf);
                
                %% global rich-club coefficient and normalized values
                for i = 7 : -1 :7
                    m = i;
                    k = min([ length(rc_FCD_Type_II{m}.R.R1) length(rc_FCD_Type_II{m}.R.R2) length(rc_HET{m}.R.R2) length(rc_PMG{m}.R.R2) ]);
                    
                    figure; hold on;
                    plot(1:k, max([rc_FCD_Type_II{m}.R.R1(1:k); rc_HET{m}.R.R1(1:k); rc_PMG{m}.R.R1(1:k)]), '--', 'LineWidth', 2);
                    plot(1:k,rc_FCD_Type_II{m}.R.R2(1:k), '--', 'LineWidth', 2);
                    plot(1:k,rc_HET{m}.R.R2(1:k), '--', 'LineWidth', 2);
                    plot(1:k,rc_PMG{m}.R.R2(1:k), '--', 'LineWidth', 2);
                    legend('control', 'FCD Type-II', 'HET', 'PMG', 'Location', 'best');
                    
                    color_lines = get(gca, 'ColorOrder');
                    templine = get(gca, 'Children');
                    set(templine(3), 'Color', color_lines(3, :));
                    set(templine(2), 'Color', color_lines(4, :));
                    set(templine(1), 'Color', color_lines(5, :));
                    plot(1:k, max([ rc_FCD_Type_II{m}.R.normR1(1:k); rc_HET{m}.R.normR1(1:k); rc_PMG{m}.R.normR1(1:k)]), 'Color', color_lines(1, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(1, :), 'LineWidth', 2);
                    plot(1:k,rc_FCD_Type_II{m}.R.normR2(1:k), 'Color', color_lines(3, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(3, :), 'LineWidth', 2);
                    plot(1:k,rc_HET{m}.R.normR2(1:k), 'Color', color_lines(4, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(4, :), 'LineWidth', 2);
                    plot(1:k,rc_PMG{m}.R.normR2(1:k), 'Color', color_lines(5, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(5, :), 'LineWidth', 2);
                    plot(0:k+1, ones(1, k+2), ':', 'Color', [0.6 0.6 0.6]);
                end
                exportfigbo(gcf,[OUTPATH 'rich_club_coefficient1.png'], 'png', 13); close(gcf);
             
                %% Connection density: rich club, feeder, local
                m = 1;
                k_sig=11; % signficant differences between patients and controls
                groups   = { 'cont', 'FCD_Type_II', 'HET', 'PMG', 'FCD_Type_I' };
                groupstr = { '', 'Control', 'FCD-II', 'HET', 'PMG', 'FCD-I', ''};
                
                rc_conn             = [];
                feeder_conn         = [];
                local_conn          = [];
                rc_feeder_conn      = [];
                feeder_local_conn   = [];
                for i = 1 : length(groups)
                    
                    eval(['R_matrix = corrcoef(asso_mat_' groups{i} '_res_val);' ]);
                    
                    wmatrix = threshold_proportional(R_matrix, myspar(m));
                    corrmatrix = double(wmatrix>0);
                    [deg] = degrees_und(corrmatrix);
                    rich_club_nodes = find(deg>=k_sig);
                    non_rich_club_nodes = find(deg<k_sig);
                    
                    % - rich club connections
                    rc_conn = [ rc_conn sum(sum(corrmatrix(rich_club_nodes, rich_club_nodes)))/(length(rich_club_nodes)*(length(rich_club_nodes)-1)) ];
                    % - feeder connections
                    feeder_conn = [ feeder_conn sum(sum(corrmatrix(rich_club_nodes, non_rich_club_nodes)))/(length(rich_club_nodes)*(length(non_rich_club_nodes))) ];
                    % - local connections
                    local_conn = [ local_conn sum(sum(corrmatrix(non_rich_club_nodes, non_rich_club_nodes)))/(length(non_rich_club_nodes)*(length(non_rich_club_nodes)-1)) ];
                    % - rich + feeder connections
                    rc_feeder_conn = [ rc_feeder_conn sum(sum(corrmatrix(rich_club_nodes, [ rich_club_nodes non_rich_club_nodes ])))/(length(rich_club_nodes)*(length([ rich_club_nodes non_rich_club_nodes ])-1)) ];
                    % - feeder + local connections
                    feeder_local_conn = [ feeder_local_conn sum(sum(corrmatrix(non_rich_club_nodes, [ rich_club_nodes non_rich_club_nodes ])))/(length(non_rich_club_nodes)*(length([ rich_club_nodes non_rich_club_nodes ])-1)) ];
                    
                end
                org_conn = [ rc_conn; feeder_conn; local_conn; rc_feeder_conn; feeder_local_conn ]';
                
                p_set = zeros(length(groups)-1, 5);
                for i = 2 : length(groups)
                    
                    eval(['N_sub = [ size(asso_mat_cont_res_val, 1) size(asso_mat_' groups{i} '_res_val, 1) ];']);
                    N1 = N_sub(1); N2 = N_sub(2);
                    eval(['rand_rc = rand_rc_' groups{i} ';']);
                    eval(['tgroup_Reg = [asso_mat_cont_res_val;asso_mat_' groups{i} '_res_val];']);
                    
                    rc_conn_rand1               = zeros(randomization, 1);
                    feeder_conn_rand1           = zeros(randomization, 1);
                    local_conn_rand1            = zeros(randomization, 1);
                    rc_feeder_conn_rand1        = zeros(randomization, 1);
                    feeder_local_conn_rand1     = zeros(randomization, 1);
                    
                    rc_conn_rand2               = zeros(randomization, 1);
                    feeder_conn_rand2           = zeros(randomization, 1);
                    local_conn_rand2            = zeros(randomization, 1);
                    rc_feeder_conn_rand2        = zeros(randomization, 1);
                    feeder_local_conn_rand2     = zeros(randomization, 1);
                    
                    for r = 1 : randomization
                        
                        rp = rand_rc.rp_set(r, :);
                        rerand_corr_1 = corrcoef(tgroup_Reg(rp(1:N1),:));
                        rerand_corr_2 = corrcoef(tgroup_Reg(rp((N1+1):end),:));
                        rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                        rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                        
                        % group 1
                        wmatrix                     = threshold_proportional(rerand_corr_1, myspar(m));
                        corrmatrix                  = double(wmatrix>0);
                        [deg]                       = degrees_und(corrmatrix);
                        rich_club_nodes_rand        = find(deg>=k_sig);
                        non_rich_club_nodes_rand    = find(deg<k_sig);
                        
                        % - rich club connections
                        rc_conn_rand1(r)            = sum(sum(corrmatrix(rich_club_nodes_rand, rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(rich_club_nodes_rand)-1));
                        % - feeder connections
                        feeder_conn_rand1(r)        = sum(sum(corrmatrix(rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)));
                        % - local connections
                        local_conn_rand1(r)         = sum(sum(corrmatrix(non_rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(non_rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)-1));
                        % - rc + feeder connections
                        rc_feeder_conn_rand1(r)     = sum(sum(corrmatrix(rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        % - feeder + local connections
                        feeder_local_conn_rand1(r)  = sum(sum(corrmatrix(non_rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(non_rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        
                        % group 2
                        wmatrix                     = threshold_proportional(rerand_corr_2, myspar(m));
                        corrmatrix                  = double(wmatrix>0);
                        [deg]                       = degrees_und(corrmatrix);
                        rich_club_nodes_rand        = find(deg>=k_sig);
                        non_rich_club_nodes_rand    = find(deg<k_sig);
                        
                        % - rich club connections
                        rc_conn_rand2(r)            = sum(sum(corrmatrix(rich_club_nodes_rand, rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(rich_club_nodes_rand)-1));
                        % - feeder connections
                        feeder_conn_rand2(r)        = sum(sum(corrmatrix(rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)));
                        % - local connections
                        local_conn_rand2(r)         = sum(sum(corrmatrix(non_rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(non_rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)-1));
                        % - rc + feeder connections
                        rc_feeder_conn_rand2(r)     = sum(sum(corrmatrix(rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        % - feeder + local connections
                        feeder_local_conn_rand2(r)  = sum(sum(corrmatrix(non_rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(non_rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        
                    end
                    
                    rand1_features = [ rc_conn_rand1 feeder_conn_rand1 local_conn_rand1 rc_feeder_conn_rand1 feeder_local_conn_rand1 ];
                    rand2_features = [ rc_conn_rand2 feeder_conn_rand2 local_conn_rand2 rc_feeder_conn_rand2 feeder_local_conn_rand2 ];
                    
                    for c = 1 : 5
                        
                        org_diff = org_conn(1, c) - org_conn(i, c);
                        if(org_diff>=0)
                            p_set(i-1, c) = sum(org_diff<(rand1_features(:, c)-rand2_features(:, c)))/randomization;
                        elseif(org_diff<0)
                            p_set(i-1, c) = sum(org_diff>(rand1_features(:, c)-rand2_features(:, c)))/randomization;
                        end
                        
                    end
                    
                end
                
                % bar graph
                for c = 1 : 5
                    
                    figure; hold on;
                    bar(1, org_conn(1, c), 'FaceColor', color_lines(1, :), 'EdgeColor', color_lines(1, :));
                    for i = 2 : length(groups)
                        
                        bar(i, org_conn(i, c), 'FaceColor', color_lines(i, :), 'EdgeColor', color_lines(i, :));
                        
                    end
                    
                    set(gca,'ygrid','on');
                    temp1 = get(gca, 'YLim');
                    temp2 = get(gca, 'YTick');
                    interval = (temp2(2) - temp2(1))*0.5;
                    
                    for i = 2 : length(groups)
                        
                        if(p_set(i-1, c)<=0.0125)
                            
                            line([ 1-0.2, i+0.2 ], [ temp1(2)+(i-1)*interval temp1(2)+(i-1)*interval ], 'Color', 'k');
                            text((1+i)/2, temp1(2)+(i-1)*interval+interval/5, '**')
                            
                        elseif(p_set(i-1, c)>0.0125 & p_set(i-1, c)<0.05)
                            
                            line([ 1-0.2, i+0.2 ], [ temp1(2)+(i-1)*interval temp1(2)+(i-1)*interval ], 'Color', 'k');
                            text((1+i)/2, temp1(2)+(i-1)*interval+interval/5, '*')
                            
                        end
                        
                    end
                    
                    set(gca, 'XTickLabel', groupstr);
                    exportfigbo(gcf,[OUTPATH 'rich_club_connectivity' num2str(c) '.png'], 'png', 13); close(gcf);
                    
                end
                
                %% rich club node representation on the surface maps
                for node_rep = 1
                    
                    % node position of each parcel in AAL
                    for node_position = 1
                        
                        centroid = zeros(1, length(Anatomical_Label_structs_idx));
                        for i = 1 : length(Anatomical_Label_structs_idx)
                            
                            vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                            vert_coord = S.coord(:, vert_idx);
                            [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                            centroid(i) = vert_idx(b);
                            
                        end
                        nodePos = S.coord(:, centroid);
                        
                    end
                    
                    for variable_set = 1
                        
                        m = 1;
                        k_sig=11; % signficant differences between patients and controls
                        colormap_curr = hsv(255);
                        colormap_curr(1, :) = [ 0.8 0.8 0.8 ];
                        link_width = 2;
                        links_visible = 0.85;
                        link_colormap = [ repmat([1 1 0], 85, 1); repmat([1 0.5 0], 85, 1); repmat([1 0 0], 85, 1); ]
                        
                    end
                    
                    % - control
                    for control = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_cont_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        rc_nodes(18) = 1;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        matrix_reorg(18, :) = matrix_reorg(rc_nodes(1), :); matrix_reorg(:, 18) = matrix_reorg(:, rc_nodes(1));
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_control.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - FCD Type-II
                    for FCD_Type_II = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_FCD_Type_II_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_FCD_Type-II.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - HET
                    for HET = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_HET_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_HET.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - PMG
                    for PMG = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_PMG_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_PMG.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - FCD Type-I
                    for FCD_Type_I = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_FCD_Type_I_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_FCD_Type-I.png'], 'png', 13); close(gcf);
                        
                    end
                    
                end
                
            end
            
        end
        
        for micro_AAL = 1
            
            num_of_cont         = length(Codes_cont);
            num_of_FCD_Type_II  = length(Codes_FCD_Type_II);
            num_of_HET          = length(Codes_HET);
            num_of_PMG          = length(Codes_PMG);
            type_mat            = 2;
            
            if(type_mat == 1)
                mat_cont            = zadjacent_matrix_control1_micro;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II1_micro;
                mat_HET             = zadjacent_matrix_HET1_micro;
                mat_PMG             = zadjacent_matrix_PMG1_micro;
            elseif(type_mat == 2)
                mat_cont            = zadjacent_matrix_control2_micro;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II2_micro;
                mat_HET             = zadjacent_matrix_HET2_micro;
                mat_PMG             = zadjacent_matrix_PMG2_micro;
            elseif(type_mat == 3)
                mat_cont            = zadjacent_matrix_control3_micro;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II3_micro;
                mat_HET             = zadjacent_matrix_HET3_micro;
                mat_PMG             = zadjacent_matrix_PMG3_micro;
            elseif(type_mat == 4)
                mat_cont            = zadjacent_matrix_control_micro;
                mat_FCD_Type_II     = zadjacent_matrix_FCD_Type_II_micro;
                mat_HET             = zadjacent_matrix_HET_micro;
                mat_PMG             = zadjacent_matrix_PMG_micro;
            end
            
            parpool(24);
            iternation_num = 1000;
            interval_thres = 0.01;
            randwiring = 10;
            myspar = 0.05:interval_thres:0.4;
            AAL_index_iter = length(Anatomical_Label_structs_idx);
            mat_thres = 0.2;
            
            % 0) binarized averaged matrix construction
            averaged_mat_cont        = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            averaged_mat_FCD_Type_II = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            averaged_mat_HET         = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            averaged_mat_PMG         = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
            
            for i = 1 : length(myspar)
                
                %% control
                mat_cont_temp = mat_cont*0;
                for j = 1 : num_of_cont
                    temp = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    mat_cont_temp(:, :, j) = temp > 0;
%                     mat_cont_temp(:, :, j) = temp;
                end
                averaged_mat_cont(:, :, i) = (sum(mat_cont_temp, 3)/num_of_cont)>mat_thres;
%                 averaged_mat_cont(:, :, i) = mean(mat_cont_temp, 3)>mat_thres;
                
                %% FCD Type II
                mat_FCD_Type_II_temp = mat_FCD_Type_II*0;
                for j = 1 : num_of_FCD_Type_II
                    temp = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    mat_FCD_Type_II_temp(:, :, j) = temp > 0;
%                     mat_FCD_Type_II_temp(:, :, j) = temp;
                end
                averaged_mat_FCD_Type_II(:, :, i) = (sum(mat_FCD_Type_II_temp, 3)/num_of_FCD_Type_II)>mat_thres;
%                 averaged_mat_FCD_Type_II(:, :, i) = mean(mat_FCD_Type_II_temp, 3)>mat_thres;
                
                %% HET
                mat_HET_temp = mat_HET*0;
                for j = 1 : num_of_HET
                    temp = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    mat_HET_temp(:, :, j) = temp > 0;
%                     mat_HET_temp(:, :, j) = temp;
                end
                averaged_mat_HET(:, :, i) = (sum(mat_HET_temp, 3)/num_of_HET)>mat_thres;
%                 averaged_mat_HET(:, :, i) = mean(mat_HET_temp, 3)>mat_thres;
                
                %% PMG
                mat_PMG_temp = mat_PMG*0;
                for j = 1 : num_of_PMG
                    temp = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    mat_PMG_temp(:, :, j) = temp > 0;
%                     mat_PMG_temp(:, :, j) = temp;
                end
                averaged_mat_PMG(:, :, i) = (sum(mat_PMG_temp, 3)/num_of_PMG)>mat_thres;
%                 averaged_mat_PMG(:, :, i) = mean(mat_PMG_temp, 3)>mat_thres;
                
            end
            
            for rc_comp_individual = 1
                
                rc_uw_pos_cont         = cell(length(myspar), num_of_cont);
                rc_uw_pos_FCD_Type_II  = cell(length(myspar), num_of_FCD_Type_II);
                rc_uw_pos_HET          = cell(length(myspar), num_of_HET);
                rc_uw_pos_PMG          = cell(length(myspar), num_of_PMG);
                
                for i = 1 : length(myspar)
                    
                    %% control
                    parfor j = 1 : num_of_cont
                        wmatrix                          = threshold_proportional(mat_cont(:, :, j), myspar(i));
                        rc_uw_pos_cont{i, j}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_cont{i, j}             = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                    %% FCD Type II
                    parfor j = 1 : num_of_FCD_Type_II
                        wmatrix                          = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                        rc_uw_pos_FCD_Type_II{i, j}      = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_FCD_Type_II{i, j}      = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                    %% HET
                    parfor j = 1 : num_of_HET
                        wmatrix                         = threshold_proportional(mat_HET(:, :, j), myspar(i));
                        rc_uw_pos_HET{i, j}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_HET{i, j}             = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                    %% PMG
                    parfor j = 1 : num_of_PMG
                        wmatrix                          = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                        rc_uw_pos_PMG{i, j}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                        %                     rc_uw_pos_PMG{i, j}             = hrich_club_wu(wmatrix, 1, randwiring);
                    end
                    
                end
                
                min_cont        = zeros(length(myspar), 1);
                min_FCD_Type_II = zeros(length(myspar), 1);
                min_HET         = zeros(length(myspar), 1);
                min_PMG         = zeros(length(myspar), 1);
                for i = 1 : length(myspar)
                    
                    min_val = 100;
                    for j = 1 : num_of_cont
                        temp = rc_uw_pos_cont{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_cont(i) = min_val;
                    
                    min_val = 100;
                    for j = 1 : num_of_FCD_Type_II
                        temp = rc_uw_pos_FCD_Type_II{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_FCD_Type_II(i) = min_val;
                    
                    min_val = 100;
                    for j = 1 : num_of_HET
                        temp = rc_uw_pos_HET{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_HET(i) = min_val;
                    
                    min_val = 100;
                    for j = 1 : num_of_PMG
                        temp = rc_uw_pos_PMG{i, j}.normR;
                        if(min_val > sum(temp~=0))
                            min_val = sum(temp~=0);
                        end
                    end
                    min_PMG(i) = min_val;
                    
                end
                
                min_k = min([ min_cont min_FCD_Type_II min_HET min_PMG ], [], 2);
                
                rc_uw_pos_cont_final        = cell(length(myspar), 1);
                rc_uw_pos_FCD_Type_II_final = cell(length(myspar), 1);
                rc_uw_pos_HET_final         = cell(length(myspar), 1);
                rc_uw_pos_PMG_final         = cell(length(myspar), 1);
                
                for i = 1 : length(myspar)
                    
                    for j = 1 : num_of_cont
                        rc_uw_pos_cont_final{i} = [ rc_uw_pos_cont_final{i}; rc_uw_pos_cont{i, j}.normR(1:min_k(i)); ];
                    end
                    for j = 1 : num_of_FCD_Type_II
                        rc_uw_pos_FCD_Type_II_final{i} = [ rc_uw_pos_FCD_Type_II_final{i}; rc_uw_pos_FCD_Type_II{i, j}.normR(1:min_k(i)); ];
                    end
                    for j = 1 : num_of_HET
                        rc_uw_pos_HET_final{i} = [ rc_uw_pos_HET_final{i}; rc_uw_pos_HET{i, j}.normR(1:min_k(i)); ];
                    end
                    for j = 1 : num_of_PMG
                        rc_uw_pos_PMG_final{i} = [ rc_uw_pos_PMG_final{i}; rc_uw_pos_PMG{i, j}.normR(1:min_k(i)); ];
                    end
                    
                end
                
                rc_uw_pos_group = cell(length(myspar), 1);
                
                for i = 1 : length(myspar)
                    
                    rc_uw_pos_group{i} = [ mean(rc_uw_pos_cont_final{i}, 1); mean(rc_uw_pos_FCD_Type_II_final{i}, 1); mean(rc_uw_pos_HET_final{i}, 1); mean(rc_uw_pos_PMG_final{i}, 1) ];
                    
                end
            end
            
            for rc_comp_group_mat = 1
                
                rc_uw_pos_cont         = cell(length(myspar), 1);
                rc_uw_pos_FCD_Type_II  = cell(length(myspar), 1);
                rc_uw_pos_HET          = cell(length(myspar), 1);
                rc_uw_pos_PMG          = cell(length(myspar), 1);
                
                parfor i = 1 : length(myspar)
                    
                    %% control
                    wmatrix                       = threshold_proportional(mean(mat_cont, 3), myspar(i));
                    rc_uw_pos_cont{i}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_cont{i}           = hrich_club_wu(wmatrix, 1, randwiring);

                    %% FCD Type II
                    wmatrix                        = threshold_proportional(mean(mat_FCD_Type_II, 3), myspar(i));
                    rc_uw_pos_FCD_Type_II{i}       = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_FCD_Type_II{i}     = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% HET
                    wmatrix                     = threshold_proportional(mean(mat_HET, 3), myspar(i));
                    rc_uw_pos_HET{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_HET{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% PMG
                    wmatrix                     = threshold_proportional(mean(mat_PMG, 3), myspar(i));
                    rc_uw_pos_PMG{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_PMG{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                end
                
            end
            
            for rc_comp_group_mat_threshold = 1
                
                rc_uw_pos_cont         = cell(length(myspar), 1);
                rc_uw_pos_FCD_Type_II  = cell(length(myspar), 1);
                rc_uw_pos_HET          = cell(length(myspar), 1);
                rc_uw_pos_PMG          = cell(length(myspar), 1);
                
                parfor i = 1 : length(myspar)
                    
                    %% control
                    wmatrix                       = averaged_mat_cont(:, :, i);
                    rc_uw_pos_cont{i}             = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_cont{i}           = hrich_club_wu(wmatrix, 1, randwiring);

                    %% FCD Type II
                    wmatrix                        = averaged_mat_FCD_Type_II(:, :, i);
                    rc_uw_pos_FCD_Type_II{i}       = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_FCD_Type_II{i}     = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% HET
                    wmatrix                        = averaged_mat_HET(:, :, i);
                    rc_uw_pos_HET{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_HET{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                    %% PMG
                    wmatrix                        = averaged_mat_PMG(:, :, i);
                    rc_uw_pos_PMG{i}            = hrich_club_bu(double(wmatrix>0), 1, randwiring);
                    % rc_uw_pos_PMG{i, j}          = hrich_club_wu(wmatrix, 1, randwiring);
                    
                end
                
                randomization = 100;
                brandnetwork = 1;
                interval_thres = 0.01;
                myspar = 0.05:interval_thres:0.36;
                option.bagging = 1;
                
                [ rc_FCD_Type_II rand_rc_FCD_Type_II ] = rcdiff_sparsity_parfor_rsfMRI(mat_cont, mat_FCD_Type_II, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);  figure; pause(5);
                [ rc_HET         rand_rc_HET         ] = rcdiff_sparsity_parfor_rsfMRI(mat_cont, mat_HET, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);          figure; pause(5);
                [ rc_PMG         rand_rc_PMG         ] = rcdiff_sparsity_parfor_rsfMRI(mat_cont, mat_PMG, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);          figure; pause(5);
                
            end
            
            if(visualization == 1)
                
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/03_rich_club_analysis/';
                
                for i = 9: -1 : 9
                    
                    k = length(rc_uw_pos_cont{i}.normR);
                    
                    temp_cont = rc_uw_pos_cont{i}.R;
                    if(length(rc_uw_pos_FCD_Type_II{i}.R) >= k)
                        temp_FCD_Type_II = rc_uw_pos_FCD_Type_II{i}.R(1:k);
                    else
                        temp_FCD_Type_II = [ rc_uw_pos_FCD_Type_II{i}.R zeros(1, k-length(rc_uw_pos_FCD_Type_II{i}.R)) ];
                    end                        
                    if(length(rc_uw_pos_HET{i}.R) >= k)
                        temp_HET = rc_uw_pos_HET{i}.R(1:k);
                    else
                        temp_HET = [ rc_uw_pos_HET{i}.R zeros(1, k-length(rc_uw_pos_HET{i}.R)) ];
                    end
                    if(length(rc_uw_pos_PMG{i}.R) >= k)
                        temp_PMG = rc_uw_pos_PMG{i}.R(1:k);
                    else
                        temp_PMG = [ rc_uw_pos_PMG{i}.R zeros(1, k-length(rc_uw_pos_PMG{i}.R)) ];
                    end                    
                    
                    figure; hold on;
                    plot(1:k, temp_cont, '--', 'LineWidth', 2);
                    plot(1:k, temp_FCD_Type_II, '--', 'LineWidth', 2);
                    plot(1:k, temp_HET, '--', 'LineWidth', 2);
                    plot(1:k, temp_PMG, '--', 'LineWidth', 2);
                    legend('control', 'FCD Type-II', 'HET', 'PMG', 'Location', 'best');
                    
                    color_lines = get(gca, 'ColorOrder');
                    templine = get(gca, 'Children');
                    set(templine(3), 'Color', color_lines(3, :));
                    set(templine(2), 'Color', color_lines(4, :));
                    set(templine(1), 'Color', color_lines(5, :));
                    
                    temp_cont = rc_uw_pos_cont{i}.normR;
                    if(length(rc_uw_pos_FCD_Type_II{i}.normR) >= k)
                        temp_FCD_Type_II = rc_uw_pos_FCD_Type_II{i}.normR(1:k);
                    else
                        temp_FCD_Type_II = [ rc_uw_pos_FCD_Type_II{i}.normR zeros(1, k-length(rc_uw_pos_FCD_Type_II{i}.normR)) ];
                    end
                    if(length(rc_uw_pos_HET{i}.normR) >= k)
                        temp_HET = rc_uw_pos_HET{i}.normR(1:k);
                    else
                        temp_HET = [ rc_uw_pos_HET{i}.normR zeros(1, k-length(rc_uw_pos_HET{i}.normR)) ];
                    end
                    if(length(rc_uw_pos_PMG{i}.normR) >= k)
                        temp_PMG = rc_uw_pos_PMG{i}.normR(1:k);
                    else
                        temp_PMG = [ rc_uw_pos_PMG{i}.normR zeros(1, k-length(rc_uw_pos_PMG{i}.normR)) ];
                    end
                    
                    plot(1:k, temp_cont, 'Color', color_lines(1, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(1, :), 'LineWidth', 2);
                    plot(1:k, temp_FCD_Type_II, 'Color', color_lines(3, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(3, :), 'LineWidth', 2);
                    plot(1:k, temp_HET, 'Color', color_lines(4, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(4, :), 'LineWidth', 2);
                    plot(1:k, temp_PMG, 'Color', color_lines(5, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(5, :), 'LineWidth', 2);
                    plot(0:k+1, ones(1, k+2), ':', 'Color', [0.6 0.6 0.6]);
                    
                end
                xlim([0 21]); ylim([0 2.5]);
                exportfigbo(gcf,[OUTPATH 'rich_club_coefficient2_rsfMRI.png'], 'png', 13); close(gcf);
                
                %% global rich-club coefficient and normalized values
                for i = 4 : 1 : 4
                    m = i;
                    k = min([ length(rc_FCD_Type_II{m}.R.R1) length(rc_FCD_Type_II{m}.R.R2) length(rc_HET{m}.R.R2) length(rc_PMG{m}.R.R2) ]);
                    
                    figure; hold on;
                    plot(1:k, max([rc_FCD_Type_II{m}.R.R1(1:k); rc_HET{m}.R.R1(1:k); rc_PMG{m}.R.R1(1:k)]), '--', 'LineWidth', 2);
                    plot(1:k,rc_FCD_Type_II{m}.R.R2(1:k), '--', 'LineWidth', 2);
                    plot(1:k,rc_HET{m}.R.R2(1:k), '--', 'LineWidth', 2);
                    plot(1:k,rc_PMG{m}.R.R2(1:k), '--', 'LineWidth', 2);
                    legend('control', 'FCD Type-II', 'HET', 'PMG', 'Location', 'best');
                    
                    color_lines = get(gca, 'ColorOrder');
                    templine = get(gca, 'Children');
                    set(templine(3), 'Color', color_lines(3, :));
                    set(templine(2), 'Color', color_lines(4, :));
                    set(templine(1), 'Color', color_lines(5, :));
                    plot(1:k,rc_HET{m}.R.normR2(1:k), 'Color', color_lines(4, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(4, :), 'LineWidth', 2);
                    plot(1:k,rc_PMG{m}.R.normR2(1:k), 'Color', color_lines(5, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(5, :), 'LineWidth', 2);
                    plot(1:k,rc_FCD_Type_II{m}.R.normR2(1:k), 'Color', color_lines(3, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(3, :), 'LineWidth', 2);
                    plot(1:k, max([ rc_FCD_Type_II{m}.R.normR1(1:k); rc_HET{m}.R.normR1(1:k); rc_PMG{m}.R.normR1(1:k)]), 'Color', color_lines(1, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(1, :), 'LineWidth', 2);
                    plot(0:k+1, ones(1, k+2), ':', 'Color', [0.6 0.6 0.6]);
                    
                    FDR_FCD_TypeII  = FDR(rc_FCD_Type_II{m}.R.p_delta_normR, 0.05);
                    FDR_HET         = FDR(rc_HET{m}.R.p_delta_normR, 0.05);
                    FDR_PMG         = FDR(rc_PMG{m}.R.p_delta_normR, 0.05);
                    
                    plot(find(rc_FCD_Type_II{m}.R.p_delta_normR<=FDR_FCD_TypeII & (rc_FCD_Type_II{m}.R.normR1-rc_FCD_Type_II{m}.R.normR2)>0.1), 4.5, '*', 'color', color_lines(3, :));
                    plot(find(rc_FCD_Type_II{m}.R.p_delta_normR<=0.05 & rc_FCD_Type_II{m}.R.p_delta_normR>FDR_FCD_TypeII & rc_FCD_Type_II{m}.R.normR1>rc_FCD_Type_II{m}.R.normR2), 4.5, 'o', 'color', color_lines(3, :));
                    
                    plot(find(rc_HET{m}.R.p_delta_normR<=FDR_HET & (rc_HET{m}.R.normR1-rc_HET{m}.R.normR2)>0.1), 4.7, '*', 'color', color_lines(4, :));
                    plot(find(rc_HET{m}.R.p_delta_normR<=0.05 & rc_HET{m}.R.p_delta_normR>FDR_HET & rc_HET{m}.R.normR1>rc_HET{m}.R.normR2), 4.7, 'o', 'color', color_lines(4, :));
                    
                    plot(find(rc_PMG{m}.R.p_delta_normR<=FDR_PMG & (rc_PMG{m}.R.normR1-rc_PMG{m}.R.normR2)>0.1), 4.9, '*', 'color', color_lines(5, :));
                    plot(find(rc_PMG{m}.R.p_delta_normR<=0.09 & rc_PMG{m}.R.p_delta_normR>FDR_PMG & rc_PMG{m}.R.normR1>rc_PMG{m}.R.normR2), 4.9, 'o', 'color', color_lines(5, :));
                    
                end
                exportfigbo(gcf,[OUTPATH 'rich_club_coefficient1_rsfMRI_micro.png'], 'png', 13); close(gcf);
             
                %% Connection density: rich club, feeder, local
                m = 1;
                k_sig=11; % signficant differences between patients and controls
                groups   = { 'cont', 'FCD_Type_II', 'HET', 'PMG', 'FCD_Type_I' };
                groupstr = { '', 'Control', 'FCD-II', 'HET', 'PMG', 'FCD-I', ''};
                
                rc_conn             = [];
                feeder_conn         = [];
                local_conn          = [];
                rc_feeder_conn      = [];
                feeder_local_conn   = [];
                for i = 1 : length(groups)
                    
                    eval(['R_matrix = corrcoef(asso_mat_' groups{i} '_res_val);' ]);
                    
                    wmatrix = threshold_proportional(R_matrix, myspar(m));
                    corrmatrix = double(wmatrix>0);
                    [deg] = degrees_und(corrmatrix);
                    rich_club_nodes = find(deg>=k_sig);
                    non_rich_club_nodes = find(deg<k_sig);
                    
                    % - rich club connections
                    rc_conn = [ rc_conn sum(sum(corrmatrix(rich_club_nodes, rich_club_nodes)))/(length(rich_club_nodes)*(length(rich_club_nodes)-1)) ];
                    % - feeder connections
                    feeder_conn = [ feeder_conn sum(sum(corrmatrix(rich_club_nodes, non_rich_club_nodes)))/(length(rich_club_nodes)*(length(non_rich_club_nodes))) ];
                    % - local connections
                    local_conn = [ local_conn sum(sum(corrmatrix(non_rich_club_nodes, non_rich_club_nodes)))/(length(non_rich_club_nodes)*(length(non_rich_club_nodes)-1)) ];
                    % - rich + feeder connections
                    rc_feeder_conn = [ rc_feeder_conn sum(sum(corrmatrix(rich_club_nodes, [ rich_club_nodes non_rich_club_nodes ])))/(length(rich_club_nodes)*(length([ rich_club_nodes non_rich_club_nodes ])-1)) ];
                    % - feeder + local connections
                    feeder_local_conn = [ feeder_local_conn sum(sum(corrmatrix(non_rich_club_nodes, [ rich_club_nodes non_rich_club_nodes ])))/(length(non_rich_club_nodes)*(length([ rich_club_nodes non_rich_club_nodes ])-1)) ];
                    
                end
                org_conn = [ rc_conn; feeder_conn; local_conn; rc_feeder_conn; feeder_local_conn ]';
                
                p_set = zeros(length(groups)-1, 5);
                for i = 2 : length(groups)
                    
                    eval(['N_sub = [ size(asso_mat_cont_res_val, 1) size(asso_mat_' groups{i} '_res_val, 1) ];']);
                    N1 = N_sub(1); N2 = N_sub(2);
                    eval(['rand_rc = rand_rc_' groups{i} ';']);
                    eval(['tgroup_Reg = [asso_mat_cont_res_val;asso_mat_' groups{i} '_res_val];']);
                    
                    rc_conn_rand1               = zeros(randomization, 1);
                    feeder_conn_rand1           = zeros(randomization, 1);
                    local_conn_rand1            = zeros(randomization, 1);
                    rc_feeder_conn_rand1        = zeros(randomization, 1);
                    feeder_local_conn_rand1     = zeros(randomization, 1);
                    
                    rc_conn_rand2               = zeros(randomization, 1);
                    feeder_conn_rand2           = zeros(randomization, 1);
                    local_conn_rand2            = zeros(randomization, 1);
                    rc_feeder_conn_rand2        = zeros(randomization, 1);
                    feeder_local_conn_rand2     = zeros(randomization, 1);
                    
                    for r = 1 : randomization
                        
                        rp = rand_rc.rp_set(r, :);
                        rerand_corr_1 = corrcoef(tgroup_Reg(rp(1:N1),:));
                        rerand_corr_2 = corrcoef(tgroup_Reg(rp((N1+1):end),:));
                        rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                        rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                        
                        % group 1
                        wmatrix                     = threshold_proportional(rerand_corr_1, myspar(m));
                        corrmatrix                  = double(wmatrix>0);
                        [deg]                       = degrees_und(corrmatrix);
                        rich_club_nodes_rand        = find(deg>=k_sig);
                        non_rich_club_nodes_rand    = find(deg<k_sig);
                        
                        % - rich club connections
                        rc_conn_rand1(r)            = sum(sum(corrmatrix(rich_club_nodes_rand, rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(rich_club_nodes_rand)-1));
                        % - feeder connections
                        feeder_conn_rand1(r)        = sum(sum(corrmatrix(rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)));
                        % - local connections
                        local_conn_rand1(r)         = sum(sum(corrmatrix(non_rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(non_rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)-1));
                        % - rc + feeder connections
                        rc_feeder_conn_rand1(r)     = sum(sum(corrmatrix(rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        % - feeder + local connections
                        feeder_local_conn_rand1(r)  = sum(sum(corrmatrix(non_rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(non_rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        
                        % group 2
                        wmatrix                     = threshold_proportional(rerand_corr_2, myspar(m));
                        corrmatrix                  = double(wmatrix>0);
                        [deg]                       = degrees_und(corrmatrix);
                        rich_club_nodes_rand        = find(deg>=k_sig);
                        non_rich_club_nodes_rand    = find(deg<k_sig);
                        
                        % - rich club connections
                        rc_conn_rand2(r)            = sum(sum(corrmatrix(rich_club_nodes_rand, rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(rich_club_nodes_rand)-1));
                        % - feeder connections
                        feeder_conn_rand2(r)        = sum(sum(corrmatrix(rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)));
                        % - local connections
                        local_conn_rand2(r)         = sum(sum(corrmatrix(non_rich_club_nodes_rand, non_rich_club_nodes_rand)))/(length(non_rich_club_nodes_rand)*(length(non_rich_club_nodes_rand)-1));
                        % - rc + feeder connections
                        rc_feeder_conn_rand2(r)     = sum(sum(corrmatrix(rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        % - feeder + local connections
                        feeder_local_conn_rand2(r)  = sum(sum(corrmatrix(non_rich_club_nodes_rand, [ rich_club_nodes_rand non_rich_club_nodes_rand ])))/(length(non_rich_club_nodes_rand)*(length([ rich_club_nodes_rand non_rich_club_nodes_rand ])-1));
                        
                    end
                    
                    rand1_features = [ rc_conn_rand1 feeder_conn_rand1 local_conn_rand1 rc_feeder_conn_rand1 feeder_local_conn_rand1 ];
                    rand2_features = [ rc_conn_rand2 feeder_conn_rand2 local_conn_rand2 rc_feeder_conn_rand2 feeder_local_conn_rand2 ];
                    
                    for c = 1 : 5
                        
                        org_diff = org_conn(1, c) - org_conn(i, c);
                        if(org_diff>=0)
                            p_set(i-1, c) = sum(org_diff<(rand1_features(:, c)-rand2_features(:, c)))/randomization;
                        elseif(org_diff<0)
                            p_set(i-1, c) = sum(org_diff>(rand1_features(:, c)-rand2_features(:, c)))/randomization;
                        end
                        
                    end
                    
                end
                
                % bar graph
                for c = 1 : 5
                    
                    figure; hold on;
                    bar(1, org_conn(1, c), 'FaceColor', color_lines(1, :), 'EdgeColor', color_lines(1, :));
                    for i = 2 : length(groups)
                        
                        bar(i, org_conn(i, c), 'FaceColor', color_lines(i, :), 'EdgeColor', color_lines(i, :));
                        
                    end
                    
                    set(gca,'ygrid','on');
                    temp1 = get(gca, 'YLim');
                    temp2 = get(gca, 'YTick');
                    interval = (temp2(2) - temp2(1))*0.5;
                    
                    for i = 2 : length(groups)
                        
                        if(p_set(i-1, c)<=0.0125)
                            
                            line([ 1-0.2, i+0.2 ], [ temp1(2)+(i-1)*interval temp1(2)+(i-1)*interval ], 'Color', 'k');
                            text((1+i)/2, temp1(2)+(i-1)*interval+interval/5, '**')
                            
                        elseif(p_set(i-1, c)>0.0125 & p_set(i-1, c)<0.05)
                            
                            line([ 1-0.2, i+0.2 ], [ temp1(2)+(i-1)*interval temp1(2)+(i-1)*interval ], 'Color', 'k');
                            text((1+i)/2, temp1(2)+(i-1)*interval+interval/5, '*')
                            
                        end
                        
                    end
                    
                    set(gca, 'XTickLabel', groupstr);
                    exportfigbo(gcf,[OUTPATH 'rich_club_connectivity' num2str(c) '.png'], 'png', 13); close(gcf);
                    
                end
                
                %% rich club node representation on the surface maps
                for node_rep = 1
                    
                    % node position of each parcel in AAL
                    for node_position = 1
                        
                        centroid = zeros(1, length(Anatomical_Label_structs_idx));
                        for i = 1 : length(Anatomical_Label_structs_idx)
                            
                            vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                            vert_coord = S.coord(:, vert_idx);
                            [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                            centroid(i) = vert_idx(b);
                            
                        end
                        nodePos = S.coord(:, centroid);
                        
                    end
                    
                    for variable_set = 1
                        
                        m = 1;
                        k_sig=11; % signficant differences between patients and controls
                        colormap_curr = hsv(255);
                        colormap_curr(1, :) = [ 0.8 0.8 0.8 ];
                        link_width = 2;
                        links_visible = 0.85;
                        link_colormap = [ repmat([1 1 0], 85, 1); repmat([1 0.5 0], 85, 1); repmat([1 0 0], 85, 1); ]
                        
                    end
                    
                    % - control
                    for control = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_cont_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        rc_nodes(18) = 1;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        matrix_reorg(18, :) = matrix_reorg(rc_nodes(1), :); matrix_reorg(:, 18) = matrix_reorg(:, rc_nodes(1));
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_control.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - FCD Type-II
                    for FCD_Type_II = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_FCD_Type_II_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_FCD_Type-II.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - HET
                    for HET = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_HET_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_HET.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - PMG
                    for PMG = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_PMG_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_PMG.png'], 'png', 13); close(gcf);
                        
                    end
                    
                    % - FCD Type-I
                    for FCD_Type_I = 1
                        
                        R_matrix   = corrcoef(temp.asso_mat_FCD_Type_I_res_val);
                        wmatrix = threshold_proportional(R_matrix, myspar(m));
                        corrmatrix = double(wmatrix>0);
                        [deg] = degrees_und(corrmatrix);
                        rc_nodes = deg>=k_sig;
                        Anatomical_Label_structs_idx(find(rc_nodes))
                        NodeSize = 3*(rc_nodes) + (rc_nodes==0);
                        
                        matrix_reorg = zeros(size(corrmatrix));
                        % connectivity reorganization
                        for i = 1 : num_of_ROIs
                            
                            for j = 1 : num_of_ROIs
                                
                                if(corrmatrix(i, j) == 1)
                                    if(rc_nodes(i) == 1)
                                        if(rc_nodes(j) == 1)
                                            % rich club connections
                                            matrix_reorg(i, j) = 3;
                                        else
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        end
                                    else
                                        if(rc_nodes(j) == 1)
                                            % feeder connections
                                            matrix_reorg(i, j) = 2;
                                        else
                                            % peripheral connections
                                            matrix_reorg(i, j) = 1;
                                        end
                                    end
                                end
                                
                            end
                            
                        end
                        
                        figure;
                        [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf_RichClub(NodeSize, corrmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'node_color_by', 'nodalMeasure', 'node_radius', 3.5, 'Node_colormap', colormap_curr);
                        NOELNetVisuNetConnections_RichClub(aa1, cb1, matrix_reorg, (1:78), nodePos, 'link_colormap', link_colormap, 'link_width', link_width, 'link_width_weight', 'yes', 'links_visible', links_visible)
                        exportfigbo(gcf,[OUTPATH 'rich_club_FCD_Type-I.png'], 'png', 13); close(gcf);
                        
                    end
                    
                end
                
            end
            
        end
        
    end  
    
end