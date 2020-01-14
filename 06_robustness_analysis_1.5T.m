clear all
close all

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath('/host/gypsy/local_raid/seokjun/03_downloads/GRETNA/NetFunctions/');
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
%        load( [ MATDIR '/01_generate_association_matrix_MCD_maskingout_lesion.mat' ] );        
%        load( [ MATDIR '/05_hierarchical_analysis_result_MCD_maskingout_lesion.mat' ] );
    else
        load( [ MATDIR '/01_generate_association_matrix_MCD1.mat' ] );
        load( [ MATDIR '/06_robustness_analysis_result_MCD.mat'   ] );
    end
    
end

%% compute hierarchy coefficient
for compute_hierarchy_coeff = 1
            
    % since the Louvain algorithm is heuristic, every times it runs, it gives slighltly different values
    % so here we ran 100 times and extracted consistent modular structures among iterations.   
    num_of_ROIs = length(Anatomical_Label_structs_idx);
            
    % The moddiff_sparsity function estimates modularity and its significance with respect to random network
    for compute_significance = 1
        
        parpool(24);
        
        outputdir = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/05_robustness_analysis';
        randomization = 1000;
        normalized = 1;
        interval_thres = 0.01;
        myspar = 0.21:interval_thres:0.40;
        option = [0.04]; % FCD Type-I: [0.04, 0.04]
        
        [ robust_FCD_Type_I   robust_rand_FCD_Typ_I  ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_FCD_Type_I_res_val,  myspar(1), myspar(end), interval_thres, randomization, 'bc', 'node', outputdir, 'FCD_Type_I', normalized);
        [ robust_FCD_Type_II  robust_rand_FCD_Typ_II ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, randomization, 'bc', 'node', outputdir, 'FCD_Type_II', normalized);
        [ robust_HET          robust_rand_HET        ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_HET_res_val,         myspar(1), myspar(end), interval_thres, randomization, 'bc', 'node', outputdir, 'HET', normalized);
        [ robust_PMG          robust_rand_PMG        ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_PMG_res_val,         myspar(1), myspar(end), interval_thres, randomization, 'bc', 'node', outputdir, 'PMG', normalized);
                       
        [ robust_FCD_Type_I_edge   robust_rand_FCD_Typ_I_edge  ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_FCD_Type_I_res_val,  myspar(1), myspar(end), interval_thres, randomization, 'bc', 'edge', outputdir, 'FCD_Type_I', normalized);
        [ robust_FCD_Type_II_edge  robust_rand_FCD_Typ_II_edge ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, randomization, 'bc', 'edge', outputdir, 'FCD_Type_II', normalized);
        [ robust_HET_edge          robust_rand_HET_edge        ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_HET_res_val,         myspar(1), myspar(end), interval_thres, randomization, 'bc', 'edge', outputdir, 'HET', normalized);
        [ robust_PMG_edge          robust_rand_PMG_edge        ] = robustdiff_sparsity_parfor_sonic(asso_mat_cont_res_val, asso_mat_PMG_res_val,         myspar(1), myspar(end), interval_thres, randomization, 'bc', 'edge', outputdir, 'PMG', normalized);
        
        if(visualization == 1)
            
            linewidth = 2;
            
            for node = 1
                
                %% target attack
                rc_k = 11;
                m = 6;                
                variable_sets = { 'FCD_Type_I', 'FCD_Type_II', 'HET', 'PMG' };
                
                for i = 1 : length(variable_sets)
                    
                    temp = load([ outputdir '/robust_result_' num2str(myspar(m)) '_' variable_sets{i} '_node.mat']);
                    eval(['robust_' variable_sets{i} '.Ge1{m}                        = temp.robust_curr.Ge1;']);
                    eval(['robust_' variable_sets{i} '.Ge2{m}                        = temp.robust_curr.Ge2;']);
                    eval(['robust_' variable_sets{i} '.deltaGe{m}                    = temp.robust_curr.deltaGe;']);
                    eval(['robust_' variable_sets{i} '.Ge12_diff_p{m}                = temp.robust_curr.Ge12_diff_p;']);
                    eval(['robust_' variable_sets{i} '.aGe12_diff_p.targetattack(m)  = temp.robust_curr.aGe12_diff_p.targetattack;']);
                    eval(['robust_' variable_sets{i} '.aGe12_diff_p.randomerror(m)   = temp.robust_curr.aGe12_diff_p.randomerror;']);
                    
                    eval(['robust_rand_' variable_sets{i} '.rp                        = temp.rerandRobust_curr.rp;']);
                    eval(['robust_rand_' variable_sets{i} '.delta{m}                  = temp.rerandRobust_curr.deltaGe;']);
                    eval(['robust_rand_' variable_sets{i} '.lower95TA12diff{m}        = temp.rerandRobust_curr.lower95TA12diff;']);
                    eval(['robust_rand_' variable_sets{i} '.upper95TA12diff{m}        = temp.rerandRobust_curr.upper95TA12diff;']);
                    eval(['robust_rand_' variable_sets{i} '.rerand_TA12diff_mean{m}   = temp.rerandRobust_curr.rerand_TA12diff_mean;']);
                    eval(['robust_rand_' variable_sets{i} '.lower95RE12diff{m}        = temp.rerandRobust_curr.lower95RE12diff;']);
                    eval(['robust_rand_' variable_sets{i} '.upper95RE12diff{m}        = temp.rerandRobust_curr.upper95RE12diff;']);
                    eval(['robust_rand_' variable_sets{i} '.rerand_RE12diff_mean{m}   = temp.rerandRobust_curr.rerand_RE12diff_mean;']);
                    
                end
                
                k = min([ length(robust_FCD_Type_I.Ge1{m}.targetattack) length(robust_FCD_Type_II.Ge1{m}.targetattack) length(robust_HET.Ge1{m}.targetattack) length(robust_PMG.Ge1{m}.targetattack) ]);
                figure; hold on;
                plot(1:k, robust_FCD_Type_I.Ge1{m}.targetattack(1:k)/robust_FCD_Type_I.Ge1{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_I.Ge2{m}.targetattack(1:k)/robust_FCD_Type_I.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_II.Ge2{m}.targetattack(1:k)/robust_FCD_Type_II.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_HET.Ge2{m}.targetattack(1:k)/robust_HET.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_PMG.Ge2{m}.targetattack(1:k)/robust_PMG.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                temp = get(gca, 'Children');
                set(temp(4), 'Visible', 'off');
                color_lines = get(gca, 'ColorOrder');
                
%                 Q = FDR(robust_FCD_Type_I.Ge12_diff_p{m}.targetattack, 0.05);
%                 if(isempty(Q))
%                     Q = 0;
%                 end
%                 uncorr_sig = find(robust_FCD_Type_I.Ge12_diff_p{m}.targetattack<=0.05 & robust_FCD_Type_I.Ge12_diff_p{m}.targetattack>Q);
%                 FDR_sig    = find(robust_FCD_Type_I.Ge12_diff_p{m}.targetattack<=Q);
%                 scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.18, 'V', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
%                 scatter(FDR_sig, ones(1, length(FDR_sig))*1.18, '*', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
%                 scatter(find(rc_FCD_Type_I), robust_FCD_Type_I.Ge2{m}.targetattack(rc_FCD_Type_I)/robust_FCD_Type_I.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                
                Q = FDR(robust_FCD_Type_II.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_II.Ge12_diff_p{m}.targetattack<=0.05 & robust_FCD_Type_II.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_FCD_Type_II.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.03, 'V', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.03, '*', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                %             scatter(find(rc_FCD_Type_II), robust_FCD_Type_II.Ge2{m}.targetattack(rc_FCD_Type_II)/robust_FCD_Type_II.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_cont), robust_FCD_Type_I.Ge1{m}.targetattack(rc_cont)/robust_FCD_Type_I.Ge1{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                                
                Q = FDR(robust_HET.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_HET.Ge12_diff_p{m}.targetattack<=0.05 & robust_HET.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_HET.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.08, 'V', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.08, '*', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                %             scatter(find(rc_HET), robust_HET.Ge2{m}.targetattack(rc_HET)/robust_HET.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                
                Q = FDR(robust_PMG.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_PMG.Ge12_diff_p{m}.targetattack<=0.05 & robust_PMG.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_PMG.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.13, 'V', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.13, '*', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                %             scatter(find(rc_PMG), robust_PMG.Ge2{m}.targetattack(rc_PMG)/robust_PMG.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                                
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                exportfigbo(gcf,[outputdir '/node_targetattack2.png'], 'png', 13); close(gcf);
                
                %% random error
                m = 1;
                k = min([ length(robust_FCD_Type_I.Ge1{m}.randomerror) length(robust_FCD_Type_II.Ge1{m}.randomerror) length(robust_HET.Ge1{m}.randomerror) length(robust_PMG.Ge1{m}.randomerror) ]);
                
                figure; hold on;
                plot(1:k, robust_FCD_Type_I.Ge1{m}.randomerror(1:k)/robust_FCD_Type_I.Ge1{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_I.Ge2{m}.randomerror(1:k)/robust_FCD_Type_I.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_II.Ge2{m}.randomerror(1:k)/robust_FCD_Type_II.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_HET.Ge2{m}.randomerror(1:k)/robust_HET.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_PMG.Ge2{m}.randomerror(1:k)/robust_PMG.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                
                temp = get(gca, 'Children');
                set(temp(4), 'Visible', 'off');
                color_lines = get(gca, 'ColorOrder');

                Q = FDR(robust_FCD_Type_I.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_I.Ge12_diff_p{m}.randomerror<=0.05 & robust_FCD_Type_I.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_FCD_Type_I.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.18, 'V', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.18, '*', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_FCD_Type_I), robust_FCD_Type_I.Ge2{m}.randomerror(rc_FCD_Type_I)/robust_FCD_Type_I.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                
                Q = FDR(robust_FCD_Type_II.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_II.Ge12_diff_p{m}.randomerror<=0.05 & robust_FCD_Type_II.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_FCD_Type_II.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.03, 'V', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.03, '*', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                %             scatter(find(rc_FCD_Type_II), robust_FCD_Type_II.Ge2{m}.randomerror(rc_FCD_Type_II)/robust_FCD_Type_II.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_cont), robust_FCD_Type_I.Ge1{m}.randomerror(rc_cont)/robust_FCD_Type_I.Ge1{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                
                Q = FDR(robust_HET.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_HET.Ge12_diff_p{m}.randomerror<=0.05 & robust_HET.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_HET.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.08, 'V', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.08, '*', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                %             scatter(find(rc_HET), robust_HET.Ge2{m}.randomerror(rc_HET)/robust_HET.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                
                Q = FDR(robust_PMG.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_PMG.Ge12_diff_p{m}.randomerror<=0.05 & robust_PMG.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_PMG.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.13, 'V', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.13, '*', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                %             scatter(find(rc_PMG), robust_PMG.Ge2{m}.randomerror(rc_PMG)/robust_PMG.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                
                
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                exportfigbo(gcf,[outputdir '/node_randomerror2.png'], 'png', 13); close(gcf);
                
            end
            
            for edge = 1
                
                %% target attack
                rc_k = 11;
                m = 2;
                
                variable_sets = { 'FCD_Type_I', 'FCD_Type_II', 'HET', 'PMG' };
                
                for i = 1 : length(variable_sets)
                    
                    temp = load([ outputdir '/robust_result_' num2str(myspar(m)) '_' variable_sets{i} '_edge.mat']);
                    eval(['robust_' variable_sets{i} '_edge.Ge1{m}                        = temp.robust_curr.Ge1;']);
                    eval(['robust_' variable_sets{i} '_edge.Ge2{m}                        = temp.robust_curr.Ge2;']);
                    eval(['robust_' variable_sets{i} '_edge.deltaGe{m}                    = temp.robust_curr.deltaGe;']);
                    eval(['robust_' variable_sets{i} '_edge.Ge12_diff_p{m}                = temp.robust_curr.Ge12_diff_p;']);
                    eval(['robust_' variable_sets{i} '_edge.aGe12_diff_p.targetattack(m)  = temp.robust_curr.aGe12_diff_p.targetattack;']);
                    eval(['robust_' variable_sets{i} '_edge.aGe12_diff_p.randomerror(m)   = temp.robust_curr.aGe12_diff_p.randomerror;']);
                    
                    eval(['robust_rand_' variable_sets{i} '_edge.rp                        = temp.rerandRobust_curr.rp;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.delta{m}                  = temp.rerandRobust_curr.deltaGe;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.lower95TA12diff{m}        = temp.rerandRobust_curr.lower95TA12diff;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.upper95TA12diff{m}        = temp.rerandRobust_curr.upper95TA12diff;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.rerand_TA12diff_mean{m}   = temp.rerandRobust_curr.rerand_TA12diff_mean;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.lower95RE12diff{m}        = temp.rerandRobust_curr.lower95RE12diff;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.upper95RE12diff{m}        = temp.rerandRobust_curr.upper95RE12diff;']);
                    eval(['robust_rand_' variable_sets{i} '_edge.rerand_RE12diff_mean{m}   = temp.rerandRobust_curr.rerand_RE12diff_mean;']);
                    
                end
                
                k = min([ length(robust_FCD_Type_I_edge.Ge1{m}.targetattack) length(robust_FCD_Type_II_edge.Ge1{m}.targetattack) length(robust_HET_edge.Ge1{m}.targetattack) length(robust_PMG_edge.Ge1{m}.targetattack) ]);
                                
                figure; hold on;
                plot(1:k, robust_FCD_Type_I_edge.Ge1{m}.targetattack(1:k)/robust_FCD_Type_I_edge.Ge1{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_I_edge.Ge2{m}.targetattack(1:k)/robust_FCD_Type_I_edge.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_II_edge.Ge2{m}.targetattack(1:k)/robust_FCD_Type_II_edge.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_HET_edge.Ge2{m}.targetattack(1:k)/robust_HET_edge.Ge2{m}.targetattack(1), 'LineWidth', linewidth);
                plot(1:k, robust_PMG_edge.Ge2{m}.targetattack(1:k)/robust_PMG_edge.Ge2{m}.targetattack(1), 'LineWidth', linewidth);                
                temp = get(gca, 'Children');
                set(temp(4), 'Visible', 'off');
                
%                 figure; hold on;
%                 plot(1:k, robust_FCD_Type_I_edge.Ge1{m}.targetattack(1:k));
%                 plot(1:k, robust_FCD_Type_II_edge.Ge2{m}.targetattack(1:k));
%                 plot(1:k, robust_HET_edge.Ge2{m}.targetattack(1:k));
%                 plot(1:k, robust_PMG_edge.Ge2{m}.targetattack(1:k));
%                 plot(1:k, robust_FCD_Type_I_edge.Ge2{m}.targetattack(1:k));
                
                color_lines = get(gca, 'ColorOrder');
                
                Q = FDR(robust_FCD_Type_I_edge.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_I_edge.Ge12_diff_p{m}.targetattack<=0.05 & robust_FCD_Type_I_edge.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_FCD_Type_I_edge.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.18, 'V', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.18, '*', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_FCD_Type_I), robust_FCD_Type_I.Ge2{m}.targetattack(rc_FCD_Type_I)/robust_FCD_Type_I.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                
                Q = FDR(robust_FCD_Type_II_edge.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_II_edge.Ge12_diff_p{m}.targetattack<=0.05 & robust_FCD_Type_II_edge.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_FCD_Type_II_edge.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.03, 'V', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.03, '*', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                %             scatter(find(rc_FCD_Type_II), robust_FCD_Type_II.Ge2{m}.targetattack(rc_FCD_Type_II)/robust_FCD_Type_II.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_cont), robust_FCD_Type_I.Ge1{m}.targetattack(rc_cont)/robust_FCD_Type_I.Ge1{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                                
                Q = FDR(robust_HET_edge.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_HET_edge.Ge12_diff_p{m}.targetattack<=0.05 & robust_HET_edge.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_HET_edge.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.08, 'V', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.08, '*', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                %             scatter(find(rc_HET), robust_HET.Ge2{m}.targetattack(rc_HET)/robust_HET.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                
                Q = FDR(robust_PMG_edge.Ge12_diff_p{m}.targetattack, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_PMG_edge.Ge12_diff_p{m}.targetattack<=0.05 & robust_PMG_edge.Ge12_diff_p{m}.targetattack>Q);
                FDR_sig    = find(robust_PMG_edge.Ge12_diff_p{m}.targetattack<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.13, 'V', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.13, '*', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                %             scatter(find(rc_PMG), robust_PMG.Ge2{m}.targetattack(rc_PMG)/robust_PMG.Ge2{m}.targetattack(1), 'o', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                                
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                exportfigbo(gcf,[outputdir '/edge_targetattack2.png'], 'png', 13); close(gcf);
                
                %% random error
                m = 1;
                k = min([ length(robust_FCD_Type_I_edge.Ge1{m}.randomerror) length(robust_FCD_Type_II_edge.Ge1{m}.randomerror) length(robust_HET_edge.Ge1{m}.randomerror) length(robust_PMG_edge.Ge1{m}.randomerror) ]);
                
                figure; hold on;
                plot(1:k, robust_FCD_Type_I_edge.Ge1{m}.randomerror(1:k)/robust_FCD_Type_I.Ge1{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_I_edge.Ge2{m}.randomerror(1:k)/robust_FCD_Type_I.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_FCD_Type_II_edge.Ge2{m}.randomerror(1:k)/robust_FCD_Type_II.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_HET_edge.Ge2{m}.randomerror(1:k)/robust_HET.Ge2{m}.randomerror(1), 'LineWidth', linewidth);
                plot(1:k, robust_PMG_edge.Ge2{m}.randomerror(1:k)/robust_PMG.Ge2{m}.randomerror(1), 'LineWidth', linewidth);  
                temp = get(gca, 'Children');
                set(temp(4), 'Visible', 'off');
                
                color_lines = get(gca, 'ColorOrder');
                
                Q = FDR(robust_FCD_Type_I_edge.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_I_edge.Ge12_diff_p{m}.randomerror<=0.05 & robust_FCD_Type_I_edge.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_FCD_Type_I_edge.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.18, 'V', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.18, '*', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_FCD_Type_I), robust_FCD_Type_I.Ge2{m}.randomerror(rc_FCD_Type_I)/robust_FCD_Type_I.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                
                Q = FDR(robust_FCD_Type_II_edge.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_FCD_Type_II_edge.Ge12_diff_p{m}.randomerror<=0.05 & robust_FCD_Type_II_edge.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_FCD_Type_II_edge.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.03, 'V', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.03, '*', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                %             scatter(find(rc_FCD_Type_II), robust_FCD_Type_II.Ge2{m}.randomerror(rc_FCD_Type_II)/robust_FCD_Type_II.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                %             scatter(find(rc_cont), robust_FCD_Type_I.Ge1{m}.randomerror(rc_cont)/robust_FCD_Type_I.Ge1{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                
                
                Q = FDR(robust_HET_edge.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_HET_edge.Ge12_diff_p{m}.randomerror<=0.05 & robust_HET_edge.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_HET_edge.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.08, 'V', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.08, '*', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                %             scatter(find(rc_HET), robust_HET.Ge2{m}.randomerror(rc_HET)/robust_HET.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                
                Q = FDR(robust_PMG_edge.Ge12_diff_p{m}.randomerror, 0.05);
                if(isempty(Q))
                    Q = 0;
                end
                uncorr_sig = find(robust_PMG_edge.Ge12_diff_p{m}.randomerror<=0.05 & robust_PMG_edge.Ge12_diff_p{m}.randomerror>Q);
                FDR_sig    = find(robust_PMG_edge.Ge12_diff_p{m}.randomerror<=Q);
                scatter(uncorr_sig, ones(1, length(uncorr_sig))*1.13, 'V', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                scatter(FDR_sig, ones(1, length(FDR_sig))*1.13, '*', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                %             scatter(find(rc_PMG), robust_PMG.Ge2{m}.randomerror(rc_PMG)/robust_PMG.Ge2{m}.randomerror(1), 'o', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                

                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                exportfigbo(gcf,[outputdir '/edge_randomerror2.png'], 'png', 13); close(gcf);
                
            end            
            
        end
                
        if(save_file)
            if(lesion_exclusion)
                save([ MATDIR '04_rich_club_analysis_result_MCD_maskingout_lesion.mat' ], ...
                    '');
            else
                save([ MATDIR '04_rich_club_analysis_result_MCD.mat' ], ...
                    '');
            end
        end
        
    end  
    
end