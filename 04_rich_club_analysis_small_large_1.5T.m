clear all
close all

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/'));
P='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/';
RDIR='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_FCD_Type_I_HET/02_modularity/'; 
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
    
    AAL_surf_data_both = SurfStatReadData1('aal_both_rsl_final.txt');
    AAL_surf_data_left = AAL_surf_data_both(1, 1:40962);
    AAL_surf_data_right = AAL_surf_data_both(1, 40963:81924);
    AAL_surf_data_both = [ AAL_surf_data_left AAL_surf_data_left + 1];
    AAL_surf_data_both(find(AAL_surf_data_right==0)+40962) = 0;
    
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
    
    fid = fopen('DemographicNGrouping_for_PMG_Data_small_large.csv');
    C = textscan(fid, '%s%f%s%d%d%s%s%s%d','Delimiter',',','CollectOutput', 1);
    
    Inclusion                       = C{4}(:, 1);
    Codes_PMG                       = C{1}(Inclusion>0);
    Age_PMG                         = C{2}(Inclusion>0);     % mean/sd: 28.8/11.2
    Sex_PMG                         = C{3}(Inclusion>0);     % male/female: 12/9      
    Left_PMG                        = C{4}(Inclusion>0, 2) == 1; 
    Right_PMG                       = C{4}(Inclusion>0, 2) == 2;
    Bilateral_PMG                   = C{4}(Inclusion>0, 2) == 3;  % Left/Right/Bilateral: 6/5/10  
    Path_PMG                        = C{5}(Inclusion>0, 1);
    Size_PMG                        = C{6}(Inclusion>0, 1);
    
    %% to subgroup into small and large lesion cases
    %% condition1   : If the PMG is bilateral and forms a relatively large lesion across whole brain such as a fronto-parietal region: large
    %% condition2-1 : If the PMG is bilateral but the lesion covers only limited extend around the peri-sylvian area and the other areas are relatively intact: small
    %% condition2-2 : If the PMG is bilateral and although the lesion covers only limited extend around the peri-sylvian area but lesional size is big: large
    %% condition3   : If the PMG is unilateral and cover wide cortical areas (i.e., >3 multiple lobs): large
    %% condition4   : If the PMG is unilateral and limited to only to a single or two lobes or mildly encroached around peri-sylvian area : small
    
    Codes_PMG_s                     = Codes_PMG(Size_PMG == 1);
    Codes_PMG_l                     = Codes_PMG(Size_PMG == 2);
    Age_PMG_s                       = Age_PMG(Size_PMG == 1);
    Age_PMG_l                       = Age_PMG(Size_PMG == 2);
    Sex_PMG_s                       = Sex_PMG(Size_PMG == 1);
    Sex_PMG_l                       = Sex_PMG(Size_PMG == 2);
    Left_PMG_s                      = Left_PMG(Size_PMG == 1);
    Left_PMG_l                      = Left_PMG(Size_PMG == 2);
    Right_PMG_s                     = Right_PMG(Size_PMG == 1);
    Right_PMG_l                     = Right_PMG(Size_PMG == 2);
    Bilateral_PMG_s                 = Bilateral_PMG(Size_PMG == 1);
    Bilateral_PMG_l                 = Bilateral_PMG(Size_PMG == 2);
    Path_PMG_s                      = Path_PMG(Size_PMG == 1);
    Path_PMG_l                      = Path_PMG(Size_PMG == 2);
    
    fclose(fid);
    
end

%% load matfiles
for load_files = 1
    
    if(lesion_exclusion)
    else
        load( [ MATDIR '/01_generate_association_matrix_PMG_small_vs_large.mat' ] );
        load( [ MATDIR '/04_rich_club_analysis_result_PMG_small_vs_large.mat' ] );
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
            
            parpool(24);
            
            randomization = 1000;
            brandnetwork = 1;
            interval_thres = 0.01;
            myspar = 0.11:interval_thres:0.4;
            option = [0.04]; % FCD Type-I: [0.04, 0.04]
            
            [ rc_PMG         rand_rc_PMG         ] = rcdiff_sparsity_parfor(asso_mat_PMG_res_val(Size_PMG==1, :), asso_mat_PMG_res_val(Size_PMG==2, :), myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);          figure;
            
            if(visualization == 1)
                
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/03_rich_club_analysis/';
                
                %% global rich-club coefficient and normalized values
                m = 3;
                k = min([ length(rc_PMG{m}.R.R1) length(rc_PMG{m}.R.R2) ]);
                
                figure; hold on;                
                plot(1:k,rc_PMG{m}.R.R1(1:k), '--', 'LineWidth', 2);
                plot(1:k,rc_PMG{m}.R.R2(1:k), '--', 'LineWidth', 2);               
                legend('PMG small', 'PMG large', 'Location', 'best');
                
                color_lines = get(gca, 'ColorOrder');                
                plot(1:k,rc_PMG{m}.R.normR1(1:k), 'Color', color_lines(1, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(1, :), 'LineWidth', 2);
                plot(1:k,rc_PMG{m}.R.normR2(1:k), 'Color', color_lines(2, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(2, :), 'LineWidth', 2);
                plot(0:k+1, ones(1, k+2), ':', 'Color', [0.6 0.6 0.6]);
                exportfigbo(gcf,[OUTPATH 'rich_club_coefficient1.png'], 'png', 13); close(gcf);
                
                figure; hold on;
                plot(1:k, max([rc_FCD_Type_I{m}.R.R1(1:k); rc_FCD_Type_II{m}.R.R1(1:k); rc_HET{m}.R.R1(1:k); rc_PMG{m}.R.R1(1:k)]), '--', 'LineWidth', 2, 'Color', color_lines(1, :));
                plot(1:k,rc_FCD_Type_II{m}.R.R2(1:k), '--', 'LineWidth', 2, 'Color', color_lines(3, :));
                plot(1:k,rc_HET{m}.R.R2(1:k), '--', 'LineWidth', 2, 'Color', color_lines(4, :));
                plot(1:k,rc_PMG{m}.R.R2(1:k), '--', 'LineWidth', 2, 'Color', color_lines(5, :));
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'Location', 'best');
                
                plot(1:k, max([rc_FCD_Type_I{m}.R.normR1(1:k); rc_FCD_Type_II{m}.R.normR1(1:k); rc_HET{m}.R.normR1(1:k); rc_PMG{m}.R.normR1(1:k)]), 'Color', color_lines(1, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(1, :), 'LineWidth', 2);
                plot(1:k,rc_FCD_Type_II{m}.R.normR2(1:k), 'Color', color_lines(3, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(3, :), 'LineWidth', 2);
                plot(1:k,rc_HET{m}.R.normR2(1:k), 'Color', color_lines(4, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(4, :), 'LineWidth', 2);
                plot(1:k,rc_PMG{m}.R.normR2(1:k), 'Color', color_lines(5, :), 'Marker', 'o', 'MarkerFaceColor', color_lines(5, :), 'LineWidth', 2);
                plot(0:k+1, ones(1, k+2), ':', 'Color', [0.6 0.6 0.6]);
                exportfigbo(gcf,[OUTPATH 'rich_club_coefficient2.png'], 'png', 13); close(gcf);
                
            end
            
            if(save_file)
                if(lesion_exclusion)
                    save([ MATDIR '04_rich_club_analysis_result_MCD_maskingout_lesion.mat' ], ...
                        '');
                else
                    save([ MATDIR '04_rich_club_analysis_result_PMG_small_large.mat' ], ...
                        'rc_PMG', 'rand_rc_PMG');
                end
            end
            
        end        
        
    end  
    
end

%% assess the betweenness centrality with respect to rich club
for betweeness_centrality = 1
    
    interval_thres = 0.01;
    myspar = 0.11:interval_thres:0.4;
    AAL_index_iter = length(Anatomical_Label_structs_idx);
    
    %% Node and edge centrality
    betweenness_centrality_cont         = zeros(length(myspar), AAL_index_iter);
    betweenness_centrality_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
    betweenness_centrality_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
    betweenness_centrality_HET          = zeros(length(myspar), AAL_index_iter);
    betweenness_centrality_PMG          = zeros(length(myspar), AAL_index_iter);
    
    betweenness_centrality_edge_cont         = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
    betweenness_centrality_edge_FCD_Type_II  = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
    betweenness_centrality_edge_FCD_Type_II  = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
    betweenness_centrality_edge_HET          = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
    betweenness_centrality_edge_PMG          = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
    
    for i = 1 : length(myspar)
        
        %% controls
        wmatrix = threshold_proportional(adjacent_matrix_control, myspar(i));
        [ betweenness_centrality_edge_cont(:, :, i) betweenness_centrality_cont(i, :) ] = edge_betweenness_bin(double(wmatrix>0));
%         bc_edge  = betweenness_centrality_edge_cont(:, :, i);                
%         max_edge = max(bc_edge((wmatrix>0)));
%         min_edge = min(bc_edge((wmatrix>0)));
%         
%         bc_node  = betweenness_centrality_cont(i, :);
%         max_node = max(bc_node(bc_node>0));
%         min_node = min(bc_node(bc_node>0));
%         betweenness_centrality_edge_cont(:, :, i) = (betweenness_centrality_edge_cont(:, :, i) - min_edge)/(max_edge - min_edge);
%         betweenness_centrality_cont(i, :)         = (betweenness_centrality_cont(i, :) - min_node)/(max_node - min_node);
        betweenness_centrality_edge_cont(:, :, i) = betweenness_centrality_edge_cont(:, :, i) / mean(mean(betweenness_centrality_edge_cont(:, :, i)));
        betweenness_centrality_cont(i, :)         = betweenness_centrality_cont(i, :) / mean(betweenness_centrality_cont(i, :));
        
        %% FCD Type II
        wmatrix = threshold_proportional(adjacent_matrix_FCD_Type_II, myspar(i));
        [ betweenness_centrality_edge_FCD_Type_II(:, :, i) betweenness_centrality_FCD_Type_II(i, :) ] = edge_betweenness_bin(double(wmatrix>0));
%         bc_edge  = betweenness_centrality_edge_FCD_Type_II(:, :, i);                
%         max_edge = max(bc_edge((wmatrix>0)));
%         min_edge = min(bc_edge((wmatrix>0)));
%         
%         bc_node  = betweenness_centrality_FCD_Type_II(i, :);
%         max_node = max(bc_node(bc_node>0));
%         min_node = min(bc_node(bc_node>0));
%         betweenness_centrality_edge_FCD_Type_II(:, :, i) = (betweenness_centrality_edge_FCD_Type_II(:, :, i) - min_edge)/(max_edge - min_edge);
%         betweenness_centrality_FCD_Type_II(i, :)         = (betweenness_centrality_FCD_Type_II(i, :) - min_node)/(max_node - min_node);
        betweenness_centrality_edge_FCD_Type_II(:, :, i) = betweenness_centrality_edge_FCD_Type_II(:, :, i) / mean(mean(betweenness_centrality_edge_FCD_Type_II(:, :, i)));
        betweenness_centrality_FCD_Type_II(i, :)         = betweenness_centrality_FCD_Type_II(i, :) / mean(betweenness_centrality_FCD_Type_II(i, :));
        
        %% FCD Type I
        wmatrix = threshold_proportional(adjacent_matrix_FCD_Type_I, myspar(i));
        [ betweenness_centrality_edge_FCD_Type_I(:, :, i) betweenness_centrality_FCD_Type_I(i, :) ] = edge_betweenness_bin(double(wmatrix>0));
%         bc_edge  = betweenness_centrality_edge_FCD_Type_I(:, :, i);                
%         max_edge = max(bc_edge((wmatrix>0)));
%         min_edge = min(bc_edge((wmatrix>0)));
%         
%         bc_node  = betweenness_centrality_FCD_Type_I(i, :);
%         max_node = max(bc_node(bc_node>0));
%         min_node = min(bc_node(bc_node>0));
%         betweenness_centrality_edge_FCD_Type_I(:, :, i) = (betweenness_centrality_edge_FCD_Type_I(:, :, i) - min_edge)/(max_edge - min_edge);
%         betweenness_centrality_FCD_Type_I(i, :)         = (betweenness_centrality_FCD_Type_I(i, :) - min_node)/(max_node - min_node);
        betweenness_centrality_edge_FCD_Type_I(:, :, i) = betweenness_centrality_edge_FCD_Type_I(:, :, i) / mean(mean(betweenness_centrality_edge_FCD_Type_I(:, :, i)));
        betweenness_centrality_FCD_Type_I(i, :)         = betweenness_centrality_FCD_Type_I(i, :) / mean(betweenness_centrality_FCD_Type_I(i, :));
        
        %% HET
        wmatrix = threshold_proportional(adjacent_matrix_HET, myspar(i));
        [ betweenness_centrality_edge_HET(:, :, i) betweenness_centrality_HET(i, :) ] = edge_betweenness_bin(double(wmatrix>0));
%         bc_edge  = betweenness_centrality_edge_HET(:, :, i);                
%         max_edge = max(bc_edge((wmatrix>0)));
%         min_edge = min(bc_edge((wmatrix>0)));
%         
%         bc_node  = betweenness_centrality_HET(i, :);
%         max_node = max(bc_node(bc_node>0));
%         min_node = min(bc_node(bc_node>0));
%         betweenness_centrality_edge_HET(:, :, i) = (betweenness_centrality_edge_HET(:, :, i) - min_edge)/(max_edge - min_edge);
%         betweenness_centrality_HET(i, :)         = (betweenness_centrality_HET(i, :) - min_node)/(max_node - min_node);
        betweenness_centrality_edge_HET(:, :, i) = betweenness_centrality_edge_HET(:, :, i) / mean(mean(betweenness_centrality_edge_HET(:, :, i)));
        betweenness_centrality_HET(i, :)         = betweenness_centrality_HET(i, :) / mean(betweenness_centrality_HET(i, :));
        
        %% PMG
        wmatrix = threshold_proportional(adjacent_matrix_PMG, myspar(i));
        [ betweenness_centrality_edge_PMG(:, :, i) betweenness_centrality_PMG(i, :) ] = edge_betweenness_bin(double(wmatrix>0));
%         bc_edge  = betweenness_centrality_edge_PMG(:, :, i);                
%         max_edge = max(bc_edge((wmatrix>0)));
%         min_edge = min(bc_edge((wmatrix>0)));
%         
%         bc_node  = betweenness_centrality_PMG(i, :);
%         max_node = max(bc_node(bc_node>0));
%         min_node = min(bc_node(bc_node>0));
%         betweenness_centrality_edge_PMG(:, :, i) = (betweenness_centrality_edge_PMG(:, :, i) - min_edge)/(max_edge - min_edge);
%         betweenness_centrality_PMG(i, :)         = (betweenness_centrality_PMG(i, :) - min_node)/(max_node - min_node);
        betweenness_centrality_edge_PMG(:, :, i) = betweenness_centrality_edge_PMG(:, :, i) / mean(mean(betweenness_centrality_edge_PMG(:, :, i)));
        betweenness_centrality_PMG(i, :)         = betweenness_centrality_PMG(i, :) / mean(betweenness_centrality_PMG(i, :));
        
    end
    
    %% Combined analysis of centrality and rich-club
    for rich_club_based_centrality_analysis = 1
        
        m = 1;
        k_sig = 11;
        
        % - control
        for control = 1
            
            R_matrix   = corrcoef(asso_mat_cont_res_val);
            wmatrix = threshold_proportional(R_matrix, myspar(m));
            corrmatrix = double(wmatrix>0);
            [deg] = degrees_und(corrmatrix);
            rc_nodes = deg>=k_sig;
            num_rc_node_cont     = sum(rc_nodes);
            num_nonrc_node_cont  = sum(~rc_nodes);
            num_rc_link_cont     = sum(sum(corrmatrix(rc_nodes, rc_nodes)))/2;
            num_feeder_node_cont = sum(sum(corrmatrix(rc_nodes, ~rc_nodes)));
            num_local_node_cont  = sum(sum(corrmatrix(~rc_nodes, ~rc_nodes)))/2;
            
            % Rich club nodes vs. non rich club nodes 
            rc_bc_cont      = [ mean(betweenness_centrality_cont(m, rc_nodes))  std(betweenness_centrality_cont(m, rc_nodes))  ];
            nonrc_bc_cont   = [ mean(betweenness_centrality_cont(m, ~rc_nodes)) std(betweenness_centrality_cont(m, ~rc_nodes)) ];

            % Rich club links vs. feeder links vs. peripheral links            
            bc_curr = betweenness_centrality_edge_cont(rc_nodes, rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, rc_nodes))); 
            rc_bc_link_cont      = [ mean(bc_curr)  std(bc_curr)  ];
            
            bc_curr = betweenness_centrality_edge_cont(rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, ~rc_nodes)));             
            feeder_bc_link_cont   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_cont(~rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, ~rc_nodes)));
            peripheral_bc_link_cont   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_cont(rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, :)));
            rc_feeder_bc_link_cont   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_cont(~rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, :)));
            feeder_peripheral_bc_link_cont   = [ mean(bc_curr) std(bc_curr) ];
            
        end
        
        % - FCD Type-II
        for FCD_Type_II = 1
            
            R_matrix   = corrcoef(asso_mat_FCD_Type_II_res_val);
            wmatrix = threshold_proportional(R_matrix, myspar(m));
            corrmatrix = double(wmatrix>0);
            [deg] = degrees_und(corrmatrix);
            rc_nodes = deg>=k_sig;
            num_rc_node_FCD_Type_II     = sum(rc_nodes);
            num_nonrc_node_FCD_Type_II  = sum(~rc_nodes);
            num_rc_link_FCD_Type_II     = sum(sum(corrmatrix(rc_nodes, rc_nodes)))/2;
            num_feeder_node_FCD_Type_II = sum(sum(corrmatrix(rc_nodes, ~rc_nodes)));
            num_local_node_FCD_Type_II  = sum(sum(corrmatrix(~rc_nodes, ~rc_nodes)))/2;
            
            % Rich club nodes vs. non rich club nodes
            rc_bc_FCD_Type_II      = [ mean(betweenness_centrality_FCD_Type_II(m, rc_nodes))  std(betweenness_centrality_FCD_Type_II(m, rc_nodes))  ];
            nonrc_bc_FCD_Type_II   = [ mean(betweenness_centrality_FCD_Type_II(m, ~rc_nodes)) std(betweenness_centrality_FCD_Type_II(m, ~rc_nodes)) ];

            % Rich club links vs. feeder links vs. peripheral links            
            bc_curr = betweenness_centrality_edge_FCD_Type_II(rc_nodes, rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, rc_nodes))); 
            rc_bc_link_FCD_Type_II      = [ mean(bc_curr)  std(bc_curr)  ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_II(rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, ~rc_nodes)));             
            feeder_bc_link_FCD_Type_II   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_II(~rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, ~rc_nodes)));             
            peripheral_bc_link_FCD_Type_II   = [ mean(bc_curr) std(bc_curr) ];
                        
            bc_curr = betweenness_centrality_edge_FCD_Type_II(rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, :)));
            rc_feeder_bc_link_FCD_Type_II   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_II(~rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, :)));
            feeder_peripheral_bc_link_FCD_Type_II   = [ mean(bc_curr) std(bc_curr) ];
            
        end
        
        % - HET
        for HET = 1
            
            R_matrix   = corrcoef(asso_mat_HET_res_val);
            wmatrix = threshold_proportional(R_matrix, myspar(m));
            corrmatrix = double(wmatrix>0);
            [deg] = degrees_und(corrmatrix);
            rc_nodes = deg>=k_sig;
            num_rc_node_HET     = sum(rc_nodes);
            num_nonrc_node_HET  = sum(~rc_nodes);
            num_rc_link_HET     = sum(sum(corrmatrix(rc_nodes, rc_nodes)))/2;
            num_feeder_node_HET = sum(sum(corrmatrix(rc_nodes, ~rc_nodes)));
            num_local_node_HET  = sum(sum(corrmatrix(~rc_nodes, ~rc_nodes)))/2;
            
            % Rich club nodes vs. non rich club nodes
            rc_bc_HET      = [ mean(betweenness_centrality_HET(m, rc_nodes))  std(betweenness_centrality_HET(m, rc_nodes))  ];
            nonrc_bc_HET   = [ mean(betweenness_centrality_HET(m, ~rc_nodes)) std(betweenness_centrality_HET(m, ~rc_nodes)) ];

            % Rich club links vs. feeder links vs. peripheral links            
            bc_curr = betweenness_centrality_edge_HET(rc_nodes, rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, rc_nodes))); 
            rc_bc_link_HET      = [ mean(bc_curr)  std(bc_curr)  ];
            
            bc_curr = betweenness_centrality_edge_HET(rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, ~rc_nodes)));             
            feeder_bc_link_HET   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_HET(~rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, ~rc_nodes)));             
            peripheral_bc_link_HET   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_HET(rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, :)));
            rc_feeder_bc_link_HET   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_HET(~rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, :)));
            feeder_peripheral_bc_link_HET   = [ mean(bc_curr) std(bc_curr) ];            
            
        end
        
        % - PMG
        for PMG = 1
            
            R_matrix   = corrcoef(asso_mat_PMG_res_val);
            wmatrix = threshold_proportional(R_matrix, myspar(m));
            corrmatrix = double(wmatrix>0);
            [deg] = degrees_und(corrmatrix);
            rc_nodes = deg>=k_sig;
            num_rc_node_PMG     = sum(rc_nodes);
            num_nonrc_node_PMG  = sum(~rc_nodes);
            num_rc_link_PMG     = sum(sum(corrmatrix(rc_nodes, rc_nodes)))/2;
            num_feeder_node_PMG = sum(sum(corrmatrix(rc_nodes, ~rc_nodes)));
            num_local_node_PMG  = sum(sum(corrmatrix(~rc_nodes, ~rc_nodes)))/2;            
            
            % Rich club nodes vs. non rich club nodes
            rc_bc_PMG      = [ mean(betweenness_centrality_PMG(m, rc_nodes))  std(betweenness_centrality_PMG(m, rc_nodes))  ];
            nonrc_bc_PMG   = [ mean(betweenness_centrality_PMG(m, ~rc_nodes)) std(betweenness_centrality_PMG(m, ~rc_nodes)) ];

            % Rich club links vs. feeder links vs. peripheral links            
            bc_curr = betweenness_centrality_edge_PMG(rc_nodes, rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, rc_nodes))); 
            rc_bc_link_PMG      = [ mean(bc_curr)  std(bc_curr)  ];
            
            bc_curr = betweenness_centrality_edge_PMG(rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, ~rc_nodes)));             
            feeder_bc_link_PMG   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_PMG(~rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, ~rc_nodes)));             
            peripheral_bc_link_PMG   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_PMG(rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, :)));
            rc_feeder_bc_link_PMG   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_PMG(~rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, :)));
            feeder_peripheral_bc_link_PMG   = [ mean(bc_curr) std(bc_curr) ];            
            
        end
        
        % - FCD Type-I
        for FCD_Type_I = 1
            
            R_matrix   = corrcoef(asso_mat_FCD_Type_I_res_val);
            wmatrix = threshold_proportional(R_matrix, myspar(m));
            corrmatrix = double(wmatrix>0);
            [deg] = degrees_und(corrmatrix);
            rc_nodes = deg>=k_sig;
            num_rc_node_FCD_Type_I     = sum(rc_nodes);
            num_nonrc_node_FCD_Type_I  = sum(~rc_nodes);
            num_rc_link_FCD_Type_I     = sum(sum(corrmatrix(rc_nodes, rc_nodes)))/2;
            num_feeder_node_FCD_Type_I = sum(sum(corrmatrix(rc_nodes, ~rc_nodes)));
            num_local_node_FCD_Type_I  = sum(sum(corrmatrix(~rc_nodes, ~rc_nodes)))/2;            
            
            % Rich club nodes vs. non rich club nodes
            rc_bc_FCD_Type_I      = [ mean(betweenness_centrality_FCD_Type_I(m, rc_nodes))  std(betweenness_centrality_FCD_Type_I(m, rc_nodes))  ];
            nonrc_bc_FCD_Type_I   = [ mean(betweenness_centrality_FCD_Type_I(m, ~rc_nodes)) std(betweenness_centrality_FCD_Type_I(m, ~rc_nodes)) ];

            % Rich club links vs. feeder links vs. peripheral links            
            bc_curr = betweenness_centrality_edge_FCD_Type_I(rc_nodes, rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, rc_nodes))); 
            rc_bc_link_FCD_Type_I      = [ mean(bc_curr)  std(bc_curr)  ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_I(rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, ~rc_nodes)));             
            feeder_bc_link_FCD_Type_I   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_I(~rc_nodes, ~rc_nodes, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, ~rc_nodes)));             
            peripheral_bc_link_FCD_Type_I   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_I(rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(rc_nodes, :)));
            rc_feeder_bc_link_FCD_Type_I   = [ mean(bc_curr) std(bc_curr) ];
            
            bc_curr = betweenness_centrality_edge_FCD_Type_I(~rc_nodes, :, m);
            bc_curr = bc_curr(logical(corrmatrix(~rc_nodes, :)));
            feeder_peripheral_bc_link_FCD_Type_I   = [ mean(bc_curr) std(bc_curr) ];
            
        end
        
        % final result collections
        for final_result = 1
            
            a = [ rc_bc_cont nonrc_bc_cont; ...
                  rc_bc_FCD_Type_II nonrc_bc_FCD_Type_II; ...
                  rc_bc_HET nonrc_bc_HET; ...
                  rc_bc_PMG nonrc_bc_PMG; ...
                  rc_bc_FCD_Type_I nonrc_bc_FCD_Type_I ]
            b = [ rc_bc_link_cont        feeder_bc_link_cont        peripheral_bc_link_cont        rc_feeder_bc_link_cont        feeder_peripheral_bc_link_cont; ...
                  rc_bc_link_FCD_Type_II feeder_bc_link_FCD_Type_II peripheral_bc_link_FCD_Type_II rc_feeder_bc_link_FCD_Type_II feeder_peripheral_bc_link_FCD_Type_II; ...
                  rc_bc_link_HET         feeder_bc_link_HET         peripheral_bc_link_HET         rc_feeder_bc_link_HET         feeder_peripheral_bc_link_HET; ...
                  rc_bc_link_PMG         feeder_bc_link_PMG         peripheral_bc_link_PMG         rc_feeder_bc_link_PMG         feeder_peripheral_bc_link_PMG; ...
                  rc_bc_link_FCD_Type_I  feeder_bc_link_FCD_Type_I  peripheral_bc_link_FCD_Type_I  rc_feeder_bc_link_FCD_Type_I  feeder_peripheral_bc_link_FCD_Type_I ]
            
            org_diff_final = [ a b ];
            
        end
        
        % significance
        for significance_comp = 1
            
            parpool(24);
            
            m = 1;
            k_sig=11; % signficant differences between patients and controls
            groups   = { 'cont', 'FCD_Type_II', 'HET', 'PMG', 'FCD_Type_I' };
            groupstr = { '', 'Control', 'FCD-II', 'HET', 'PMG', 'FCD-I', ''};            
                                    
            p_set = zeros(length(groups)-1, 5);
            for i = 2 : length(groups)
                
                eval(['N_sub = [ size(asso_mat_cont_res_val, 1) size(asso_mat_' groups{i} '_res_val, 1) ];']);
                N1 = N_sub(1); N2 = N_sub(2);
                eval(['rand_rc = rand_rc_' groups{i} ';']);
                eval(['tgroup_Reg = [asso_mat_cont_res_val;asso_mat_' groups{i} '_res_val];']);
                
                rc_bc_set                           = cell(1, randomization);
                delta_rc_bc_rand                    = cell(1, randomization);
                delta_nonrc_bc_rand                 = cell(1, randomization);
                delta_rc_bc_link_rand               = cell(1, randomization);
                delta_feeder_bc_link_rand           = cell(1, randomization);
                delta_peripheral_bc_link_rand       = cell(1, randomization);
                delta_nonrc_bc_link_rand            = cell(1, randomization);
                
                parfor r = 1 : randomization
                    
                    r
                    
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
                
                    % - betweenness centrality
                    [ edge_bc node_bc ] = edge_betweenness_bin(corrmatrix);
                    edge_bc = edge_bc / mean(mean(edge_bc));
                    node_bc = node_bc / mean(node_bc);
                    
                    rc_bc_group1_rand               =  [ mean(node_bc(rich_club_nodes_rand))      std(node_bc(rich_club_nodes_rand))      ];
                    nonrc_bc_group1_rand            =  [ mean(node_bc(non_rich_club_nodes_rand))  std(node_bc(non_rich_club_nodes_rand))  ];
                    
                    bc_curr = edge_bc(rich_club_nodes_rand, rich_club_nodes_rand);
                    bc_curr = bc_curr(logical(corrmatrix(rich_club_nodes_rand, rich_club_nodes_rand)));
                    rc_bc_link_group1_rand          = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    bc_curr = edge_bc(rich_club_nodes_rand, non_rich_club_nodes_rand);
                    bc_curr = bc_curr(logical(corrmatrix(rich_club_nodes_rand, non_rich_club_nodes_rand)));
                    feeder_bc_link_group1_rand      = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    bc_curr = edge_bc(non_rich_club_nodes_rand, non_rich_club_nodes_rand);
                    bc_curr = bc_curr(logical(corrmatrix(non_rich_club_nodes_rand, non_rich_club_nodes_rand)));
                    peripheral_bc_link_group1_rand  = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    bc_curr = edge_bc(non_rich_club_nodes_rand, :);
                    bc_curr = bc_curr(logical(corrmatrix(non_rich_club_nodes_rand, :)));
                    nonrc_bc_link_group1_rand  = [ mean(bc_curr)  std(bc_curr)  ];
                                                            
                    % group 2
                    wmatrix                     = threshold_proportional(rerand_corr_2, myspar(m));
                    corrmatrix                  = double(wmatrix>0);
                    [deg]                       = degrees_und(corrmatrix);
                    rich_club_nodes_rand        = find(deg>=k_sig);
                    non_rich_club_nodes_rand    = find(deg<k_sig);
                
                    % - betweenness centrality
                    [ edge_bc node_bc ] = edge_betweenness_bin(corrmatrix);
                    edge_bc = edge_bc / mean(mean(edge_bc));
                    node_bc = node_bc / mean(node_bc);
                    
                    rc_bc_group2_rand               =  [ mean(node_bc(rich_club_nodes_rand))      std(node_bc(rich_club_nodes_rand))      ];
                    nonrc_bc_group2_rand            =  [ mean(node_bc(non_rich_club_nodes_rand))  std(node_bc(non_rich_club_nodes_rand))  ];
                    
                    bc_curr = edge_bc(rich_club_nodes_rand, rich_club_nodes_rand);
                    bc_curr = bc_curr(logical(corrmatrix(rich_club_nodes_rand, rich_club_nodes_rand)));
                    rc_bc_link_group2_rand          = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    bc_curr = edge_bc(rich_club_nodes_rand, non_rich_club_nodes_rand);
                    bc_curr = bc_curr(logical(corrmatrix(rich_club_nodes_rand, non_rich_club_nodes_rand)));
                    feeder_bc_link_group2_rand      = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    bc_curr = edge_bc(non_rich_club_nodes_rand, non_rich_club_nodes_rand);
                    bc_curr = bc_curr(logical(corrmatrix(non_rich_club_nodes_rand, non_rich_club_nodes_rand)));
                    peripheral_bc_link_group2_rand  = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    bc_curr = edge_bc(non_rich_club_nodes_rand, :);
                    bc_curr = bc_curr(logical(corrmatrix(non_rich_club_nodes_rand, :)));
                    nonrc_bc_link_group2_rand  = [ mean(bc_curr)  std(bc_curr)  ];
                    
                    rc_bc_set{r}.rc_bc_group1_rand              = rc_bc_group1_rand;
                    rc_bc_set{r}.nonrc_bc_group1_rand           = nonrc_bc_group1_rand;
                    rc_bc_set{r}.rc_bc_link_group1_rand         = rc_bc_link_group1_rand;
                    rc_bc_set{r}.feeder_bc_link_group1_rand     = feeder_bc_link_group1_rand;
                    rc_bc_set{r}.peripheral_bc_link_group1_rand = peripheral_bc_link_group1_rand;
                    rc_bc_set{r}.nonrc_bc_link_group1_rand      = nonrc_bc_link_group1_rand;
                    
                    rc_bc_set{r}.rc_bc_group2_rand              = rc_bc_group2_rand;
                    rc_bc_set{r}.nonrc_bc_group2_rand           = nonrc_bc_group2_rand;
                    rc_bc_set{r}.rc_bc_link_group2_rand         = rc_bc_link_group2_rand;
                    rc_bc_set{r}.feeder_bc_link_group2_rand     = feeder_bc_link_group2_rand;
                    rc_bc_set{r}.peripheral_bc_link_group2_rand = peripheral_bc_link_group2_rand;
                    rc_bc_set{r}.nonrc_bc_link_group2_rand      = nonrc_bc_link_group2_rand;
                    
                    delta_rc_bc_rand{r}                    = rc_bc_group1_rand(1)              - rc_bc_group2_rand(1);
                    delta_nonrc_bc_rand{r}                 = nonrc_bc_group1_rand(1)           - nonrc_bc_group2_rand(1);
                    delta_rc_bc_link_rand{r}               = rc_bc_link_group1_rand(1)         - rc_bc_link_group2_rand(1);
                    delta_feeder_bc_link_rand{r}           = feeder_bc_link_group1_rand(1)     - feeder_bc_link_group2_rand(1);
                    delta_peripheral_bc_link_rand{r}       = peripheral_bc_link_group1_rand(1) - peripheral_bc_link_group2_rand(1);
                    delta_nonrc_bc_link_rand{r}            = nonrc_bc_link_group1_rand(1)      - nonrc_bc_link_group2_rand(1);
                    
                end
                
                delta_rc_bc_rand                    = cell2mat(delta_rc_bc_rand);
                delta_nonrc_bc_rand                 = cell2mat(delta_nonrc_bc_rand);
                delta_rc_bc_link_rand               = cell2mat(delta_rc_bc_link_rand);
                delta_feeder_bc_link_rand           = cell2mat(delta_feeder_bc_link_rand);
                delta_peripheral_bc_link_rand       = cell2mat(delta_peripheral_bc_link_rand);
                delta_nonrc_bc_link_rand            = cell2mat(delta_nonrc_bc_link_rand);
                
                delta_rand = [ delta_rc_bc_rand' delta_nonrc_bc_rand' delta_rc_bc_link_rand' delta_feeder_bc_link_rand' delta_peripheral_bc_link_rand' ];
                
                for c = 1 : 5
                    
                    org_diff = org_diff_final(1, (c-1)*2+1) - org_diff_final(i, (c-1)*2+1);
                    if(org_diff>=0)
                        p_set(i-1, c) = sum(org_diff<delta_rand(:, c))/randomization;
                    elseif(org_diff<0)
                        p_set(i-1, c) = sum(org_diff>delta_rand(:, c))/randomization;
                    end
                    
                end
                    
            end
        end
          
    end
        
end