clear all
close all

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath('/local_raid/seokjun/03_downloads/GRETNA/NetFunctions/');
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
        load( [ MATDIR '/05_hierarchical_analysis_result_MCD.mat' ] );
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
        
        randomization = 1000;
        brandnetwork = 1;
        interval_thres = 0.01;
        myspar = 0.11:interval_thres:0.4;
        option = [0.20, 0.25]; % FCD Type-I: [0.04, 0.04]
        
        [ beta_FCD_Type_I   beta_FCD_Type_I_rand  ] = hierdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_FCD_Type_I_res_val, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);
        [ beta_FCD_Type_II  beta_FCD_Type_II_rand ] = hierdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);
        [ beta_HET          beta_HET_rand         ] = hierdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);
        [ beta_PMG          beta_PMG_rand         ] = hierdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, randomization, brandnetwork, option);
        
        
        if(visualization == 1)
            
            OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/04_hierarchical_analysis/';
            MarkerSize = 20;
                        
            for beta = 1
                
                %% 1) presentation across groups
                figure; hold on; xlim([10 42]);
                plot(myspar*100, [ beta_FCD_Type_I.beta1; beta_FCD_Type_I.beta2; beta_FCD_Type_II.beta2; beta_HET.beta2; beta_PMG.beta2 ], 'LineWidth', 2);
                color_lines = get(gca, 'ColorOrder');
                scatter(myspar*100', [ beta_FCD_Type_I.beta1' ], MarkerSize, color_lines(1, :), 'filled', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                scatter(myspar*100', [ beta_FCD_Type_I.beta2' ], MarkerSize, color_lines(2, :), 'filled', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                scatter(myspar*100', [ beta_FCD_Type_II.beta2' ], MarkerSize, color_lines(3, :), 'filled', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                scatter(myspar*100', [ beta_HET.beta2' ], MarkerSize, color_lines(4, :), 'filled', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                scatter(myspar*100', [ beta_PMG.beta2' ], MarkerSize, color_lines(5, :), 'filled', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'hierarchical_global_across_groups.png'], 'png', 13); close(gcf);
                end
                
                %% 2) delta, FCD Type I
                colorline_idx = 2;
                sig_pos = -0.8;
                feature1 = beta_FCD_Type_I;
                feature2 = beta_FCD_Type_I_rand;
                FDRp = FDR(feature1.beta12_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1-feature2.beta2, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95Beta12diff; feature1.upper95Beta12diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1-feature1.beta2, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1-feature1.beta2, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta12_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta12_diff_p(i)>FDRp & feature1.beta12_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_FCD_TypeI.png'], 'png', 13); close(gcf);
                end
                
                %% 3) delta, FCD Type II
                colorline_idx = 3;
                sig_pos = -0.8;
                feature1 = beta_FCD_Type_II;
                feature2 = beta_FCD_Type_II_rand;
                FDRp = FDR(feature1.beta12_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1-feature2.beta2, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95Beta12diff; feature1.upper95Beta12diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1-feature1.beta2, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1-feature1.beta2, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta12_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta12_diff_p(i)>FDRp & feature1.beta12_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_FCD_TypeII.png'], 'png', 13); close(gcf);
                end
                
                %% 4) delta, HET
                colorline_idx = 4;
                sig_pos = -0.8;
                feature1 = beta_HET;
                feature2 = beta_HET_rand;
                FDRp = FDR(feature1.beta12_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1-feature2.beta2, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95Beta12diff; feature1.upper95Beta12diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1-feature1.beta2, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1-feature1.beta2, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta12_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta12_diff_p(i)>FDRp & feature1.beta12_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_HET.png'], 'png', 13); close(gcf);
                end
                
                %% 5) delta, PMG
                colorline_idx = 5;
                sig_pos = -0.8;
                feature1 = beta_PMG;
                feature2 = beta_PMG_rand;
                FDRp = FDR(feature1.beta12_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1-feature2.beta2, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95Beta12diff; feature1.upper95Beta12diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1-feature1.beta2, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1-feature1.beta2, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta12_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta12_diff_p(i)>FDRp & feature1.beta12_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_PMG.png'], 'png', 13); close(gcf);
                end
                
            end
            
            for betaz = 1
                
                %% 1) presentation across groups
                figure; hold on; xlim([10 42]);
                plot(myspar*100, [ beta_FCD_Type_I.beta1z; beta_FCD_Type_I.beta2z; beta_FCD_Type_II.beta2z; beta_HET.beta2z; beta_PMG.beta2z ], 'LineWidth', 2);
                color_lines = get(gca, 'ColorOrder');
                scatter(myspar*100', [ beta_FCD_Type_I.beta1z' ], MarkerSize, color_lines(1, :), 'filled', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                scatter(myspar*100', [ beta_FCD_Type_I.beta2z' ], MarkerSize, color_lines(2, :), 'filled', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                scatter(myspar*100', [ beta_FCD_Type_II.beta2z' ], MarkerSize, color_lines(3, :), 'filled', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                scatter(myspar*100', [ beta_HET.beta2z' ], MarkerSize, color_lines(4, :), 'filled', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                scatter(myspar*100', [ beta_PMG.beta2z' ], MarkerSize, color_lines(5, :), 'filled', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'hierarchical_global_across_groups.png'], 'png', 13); close(gcf);
                end
                
                %% 2) delta, FCD Type I
                colorline_idx = 2;
                sig_pos = -0.8;
                feature1 = beta_FCD_Type_I;
                feature2 = beta_FCD_Type_I_rand;
                FDRp = FDR(feature1.beta1z2_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1z-feature2.beta2z, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95beta1z2diff; feature1.upper95beta1z2diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1z-feature1.beta2z, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1z-feature1.beta2z, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta1z2_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta1z2_diff_p(i)>FDRp & feature1.beta1z2_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_FCD_TypeI.png'], 'png', 13); close(gcf);
                end
                
                %% 3) delta, FCD Type II
                colorline_idx = 3;
                sig_pos = -0.8;
                feature1 = beta_FCD_Type_II;
                feature2 = beta_FCD_Type_II_rand;
                FDRp = FDR(feature1.beta1z2_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1z-feature2.beta2z, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95beta1z2diff; feature1.upper95beta1z2diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1z-feature1.beta2z, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1z-feature1.beta2z, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta1z2_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta1z2_diff_p(i)>FDRp & feature1.beta1z2_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_FCD_TypeII.png'], 'png', 13); close(gcf);
                end
                
                %% 4) delta, HET
                colorline_idx = 4;
                sig_pos = -0.8;
                feature1 = beta_HET;
                feature2 = beta_HET_rand;
                FDRp = FDR(feature1.beta1z2_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1z-feature2.beta2z, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95beta1z2diff; feature1.upper95beta1z2diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1z-feature1.beta2z, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1z-feature1.beta2z, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta1z2_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta1z2_diff_p(i)>FDRp & feature1.beta1z2_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_HET.png'], 'png', 13); close(gcf);
                end
                
                %% 5) delta, PMG
                colorline_idx = 5;
                sig_pos = -0.8;
                feature1 = beta_PMG;
                feature2 = beta_PMG_rand;
                FDRp = FDR(feature1.beta1z2_diff_p, 0.05);
                
                if(isempty(FDRp))
                    FDRp = 0;
                end
                
                figure; hold on;
                plot(myspar*100, feature2.beta1z-feature2.beta2z, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                plot(myspar*100, [ feature1.lower95beta1z2diff; feature1.upper95beta1z2diff ], '--k', 'LineWidth', 1);
                plot(myspar*100, feature1.rerand_betaz12diff_mean, 'k', 'LineWidth', 2);
                scatter((myspar*100)', feature1.rerand_betaz12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                plot(myspar*100, feature1.beta1z-feature1.beta2z, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                scatter(myspar*100', feature1.beta1z-feature1.beta2z, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                for i = 1: length(myspar)
                    if(feature1.beta1z2_diff_p(i)<=FDRp)
                        text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    elseif(feature1.beta1z2_diff_p(i)>FDRp & feature1.beta1z2_diff_p(i)<=0.055)
                        text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                    end
                end
                ylim([-1 1]);
                
                if(save_file)
                    exportfigbo(gcf,[OUTPATH 'delta_hierarchical_global_diff_cont_vs_PMG.png'], 'png', 13); close(gcf);
                end
                
            end
            
        end
                
        if(save_file)
            if(lesion_exclusion)
                save([ MATDIR '05_hierarchical_analysis_result_MCD_maskingout_lesion.mat' ], ...
                    '');
            else
                save([ MATDIR '05_hierarchical_analysis_result_MCD.mat' ], ...
                    '');
            end
        end
        
    end  
    
end