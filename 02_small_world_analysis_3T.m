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
CIVET_ver2 = 1; % 1 -> HET data processed after the removal of nodules | 2 -> HET data processed before the removal of nodules

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
    
    fid = fopen('DemographicNGrouping_for_controls_Data_3T.csv');
    C = textscan(fid, '%s%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    Codes_cont                      = C{1}(:, 1);    
    Path_cont                       = strcat(C{1}(:, 2), repmat(['ver' num2str(CIVET_ver1) '/' ], length(Codes_cont), 1));
    Age_cont                        = C{2};                 % mean/sd: 30.9/11.0
    Sex_cont                        = C{3}(:, 1);           % male/female: 16/25    
    fclose(fid);
    
end

%% read FCD Type-II demographics
for read_FCD_Type_II = 1
    
    fid = fopen('DemographicNGrouping_for_FCD_Type-II_Data_3T.csv');
    C = textscan(fid, '%s%f%s%d%d%d%d%s%d%d%d%d','Delimiter',',','CollectOutput', 1);
    Codes_FCD_Type_II               = C{1};
    if(CIVET_ver1 == 2)
        for i = 1 : length(Codes_FCD_Type_II)
            Codes_FCD_Type_II{i} = Codes_FCD_Type_II{i}(1:3);
        end
    end
    Age_FCD_Type_II                 = C{2};                 % mean/sd: 27.8/10.1
    Sex_FCD_Type_II                 = C{3};                 % male/female: 16/14        
    Left_FCD_Type_II                = C{4}(:, 1);           % 15
    Right_FCD_Type_II               = C{4}(:, 2);           % 15
    Pre_Central_FCD_Type_II         = C{4}(:, 3);           % 16
    Post_Central_FCD_Type_II        = C{4}(:, 4);           % 14    
    Path_FCD_Type_II                = strcat(C{5}, repmat(['ver' num2str(CIVET_ver1) '/' ], length(Codes_FCD_Type_II), 1));    
    Duration_FCD_Type_II            = double(C{6}(:, 1));   % mean/sd: 20.9/11.9
    Operation_FCD_Type_II           = double(C{6}(:, 2));   % 30
    Scanner_FCD_Type_II             = double(C{6}(:, 3));
    Lesion_Volume                   = double(C{6}(:, 4));    
    
    LabelPath_FCD = '/host/gypsy/local_raid/seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label_3T/';
    fclose(fid); 
    
end

%% read HET demographics
for read_HET = 1
    
    fid = fopen('DemographicNGrouping_for_HET_Data_3T.csv');
    C = textscan(fid, '%s%f%s%d%s%s%d','Delimiter',',','CollectOutput', 1);    
    Codes_HET                       = C{1};
    Age_HET                         = C{2};    % mean/sd: 28.1/11.5
    Sex_HET                         = C{3};    % male/female: 11/16
    Left_HET                        = C{4} == 1; 
    Right_HET                       = C{4} == 2;
    Bilateral_HET                   = C{4} == 3;  % Left/Right/Bilateral: 6/10/11  
    Path_HET                        = strcat(C{5}(:, 1), repmat(['ver' num2str(CIVET_ver2) '/' ], length(Codes_HET), 1));    
    Subtype_HET                     = C{5}(:, 2); % PVNH/SUBC: 25/2
    Lesion_label_HET                = C{6};
    
    LabelPath_HET = '/host/gypsy/local_raid/seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label/';
    fclose(fid);
    
end

%% read PMG demographics
for read_PMG = 1
    
    fid = fopen('DemographicNGrouping_for_PMG_Data_3T.csv');
    C = textscan(fid, '%s%f%s%d%s%s','Delimiter',',','CollectOutput', 1);   
    Codes_PMG                       = C{1};
    Age_PMG                         = C{2};     % mean/sd: 28.8/11.2
    Sex_PMG                         = C{3};     % male/female: 12/9      
    Left_PMG                        = C{4} == 1; 
    Right_PMG                       = C{4} == 2;
    Bilateral_PMG                   = C{4} == 3;  % Left/Right/Bilateral: 6/6/9  
    Path_PMG                        = C{5}(:, 1);
    
    fclose(fid); 
    
end

%% configure AAL parcellation
for AAL_configuration = 1
    
    %% Tzourui-Mazoyer et al. NIMG 2002
    %% The original automated anatomical labeling (AAL) system consists of 90 ROIs including cerebral and subcortical areas.
    %% In this study, we excluded 6 subcortical areas (i.e., Hippocampus, Amygdala, Caudate nucleus, Putamen, Pallidum and Thalamus).
    %% The final number of parcels in this study is 78 (=90-2*6).
    
    %% Loading the AAL parcellation scheme ...
    AAL_surf_data_both = SurfStatReadData1('aal_both_rsl_final.txt');
    AAL_surf_data_left = AAL_surf_data_both(1, 1:40962);
    AAL_surf_data_right = AAL_surf_data_both(1, 40963:81924);
    AAL_surf_data_both = [ AAL_surf_data_left AAL_surf_data_left + 1];
    AAL_surf_data_both(find(AAL_surf_data_right==0)+40962) = 0;
    
    %% anatomical description and labels ...
    
    %% frontal lateral lobe
    %                                                     1                                      2                                3                                        4
    %                                                     40                                     41                               42                                       43
    frontal_lobe_lateral_surface_labels = {'superior_frontal_gyrus_dorsolateral_part', 'middle_frontal_gyrus', 'inferior_frontal_gyrus_opercular_part', 'inrerior_frontal_gyrus_triangular_part'};
    frontal_lobe_lateral_surface_L_index = [              3                          ,           7           ,                    11                  ,                   13                    ];
    frontal_lobe_lateral_surface_R_index = [              4                          ,           8           ,                    12                  ,                   14                    ];
    
    %% frontal medial lobe
    %                                                     5                                 6                        7
    %                                                     44                                45                       46
    frontal_lobe_medial_surface_labels = {'superior_frontal_gyrus_medial_part', 'supplementary_motor_area', 'paracentral_lobule'};
    frontal_lobe_medial_surface_L_index = [               23                  ,             19            ,         69          ];
    frontal_lobe_medial_surface_R_index = [               24                  ,             20            ,         70          ];
    
    %% frontal orbital lobe
    %                                                     8                                      9                                10                                       11                                   12                  13
    %                                                     47                                     48                               49                                       50                                   51                  52
    frontal_lobe_orbital_surace_labels = {'superior_frontal_gyrus_orbital_part', 'superior_frontal_gyrus_medial_orbital_part', 'middle_frontal_gyrus_orbital_part', 'inferior_frontal_gyrus_orbital_part', 'gyrus_rectus', 'olfactory_cortex'};
    frontal_lobe_orbital_surace_L_index = [               5                    ,                    25                       ,                  9                 ,                  15                  ,       27      ,          21       ];
    frontal_lobe_orbital_surace_R_index = [               6                    ,                    26                       ,                  10                ,                  16                  ,       28      ,          22       ];
    
    %% central 
    %                                14                  15                    16
    %                                53                  54                    55
    central_region_labels = {'precentral_gyrus', 'postcentral_gyrus', 'rolandic_operculum'};
    central_region_L_index = [       1         ,        57          ,         17          ];
    central_region_R_index = [       2         ,        58          ,         18          ];
    
    %% parietal lateral lobe
    %                                                   17                                   18                                   19                 20
    %                                                   56                                   57                                   58                 59
    parietal_lobe_lateral_surface_labels = {'superior_parietal_gyrus', 'inferior_parietal_supramarginal_and_angular_gyri', 'angular_gyrus', 'supramarginal_gyrus'};
    parietal_lobe_lateral_surface_L_index = [          59            ,                       61                          ,       65       ,          63          ];
    parietal_lobe_lateral_surface_R_index = [          60            ,                       62                          ,       66       ,          64          ];
    
    %% parietal medial lobe
    %                                           21
    %                                           60
    parietal_lobe_medial_surface_labels = {'precuneuos'};
    parietal_lobe_medial_surface_L_index = [    67     ];
    parietal_lobe_medial_surface_R_index = [    68     ];
    
    %% occipital lateral lobe
    %                                                    22                        23                         24
    %                                                    61                        62                         63
    occipital_lobe_lateral_surface_labels = {'superior_occipital_gyrus', 'middle_occipital_gyrus', 'inferior_occipital_gyrus'};
    occipital_lobe_lateral_surface_L_index = [          49             ,           51            ,            53             ];
    occipital_lobe_lateral_surface_R_index = [          50             ,           52            ,            54             ];
    
    %% occipital medial lobe
    %                                                   25                         26                             27                28
    %                                                   64                         65                             66                67
    occipital_lobe_medial_inferior_surface_labels = {'cuneus', 'calcarine_fissure_and_surrounding_cortex', 'lingual_gyrus', 'fusiform_gyrus'};
    occipital_lobe_medial_inferior_surface_L_index = [  45   ,                     43                    ,        47      ,         55      ];
    occipital_lobe_medial_inferior_surface_R_index = [  46   ,                     44                    ,        48      ,         56      ];
    
    %% temporal lobe
    %                                                   29                   30                  31                         32
    %                                                   68                   69                  70                         71
    temporal_lobe_lateral_surface_labels = {'superior_temporal_gyrus', 'heschl_gyrus', 'middle_temporal_gyrus', 'inferior_temporal_gyrus'};
    temporal_lobe_lateral_surface_L_index = [           81           ,       79      ,           85           ,              89          ];
    temporal_lobe_lateral_surface_R_index = [           82           ,       80      ,           86           ,              90          ];
    
    %% limbic
    %                                     33                                        34                                       35                                         36                                     37                    38
    %                                     72                                        73                                       74                                         75                                     76                    77
    limbic_lobe_labels = {'temporal_pole_superior_temporal_gyrus', 'temporal_pole_middle_temporal_gyrus', 'anterior_cingulate_and_paracingulate_gyri', 'median_cingulate_and_paracingulate_gyri', 'posterior_cingulate_gyrus', 'parahippocampus'};
    limbic_lobe_L_index = [                 83                   ,                   87                 ,                    31                      ,                    33                    ,              35            ,         39       ];
    limbic_lobe_R_index = [                 84                   ,                   88                 ,                    32                      ,                    34                    ,              36            ,         40       ];
    
    %% insula
    %                   39
    %                   78
    insula_labels = {'insula'};
    insula_L_index = [  29   ];
    insula_R_index = [  30   ];    
    
    %% excluded structures
    sub_cortical_gray_nuclei_labels  = {'hippocampus', 'amygdala', 'caudate_nucleus', 'lenticular_nucleus_putamen', 'lenticular_nucleus_pallidum', 'thalamus'};
    sub_cortical_gray_nuclei_L_index = [   37,            41    ,        71        ,             73              ,               75             ,     77    ];
    sub_cortical_gray_nuclei_R_index = [   38,            42    ,        72        ,             74              ,               76             ,     78    ];
    
    %% parcel structure
    Anatomical_Label_structs_temp = {'frontal_lobe_lateral_surface',           frontal_lobe_lateral_surface_labels, frontal_lobe_lateral_surface_L_index, frontal_lobe_lateral_surface_R_index; ...
                                     'frontal_lobe_medial_surface',            frontal_lobe_medial_surface_labels, frontal_lobe_medial_surface_L_index, frontal_lobe_medial_surface_R_index; ...
                                     'frontal_lobe_orbital_surace',            frontal_lobe_orbital_surace_labels, frontal_lobe_orbital_surace_L_index, frontal_lobe_orbital_surace_R_index; ...
                                     'insula',                                 insula_labels, insula_L_index, insula_R_index; ...
                                     'central_region',                         central_region_labels, central_region_L_index, central_region_R_index; ...                                     
                                     'parietal_lobe_lateral_surface',          parietal_lobe_lateral_surface_labels, parietal_lobe_lateral_surface_L_index, parietal_lobe_lateral_surface_R_index; ...
                                     'parietal_lobe_medial_surface',           parietal_lobe_medial_surface_labels, parietal_lobe_medial_surface_L_index, parietal_lobe_medial_surface_R_index; ...
                                     'occipital_lobe_lateral_surface',         occipital_lobe_lateral_surface_labels, occipital_lobe_lateral_surface_L_index, occipital_lobe_lateral_surface_R_index; ...
                                     'occipital_lobe_medial_inferior_surface', occipital_lobe_medial_inferior_surface_labels, occipital_lobe_medial_inferior_surface_L_index, occipital_lobe_medial_inferior_surface_R_index; ...
                                     'temporal_lobe_lateral_surface',          temporal_lobe_lateral_surface_labels, temporal_lobe_lateral_surface_L_index, temporal_lobe_lateral_surface_R_index; ...
                                     'limbic_lobe',                            limbic_lobe_labels, limbic_lobe_L_index, limbic_lobe_R_index; ...                                     
                                     'sub_cortical_gray_nuclei',               sub_cortical_gray_nuclei_labels, sub_cortical_gray_nuclei_L_index, sub_cortical_gray_nuclei_R_index;};
                                 
                                 
    Anatomical_Label_structs_idx = [ frontal_lobe_lateral_surface_L_index, ... 
                                     frontal_lobe_medial_surface_L_index, ... 
                                     frontal_lobe_orbital_surace_L_index, ... 
                                     insula_L_index, ...
                                     central_region_L_index, ... 
                                     parietal_lobe_lateral_surface_L_index, ... 
                                     parietal_lobe_medial_surface_L_index, ... 
                                     occipital_lobe_lateral_surface_L_index, ... 
                                     occipital_lobe_medial_inferior_surface_L_index, ... 
                                     temporal_lobe_lateral_surface_L_index, ... 
                                     limbic_lobe_L_index, ...
                                       ];
                                 
    Anatomical_Label_structs_idx = [ Anatomical_Label_structs_idx ...
                                   [ frontal_lobe_lateral_surface_R_index, ... 
                                     frontal_lobe_medial_surface_R_index, ... 
                                     frontal_lobe_orbital_surace_R_index, ... 
                                     insula_R_index, ...
                                     central_region_R_index, ... 
                                     parietal_lobe_lateral_surface_R_index, ... 
                                     parietal_lobe_medial_surface_R_index, ... 
                                     occipital_lobe_lateral_surface_R_index, ... 
                                     occipital_lobe_medial_inferior_surface_R_index, ... 
                                     temporal_lobe_lateral_surface_R_index, ... 
                                     limbic_lobe_R_index, ...
                                      ]];
    
    
    field = {'name', 'subarea', 'AAL_L_index', 'AAL_R_index'};
    Anatomical_Label_structs = cell2struct(Anatomical_Label_structs_temp, field, 2);
    
    %% Loading microstructural parcellation scheme ...
    %% Refer /local_raid/seokjun/01_project/Parcellation_Network_Analysis/
    using_whichparcellation = 3;
    if(using_whichparcellation == 1)
        filename = 'surf_reg_model_left_parcellation_without_sulcal_line.txt';
        parcellation_map_without_sulcal_line = SurfStatReadData1(filename);
        microAAL_surf_data_both = parcellation_map_with_sulcal_line;
    elseif(using_whichparcellation == 2)
        filename = 'surf_reg_model_left_parcellation_with_sulcal_line.txt';
        parcellation_map_with_sulcal_line = SurfStatReadData1(filename);
        microAAL_surf_data_both = parcellation_map_without_sulcal_line;
    elseif(using_whichparcellation == 3)
        filename = 'surf_reg_model_left_parcellation_kmeans.txt';
        parcellation_map_kmeans = SurfStatReadData1(filename);
        microAAL_surf_data_both = [ parcellation_map_kmeans max(parcellation_map_kmeans) + parcellation_map_kmeans ];
        microAAL_surf_data_both(find(microAAL_surf_data_both==0)+40962) = 0;
    end
    
end

%% load matfiles
for load_files = 1
    
    if(lesion_exclusion)
%         load( [ MATDIR '/01_generate_association_matrix_MCD_maskingout_lesion.mat' ] );
%         load( [ MATDIR '/02_small_world_analysis_result_MCD_maskingout_lesion.mat' ] );
    else
        load( [ MATDIR '/07_generate_association_matrix_MCD1_3T_HET_with_NODULE.mat' ] );
        %load( [ MATDIR '/01_generate_micro_association_matrix_MCD2_3T.mat' ] );
        load( [ MATDIR '/02_small_world_analysis_result_MCD_3T_16idx.mat' ] ); 
        %load( [ MATDIR '/02_small_world_analysis_result_MCD_micro.mat' ] ); 
    end
    
end

%% compute nodal degree (i.e., degree, strength, CC, Lp, global efficiency, local efficiency and small-worldness)
for compute_nodal_measure = 1
    
    for original_AAL = 1
        
        iternation_num = 1000;
        interval_thres = 0.01;
        myspar = 0.05:interval_thres:0.4;
        AAL_index_iter = length(Anatomical_Label_structs_idx);
        
        % 1) nodal degrees and strengths:   unweighted positive network
        for degree_strength = 1
            
            node_degree_uw_pos_cont             = zeros(length(myspar), AAL_index_iter);            
            node_degree_uw_pos_FCD_Type_II      = zeros(length(myspar), AAL_index_iter);
            node_degree_uw_pos_HET              = zeros(length(myspar), AAL_index_iter);
            node_degree_uw_pos_PMG              = zeros(length(myspar), AAL_index_iter);
            
            node_strength_uw_pos_cont             = zeros(length(myspar), AAL_index_iter);            
            node_strength_uw_pos_FCD_Type_II      = zeros(length(myspar), AAL_index_iter);
            node_strength_uw_pos_HET            = zeros(length(myspar), AAL_index_iter);
            node_strength_uw_pos_PMG            = zeros(length(myspar), AAL_index_iter);
            
            node_fully_connected                = ones(1, length(myspar));
            
            for i = 1 : length(myspar)
                
                %% controls
                wmatrix                                  = threshold_proportional(zadjacent_matrix_control, myspar(i));
                node_degree_uw_pos_cont(i, :)            = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_cont(i, :)          = strengths_und(double(wmatrix>0));               
                
                %% FCD Type-II
                wmatrix                                  = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
                node_degree_uw_pos_FCD_Type_II(i, :)     = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_FCD_Type_II(i, :)   = strengths_und(double(wmatrix>0));
                
                %% HET
                wmatrix                                  = threshold_proportional(zadjacent_matrix_HET, myspar(i));
                node_degree_uw_pos_HET(i, :)             = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_HET(i, :)           = strengths_und(double(wmatrix>0));
                
                %% PMG
                wmatrix                                  = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
                node_degree_uw_pos_PMG(i, :)             = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_PMG(i, :)           = strengths_und(double(wmatrix>0));
                
                %% check if fully connected
                if((sum(node_degree_uw_pos_cont(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_II(i, :)==0) + sum(node_degree_uw_pos_HET(i, :)==0) + sum(node_degree_uw_pos_PMG(i, :)==0)) > 0)
                    node_fully_connected(i) = 0;
                end
                
            end
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_degree_uw_pos_cont, 2) mean(node_degree_uw_pos_FCD_Type_II, 2) mean(node_degree_uw_pos_HET, 2) mean(node_degree_uw_pos_PMG, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
                
                figure; plot(1:length(myspar), [ mean(node_strength_uw_pos_cont, 2) mean(node_strength_uw_pos_FCD_Type_II, 2) mean(node_strength_uw_pos_HET, 2) mean(node_strength_uw_pos_PMG, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 2) nodal clustering coefficient:	unweighted positive network
        for CC = 1
            
            node_cc_uw_pos_cont        = zeros(length(myspar), AAL_index_iter);            
            node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_HET         = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_PMG         = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                              = threshold_proportional(zadjacent_matrix_control, myspar(i));
                node_cc_uw_pos_cont(i, :)            = clustering_coef_bu(double(wmatrix>0))';
                                
                %% FCD Type II
                wmatrix                              = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
                node_cc_uw_pos_FCD_Type_II(i, :)     = clustering_coef_bu(double(wmatrix>0))';
                
                %% HET
                wmatrix                              = threshold_proportional(zadjacent_matrix_HET, myspar(i));
                node_cc_uw_pos_HET(i, :)             = clustering_coef_bu(double(wmatrix>0))';
                
                %% PMG
                wmatrix                              = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
                node_cc_uw_pos_PMG(i, :)             = clustering_coef_bu(double(wmatrix>0))';
                
            end
            
        end
        
        % 3) nodal shortest path length:    unweighted positive network
        for Lp = 1
            
            % Characteristic path length using BCT
            node_spl_uw_pos_cont            = zeros(length(myspar), AAL_index_iter);            
            node_spl_uw_pos_FCD_Type_II     = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_HET             = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_PMG             = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                             = threshold_proportional(zadjacent_matrix_control, myspar(i));
                node_spl_uw_pos_cont(i, :)          = mean(distance_bin(double(wmatrix>0)))';
               
                %% FCD Type II
                wmatrix                             = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
                node_spl_uw_pos_FCD_Type_II(i, :)   = mean(distance_bin(double(wmatrix>0)))';
                
                %% HET
                wmatrix                             = threshold_proportional(zadjacent_matrix_HET, myspar(i));
                node_spl_uw_pos_HET(i, :)           = mean(distance_bin(double(wmatrix>0)))';
                
                %% PMG
                wmatrix                             = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
                node_spl_uw_pos_PMG(i, :)           = mean(distance_bin(double(wmatrix>0)))';
                
            end
            
        end
        
        % 4) global efficiency:             unweighted positive network
        for glob_eff = 1
            
            node_cp_uw_pos_cont         = zeros(length(myspar), 2);           
            node_cp_uw_pos_FCD_Type_II  = zeros(length(myspar), 2);
            node_cp_uw_pos_HET          = zeros(length(myspar), 2);
            node_cp_uw_pos_PMG          = zeros(length(myspar), 2);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix = threshold_proportional(zadjacent_matrix_control, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                node_cp_uw_pos_cont(i, :) = [lamda efficiency];
               
                %% FCD Type II
                wmatrix = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                node_cp_uw_pos_FCD_Type_II(i, :)  = [lamda efficiency];
                
                %% HET
                wmatrix = threshold_proportional(zadjacent_matrix_HET, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                node_cp_uw_pos_HET(i, :) = [lamda efficiency];
                
                %% FCD Type II
                wmatrix = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                node_cp_uw_pos_PMG(i, :)  = [lamda efficiency];
                
            end
            
            clear('efficiency');
            
            if(printfigs)
                figure; plot(1:length(myspar), [ node_cp_uw_pos_cont(:, 2) node_cp_uw_pos_FCD_Type_II(:, 2) node_cp_uw_pos_HET(:, 2) node_cp_uw_pos_PMG(:, 2) ]');
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 5) local efficiency:              unweighted positive network
        for local_eff = 1
            
            node_le_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);          
            node_le_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_le_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_le_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                              = threshold_proportional(zadjacent_matrix_control, myspar(i));
                node_le_uw_pos_cont(i, :)            = efficiency(double(wmatrix>0), 1)';
                               
                %% FCD Type II
                wmatrix                              = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
                node_le_uw_pos_FCD_Type_II(i, :)     = efficiency(double(wmatrix>0), 1)';
                
                %% HET
                wmatrix                              = threshold_proportional(zadjacent_matrix_HET, myspar(i));
                node_le_uw_pos_HET(i, :)             = efficiency(double(wmatrix>0), 1)';
                
                %% PMG
                wmatrix                              = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
                node_le_uw_pos_PMG(i, :)             = efficiency(double(wmatrix>0), 1)';
                
            end
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_le_uw_pos_cont, 2) mean(node_le_uw_pos_FCD_Type_II, 2) mean(node_le_uw_pos_HET, 2) mean(node_le_uw_pos_PMG, 2) ]);
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 6) small world parameters without random network comparison (quick running)
        for small_world_ver1 = 1
            
            parpool(24); 
            
            n = 1000;
            sw_cont_set = cell(1, length(myspar));            
            sw_FCD_Type_II_set = cell(1, length(myspar));
            sw_HET_set = cell(1, length(myspar));
            sw_PMG_set = cell(1, length(myspar));
            
            parfor i = 1 : length(myspar)
                
                i
                
                %% control
                wmatrix          = threshold_proportional(zadjacent_matrix_control, myspar(i));
                sw_cont          = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% FCD Type II
                wmatrix          = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
                sw_FCD_Type_II   = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% HET
                wmatrix          = threshold_proportional(zadjacent_matrix_HET, myspar(i));
                sw_HET           = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% PMG
                wmatrix          = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
                sw_PMG           = smallworldness_bu(double(wmatrix>0), 1, n);
                
                sw_cont_set{i}          = sw_cont;                
                sw_FCD_Type_II_set{i}   = sw_FCD_Type_II;
                sw_HET_set{i}           = sw_HET;
                sw_PMG_set{i}           = sw_PMG;
                
            end
            
            sw_sigma  = [];
            sw_Lp     = [];
            sw_Cp     = [];
            sw_lambda = [];
            sw_gamma  = [];
            
            for i = 1 : length(myspar)
                sw_sigma  = [ sw_sigma;  sw_cont_set{i}.Sigma  sw_FCD_Type_II_set{i}.Sigma sw_HET_set{i}.Sigma sw_PMG_set{i}.Sigma     ];
                sw_Lp     = [ sw_Lp;     sw_cont_set{i}.Lp     sw_FCD_Type_II_set{i}.Lp    sw_HET_set{i}.Lp    sw_PMG_set{i}.Lp        ];
                sw_Cp     = [ sw_Cp;  sw_cont_set{i}.Cp     sw_FCD_Type_II_set{i}.Cp    sw_HET_set{i}.Cp    sw_PMG_set{i}.Cp           ];
                sw_lambda = [ sw_lambda; sw_cont_set{i}.Lambda sw_FCD_Type_II_set{i}.Lambda sw_HET_set{i}.Lambda sw_PMG_set{i}.Lambda ];
                sw_gamma  = [ sw_gamma;  sw_cont_set{i}.Gamma  sw_FCD_Type_II_set{i}.Gamma sw_HET_set{i}.Gamma sw_PMG_set{i}.Gamma     ];
            end
            
            if(printfigs)
                graph_metrics = sw_gamma;
                figure; plot(1:length(myspar), graph_metrics');
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 7) small world parameters with random network comparison (complete running)
        for small_world_ver2 = 1
            
            % density = 0.05 - 0.4 with 0.01 interval and 1000 randomization
            parpool(24);
            option.homogeneity_thres = [0.3, 0.3];
            option.bagging = 1;

            %% control vs. FCD Type II with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_FCD_Type_II, ClustCoff_corr_FCD_Type_II, Sigma_corr_FCD_Type_II, rerandPathLen_FCD_Type_II, rerandClustCoff_FCD_Type_II ] = ...
                swdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);            
            
            %% control vs. Heterotopia with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_HET, ClustCoff_corr_HET, Sigma_corr_HET, rerandPathLen_HET, rerandClustCoff_HET ] = ...
                swdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% control vs. Polymicrogria with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_PMG, ClustCoff_corr_PMG, Sigma_corr_PMG, rerandPathLen_PMG, rerandClustCoff_PMG ] = ...
                swdiff_sparsity_parfor(asso_mat_cont_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% FCD Type II vs. HET with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_patient1, ClustCoff_corr_patient1, Sigma_corr_patient1, rerandPathLen_patien1, rerandClustCoff_patient1 ] = ...
                swdiff_sparsity_parfor(asso_mat_FCD_Type_II_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% FCD Type II vs. PMG with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_patient2, ClustCoff_corr_patient2, Sigma_corr_patient2, rerandPathLen_patient2, rerandClustCoff_patient2 ] = ...
                swdiff_sparsity_parfor(asso_mat_FCD_Type_II_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% visualization
            if(visualization)
                
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                MarkerSize = 20;
                
                %% Path length
                for path_legnth = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([4 42]); ylim([1 6]);
                    plot(myspar*100, [ mean([ PathLen_corr_FCD_Type_II.Lp.Lp1; PathLen_corr_HET.Lp.Lp1; PathLen_corr_PMG.Lp.Lp1; ], 1);
                        PathLen_corr_FCD_Type_II.Lp.Lp2;
                        PathLen_corr_HET.Lp.Lp2;
                        PathLen_corr_PMG.Lp.Lp2 ], 'LineWidth', 2);
                    color_lines = get(gca, 'ColorOrder');
                    templine = get(gca, 'Children');
                    set(templine(3), 'Color', color_lines(3, :));
                    set(templine(2), 'Color', color_lines(4, :));
                    set(templine(1), 'Color', color_lines(5, :));
                    
                    scatter(myspar*100', [ mean([PathLen_corr_FCD_Type_II.Lp.Lp1; PathLen_corr_HET.Lp.Lp1; PathLen_corr_PMG.Lp.Lp1; ], 1)' ], MarkerSize, color_lines(1, :), 'filled', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));                    
                    scatter(myspar*100', [ PathLen_corr_FCD_Type_II.Lp.Lp2' ], MarkerSize, color_lines(3, :), 'filled', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                    scatter(myspar*100', [ PathLen_corr_HET.Lp.Lp2' ], MarkerSize, color_lines(4, :), 'filled', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                    scatter(myspar*100', [ PathLen_corr_PMG.Lp.Lp2' ], MarkerSize, color_lines(5, :), 'filled', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                    legend('control', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                    
                    exportfigbo(gcf,[OUTPATH 'path_length_global_diff_cont_vs_pats_3T_16idx.png'], 'png', 13); close(gcf);
                    
                    %% 2) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -2.2;
                    feature1 = PathLen_corr_FCD_Type_II;
                    feature2 = rerandPathLen_FCD_Type_II;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, feature2.Lp.Lp2-feature2.Lp.Lp1, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    ylim([-2.5 2.5]);
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_FCD_TypeII_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% 3) delta, HET
                    colorline_idx = 4;
                    sig_pos = -0.4;
                    feature1 = PathLen_corr_HET;
                    feature2 = rerandPathLen_HET;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, feature2.Lp.Lp2-feature2.Lp.Lp1, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    ylim([-2 2]);
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_HET_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -1;
                    feature1 = PathLen_corr_PMG;
                    feature2 = rerandPathLen_PMG;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, feature2.Lp.Lp2-feature2.Lp.Lp1, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    ylim([-2.5 2.5]);
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',20);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_PMG_3T_20idx.png'], 'png', 13); close(gcf);
                    
                end
                
                %% Clustering coefficient
                for clust_coeff = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([4 42]); ylim([0.25 0.7]);
                    plot(myspar*100, [ mean([ ClustCoff_corr_FCD_Type_II.Cp.Cp1; ClustCoff_corr_HET.Cp.Cp1; ClustCoff_corr_PMG.Cp.Cp1; ], 1);                        
                        ClustCoff_corr_FCD_Type_II.Cp.Cp2;
                        ClustCoff_corr_HET.Cp.Cp2;
                        ClustCoff_corr_PMG.Cp.Cp2 ], 'LineWidth', 2);
                    color_lines = get(gca, 'ColorOrder');
                    templine = get(gca, 'Children');
                    set(templine(3), 'Color', color_lines(3, :));
                    set(templine(2), 'Color', color_lines(4, :));
                    set(templine(1), 'Color', color_lines(5, :));
                    
                    scatter(myspar*100', [ mean([ ClustCoff_corr_FCD_Type_II.Cp.Cp1; ClustCoff_corr_HET.Cp.Cp1; ClustCoff_corr_PMG.Cp.Cp1; ], 1)' ], MarkerSize, color_lines(1, :), 'filled', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));                    
                    scatter(myspar*100', [ ClustCoff_corr_FCD_Type_II.Cp.Cp2' ], MarkerSize, color_lines(3, :), 'filled', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                    scatter(myspar*100', [ ClustCoff_corr_HET.Cp.Cp2' ], MarkerSize, color_lines(4, :), 'filled', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                    scatter(myspar*100', [ ClustCoff_corr_PMG.Cp.Cp2' ], MarkerSize, color_lines(5, :), 'filled', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                    legend('control', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                    
                    exportfigbo(gcf,[OUTPATH 'clustering_coeff_global_diff_cont_vs_pats_3T_16idx.png'], 'png', 13); close(gcf);
                           
                    %% 2) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -0.15;
                    feature1 = ClustCoff_corr_FCD_Type_II;
                    feature2 = rerandClustCoff_FCD_Type_II;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, feature2.Cp.Cp2-feature2.Cp.Cp1, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.25 0.25]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_FCD_TypeII_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, HET
                    colorline_idx = 4;
                    sig_pos = -0.15;
                    feature1 = ClustCoff_corr_HET;
                    feature2 = rerandClustCoff_HET;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, feature2.Cp.Cp2-feature2.Cp.Cp1, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.2 0.2]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_HET_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% 5) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -0.17;
                    feature1 = ClustCoff_corr_PMG;
                    feature2 = rerandClustCoff_PMG;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, feature2.Cp.Cp2-feature2.Cp.Cp1, 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.2 0.25]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_PMG_3T_20idx.png'], 'png', 13); close(gcf);
                    
                end
                
            end
            
        end
        
        % 8) Lp and Cp at nodal level
        for Lp_Cp_nodal = 1
            
            parpool(24);
            
            %% control vs. FCD Type II
            [ myspar, PathLen_corr_FCD_Type_II_nodal, ClustCoff_corr_FCD_Type_II_nodal, rerandPathLen_FCD_Type_II_nodal, rerandClustCoff_FCD_Type_II_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Heterotopia
            [myspar, PathLen_corr_HET_nodal, ClustCoff_corr_HET_nodal, rerandPathLen_HET_nodal, rerandClustCoff_HET_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Polymicrogria
            [myspar, PathLen_corr_PMG_nodal, ClustCoff_corr_PMG_nodal, rerandPathLen_PMG_nodal, rerandClustCoff_PMG_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% FCD Type II vs. HET
            [myspar, PathLen_corr_patient1_nodal, ClustCoff_corr_patient1_nodal, rerandPathLen_patien1_nodal, rerandClustCoff_patient1_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_FCD_Type_II_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% FCD Type II vs. PMG
            [myspar, PathLen_corr_patient2_nodal, ClustCoff_corr_patient2_nodal, rerandPathLen_patient2_nodal, rerandClustCoff_patient2_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_FCD_Type_II_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% visualization
            if(visualization)
                
                %% variable setup
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                
                spar            = min(find(node_fully_connected));
                uncorrp         = 0.05;
                FDRQ            = 0.05;
                diff_measure_all = [ (PathLen_corr_FCD_Type_II_nodal.Lp.Lp2(:, spar)' - PathLen_corr_FCD_Type_II_nodal.Lp.Lp1(:, spar)')  ...
                                     (PathLen_corr_HET_nodal.Lp.Lp2(:, spar)' - PathLen_corr_HET_nodal.Lp.Lp1(:, spar)')  ...
                                     (PathLen_corr_PMG_nodal.Lp.Lp2(:, spar)' - PathLen_corr_PMG_nodal.Lp.Lp1(:, spar)')  ...
                                   ];
                
                %% node position of each parcel in AAL
                centroid = zeros(1, length(Anatomical_Label_structs_idx));
                for i = 1 : length(Anatomical_Label_structs_idx)
                    
                    vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                    vert_coord = S.coord(:, vert_idx);
                    [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                    centroid(i) = vert_idx(b);
                    
                end
                
                nodePos = S.coord(:, centroid);
                
                %% path length
                for path_length = 1
                    
                    %% control vs. FCD Type II
                    graph_feature     = PathLen_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = zeros(78,78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_FCD_Type_II_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = PathLen_corr_HET_nodal;
                    diff_real_wmatrix = zeros(78,78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                                        
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_HET_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = PathLen_corr_PMG_nodal;
                    diff_real_wmatrix = zeros(78,78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_PMG_3T_20idx.png'], 'png', 13); close(gcf);
                    
                end
                
                %% clustering coefficient
                for clustering_coeff = 1
                    
                    %% control vs. FCD Type II
                    graph_feature     = ClustCoff_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = zeros(78,78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_FCD_Type_II_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = ClustCoff_corr_HET_nodal;
                    diff_real_wmatrix = zeros(78,78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_HET_3T_20idx.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = ClustCoff_corr_PMG_nodal;
                    diff_real_wmatrix = zeros(78,78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_PMG_3T_20idx.png'], 'png', 13); close(gcf);
                    
                end
                
            end
            
        end
        
        % 9) save all resultant variables
        for save_files = 1
            
            %         if(lesion_exclusion)
            %             filename = [ MATDIR '/01_generate_association_matrix_MCD_maskingout_lesion.mat' ];
            %
            %         else
            %             filename = [ MATDIR '/01_generate_association_matrix_MCD.mat' ];
            %         end
            %
            %         save(filename, ...
            %             'node_cc_uw_pos_cont', 'node_cc_uw_pos_FCD_Type_I', 'node_cc_uw_pos_FCD_Type_II', 'node_cc_uw_pos_HET', 'node_cc_uw_pos_PMG', ...
            %             'node_cp_uw_pos_cont', 'node_cp_uw_pos_FCD_Type_I', 'node_cp_uw_pos_FCD_Type_II', 'node_cp_uw_pos_HET', 'node_cp_uw_pos_PMG', ...
            %             'node_degree_uw_pos_cont', 'node_degree_uw_pos_FCD_Type_I', 'node_degree_uw_pos_FCD_Type_II', 'node_degree_uw_pos_HET', 'node_degree_uw_pos_PMG', ...
            %             'node_le_uw_pos_cont', 'node_le_uw_pos_FCD_Type_I', 'node_le_uw_pos_FCD_Type_II', 'node_le_uw_pos_HET', 'node_le_uw_pos_PMG', ...
            %             'node_spl_uw_pos_cont', 'node_spl_uw_pos_FCD_Type_'I, 'node_spl_uw_pos_FCD_Type_II', 'node_spl_uw_pos_HET', 'node_spl_uw_pos_PMG', ...
            %             'node_strength_uw_pos_cont', 'node_strength_uw_pos_FCD_Type_I', 'node_strength_uw_pos_FCD_Type_II', 'node_strength_uw_pos_HET', 'node_strength_uw_pos_PMG', ...
            %             'sw_cont_set', 'sw_FCD_Type_I_set', 'sw_FCD_Type_II_set', 'sw_HET_set', 'sw_PMG_set', ...
            %             'PathLen_corr_FCD_Type_II', 'PathLen_corr_FCD_Type_I', 'PathLen_corr_HET', 'PathLen_corr_PMG', 'PathLen_corr_patient1', 'PathLen_corr_patient2', ...
            %             'ClustCoff_corr_FCD_Type_I', 'ClustCoff_corr_FCD_Type_II', 'ClustCoff_corr_HET', 'ClustCoff_corr_patient1', 'ClustCoff_corr_patient2', 'ClustCoff_corr_PMG', ...
            %             'Sigma_corr_FCD_Type_I', 'Sigma_corr_FCD_Type_II', 'Sigma_corr_HET', 'Sigma_corr_patient1', 'Sigma_corr_patient2', 'Sigma_corr_PMG', ...
            %             'rerandClustCoff_FCD_Type_II', 'rerandClustCoff_FCD_Type_I', 'rerandClustCoff_HET', 'rerandClustCoff_patient1', 'rerandClustCoff_patient2', 'rerandClustCoff_PMG', ...
            %             'rerandPathLen_FCD_Type_I', 'rerandPathLen_FCD_Type_II', 'rerandPathLen_HET', 'rerandPathLen_patien1', 'rerandPathLen_patient2', 'rerandPathLen_PMG');
            
            
        end
        
    end
    
    for micro_AAL = 1
        
        iteration_num = 100;
        interval_thres = 0.01;
        myspar = 0.05:interval_thres:0.40;
        AAL_index_iter = max(microAAL_surf_data_both);
        
        % 1) nodal degrees and strengths:   unweighted positive network
        for degree_strength = 1
            
            micro_node_degree_uw_pos_cont             = zeros(length(myspar), AAL_index_iter);
            micro_node_degree_uw_pos_FCD_Type_I       = zeros(length(myspar), AAL_index_iter);
            micro_node_degree_uw_pos_FCD_Type_II      = zeros(length(myspar), AAL_index_iter);
            micro_node_degree_uw_pos_HET              = zeros(length(myspar), AAL_index_iter);
            micro_node_degree_uw_pos_PMG              = zeros(length(myspar), AAL_index_iter);
            
            micro_node_strength_uw_pos_cont           = zeros(length(myspar), AAL_index_iter);
            micro_node_strength_uw_pos_FCD_Type_I     = zeros(length(myspar), AAL_index_iter);
            micro_node_strength_uw_pos_FCD_Type_II    = zeros(length(myspar), AAL_index_iter);
            micro_node_strength_uw_pos_HET            = zeros(length(myspar), AAL_index_iter);
            micro_node_strength_uw_pos_PMG            = zeros(length(myspar), AAL_index_iter);
            
            micro_node_fully_connected                = ones(1, length(myspar));
            
            for i = 1 : length(myspar)
                
                %% controls
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
                micro_node_degree_uw_pos_cont(i, :)      = degrees_und(double(wmatrix>0));
                micro_node_strength_uw_pos_cont(i, :)    = strengths_und(double(wmatrix>0));
                
                %% FCD Type-I
                wmatrix                                         = threshold_proportional(micro_zadjacent_matrix_FCD_Type_I, myspar(i));
                micro_node_degree_uw_pos_FCD_Type_I(i, :)       = degrees_und(double(wmatrix>0));
                micro_node_strength_uw_pos_FCD_Type_I(i, :)     = strengths_und(double(wmatrix>0));
                
                %% FCD Type-II
                wmatrix                                         = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
                micro_node_degree_uw_pos_FCD_Type_II(i, :)      = degrees_und(double(wmatrix>0));
                micro_node_strength_uw_pos_FCD_Type_II(i, :)    = strengths_und(double(wmatrix>0));
                
                %% HET
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
                micro_node_degree_uw_pos_HET(i, :)       = degrees_und(double(wmatrix>0));
                micro_node_strength_uw_pos_HET(i, :)     = strengths_und(double(wmatrix>0));
                
                %% PMG
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
                micro_node_degree_uw_pos_PMG(i, :)       = degrees_und(double(wmatrix>0));
                micro_node_strength_uw_pos_PMG(i, :)     = strengths_und(double(wmatrix>0));
                
                %% check if fully connected
                if((sum(node_degree_uw_pos_cont(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_I(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_II(i, :)==0) + sum(node_degree_uw_pos_HET(i, :)==0) + sum(node_degree_uw_pos_PMG(i, :)==0)) > 0)
                    micro_node_fully_connected(i) = 0;
                end
                
            end
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(micro_node_degree_uw_pos_cont, 2) mean(micro_node_degree_uw_pos_FCD_Type_I, 2) mean(micro_node_degree_uw_pos_FCD_Type_II, 2) mean(micro_node_degree_uw_pos_HET, 2) mean(micro_node_degree_uw_pos_PMG, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
                
                figure; plot(1:length(myspar), [ mean(micro_node_strength_uw_pos_cont, 2) mean(micro_node_strength_uw_pos_FCD_Type_I, 2) mean(micro_node_strength_uw_pos_FCD_Type_II, 2) mean(micro_node_strength_uw_pos_HET, 2) mean(micro_node_strength_uw_pos_PMG, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 2) nodal clustering coefficient:	unweighted positive network
        for CC = 1
            
            micro_node_cc_uw_pos_cont        = zeros(length(myspar), AAL_index_iter);
            micro_node_cc_uw_pos_FCD_Type_I  = zeros(length(myspar), AAL_index_iter);
            micro_node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            micro_node_cc_uw_pos_HET         = zeros(length(myspar), AAL_index_iter);
            micro_node_cc_uw_pos_PMG         = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
                micro_node_cc_uw_pos_cont(i, :)            = clustering_coef_bu(double(wmatrix>0))';
                
                %% FCD Type I
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_FCD_Type_I, myspar(i));
                micro_node_cc_uw_pos_FCD_Type_I(i, :)      = clustering_coef_bu(double(wmatrix>0))';
                
                %% FCD Type II
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
                micro_node_cc_uw_pos_FCD_Type_II(i, :)     = clustering_coef_bu(double(wmatrix>0))';
                
                %% HET
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
                micro_node_cc_uw_pos_HET(i, :)             = clustering_coef_bu(double(wmatrix>0))';
                
                %% PMG
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
                micro_node_cc_uw_pos_PMG(i, :)             = clustering_coef_bu(double(wmatrix>0))';
                
            end
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(micro_node_cc_uw_pos_cont, 2) mean(micro_node_cc_uw_pos_FCD_Type_I, 2) mean(micro_node_cc_uw_pos_FCD_Type_II, 2) mean(micro_node_cc_uw_pos_HET, 2) mean(micro_node_cc_uw_pos_PMG, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 3) nodal shortest path length:    unweighted positive network
        for Lp = 1
            
            % Characteristic path length using BCT
            micro_node_spl_uw_pos_cont            = zeros(length(myspar), AAL_index_iter);
            micro_node_spl_uw_pos_FCD_Type_I      = zeros(length(myspar), AAL_index_iter);
            micro_node_spl_uw_pos_FCD_Type_II     = zeros(length(myspar), AAL_index_iter);
            micro_node_spl_uw_pos_HET             = zeros(length(myspar), AAL_index_iter);
            micro_node_spl_uw_pos_PMG             = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                             = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
                micro_node_spl_uw_pos_cont(i, :)          = mean(distance_bin(double(wmatrix>0)))';
                
                %% FCD Type I
                wmatrix                             = threshold_proportional(micro_zadjacent_matrix_FCD_Type_I, myspar(i));
                micro_node_spl_uw_pos_FCD_Type_I(i, :)    = mean(distance_bin(double(wmatrix>0)))';
                
                %% FCD Type II
                wmatrix                             = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
                micro_node_spl_uw_pos_FCD_Type_II(i, :)   = mean(distance_bin(double(wmatrix>0)))';
                
                %% HET
                wmatrix                             = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
                micro_node_spl_uw_pos_HET(i, :)           = mean(distance_bin(double(wmatrix>0)))';
                
                %% PMG
                wmatrix                             = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
                micro_node_spl_uw_pos_PMG(i, :)           = mean(distance_bin(double(wmatrix>0)))';
                
            end
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(micro_node_spl_uw_pos_cont, 2) mean(micro_node_spl_uw_pos_FCD_Type_I, 2) mean(micro_node_spl_uw_pos_FCD_Type_II, 2) mean(micro_node_spl_uw_pos_HET, 2) mean(micro_node_spl_uw_pos_PMG, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 4) global efficiency:             unweighted positive network
        for glob_eff = 1
            
            micro_node_cp_uw_pos_cont         = zeros(length(myspar), 2);
            micro_node_cp_uw_pos_FCD_Type_I   = zeros(length(myspar), 2);
            micro_node_cp_uw_pos_FCD_Type_II  = zeros(length(myspar), 2);
            micro_node_cp_uw_pos_HET          = zeros(length(myspar), 2);
            micro_node_cp_uw_pos_PMG          = zeros(length(myspar), 2);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                micro_node_cp_uw_pos_cont(i, :) = [lamda efficiency];
                
                %% FCD Type I
                wmatrix = threshold_proportional(micro_zadjacent_matrix_FCD_Type_I, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                micro_node_cp_uw_pos_FCD_Type_I(i, :) = [lamda efficiency];
                
                %% FCD Type II
                wmatrix = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                micro_node_cp_uw_pos_FCD_Type_II(i, :)  = [lamda efficiency];
                
                %% HET
                wmatrix = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                micro_node_cp_uw_pos_HET(i, :) = [lamda efficiency];
                
                %% FCD Type II
                wmatrix = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
                [lamda efficiency] = charpath(distance_bin(double(wmatrix>0)));
                micro_node_cp_uw_pos_PMG(i, :)  = [lamda efficiency];
                
            end
            
            clear('efficiency');
            
            if(printfigs)
                figure; plot(1:length(myspar), [ micro_node_cp_uw_pos_cont(:, 2) micro_node_cp_uw_pos_FCD_Type_I(:, 2) micro_node_cp_uw_pos_FCD_Type_II(:, 2) micro_node_cp_uw_pos_HET(:, 2) micro_node_cp_uw_pos_PMG(:, 2) ]');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 5) local efficiency:              unweighted positive network
        for local_eff = 1
            
            micro_node_le_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);
            micro_node_le_uw_pos_FCD_Type_I   = zeros(length(myspar), AAL_index_iter);
            micro_node_le_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            micro_node_le_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            micro_node_le_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
                micro_node_le_uw_pos_cont(i, :)            = efficiency(double(wmatrix>0), 1)';
                
                %% FCD Type I
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_FCD_Type_I, myspar(i));
                micro_node_le_uw_pos_FCD_Type_I(i, :)      = efficiency(double(wmatrix>0), 1)';
                
                %% FCD Type II
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
                micro_node_le_uw_pos_FCD_Type_II(i, :)     = efficiency(double(wmatrix>0), 1)';
                
                %% HET
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
                micro_node_le_uw_pos_HET(i, :)             = efficiency(double(wmatrix>0), 1)';
                
                %% PMG
                wmatrix                              = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
                micro_node_le_uw_pos_PMG(i, :)             = efficiency(double(wmatrix>0), 1)';
                
            end
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(micro_node_le_uw_pos_cont, 2) mean(micro_node_le_uw_pos_FCD_Type_I, 2) mean(micro_node_le_uw_pos_FCD_Type_II, 2) mean(micro_node_le_uw_pos_HET, 2) mean(micro_node_le_uw_pos_PMG, 2) ]);
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 6) small world parameters without random network comparison (quick running)
        for small_world_ver1 = 1
            
            n = 1000;
            micro_sw_cont_set = [];
            micro_sw_FCD_Type_I_set = [];
            micro_sw_FCD_Type_II_set = [];
            micro_sw_HET_set = [];
            micro_sw_PMG_set = [];
            
            for i = 1 : length(myspar)
                
                i
                
                %% control
                wmatrix          = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
                micro_sw_cont          = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% FCD Type I
                wmatrix          = threshold_proportional(micro_zadjacent_matrix_FCD_Type_I, myspar(i));
                micro_sw_FCD_Type_I    = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% FCD Type II
                wmatrix          = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
                micro_sw_FCD_Type_II   = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% HET
                wmatrix          = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
                micro_sw_HET           = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% PMG
                wmatrix          = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
                micro_sw_PMG           = smallworldness_bu(double(wmatrix>0), 1, n);
                
                %% output
                micro_sw_cont_set         = [ micro_sw_cont_set         micro_sw_cont ];
                micro_sw_FCD_Type_I_set   = [ micro_sw_FCD_Type_I_set   micro_sw_FCD_Type_I ];
                micro_sw_FCD_Type_II_set  = [ micro_sw_FCD_Type_II_set  micro_sw_FCD_Type_II ];
                micro_sw_HET_set          = [ micro_sw_HET_set          micro_sw_HET ];
                micro_sw_PMG_set          = [ micro_sw_PMG_set          micro_sw_PMG ];
                
            end
            
            micro_sw_sigma  = [];
            micro_sw_Lp     = [];
            micro_sw_Cp     = [];
            micro_sw_lambda = [];
            micro_sw_gamma  = [];
            
            for i = 1 : length(myspar)
                micro_sw_sigma  = [ micro_sw_sigma;  micro_sw_cont_set(i).Sigma  micro_sw_FCD_Type_I_set(i).Sigma micro_sw_FCD_Type_II_set(i).Sigma micro_sw_HET_set(i).Sigma micro_sw_PMG_set(i).Sigma     ];
                micro_sw_Lp     = [ micro_sw_Lp;     micro_sw_cont_set(i).Lp     micro_sw_FCD_Type_I_set(i).Lp    micro_sw_FCD_Type_II_set(i).Lp    micro_sw_HET_set(i).Lp    micro_sw_PMG_set(i).Lp        ];
                micro_sw_Cp     = [ micro_sw_Cp;  micro_sw_cont_set(i).Cp     micro_sw_FCD_Type_I_set(i).Cp    micro_sw_FCD_Type_II_set(i).Cp    micro_sw_HET_set(i).Cp    micro_sw_PMG_set(i).Cp           ];
                micro_sw_lambda = [ micro_sw_lambda; micro_sw_cont_set(i).Lambda micro_sw_FCD_Type_I_set(i).Lambda micro_sw_FCD_Type_II_set(i).Lambda micro_sw_HET_set(i).Lambda micro_sw_PMG_set(i).Lambda ];
                micro_sw_gamma  = [ micro_sw_gamma;  micro_sw_cont_set(i).Gamma  micro_sw_FCD_Type_I_set(i).Gamma micro_sw_FCD_Type_II_set(i).Gamma micro_sw_HET_set(i).Gamma micro_sw_PMG_set(i).Gamma     ];
            end
            
            if(printfigs)
                graph_metrics = micro_sw_Lp;
                figure; plot(1:length(myspar), graph_metrics');
                legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 7) small world parameters with random network comparison (complete running)
        for small_world_ver2 = 1
            
            parpool(24);
            
            % density = 0.05 - 0.4 with 0.01 interval and 1000 randomization
            
            %% control vs. FCD Type I with computing normalized parameters (i.e., lamda, gamma, sigma)
            [ myspar, micro_PathLen_corr_FCD_Type_I, micro_ClustCoff_corr_FCD_Type_I, micro_Sigma_corr_FCD_Type_I, micro_rerandPathLen_FCD_Type_I, micro_rerandClustCoff_FCD_Type_I ] = ...
                swdiff_sparsity_parfor_sonic(micro_asso_mat_cont_res_val, micro_asso_mat_FCD_Type_I_res_val, myspar(1), myspar(end), interval_thres, iteration_num , 1);
            
            %% control vs. FCD Type II with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, micro_PathLen_corr_FCD_Type_II, micro_ClustCoff_corr_FCD_Type_II, micro_Sigma_corr_FCD_Type_II, micro_rerandPathLen_FCD_Type_II, micro_rerandClustCoff_FCD_Type_II ] = ...
                swdiff_sparsity_parfor_sonic(micro_asso_mat_cont_res_val, micro_asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, iteration_num , 1);
            
            %% control vs. Heterotopia with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, micro_PathLen_corr_HET, micro_ClustCoff_corr_HET, micro_Sigma_corr_HET, micro_rerandPathLen_HET, micro_rerandClustCoff_HET ] = ...
                swdiff_sparsity_parfor_sonic(micro_asso_mat_cont_res_val, micro_asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iteration_num , 1);
            
            %% control vs. Polymicrogria with computing normalized parameters (i.e., lamda, gamma, sigma)
            option = [0.3, 0.3];
            [myspar, micro_PathLen_corr_PMG, micro_ClustCoff_corr_PMG, micro_Sigma_corr_PMG, micro_rerandPathLen_PMG, micro_rerandClustCoff_PMG ] = ...
                swdiff_sparsity_parfor_sonic(micro_asso_mat_cont_res_val, micro_asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iteration_num , 1, option);
            
            %% FCD Type II vs. HET with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, micro_PathLen_corr_patient1, micro_ClustCoff_corr_patient1, micro_Sigma_corr_patient1, micro_rerandPathLen_patien1, micro_rerandClustCoff_patient1 ] = ...
                swdiff_sparsity_parfor_sonic(micro_asso_mat_FCD_Type_II_res_val, micro_asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iteration_num , 1);
            
            %% FCD Type II vs. PMG with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, micro_PathLen_corr_patient2, micro_ClustCoff_corr_patient2, micro_Sigma_corr_patient2, micro_rerandPathLen_patient2, micro_rerandClustCoff_patient2 ] = ...
                swdiff_sparsity_parfor_sonic(micro_asso_mat_FCD_Type_II_res_val, micro_asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iteration_num , 1);
            
            %% visualization
            if(visualization)
                
                OUTPATH = '/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/micro_AAL/';
                MarkerSize = 20;
                
                %% Path length
                for path_legnth = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([4 42]);
                    plot(myspar*100, [ mean([ micro_PathLen_corr_FCD_Type_I.Lp.Lp1; micro_PathLen_corr_FCD_Type_II.Lp.Lp1; micro_PathLen_corr_HET.Lp.Lp1; micro_PathLen_corr_PMG.Lp.Lp1; ], 1);
                        micro_PathLen_corr_FCD_Type_I.Lp.Lp2;
                        micro_PathLen_corr_FCD_Type_II.Lp.Lp2;
                        micro_PathLen_corr_HET.Lp.Lp2;
                        micro_PathLen_corr_PMG.Lp.Lp2 ], 'LineWidth', 2);
                    color_lines = get(gca, 'ColorOrder');
                    scatter(myspar*100', [ mean([ micro_PathLen_corr_FCD_Type_I.Lp.Lp1; micro_PathLen_corr_FCD_Type_II.Lp.Lp1; micro_PathLen_corr_HET.Lp.Lp1; micro_PathLen_corr_PMG.Lp.Lp1; ], 1)' ], ...
                                                                                        MarkerSize, color_lines(1, :), 'filled', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                    scatter(myspar*100', [ micro_PathLen_corr_FCD_Type_I.Lp.Lp2' ],     MarkerSize, color_lines(2, :), 'filled', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                    scatter(myspar*100', [ micro_PathLen_corr_FCD_Type_II.Lp.Lp2' ],    MarkerSize, color_lines(3, :), 'filled', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                    scatter(myspar*100', [ micro_PathLen_corr_HET.Lp.Lp2' ],            MarkerSize, color_lines(4, :), 'filled', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                    scatter(myspar*100', [ micro_PathLen_corr_PMG.Lp.Lp2' ],            MarkerSize, color_lines(5, :), 'filled', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                    legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                    
                    exportfigbo(gcf,[OUTPATH 'path_length_global_diff_cont_vs_pats_micro.png'], 'png', 13); close(gcf);
                    
                    %% 2) delta, FCD Type I
                    colorline_idx = 2;
                    sig_pos = -0.2;
                    feature1 = micro_PathLen_corr_FCD_Type_I;
                    feature2 = micro_rerandPathLen_FCD_Type_I;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Lp.Lp2')'-cell2mat(feature2.Lp.Lp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.35 0.35]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_FCD_TypeI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 3) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -0.17;
                    feature1 = micro_PathLen_corr_FCD_Type_II;
                    feature2 = micro_rerandPathLen_FCD_Type_II;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Lp.Lp2')'-cell2mat(feature2.Lp.Lp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    ylim([-0.21 0.21]);
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_FCD_TypeII_micro.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, HET
                    colorline_idx = 4;
                    sig_pos = -0.15;
                    feature1 = micro_PathLen_corr_HET;
                    feature2 = micro_rerandPathLen_HET;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Lp.Lp2')'-cell2mat(feature2.Lp.Lp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    ylim([-0.2 0.2]);
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_HET_micro.png'], 'png', 13); close(gcf);
                    
                    %% 5) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -0.6;
                    feature1 = micro_PathLen_corr_PMG;
                    feature2 = micro_rerandPathLen_PMG;
                    FDRp = FDR(feature1.Lp.Lp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Lp.Lp2')'-cell2mat(feature2.Lp.Lp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Lp.upper95Lp12diff; -feature1.Lp.lower95Lp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Lp.rerand_Lp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Lp.rerand_Lp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Lp.Lp2-feature1.Lp.Lp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Lp.Lp2-feature1.Lp.Lp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    ylim([-0.7 0.7]);
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_PMG_micro.png'], 'png', 13); close(gcf);
                    
                end
                
                %% Clustering coefficient
                for clust_coeff = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([10 42]); ylim([0.25 0.7]);
                    plot(myspar*100, [ mean([ micro_ClustCoff_corr_FCD_Type_I.Cp.Cp1; micro_ClustCoff_corr_FCD_Type_II.Cp.Cp1; micro_ClustCoff_corr_HET.Cp.Cp1; micro_ClustCoff_corr_PMG.Cp.Cp1; ], 1);
                        micro_ClustCoff_corr_FCD_Type_I.Cp.Cp2;
                        micro_ClustCoff_corr_FCD_Type_II.Cp.Cp2;
                        micro_ClustCoff_corr_HET.Cp.Cp2;
                        micro_ClustCoff_corr_PMG.Cp.Cp2 ], 'LineWidth', 2);
                    color_lines = get(gca, 'ColorOrder');
                    scatter(myspar*100', [ mean([ micro_ClustCoff_corr_FCD_Type_I.Cp.Cp1; micro_ClustCoff_corr_FCD_Type_II.Cp.Cp1; micro_ClustCoff_corr_HET.Cp.Cp1; micro_ClustCoff_corr_PMG.Cp.Cp1; ], 1)' ], MarkerSize, color_lines(1, :), 'filled', 'MarkerFaceColor', color_lines(1, :), 'MarkerEdgeColor', color_lines(1, :));
                    scatter(myspar*100', [ micro_ClustCoff_corr_FCD_Type_I.Cp.Cp2' ], MarkerSize, color_lines(2, :), 'filled', 'MarkerFaceColor', color_lines(2, :), 'MarkerEdgeColor', color_lines(2, :));
                    scatter(myspar*100', [ micro_ClustCoff_corr_FCD_Type_II.Cp.Cp2' ], MarkerSize, color_lines(3, :), 'filled', 'MarkerFaceColor', color_lines(3, :), 'MarkerEdgeColor', color_lines(3, :));
                    scatter(myspar*100', [ micro_ClustCoff_corr_HET.Cp.Cp2' ], MarkerSize, color_lines(4, :), 'filled', 'MarkerFaceColor', color_lines(4, :), 'MarkerEdgeColor', color_lines(4, :));
                    scatter(myspar*100', [ micro_ClustCoff_corr_PMG.Cp.Cp2' ], MarkerSize, color_lines(5, :), 'filled', 'MarkerFaceColor', color_lines(5, :), 'MarkerEdgeColor', color_lines(5, :));
                    legend('control', 'FCD Type-I', 'FCD Type-II', 'HET', 'PMG', 'location', 'best');
                    
                    exportfigbo(gcf,[OUTPATH 'clustering_coeff_global_diff_cont_vs_pats_micro.png'], 'png', 13); close(gcf);
                    
                    %% 2) delta, FCD Type I
                    colorline_idx = 2;
                    sig_pos = -0.07;
                    feature1 = micro_ClustCoff_corr_FCD_Type_I;
                    feature2 = micro_rerandClustCoff_FCD_Type_I;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Cp.Cp2')'-cell2mat(feature2.Cp.Cp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.1 0.15]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_FCD_TypeI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 3) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -0.15;
                    feature1 = micro_ClustCoff_corr_FCD_Type_II;
                    feature2 = micro_rerandClustCoff_FCD_Type_II;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Cp.Cp2')'-cell2mat(feature2.Cp.Cp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.2 0.2]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_FCD_TypeII_micro.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, HET
                    colorline_idx = 4;
                    sig_pos = -0.15;
                    feature1 = micro_ClustCoff_corr_HET;
                    feature2 = micro_rerandClustCoff_HET;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Cp.Cp2')'-cell2mat(feature2.Cp.Cp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.2 0.2]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_HET_micro.png'], 'png', 13); close(gcf);
                    
                    %% 5) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -0.17;
                    feature1 = micro_ClustCoff_corr_PMG;
                    feature2 = micro_rerandClustCoff_PMG;
                    FDRp = FDR(feature1.Cp.Cp12_diff_p, 0.05);
                    if(isempty(FDRp))
                        FDRp = 0;
                    end
                    
                    figure; hold on;
                    plot(myspar*100, cell2mat(feature2.Cp.Cp2')'-cell2mat(feature2.Cp.Cp1')', 'Color',[0.8 0.8 0.8],'LineWidth',0.02);
                    plot(myspar*100, [ -feature1.Cp.upper95Cp12diff; -feature1.Cp.lower95Cp12diff], '--k', 'LineWidth', 1);
                    plot(myspar*100, -feature1.Cp.rerand_Cp12diff_mean, 'k', 'LineWidth', 2);
                    scatter((myspar*100)', -feature1.Cp.rerand_Cp12diff_mean', MarkerSize, 'k', 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
                    plot(myspar*100, feature1.Cp.Cp2-feature1.Cp.Cp1, 'Color', color_lines(colorline_idx, :), 'LineWidth', 2);
                    scatter(myspar*100', feature1.Cp.Cp2-feature1.Cp.Cp1, MarkerSize, color_lines(colorline_idx, :), 'filled', 'MarkerFaceColor', color_lines(colorline_idx, :), 'MarkerEdgeColor', color_lines(colorline_idx, :));
                    for i = 1: length(myspar)
                        if(feature1.Cp.Cp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Cp.Cp12_diff_p(i)>FDRp & feature1.Cp.Cp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-0.2 0.25]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_PMG_micro.png'], 'png', 13); close(gcf);
                    
                end
                
            end
            
        end
        
        % 8) Lp and Cp at nodal level
        for Lp_Cp_nodal = 1
            
            parpool(24);
            
            %% control vs. FCD Type I
            [ myspar, PathLen_corr_FCD_Type_I_nodal, ClustCoff_corr_FCD_Type_I_nodal, rerandPathLen_FCD_Type_I_nodal, rerandClustCoff_FCD_Type_I_nodal ] =  ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_FCD_Type_I_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. FCD Type II
            [ myspar, PathLen_corr_FCD_Type_II_nodal, ClustCoff_corr_FCD_Type_II_nodal, rerandPathLen_FCD_Type_II_nodal, rerandClustCoff_FCD_Type_II_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_FCD_Type_II_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Heterotopia
            [myspar, PathLen_corr_HET_nodal, ClustCoff_corr_HET_nodal, rerandPathLen_HET_nodal, rerandClustCoff_HET_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Polymicrogria
            [myspar, PathLen_corr_PMG_nodal, ClustCoff_corr_PMG_nodal, rerandPathLen_PMG_nodal, rerandClustCoff_PMG_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_cont_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% FCD Type II vs. HET
            [myspar, PathLen_corr_patient1_nodal, ClustCoff_corr_patient1_nodal, rerandPathLen_patien1_nodal, rerandClustCoff_patient1_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_FCD_Type_II_res_val, asso_mat_HET_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% FCD Type II vs. PMG
            [myspar, PathLen_corr_patient2_nodal, ClustCoff_corr_patient2_nodal, rerandPathLen_patient2_nodal, rerandClustCoff_patient2_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor(asso_mat_FCD_Type_II_res_val, asso_mat_PMG_res_val, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% visualization
            if(visualization)
                
                %% variable setup
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                
                spar            = 1; % min(find(node_fully_connected));
                uncorrp         = 0.05;
                FDRQ            = 0.05;
                diff_measure_all = [ (PathLen_corr_FCD_Type_I_nodal.Lp.Lp2(:, spar)' - PathLen_corr_FCD_Type_I_nodal.Lp.Lp1(:, spar)')  ...
                    (PathLen_corr_FCD_Type_II_nodal.Lp.Lp2(:, spar)' - PathLen_corr_FCD_Type_II_nodal.Lp.Lp1(:, spar)')  ...
                    (PathLen_corr_HET_nodal.Lp.Lp2(:, spar)' - PathLen_corr_HET_nodal.Lp.Lp1(:, spar)')  ...
                    (PathLen_corr_PMG_nodal.Lp.Lp2(:, spar)' - PathLen_corr_PMG_nodal.Lp.Lp1(:, spar)')  ...
                    ];
                
                %% node position of each parcel in AAL
                centroid = zeros(1, length(Anatomical_Label_structs_idx));
                for i = 1 : length(Anatomical_Label_structs_idx)
                    
                    vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                    vert_coord = S.coord(:, vert_idx);
                    [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                    centroid(i) = vert_idx(b);
                    
                end
                
                nodePos = S.coord(:, centroid);
                
                %% path length
                for path_length = 1
                    
                    %% control vs. FCD Type I
                    graph_feature     = PathLen_corr_FCD_Type_I_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_FCD_TypeI;
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_FCD_Type_I.png'], 'png', 13); close(gcf);
                    
                    %% control vs. FCD Type II
                    graph_feature     = PathLen_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_FCD_TypeII;
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_FCD_Type_II.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = PathLen_corr_HET_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_HET;
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    diff_measure(diff_measure==max(diff_measure)) = 0.9;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_HET.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = PathLen_corr_PMG_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_PMG;
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_PMG.png'], 'png', 13); close(gcf);
                    
                end
                
                %% clustering coefficient
                for clustering_coeff = 1
                    
                    %% control vs. FCD Type I
                    graph_feature     = ClustCoff_corr_FCD_Type_I_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_FCD_TypeI;
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_FCD_Type_I.png'], 'png', 13); close(gcf);
                    
                    %% control vs. FCD Type II
                    graph_feature     = ClustCoff_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_FCD_TypeII;
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_FCD_Type_II.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = ClustCoff_corr_HET_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_HET;
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_HET.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = ClustCoff_corr_PMG_nodal;
                    diff_real_wmatrix = diff_real_wmatrix_cont_vs_PMG;
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_PMG.png'], 'png', 13); close(gcf);
                    
                end
                
            end
            
        end
        
        % 9) save all resultant variables
        for save_files = 1
            
            %         if(lesion_exclusion)
            %             filename = [ MATDIR '/01_generate_association_matrix_MCD_maskingout_lesion.mat' ];
            %
            %         else
            %             filename = [ MATDIR '/01_generate_association_matrix_MCD.mat' ];
            %         end
            %
            %         save(filename, ...
            %             'node_cc_uw_pos_cont', 'node_cc_uw_pos_FCD_Type_I', 'node_cc_uw_pos_FCD_Type_II', 'node_cc_uw_pos_HET', 'node_cc_uw_pos_PMG', ...
            %             'node_cp_uw_pos_cont', 'node_cp_uw_pos_FCD_Type_I', 'node_cp_uw_pos_FCD_Type_II', 'node_cp_uw_pos_HET', 'node_cp_uw_pos_PMG', ...
            %             'node_degree_uw_pos_cont', 'node_degree_uw_pos_FCD_Type_I', 'node_degree_uw_pos_FCD_Type_II', 'node_degree_uw_pos_HET', 'node_degree_uw_pos_PMG', ...
            %             'node_le_uw_pos_cont', 'node_le_uw_pos_FCD_Type_I', 'node_le_uw_pos_FCD_Type_II', 'node_le_uw_pos_HET', 'node_le_uw_pos_PMG', ...
            %             'node_spl_uw_pos_cont', 'node_spl_uw_pos_FCD_Type_'I, 'node_spl_uw_pos_FCD_Type_II', 'node_spl_uw_pos_HET', 'node_spl_uw_pos_PMG', ...
            %             'node_strength_uw_pos_cont', 'node_strength_uw_pos_FCD_Type_I', 'node_strength_uw_pos_FCD_Type_II', 'node_strength_uw_pos_HET', 'node_strength_uw_pos_PMG', ...
            %             'sw_cont_set', 'sw_FCD_Type_I_set', 'sw_FCD_Type_II_set', 'sw_HET_set', 'sw_PMG_set', ...
            %             'PathLen_corr_FCD_Type_II', 'PathLen_corr_FCD_Type_I', 'PathLen_corr_HET', 'PathLen_corr_PMG', 'PathLen_corr_patient1', 'PathLen_corr_patient2', ...
            %             'ClustCoff_corr_FCD_Type_I', 'ClustCoff_corr_FCD_Type_II', 'ClustCoff_corr_HET', 'ClustCoff_corr_patient1', 'ClustCoff_corr_patient2', 'ClustCoff_corr_PMG', ...
            %             'Sigma_corr_FCD_Type_I', 'Sigma_corr_FCD_Type_II', 'Sigma_corr_HET', 'Sigma_corr_patient1', 'Sigma_corr_patient2', 'Sigma_corr_PMG', ...
            %             'rerandClustCoff_FCD_Type_II', 'rerandClustCoff_FCD_Type_I', 'rerandClustCoff_HET', 'rerandClustCoff_patient1', 'rerandClustCoff_patient2', 'rerandClustCoff_PMG', ...
            %             'rerandPathLen_FCD_Type_I', 'rerandPathLen_FCD_Type_II', 'rerandPathLen_HET', 'rerandPathLen_patien1', 'rerandPathLen_patient2', 'rerandPathLen_PMG');
            
            
        end
        
    end
    
end