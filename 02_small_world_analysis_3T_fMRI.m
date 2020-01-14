clear all
close all
lesion_exclusion = 0;

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/'));
addpath(genpath('/host/gypsy/local_raid/seokjun/03_downloads/GRETNA/'));
rmpath(genpath('/host/gypsy/local_raid/seokjun/01_project/05_NoelRest/04_Script/toolboxes/spm8/'));
P='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/';
MATDIR='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/01_matfiles/';
OUTPATH='/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/00_association_matrix/';
printfigs = 0; % don't flag on unless you really want to save the figures as .png, otherwise your precious time is gone...

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
        load( [ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_Age_Sex1.mat' ] );
        load( [ MATDIR '/02_small_world_analysis_result_MCD_rsfMRI.mat' ] ); 
        load( [ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_Age_Sex1_micro.mat' ]);
        load( [ MATDIR '/02_small_world_analysis_result_MCD_rsfMRI_micro.mat' ] ); 
    end
    
end

%% compute nodal degree (i.e., degree, strength, CC, Lp, global efficiency, local efficiency and small-worldness)
for compute_nodal_measure = 1
    
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
        
        iternation_num = 1000;
        interval_thres = 0.01;
        myspar = 0.05:interval_thres:0.4;
        AAL_index_iter = length(Anatomical_Label_structs_idx);
        mat_thres = 0.75;
        
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
%                  mat_cont_temp(:, :, j) = temp;
            end
            averaged_mat_cont(:, :, i) = (sum(mat_cont_temp, 3)/num_of_cont)>mat_thres;
%             averaged_mat_cont(:, :, i) = mean(mat_cont_temp, 3)>mat_thres;
            
            %% FCD Type II
            mat_FCD_Type_II_temp = mat_FCD_Type_II*0;
            for j = 1 : num_of_FCD_Type_II
                temp = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                mat_FCD_Type_II_temp(:, :, j) = temp > 0;
%                 mat_FCD_Type_II_temp(:, :, j) = temp;
            end
            averaged_mat_FCD_Type_II(:, :, i) = (sum(mat_FCD_Type_II_temp, 3)/num_of_FCD_Type_II)>mat_thres;
%             averaged_mat_FCD_Type_II(:, :, i) = mean(mat_FCD_Type_II_temp, 3)>mat_thres;
            
            %% HET
            mat_HET_temp = mat_HET*0;
            for j = 1 : num_of_HET
                temp = threshold_proportional(mat_HET(:, :, j), myspar(i));
                mat_HET_temp(:, :, j) = temp > 0;
%                 mat_HET_temp(:, :, j) = temp;
            end
            averaged_mat_HET(:, :, i) = (sum(mat_HET_temp, 3)/num_of_HET)>mat_thres;
%             averaged_mat_HET(:, :, i) = mean(mat_HET_temp, 3)>mat_thres;
            
            %% PMG
            mat_PMG_temp = mat_PMG*0;
            for j = 1 : num_of_PMG
                temp = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                mat_PMG_temp(:, :, j) = temp > 0;
%                 mat_PMG_temp(:, :, j) = temp;
            end
            averaged_mat_PMG(:, :, i) = (sum(mat_PMG_temp, 3)/num_of_PMG)>mat_thres;
%             averaged_mat_PMG(:, :, i) = mean(mat_PMG_temp, 3)>mat_thres;
                 
        end
        
        % 1) nodal clustering coefficient:	unweighted positive network
        for CC = 1
            
            node_cc_uw_pos_cont         = zeros(length(myspar), AAL_index_iter, num_of_cont);            
            node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter, num_of_FCD_Type_II);
            node_cc_uw_pos_HET          = zeros(length(myspar), AAL_index_iter, num_of_HET);
            node_cc_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter, num_of_PMG);
            
            for i = 1 : length(myspar)
                
                %% control
                for j = 1 : num_of_cont
                    wmatrix                                  = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    node_cc_uw_pos_cont(i, :, j)             = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_cont(i, :, j)             = clustering_coef_wu(wmatrix)';
                end
                                
                %% FCD Type II
                for j = 1 : num_of_FCD_Type_II
                    wmatrix                                  = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    node_cc_uw_pos_FCD_Type_II(i, :, j)      = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_FCD_Type_II(i, :, j)      = clustering_coef_wu(wmatrix)';
                end
                
                %% HET
                for j = 1 : num_of_HET
                    wmatrix                                  = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    node_cc_uw_pos_HET(i, :, j)              = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_HET(i, :, j)              = clustering_coef_wu(wmatrix)';
                end
                
                %% PMG
                for j = 1 : num_of_PMG
                    wmatrix                                  = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    node_cc_uw_pos_PMG(i, :, j)              = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_PMG(i, :, j)              = clustering_coef_wu(wmatrix)';
                end
                
            end
            
            node_cc_uw_pos_cont_final = zeros(length(myspar), num_of_cont);
            node_cc_uw_pos_FCD_Type_II_final = zeros(length(myspar), num_of_FCD_Type_II);
            node_cc_uw_pos_HET_final = zeros(length(myspar), num_of_HET);
            node_cc_uw_pos_PMG_final = zeros(length(myspar), num_of_PMG);
            
            for i = 1 : length(myspar)
                
                for j = 1 : num_of_cont
                    temp = node_cc_uw_pos_cont(i, :, j);
                    node_cc_uw_pos_cont_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_FCD_Type_II
                    temp = node_cc_uw_pos_FCD_Type_II(i, :, j);
                    node_cc_uw_pos_FCD_Type_II_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_HET
                    temp = node_cc_uw_pos_HET(i, :, j); 
                    node_cc_uw_pos_HET_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_PMG
                    temp = node_cc_uw_pos_PMG(i, :, j);
                    node_cc_uw_pos_PMG_final(i, j) = mean(temp);
                end
                
            end
            
            [ mean(node_cc_uw_pos_cont_final, 2) mean(node_cc_uw_pos_FCD_Type_II_final, 2) mean(node_cc_uw_pos_HET_final, 2) mean(node_cc_uw_pos_PMG_final, 2) ]
            [ std(node_cc_uw_pos_cont_final, 0, 2) std(node_cc_uw_pos_FCD_Type_II_final, 0, 2) std(node_cc_uw_pos_HET_final, 0, 2) std(node_cc_uw_pos_PMG_final, 0, 2) ]
            figure; plot(1:length(myspar), [ mean(node_cc_uw_pos_cont_final, 2) mean(node_cc_uw_pos_FCD_Type_II_final, 2) mean(node_cc_uw_pos_HET_final, 2) mean(node_cc_uw_pos_PMG_final, 2) ])
            legend('control', 'FCD', 'HET', 'PMG');
        end
        
        % 2) nodal shortest path length:    unweighted positive network
        for Lp = 1
            
            % Characteristic path length using BCT
            node_spl_uw_pos_cont            = zeros(length(myspar), AAL_index_iter, num_of_cont);            
            node_spl_uw_pos_FCD_Type_II     = zeros(length(myspar), AAL_index_iter, num_of_FCD_Type_II);
            node_spl_uw_pos_HET             = zeros(length(myspar), AAL_index_iter, num_of_HET);
            node_spl_uw_pos_PMG             = zeros(length(myspar), AAL_index_iter, num_of_PMG);
                        
            for i = 1 : length(myspar)
                
                %% control
                for j = 1 : num_of_cont
                    wmatrix                                 = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter * 2;
%                     [averLp Lpi]                            = gretna_node_shortestpathlength_weight(wmatrix);  
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);
                    node_spl_uw_pos_cont(i, :, j)           = Lpi;
%                     dist = distance_bin(double(wmatrix>0)); 
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_cont(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_cont(i, :, j)          = mean(dist);
                end
                
                %% FCD Type II
                for j = 1 : num_of_FCD_Type_II
                    wmatrix                                 = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter * 2;
%                     [averLp Lpi]                            = gretna_node_shortestpathlength_weight(wmatrix); 
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);                    
                    node_spl_uw_pos_FCD_Type_II(i, :, j)    = Lpi;
%                     dist = distance_bin(double(wmatrix>0));
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_FCD_Type_II(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_FCD_Type_II(i, :, j)   = mean(dist);
                end
                
                %% HET
                for j = 1 : num_of_HET
                    wmatrix                         = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0); Lpi(isinf(Lpi)) = AAL_index_iter * 2; 
%                     [averLp Lpi]                    = gretna_node_shortestpathlength_weight(wmatrix); 
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);
                    node_spl_uw_pos_HET(i, :, j)    = Lpi;
%                     dist = distance_bin(double(wmatrix>0)); 
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_HET(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_HET(i, :, j)           = mean(dist);
                end
                
                %% PMG
                for j = 1 : num_of_PMG
                    wmatrix                         = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter * 2; 
%                     [averLp Lpi]                    = gretna_node_shortestpathlength_weight(wmatrix); 
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);
                    node_spl_uw_pos_PMG(i, :, j)    = Lpi;
%                     dist = distance_bin(double(wmatrix>0)); 
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_PMG(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_PMG(i, :, j)           = mean(dist);
                end
                
            end
            
            node_spl_uw_pos_cont_final            = zeros(length(myspar), num_of_cont);            
            node_spl_uw_pos_FCD_Type_II_final     = zeros(length(myspar), num_of_FCD_Type_II);
            node_spl_uw_pos_HET_final             = zeros(length(myspar), num_of_HET);
            node_spl_uw_pos_PMG_final             = zeros(length(myspar), num_of_PMG);
            
            for i = 1 : length(myspar)
                
                for j = 1 : num_of_cont
                    wmatrix = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    temp    = node_spl_uw_pos_cont(i, :, j);
                    node_spl_uw_pos_cont_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_FCD_Type_II
                    wmatrix = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    temp = node_spl_uw_pos_FCD_Type_II(i, :, j);
                    node_spl_uw_pos_FCD_Type_II_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_HET
                    wmatrix = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    temp = node_spl_uw_pos_HET(i, :, j); 
                    node_spl_uw_pos_HET_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_PMG
                    wmatrix = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    temp = node_spl_uw_pos_PMG(i, :, j);
                    node_spl_uw_pos_PMG_final(i, j) = mean(temp);
                end
                
            end
            
            
            [ mean(node_spl_uw_pos_cont_final, 2) mean(node_spl_uw_pos_FCD_Type_II_final, 2) mean(node_spl_uw_pos_HET_final, 2) mean(node_spl_uw_pos_PMG_final, 2) ]
            [ std(node_spl_uw_pos_cont_final, 0, 2) std(node_spl_uw_pos_FCD_Type_II_final, 0, 2) std(node_spl_uw_pos_HET_final, 0, 2) std(node_spl_uw_pos_PMG_final, 0, 2) ]
            
            figure; plot(1:length(myspar), [ mean(node_spl_uw_pos_cont_final, 2) mean(node_spl_uw_pos_FCD_Type_II_final, 2) mean(node_spl_uw_pos_HET_final, 2) mean(node_spl_uw_pos_PMG_final, 2) ])
            legend('control', 'FCD', 'HET', 'PMG');
            
        end
        
        % 3) nodal clustering coefficient based on averaged matrix:	unweighted positive network
        for CC_avg_mat = 1
            
            node_cc_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                wmatrix                            = threshold_proportional(mean_cont, myspar(i));
                node_cc_uw_pos_cont(i, :)          = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_cont(i, :)          = clustering_coef_wu(wmatrix)';
                
                wmatrix                            = threshold_proportional(mean_FCD_Type_II, myspar(i));
                node_cc_uw_pos_FCD_Type_II(i, :)   = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_FCD_Type_II(i, :)   = clustering_coef_wu(wmatrix)';
                
                wmatrix                            = threshold_proportional(mean_HET, myspar(i));
                node_cc_uw_pos_HET(i, :)           = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_HET(i, :)           = clustering_coef_wu(wmatrix)';
                
                wmatrix                            = threshold_proportional(mean_PMG, myspar(i));
                node_cc_uw_pos_PMG(i, :)           = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_PMG(i, :)           = clustering_coef_wu(wmatrix)';
            end   
            
            [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ]
            [ std(node_cc_uw_pos_cont, 0, 2) std(node_cc_uw_pos_FCD_Type_II, 0, 2) std(node_cc_uw_pos_HET, 0, 2) std(node_cc_uw_pos_PMG, 0, 2) ]
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ])
                legend('control', 'FCD', 'HET', 'PMG');
            end
            
        end
        
        % 4) nodal shortest path length based on averaged matrix:    unweighted positive network
        for Lp_avg_mat = 1
            
            node_spl_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);          
            node_spl_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                            = threshold_proportional(mean_cont, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_cont(i, :)         = averLp;
                               
                %% FCD Type II
                wmatrix                            = threshold_proportional(mean_FCD_Type_II, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_FCD_Type_II(i, :)  = averLp;
                
                %% HET
                wmatrix                            = threshold_proportional(mean_HET, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_HET(i, :)          = averLp;
                
                %% PMG
                wmatrix                            = threshold_proportional(mean_PMG, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_PMG(i, :)          = averLp;
                
            end
            
            [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]            
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]);
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 5) nodal clustering coefficient based on thresholded averaged matrix:	unweighted positive network
        for CC_avg_mat_thresholded = 1
            
            node_cc_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                wmatrix                            = averaged_mat_cont(:, :, i);
                node_cc_uw_pos_cont(i, :)          = clustering_coef_bu(double(wmatrix>0))';
                
                wmatrix                            = averaged_mat_FCD_Type_II(:, :, i);
                node_cc_uw_pos_FCD_Type_II(i, :)   = clustering_coef_bu(double(wmatrix>0))';
                
                wmatrix                            = averaged_mat_HET(:, :, i);
                node_cc_uw_pos_HET(i, :)           = clustering_coef_bu(double(wmatrix>0))';
                
                wmatrix                            = averaged_mat_PMG(:, :, i);
                node_cc_uw_pos_PMG(i, :)           = clustering_coef_bu(double(wmatrix>0))';
            end   
            
            [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ]
            [ std(node_cc_uw_pos_cont, 0, 2) std(node_cc_uw_pos_FCD_Type_II, 0, 2) std(node_cc_uw_pos_HET, 0, 2) std(node_cc_uw_pos_PMG, 0, 2) ]
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ])
                legend('control', 'FCD', 'HET', 'PMG');
            end
            
        end
        
        % 6) nodal shortest path length based on thresholded averaged matrix:    unweighted positive network
        for Lp_avg_mat_thresholded = 1
            
            node_spl_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);          
            node_spl_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                            = averaged_mat_cont(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_cont(i, :)         = mean(Lpi);
                               
                %% FCD Type II
                wmatrix                            = averaged_mat_FCD_Type_II(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_FCD_Type_II(i, :)  = mean(Lpi);
                
                %% HET
                wmatrix                            = averaged_mat_HET(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_HET(i, :)          = mean(Lpi);
                
                %% PMG
                wmatrix                            = averaged_mat_PMG(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_PMG(i, :)          = mean(Lpi);
                
            end
            
            [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]);
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
                
        % 7) small world parameters with random network comparison (complete running)
        for small_world = 1
            
            % density = 0.05 - 0.4 with 0.01 interval and 1000 randomization
            parpool(24);
            option.homogeneity_thres = [5, 0.06];
            option.bagging = 1;
            iternation_num = 1000;
            interval_thres = 0.01;
            myspar = 0.05:interval_thres:0.4;
            AAL_index_iter = length(Anatomical_Label_structs_idx);

            %% control vs. FCD Type II with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_FCD_Type_II, ClustCoff_corr_FCD_Type_II, Sigma_corr_FCD_Type_II, rerandPathLen_FCD_Type_II, rerandClustCoff_FCD_Type_II ] = ...
                swdiff_sparsity_parfor_fMRI(mat_cont, mat_FCD_Type_II, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% control vs. Heterotopia with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_HET, ClustCoff_corr_HET, Sigma_corr_HET, rerandPathLen_HET, rerandClustCoff_HET ] = ...
                swdiff_sparsity_parfor_fMRI(mat_cont, mat_HET, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% control vs. Polymicrogria with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_PMG, ClustCoff_corr_PMG, Sigma_corr_PMG, rerandPathLen_PMG, rerandClustCoff_PMG ] = ...
                swdiff_sparsity_parfor_fMRI(mat_cont, mat_PMG, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
           
            %% visualization
            if(visualization)
                
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                MarkerSize = 20;
                
                %% Path length
                for path_legnth = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([5 42]); 
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
                    ylim([1 80]);
                    
                    exportfigbo(gcf,[OUTPATH 'path_length_global_diff_cont_vs_pats_fMRI.png'], 'png', 13); close(gcf);
                    
                    %% 2) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -20;
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
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-25 25]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_FCD_TypeII_fMRI.png'], 'png', 13); close(gcf);
                    
                    %% 3) delta, HET
                    colorline_idx = 4;
                    sig_pos = -20;
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
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    ylim([-25 25]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -25;
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

                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',20);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    ylim([-30 30]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
                    
                end
                
                %% Clustering coefficient
                for clust_coeff = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([5 42]); ylim([-0.03 0.5]);
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
                    
                    exportfigbo(gcf,[OUTPATH 'clustering_coeff_global_diff_cont_vs_pats_fMRI.png'], 'png', 13); close(gcf);
                           
                    %% 3) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -0.2;
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
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_FCD_TypeII_fMRI.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, HET
                    colorline_idx = 4;
                    sig_pos = -0.25;
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
                    ylim([-0.3 0.3]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
                    
                    %% 5) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -0.3;
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
                    ylim([-0.35 0.35]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
                    
                end
                
            end
            
        end
        
        % 8) Lp and Cp at nodal level
        for Lp_Cp_nodal = 1
            
            parpool(24);
             
            %% control vs. FCD Type II
            [ myspar, PathLen_corr_FCD_Type_II_nodal, ClustCoff_corr_FCD_Type_II_nodal, rerandPathLen_FCD_Type_II_nodal, rerandClustCoff_FCD_Type_II_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor_fMRI(mat_cont, mat_FCD_Type_II, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Heterotopia
            [myspar, PathLen_corr_HET_nodal, ClustCoff_corr_HET_nodal, rerandPathLen_HET_nodal, rerandClustCoff_HET_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor_fMRI(mat_cont, mat_HET, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Polymicrogria
            [myspar, PathLen_corr_PMG_nodal, ClustCoff_corr_PMG_nodal, rerandPathLen_PMG_nodal, rerandClustCoff_PMG_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor_fMRI(mat_cont, mat_PMG, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% visualization
            if(visualization)
                
                %% variable setup
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                
                spar            = 7; % min(find(node_fully_connected));
                uncorrp         = 0.05;
                FDRQ            = 0.05;
                diff_measure_all = [ ...
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
                    
                    %% control vs. FCD Type II
                    graph_feature     = PathLen_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_FCD_Type_II_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = PathLen_corr_HET_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    diff_measure(diff_measure==max(diff_measure)) = 0.9;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_HET_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = PathLen_corr_PMG_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_PMG_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                end
                
                %% clustering coefficient
                for clustering_coeff = 1
                    
                    %% control vs. FCD Type II
                    graph_feature     = ClustCoff_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_FCD_Type_II_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = ClustCoff_corr_HET_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_HET_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = ClustCoff_corr_PMG_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_PMG_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
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
        
        parpool(24);
        num_of_cont         = length(Codes_cont);
        num_of_FCD_Type_II  = length(Codes_FCD_Type_II);
        num_of_HET          = length(Codes_HET);
        num_of_PMG          = length(Codes_PMG);
        type_mat            = 1;
        
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
        
        iternation_num = 100;
        interval_thres = 0.01;
        myspar = 0.05:interval_thres:0.4;
        AAL_index_iter = max(microAAL_surf_data_both);
        mat_thres = 0.75;
        
        % 0) binarized averaged matrix construction
        for binarized_matrix = 1
            
            averaged_mat_cont        = cell(length(myspar), 1);
            averaged_mat_FCD_Type_II = cell(length(myspar), 1);
            averaged_mat_HET         = cell(length(myspar), 1);
            averaged_mat_PMG         = cell(length(myspar), 1);
            
            parfor i = 1 : length(myspar)
                
                i
                
                %% control
                mat_cont_temp = mat_cont*0;
                for j = 1 : num_of_cont
                    temp = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    mat_cont_temp(:, :, j) = temp > 0;
                    %                  mat_cont_temp(:, :, j) = temp;
                end
                averaged_mat_cont{i} = (sum(mat_cont_temp, 3)/num_of_cont)>mat_thres;
                %             averaged_mat_cont(:, :, i) = mean(mat_cont_temp, 3)>mat_thres;
                
                %% FCD Type II
                mat_FCD_Type_II_temp = mat_FCD_Type_II*0;
                for j = 1 : num_of_FCD_Type_II
                    temp = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    mat_FCD_Type_II_temp(:, :, j) = temp > 0;
                    %                 mat_FCD_Type_II_temp(:, :, j) = temp;
                end
                averaged_mat_FCD_Type_II{i} = (sum(mat_FCD_Type_II_temp, 3)/num_of_FCD_Type_II)>mat_thres;
                %             averaged_mat_FCD_Type_II(:, :, i) = mean(mat_FCD_Type_II_temp, 3)>mat_thres;
                
                %% HET
                mat_HET_temp = mat_HET*0;
                for j = 1 : num_of_HET
                    temp = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    mat_HET_temp(:, :, j) = temp > 0;
                    %                 mat_HET_temp(:, :, j) = temp;
                end
                averaged_mat_HET{i} = (sum(mat_HET_temp, 3)/num_of_HET)>mat_thres;
                %             averaged_mat_HET(:, :, i) = mean(mat_HET_temp, 3)>mat_thres;
                
                %% PMG
                mat_PMG_temp = mat_PMG*0;
                for j = 1 : num_of_PMG
                    temp = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    mat_PMG_temp(:, :, j) = temp > 0;
                    %                 mat_PMG_temp(:, :, j) = temp;
                end
                averaged_mat_PMG{i} = (sum(mat_PMG_temp, 3)/num_of_PMG)>mat_thres;
                %             averaged_mat_PMG(:, :, i) = mean(mat_PMG_temp, 3)>mat_thres;
                
            end
            
            averaged_mat_cont_temp = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
            for i = 1 : length(myspar)
                
                averaged_mat_cont_temp(:, :, i) = averaged_mat_cont{i};
                
            end
            
            averaged_mat_FCD_Type_II_temp = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
            for i = 1 : length(myspar)
                
                averaged_mat_FCD_Type_II_temp(:, :, i) = averaged_mat_FCD_Type_II{i};
                
            end
            
            averaged_mat_HET_temp = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
            for i = 1 : length(myspar)
                
                averaged_mat_HET_temp(:, :, i) = averaged_mat_HET{i};
                
            end
            
            averaged_mat_PMG_temp = zeros(AAL_index_iter, AAL_index_iter, length(myspar));
            for i = 1 : length(myspar)
                
                averaged_mat_PMG_temp(:, :, i) = averaged_mat_PMG{i};
                
            end
            
            averaged_mat_cont        = averaged_mat_cont_temp;
            averaged_mat_FCD_Type_II = averaged_mat_FCD_Type_II_temp;
            averaged_mat_HET         = averaged_mat_HET_temp;
            averaged_mat_PMG         = averaged_mat_PMG_temp;
            
        end

        % 1) nodal clustering coefficient:	unweighted positive network
        for CC = 1
          
            node_cc_uw_pos_cont         = zeros(length(myspar), AAL_index_iter, num_of_cont);            
            node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter, num_of_FCD_Type_II);
            node_cc_uw_pos_HET          = zeros(length(myspar), AAL_index_iter, num_of_HET);
            node_cc_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter, num_of_PMG);
            
            for i = 1 : length(myspar)
                
                %% control
                for j = 1 : num_of_cont
                    wmatrix                                  = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    node_cc_uw_pos_cont(i, :, j)             = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_cont(i, :, j)             = clustering_coef_wu(wmatrix)';
                end
                                
                %% FCD Type II
                for j = 1 : num_of_FCD_Type_II
                    wmatrix                                  = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    node_cc_uw_pos_FCD_Type_II(i, :, j)      = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_FCD_Type_II(i, :, j)      = clustering_coef_wu(wmatrix)';
                end
                
                %% HET
                for j = 1 : num_of_HET
                    wmatrix                                  = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    node_cc_uw_pos_HET(i, :, j)              = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_HET(i, :, j)              = clustering_coef_wu(wmatrix)';
                end
                
                %% PMG
                for j = 1 : num_of_PMG
                    wmatrix                                  = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    node_cc_uw_pos_PMG(i, :, j)              = clustering_coef_bu(double(wmatrix>0))';
%                     node_cc_uw_pos_PMG(i, :, j)              = clustering_coef_wu(wmatrix)';
                end
                
            end
            
            node_cc_uw_pos_cont_final = zeros(length(myspar), num_of_cont);
            node_cc_uw_pos_FCD_Type_II_final = zeros(length(myspar), num_of_FCD_Type_II);
            node_cc_uw_pos_HET_final = zeros(length(myspar), num_of_HET);
            node_cc_uw_pos_PMG_final = zeros(length(myspar), num_of_PMG);
            
            for i = 1 : length(myspar)
                
                for j = 1 : num_of_cont
                    temp = node_cc_uw_pos_cont(i, :, j);
                    node_cc_uw_pos_cont_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_FCD_Type_II
                    temp = node_cc_uw_pos_FCD_Type_II(i, :, j);
                    node_cc_uw_pos_FCD_Type_II_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_HET
                    temp = node_cc_uw_pos_HET(i, :, j); 
                    node_cc_uw_pos_HET_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_PMG
                    temp = node_cc_uw_pos_PMG(i, :, j);
                    node_cc_uw_pos_PMG_final(i, j) = mean(temp);
                end
                
            end
            
            [ mean(node_cc_uw_pos_cont_final, 2) mean(node_cc_uw_pos_FCD_Type_II_final, 2) mean(node_cc_uw_pos_HET_final, 2) mean(node_cc_uw_pos_PMG_final, 2) ]
            [ std(node_cc_uw_pos_cont_final, 0, 2) std(node_cc_uw_pos_FCD_Type_II_final, 0, 2) std(node_cc_uw_pos_HET_final, 0, 2) std(node_cc_uw_pos_PMG_final, 0, 2) ]
            figure; plot(1:length(myspar), [ mean(node_cc_uw_pos_cont_final, 2) mean(node_cc_uw_pos_FCD_Type_II_final, 2) mean(node_cc_uw_pos_HET_final, 2) mean(node_cc_uw_pos_PMG_final, 2) ])
            legend('control', 'FCD', 'HET', 'PMG');
            
        end
        
        % 2) nodal shortest path length:    unweighted positive network
        for Lp = 1
            
            % Characteristic path length using BCT
            node_spl_uw_pos_cont            = zeros(length(myspar), AAL_index_iter, num_of_cont);            
            node_spl_uw_pos_FCD_Type_II     = zeros(length(myspar), AAL_index_iter, num_of_FCD_Type_II);
            node_spl_uw_pos_HET             = zeros(length(myspar), AAL_index_iter, num_of_HET);
            node_spl_uw_pos_PMG             = zeros(length(myspar), AAL_index_iter, num_of_PMG);
                        
            for i = 1 : length(myspar)
                
                %% control
                for j = 1 : num_of_cont
                    wmatrix                                 = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter * 2;
%                     [averLp Lpi]                            = gretna_node_shortestpathlength_weight(wmatrix);  
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);
                    node_spl_uw_pos_cont(i, :, j)           = Lpi;
%                     dist = distance_bin(double(wmatrix>0)); 
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_cont(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_cont(i, :, j)          = mean(dist);
                end
                
                %% FCD Type II
                for j = 1 : num_of_FCD_Type_II
                    wmatrix                                 = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter * 2;
%                     [averLp Lpi]                            = gretna_node_shortestpathlength_weight(wmatrix); 
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);                    
                    node_spl_uw_pos_FCD_Type_II(i, :, j)    = Lpi;
%                     dist = distance_bin(double(wmatrix>0));
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_FCD_Type_II(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_FCD_Type_II(i, :, j)   = mean(dist);
                end
                
                %% HET
                for j = 1 : num_of_HET
                    wmatrix                         = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0); Lpi(isinf(Lpi)) = AAL_index_iter * 2; 
%                     [averLp Lpi]                    = gretna_node_shortestpathlength_weight(wmatrix); 
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);
                    node_spl_uw_pos_HET(i, :, j)    = Lpi;
%                     dist = distance_bin(double(wmatrix>0)); 
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_HET(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_HET(i, :, j)           = mean(dist);
                end
                
                %% PMG
                for j = 1 : num_of_PMG
                    wmatrix                         = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    [averLp Lpi]                            = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter * 2; 
%                     [averLp Lpi]                    = gretna_node_shortestpathlength_weight(wmatrix); 
%                     sw = smallworldness_bu(wmatrix>0, 0, 0);
                    node_spl_uw_pos_PMG(i, :, j)    = Lpi;
%                     dist = distance_bin(double(wmatrix>0)); 
%                     for k = 1 : AAL_index_iter
%                         temp = dist(k, :); node_spl_uw_pos_PMG(i, k, j) = mean(temp(~isinf(temp.*(temp~=0))));
%                     end
%                     node_spl_uw_pos_PMG(i, :, j)           = mean(dist);
                end
                
            end
            
            node_spl_uw_pos_cont_final            = zeros(length(myspar), num_of_cont);            
            node_spl_uw_pos_FCD_Type_II_final     = zeros(length(myspar), num_of_FCD_Type_II);
            node_spl_uw_pos_HET_final             = zeros(length(myspar), num_of_HET);
            node_spl_uw_pos_PMG_final             = zeros(length(myspar), num_of_PMG);
            
            for i = 1 : length(myspar)
                
                for j = 1 : num_of_cont
                    wmatrix = threshold_proportional(mat_cont(:, :, j), myspar(i));
                    temp    = node_spl_uw_pos_cont(i, :, j);
                    node_spl_uw_pos_cont_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_FCD_Type_II
                    wmatrix = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                    temp = node_spl_uw_pos_FCD_Type_II(i, :, j);
                    node_spl_uw_pos_FCD_Type_II_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_HET
                    wmatrix = threshold_proportional(mat_HET(:, :, j), myspar(i));
                    temp = node_spl_uw_pos_HET(i, :, j); 
                    node_spl_uw_pos_HET_final(i, j) = mean(temp);
                end
                for j = 1 : num_of_PMG
                    wmatrix = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                    temp = node_spl_uw_pos_PMG(i, :, j);
                    node_spl_uw_pos_PMG_final(i, j) = mean(temp);
                end
                
            end
            
            
            [ mean(node_spl_uw_pos_cont_final, 2) mean(node_spl_uw_pos_FCD_Type_II_final, 2) mean(node_spl_uw_pos_HET_final, 2) mean(node_spl_uw_pos_PMG_final, 2) ]
            [ std(node_spl_uw_pos_cont_final, 0, 2) std(node_spl_uw_pos_FCD_Type_II_final, 0, 2) std(node_spl_uw_pos_HET_final, 0, 2) std(node_spl_uw_pos_PMG_final, 0, 2) ]
            
            figure; plot(1:length(myspar), [ mean(node_spl_uw_pos_cont_final, 2) mean(node_spl_uw_pos_FCD_Type_II_final, 2) mean(node_spl_uw_pos_HET_final, 2) mean(node_spl_uw_pos_PMG_final, 2) ])
            legend('control', 'FCD', 'HET', 'PMG');
            
        end
        
        % 3) nodal clustering coefficient based on averaged matrix:	unweighted positive network
        for CC_avg_mat = 1
            
            node_cc_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_cc_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                wmatrix                            = threshold_proportional(mean_cont, myspar(i));
                node_cc_uw_pos_cont(i, :)          = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_cont(i, :)          = clustering_coef_wu(wmatrix)';
                
                wmatrix                            = threshold_proportional(mean_FCD_Type_II, myspar(i));
                node_cc_uw_pos_FCD_Type_II(i, :)   = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_FCD_Type_II(i, :)   = clustering_coef_wu(wmatrix)';
                
                wmatrix                            = threshold_proportional(mean_HET, myspar(i));
                node_cc_uw_pos_HET(i, :)           = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_HET(i, :)           = clustering_coef_wu(wmatrix)';
                
                wmatrix                            = threshold_proportional(mean_PMG, myspar(i));
                node_cc_uw_pos_PMG(i, :)           = clustering_coef_bu(double(wmatrix>0))';
%                 node_cc_uw_pos_PMG(i, :)           = clustering_coef_wu(wmatrix)';
            end   
            
            [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ]
            [ std(node_cc_uw_pos_cont, 0, 2) std(node_cc_uw_pos_FCD_Type_II, 0, 2) std(node_cc_uw_pos_HET, 0, 2) std(node_cc_uw_pos_PMG, 0, 2) ]
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ])
                legend('control', 'FCD', 'HET', 'PMG');
            end
            
        end
        
        % 4) nodal shortest path length based on averaged matrix:    unweighted positive network
        for Lp_avg_mat = 1
            
            node_spl_uw_pos_cont         = zeros(length(myspar), AAL_index_iter);          
            node_spl_uw_pos_FCD_Type_II  = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_HET          = zeros(length(myspar), AAL_index_iter);
            node_spl_uw_pos_PMG          = zeros(length(myspar), AAL_index_iter);
            
            for i = 1 : length(myspar)
                
                %% control
                wmatrix                            = threshold_proportional(mean_cont, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_cont(i, :)         = averLp;
                               
                %% FCD Type II
                wmatrix                            = threshold_proportional(mean_FCD_Type_II, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_FCD_Type_II(i, :)  = averLp;
                
                %% HET
                wmatrix                            = threshold_proportional(mean_HET, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_HET(i, :)          = averLp;
                
                %% PMG
                wmatrix                            = threshold_proportional(mean_PMG, myspar(i));
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_PMG(i, :)          = averLp;
                
            end
            
            [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]            
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]);
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
        
        % 5) nodal clustering coefficient based on thresholded averaged matrix:	unweighted positive network
        for CC_avg_mat_thresholded = 1
            
            node_cc_uw_pos_cont         = cell(length(myspar), 1);
            node_cc_uw_pos_FCD_Type_II  = cell(length(myspar), 1);
            node_cc_uw_pos_HET          = cell(length(myspar), 1);
            node_cc_uw_pos_PMG          = cell(length(myspar), 1);
            
            parfor i = 1 : length(myspar)
                
                i
                
                wmatrix                            = averaged_mat_cont(:, :, i);
                node_cc_uw_pos_cont{i}             = clustering_coef_bu(double(wmatrix>0))';
                
                wmatrix                            = averaged_mat_FCD_Type_II(:, :, i);
                node_cc_uw_pos_FCD_Type_II{i}      = clustering_coef_bu(double(wmatrix>0))';
                
                wmatrix                            = averaged_mat_HET(:, :, i);
                node_cc_uw_pos_HET{i}              = clustering_coef_bu(double(wmatrix>0))';
                
                wmatrix                            = averaged_mat_PMG(:, :, i);
                node_cc_uw_pos_PMG{i}               = clustering_coef_bu(double(wmatrix>0))';
                
            end   
            
            node_cc_uw_pos_cont = cell2mat(node_cc_uw_pos_cont);
            node_cc_uw_pos_FCD_Type_II = cell2mat(node_cc_uw_pos_FCD_Type_II);
            node_cc_uw_pos_HET = cell2mat(node_cc_uw_pos_HET);
            node_cc_uw_pos_PMG = cell2mat(node_cc_uw_pos_PMG);
            
            [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ]
            [ std(node_cc_uw_pos_cont, 0, 2) std(node_cc_uw_pos_FCD_Type_II, 0, 2) std(node_cc_uw_pos_HET, 0, 2) std(node_cc_uw_pos_PMG, 0, 2) ]
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_cc_uw_pos_cont, 2) mean(node_cc_uw_pos_FCD_Type_II, 2) mean(node_cc_uw_pos_HET, 2) mean(node_cc_uw_pos_PMG, 2) ])
                legend('control', 'FCD', 'HET', 'PMG');
            end
            
        end
        
        % 6) nodal shortest path length based on thresholded averaged matrix:    unweighted positive network
        for Lp_avg_mat_thresholded = 1
            
            node_spl_uw_pos_cont         = cell(length(myspar), 1);
            node_spl_uw_pos_FCD_Type_II  = cell(length(myspar), 1);
            node_spl_uw_pos_HET          = cell(length(myspar), 1);
            node_spl_uw_pos_PMG          = cell(length(myspar), 1);
            
            node_spl_uw_pos_cont_final         = cell(length(myspar), 1);
            node_spl_uw_pos_FCD_Type_II_final  = cell(length(myspar), 1);
            node_spl_uw_pos_HET_final          = cell(length(myspar), 1);
            node_spl_uw_pos_PMG_final          = cell(length(myspar), 1);
            
            parfor i = 1 : length(myspar)
                
                i
                
                %% control
                wmatrix                            = averaged_mat_cont(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_cont{i}            = Lpi;
                node_spl_uw_pos_cont_final{i}      = averLp;
                               
                %% FCD Type II
                wmatrix                            = averaged_mat_FCD_Type_II(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_FCD_Type_II{i}     = Lpi;
                node_spl_uw_pos_FCD_Type_II_final{i}  = averLp;
                
                %% HET
                wmatrix                            = averaged_mat_HET(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_HET{i}             = Lpi;
                node_spl_uw_pos_HET_final{i}       = averLp;
                
                %% PMG
                wmatrix                            = averaged_mat_PMG(:, :, i);
                [averLp Lpi]                       = gretna_node_shortestpathlength(wmatrix>0);  Lpi(isinf(Lpi)) = AAL_index_iter;
%                 [averLp Lpi]                       = gretna_node_shortestpathlength_weight(wmatrix);  Lpi(isinf(Lpi)) = AAL_index_iter;
                node_spl_uw_pos_PMG{i}             = Lpi;
                node_spl_uw_pos_PMG_final{i}       = averLp;
                
            end
            
            node_spl_uw_pos_cont         = cell2mat(node_spl_uw_pos_cont);
            node_spl_uw_pos_FCD_Type_II  = cell2mat(node_spl_uw_pos_FCD_Type_II);
            node_spl_uw_pos_HET          = cell2mat(node_spl_uw_pos_HET);
            node_spl_uw_pos_PMG          = cell2mat(node_spl_uw_pos_PMG);
            
            node_spl_uw_pos_cont_final         = cell2mat(node_spl_uw_pos_cont_final);
            node_spl_uw_pos_FCD_Type_II_final  = cell2mat(node_spl_uw_pos_FCD_Type_II_final);
            node_spl_uw_pos_HET_final          = cell2mat(node_spl_uw_pos_HET_final);
            node_spl_uw_pos_PMG_final          = cell2mat(node_spl_uw_pos_PMG_final);
            
            
            [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]
            [ node_spl_uw_pos_cont_final node_spl_uw_pos_FCD_Type_II_final node_spl_uw_pos_HET_final node_spl_uw_pos_PMG_final ]
            
            if(printfigs)
                figure; plot(1:length(myspar), [ mean(node_spl_uw_pos_cont, 2) mean(node_spl_uw_pos_FCD_Type_II, 2) mean(node_spl_uw_pos_HET, 2) mean(node_spl_uw_pos_PMG, 2) ]);
                legend('control', 'FCD Type-II', 'HET', 'PMG');
                
                figure; plot(1:length(myspar), [ node_spl_uw_pos_cont_final node_spl_uw_pos_FCD_Type_II_final node_spl_uw_pos_HET_final node_spl_uw_pos_PMG_final ]);
                legend('control', 'FCD Type-II', 'HET', 'PMG');
            end
            
        end
                
        % 7) small world parameters with random network comparison (complete running)
        for small_world = 1
            
            % density = 0.05 - 0.4 with 0.01 interval and 1000 randomization
            parpool(24);
            option.homogeneity_thres = [5, 0.06];
            option.bagging = 1;
            iternation_num = 100;
            interval_thres = 0.01;
            myspar = 0.05:interval_thres:0.4;
            AAL_index_iter = length(Anatomical_Label_structs_idx);

            %% control vs. FCD Type II with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_FCD_Type_II, ClustCoff_corr_FCD_Type_II, Sigma_corr_FCD_Type_II, rerandPathLen_FCD_Type_II, rerandClustCoff_FCD_Type_II ] = ...
                swdiff_sparsity_parfor_fMRI(mat_cont, mat_FCD_Type_II, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% control vs. Heterotopia with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_HET, ClustCoff_corr_HET, Sigma_corr_HET, rerandPathLen_HET, rerandClustCoff_HET ] = ...
                swdiff_sparsity_parfor_fMRI(mat_cont, mat_HET, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
            
            %% control vs. Polymicrogria with computing normalized parameters (i.e., lamda, gamma, sigma)
            [myspar, PathLen_corr_PMG, ClustCoff_corr_PMG, Sigma_corr_PMG, rerandPathLen_PMG, rerandClustCoff_PMG ] = ...
                swdiff_sparsity_parfor_fMRI(mat_cont, mat_PMG, myspar(1), myspar(end), interval_thres, iternation_num, 1, option);
           
            %% visualization
            if(visualization)
                
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                MarkerSize = 20;
                
                %% Path length
                for path_legnth = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([4 42]); 
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
                    ylim([1 17]);
                    
                    exportfigbo(gcf,[OUTPATH 'path_length_global_diff_cont_vs_pats_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 2) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -2.1;
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
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        end
                    end
                    ylim([-2.5 2.5]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_FCD_TypeII_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 3) delta, HET
                    colorline_idx = 4;
                    sig_pos = -2.1;
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
                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',12);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    ylim([-2.5 2.5]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_HET_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -2.5;
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

                    for i = 1: length(myspar)
                        if(feature1.Lp.Lp12_diff_p(i)<=FDRp)
                            text(myspar(i)*100, sig_pos, '*', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',20);
                        elseif(feature1.Lp.Lp12_diff_p(i)>FDRp & feature1.Lp.Lp12_diff_p(i)<=0.055)
                            text(myspar(i)*100, sig_pos, 'o', 'HorizontalAlignment','center', 'color',color_lines(colorline_idx, :),'FontSize',10);
                        end
                    end
                    ylim([-3 3]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_path_length_global_diff_cont_vs_PMG_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                end
                
                %% Clustering coefficient
                for clust_coeff = 1
                    
                    %% 1) All features, raw value
                    figure; hold on; xlim([4 42]); ylim([0.2 0.55]);
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
                    
                    exportfigbo(gcf,[OUTPATH 'clustering_coeff_global_diff_cont_vs_pats_fMRI_micro.png'], 'png', 13); close(gcf);
                           
                    %% 3) delta, FCD Type II
                    colorline_idx = 3;
                    sig_pos = -0.04;
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
                    ylim([-0.05 0.05]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_FCD_TypeII_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 4) delta, HET
                    colorline_idx = 4;
                    sig_pos = -0.09;
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
                    ylim([-0.1 0.06]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_HET_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                    %% 5) delta, PMG
                    colorline_idx = 5;
                    sig_pos = -0.14;
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
                    ylim([-0.15 0.1]);
                    
                    exportfigbo(gcf,[OUTPATH 'delta_clustering_coeff_global_diff_cont_vs_PMG_fMRI_micro.png'], 'png', 13); close(gcf);
                    
                end
                
            end
            
        end
        
        % 8) Lp and Cp at nodal level
        for Lp_Cp_nodal = 1
            
            parpool(24);
             
            %% control vs. FCD Type II
            [ myspar, PathLen_corr_FCD_Type_II_nodal, ClustCoff_corr_FCD_Type_II_nodal, rerandPathLen_FCD_Type_II_nodal, rerandClustCoff_FCD_Type_II_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor_fMRI(mat_cont, mat_FCD_Type_II, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Heterotopia
            [myspar, PathLen_corr_HET_nodal, ClustCoff_corr_HET_nodal, rerandPathLen_HET_nodal, rerandClustCoff_HET_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor_fMRI(mat_cont, mat_HET, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% control vs. Polymicrogria
            [myspar, PathLen_corr_PMG_nodal, ClustCoff_corr_PMG_nodal, rerandPathLen_PMG_nodal, rerandClustCoff_PMG_nodal ] = ...
                LpCpDiff_nodal_level_sparsity_parfor_fMRI(mat_cont, mat_PMG, myspar(1), myspar(end), interval_thres, iternation_num);
            
            %% visualization
            if(visualization)
                
                %% variable setup
                OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/01_sw_analysis/';
                
                spar            = 7; % min(find(node_fully_connected));
                uncorrp         = 0.05;
                FDRQ            = 0.05;
                diff_measure_all = [ ...
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
                    
                    %% control vs. FCD Type II
                    graph_feature     = PathLen_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_FCD_Type_II_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = PathLen_corr_HET_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    diff_measure(diff_measure==max(diff_measure)) = 0.9;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_HET_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = PathLen_corr_PMG_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Lp.Lp2(:, spar)' - graph_feature.Lp.Lp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Lp.Lp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'path_length_nodal_diff_cont_vs_PMG_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                end
                
                %% clustering coefficient
                for clustering_coeff = 1
                    
                    %% control vs. FCD Type II
                    graph_feature     = ClustCoff_corr_FCD_Type_II_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_FCD_Type_II_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. HET
                    graph_feature     = ClustCoff_corr_HET_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_HET_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
                    %% control vs. PMG
                    graph_feature     = ClustCoff_corr_PMG_nodal;
                    diff_real_wmatrix = zeros(78, 78);
                    diff_measure      = (graph_feature.Cp.Cp2(:, spar)' - graph_feature.Cp.Cp1(:, spar)');
                    
                    colormap_curr     = colormap_generate_nodal_measure(graph_feature, spar, uncorrp, FDRQ);
                    diff_measure(graph_feature.Cp.Cp12_diff_p(:, spar)>0.05) = 0;
                    
                    figure;
                    [aa1, cb1]=NOELNetVisuNodalMeasuresOnSurf(diff_measure, diff_real_wmatrix, nodePos, (1:78), S, 'white', 'link_width', 0, 'Node_colormap', colormap_curr, 'node_color_by', 'nodalMeasure', 'node_radius', 5);
                    exportfigbo(gcf,[OUTPATH 'cluscoef_nodal_diff_cont_vs_PMG_3T_rsfMRI.png'], 'png', 13); close(gcf);
                    
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