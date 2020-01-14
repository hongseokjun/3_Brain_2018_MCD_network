clear all
close all

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/'));
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
%     fid = fopen('DemographicNGrouping_for_controls_Data_3T_fMRI_2.csv');
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
%     fid = fopen('DemographicNGrouping_for_FCD_Type-II_Data_3T_fMRI_2.csv');
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
%     fid = fopen('DemographicNGrouping_for_HET_Data_3T_fMRI_2.csv');
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
%     fid = fopen('DemographicNGrouping_for_PMG_Data_3T_fMRI_2.csv');
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

%% statistical test for age, duration and sex distribution
for test_age_gender = 1
    
    num_of_cont         = length(Codes_cont);
    num_of_FCD_Type_II  = length(Codes_FCD_Type_II);
    num_of_HET          = length(Codes_HET);
    num_of_PMG          = length(Codes_PMG);
    
    %              N          Age         Sex         L/R/B           Duration
    % Control      41      29.2/7.1      18/15         N/A               N/A            
    % FCD Type-II  30      27.3/8.6      14/13       15/15/0          20.9/11.9   
    % HET          27      31.4/10.6      8/6        6/10/11          Not yet   
    % PMG          21      29.5/11.4      6/5        6/5/10           Not yet   
    
    % Confirmed: No age difference between controls and patients; between FCD Type-I and FCD Type-II
    for Age = 1
                
        % control vs. FCD Type-II: p=0.3632
        [h,p,ci,stats] = ttest2(Age_cont, Age_FCD_Type_II, 0.05, 'both')
        
        [ mean(Age_cont) std(Age_cont);
          mean(Age_FCD_Type_II) std(Age_FCD_Type_II); ]
        
        % control vs. HET: p=0.3904
        [h,p,ci,stats] = ttest2(Age_cont, Age_HET, 0.05, 'both')
        
        [ mean(Age_cont) std(Age_cont);
          mean(Age_HET) std(Age_HET); ]
        
        % control vs. PMG: p=0.9171
        [h,p,ci,stats] = ttest2(Age_cont, Age_PMG, 0.05, 'both')
        
        [ mean(Age_cont) std(Age_cont);
          mean(Age_PMG) std(Age_PMG); ]
      
        % FCD Type-II vs. PMG: p=0.5291
        [h,p,ci,stats] = ttest2(Age_FCD_Type_II, Age_PMG, 0.05, 'both')
        
        [ mean(Age_FCD_Type_II) std(Age_FCD_Type_II);
          mean(Age_PMG) std(Age_PMG); ]
      
        % HET vs. PMG: p=0.6591
        [h,p,ci,stats] = ttest2(Age_HET, Age_PMG, 0.05, 'both')
        
        [ mean(Age_HET) std(Age_HET);
          mean(Age_PMG) std(Age_PMG); ]
      
        % HET vs. FCD Type-II: p=0.19
        [h,p,ci,stats] = ttest2(Age_HET, Age_FCD_Type_II, 0.05, 'both')
        
        [ mean(Age_HET) std(Age_HET);
          mean(Age_FCD_Type_II) std(Age_FCD_Type_II); ]
      
    end
    
    % Confirmed: No sex difference between controls and patients; between FCD Type-I and FCD Type-II
    for Sex = 1
                
        % control vs. FCD Type-II: p=0.52
        [p, x2] = chisquarecont( [ sum(strcmp(Sex_cont, 'male')) sum(strcmp(Sex_cont, 'female')); ...
                                   sum(strcmp(Sex_FCD_Type_II, 'male')) sum(strcmp(Sex_FCD_Type_II, 'female')) ])
        
        [ sum(strcmp(Sex_cont, 'male')) sum(strcmp(Sex_cont, 'female'));
          sum(strcmp(Sex_FCD_Type_II, 'male')) sum(strcmp(Sex_FCD_Type_II, 'female')) ]

        % control vs. HET: p=0.89
        [p, x2] = chisquarecont( [ sum(strcmp(Sex_cont, 'male')) sum(strcmp(Sex_cont, 'female')); ...
                                   sum(strcmp(Sex_HET, 'male')) sum(strcmp(Sex_HET, 'female')) ])
        
        [ sum(strcmp(Sex_cont, 'male')) sum(strcmp(Sex_cont, 'female'));
          sum(strcmp(Sex_HET, 'male')) sum(strcmp(Sex_HET, 'female')) ]
      
        % control vs. PMG: 0.77 
        [p, x2] = chisquarecont( [ sum(strcmp(Sex_cont, 'male')) sum(strcmp(Sex_cont, 'female')); ...
            sum(strcmp(Sex_PMG, 'male')) sum(strcmp(Sex_PMG, 'female')) ])
        
        [ sum(strcmp(Sex_cont, 'male')) sum(strcmp(Sex_cont, 'female'));
          sum(strcmp(Sex_PMG, 'male')) sum(strcmp(Sex_PMG, 'female')) ]
      
    end    
    
    % Confrimed: No duration effect
    for Duration = 1
        
        % MRI negative vs. MRI positive FCD epilepsy: p=0.32, t=-1.0
%         [h,p,ci,stats] = ttest2(Duration_FCD_Type_I, Duration_FCD_Type_II, 0.05, 'both')
%         
%         [ mean(Duration_FCD_Type_I) std(Duration_FCD_Type_I);
%             mean(Duration_FCD_Type_II) std(Duration_FCD_Type_II); ]
%         
%         % MRI negative vs. HET: p=0.02, t=-2.32
%         [h,p,ci,stats] = ttest2(Duration_FCD_Type_I, Duration_HET, 0.05, 'both')
%         
%         [ mean(Duration_FCD_Type_I) std(Duration_FCD_Type_I);
%             mean(Duration_HET) std(Duration_HET); ]
%         
%         % MRI negative vs. PMG: p=0.02, t=-2.34
%         [h,p,ci,stats] = ttest2(Duration_FCD_Type_I, Duration_PMG, 0.05, 'both')
%         
%         [ mean(Duration_FCD_Type_I) std(Duration_FCD_Type_I);
%             mean(Duration_PMG) std(Duration_PMG); ]
%         
%         % MRI positive vs. HET
%         [h,p,ci,stats] = ttest2(Duration_FCD_Type_II, Duration_HET, 0.05, 'both')
%         
%         [ mean(Duration_FCD_Type_II) std(Duration_FCD_Type_II);
%             mean(Duration_HET) std(Duration_HET); ]
%         
%         % MRI positive vs. PMG
%         [h,p,ci,stats] = ttest2(Duration_FCD_Type_II, Duration_PMG, 0.05, 'both')
%         
%         [ mean(Duration_FCD_Type_II) std(Duration_FCD_Type_II);
%             mean(Duration_PMG) std(Duration_PMG); ]        
%         
%         % HET vs. PMG
%         [h,p,ci,stats] = ttest2(Duration_HET, Duration_PMG, 0.05, 'both')
%         
%         [ mean(Duration_FCD_Type_II) std(Duration_FCD_Type_II);
%             mean(Duration_HET) std(Duration_HET); ]                
    end
    
end

%% load lesion label in FCD Type-II and Heterotopia
for load_lesion_label = 1
    
    for FCD_Type_II_label = 1
        
%         names_lesion_left_mcd = strcat(LabelPath_FCD, 'mcd_', Codes_FCD_Type_II, '_label_union_left_20mm_rsl.txt');
%         names_lesion_right_mcd = strcat(LabelPath_FCD, 'mcd_', Codes_FCD_Type_II, '_label_union_right_20mm_rsl.txt');
%         
%         threshold_lesion = 0.5;
%         Lesion_FCD_whole_left = SurfStatReadData(names_lesion_left_mcd);
%         Lesion_FCD_whole_right = SurfStatReadData(names_lesion_right_mcd);
%         Lesion_FCD_whole = SurfStatReadData([names_lesion_left_mcd, names_lesion_right_mcd]);
%         
%         Lesion_FCD_whole_left_temp  = zeros(size(Lesion_FCD_whole_left));
%         Lesion_FCD_whole_right_temp = zeros(size(Lesion_FCD_whole_right));
%         Lesion_FCD_whole_temp       = zeros(size(Lesion_FCD_whole));
%         
%         Lesion_FCD_whole_left_temp(Lesion_FCD_whole_left > threshold_lesion) = 1;
%         Lesion_FCD_whole_right_temp(Lesion_FCD_whole_right > threshold_lesion) = 1;
%         Lesion_FCD_whole_temp(Lesion_FCD_whole > threshold_lesion) = 1;
%         
%         Lesion_FCD_whole = Lesion_FCD_whole_temp;
%         
%         for i = 1: size(Lesion_FCD_whole, 1)
%             if(Left_FCD_Type_II(i))
%                 Lesion_FCD_whole(i, 40963:end) = 0;
%             elseif(Right_FCD_Type_II(i))
%                 Lesion_FCD_whole(i, 1:40962) = 0;
%             end
%         end
        
    end
    
    for HET_label = 1
        
%         names_lesion_left_mcd = strcat(LabelPath_HET, 'mcd_', Codes_HET, '_label_union_left_20mm_rsl.txt');
%         names_lesion_right_mcd = strcat(LabelPath_HET, 'mcd_', Codes_HET, '_label_union_right_20mm_rsl.txt');
%         
%         threshold_lesion = 0.5;
%         Lesion_HET_whole_left = SurfStatReadData(names_lesion_left_mcd);
%         Lesion_HET_whole_right = SurfStatReadData(names_lesion_right_mcd);
%         Lesion_HET_whole = SurfStatReadData([names_lesion_left_mcd, names_lesion_right_mcd]);
%         
%         Lesion_HET_whole_left_temp  = zeros(size(Lesion_HET_whole_left));
%         Lesion_HET_whole_right_temp = zeros(size(Lesion_HET_whole_right));
%         Lesion_HET_whole_temp       = zeros(size(Lesion_HET_whole));
%         
%         Lesion_HET_whole_left_temp(Lesion_HET_whole_left > threshold_lesion) = 1;
%         Lesion_HET_whole_right_temp(Lesion_HET_whole_right > threshold_lesion) = 1;
%         Lesion_HET_whole_temp(Lesion_HET_whole > threshold_lesion) = 1;
%         
%         Lesion_HET_whole = Lesion_HET_whole_temp;
%         
%         for i = 1: size(Lesion_HET_whole, 1)
%             if(Left_HET_Type_II(i))
%                 Lesion_HET_whole(i, 40963:end) = 0;
%             elseif(Right_HET_Type_II(i))
%                 Lesion_HET_whole(i, 1:40962) = 0;
%             end
%         end
        
    end

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

%% read time seires
for read_TS = 1
    
    for read_time_series_by_parcel = 1
        
        rois       = AAL_surf_data_both;
        rois_micro = microAAL_surf_data_both;
        
        % control
        for control = 1
            
            PREFIX = 'TLE';
            ts_parcel_cont = zeros(145, length(Anatomical_Label_structs_idx), num_of_cont);
            for i = 1:length(Codes_cont)
                
                Codes_cont{i}
                temp_path = Path_cont{i};
                temp_code = Codes_cont{i};
                if(strcmp(temp_code, '322_1'))
                    temp_code = '322_2';
                end
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), length(Anatomical_Label_structs_idx));
                
                for roi = 1:length(Anatomical_Label_structs_idx)
                    ts_parcel(:,roi) = mean(ts(:,rois==Anatomical_Label_structs_idx(roi)),2);
                end
                ts_parcel_cont(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            ts_parcel_cont_micro = zeros(145, max(microAAL_surf_data_both), num_of_cont);
            for i = 1:length(Codes_cont)
                
                Codes_cont{i}
                temp_path = Path_cont{i};
                temp_code = Codes_cont{i};
                if(strcmp(temp_code, '322_1'))
                    temp_code = '322_2';
                end
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), max(microAAL_surf_data_both));
                
                for roi = 1:max(microAAL_surf_data_both)
                    ts_parcel(:,roi) = mean(ts(:,rois_micro==roi),2);
                end
                ts_parcel_cont_micro(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            %         ts_parcel_cont_res = zeros(size(ts_parcel_cont));
            %         for i = 1:length(Anatomical_Label_structs_idx)
            %             for j = 1:145
            %
            %                 Y = squeeze(ts_parcel_cont(j, i, :));
            %                 gms = squeeze(mean(mean(ts_parcel_cont, 2)));
            %                 AGE = term(Age_cont);
            %                 SEX = term(Sex_cont);
            %                 GMS = term(gms);
            %                 MODEL = 1 + AGE + SEX + GMS;
            %                 slm = SurfStatLinMod(Y, MODEL);
            %                 ts_parcel_cont_res(j, i, :) = Y -  slm.X*slm.coef;
            %
            %             end
            %         end
            
        end
        
        % FCD_Type_II
        for FCD_Type_II = 1
            
            PREFIX = 'mcd';
            ts_parcel_FCD_Type_II = zeros(145, length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
            for i = 1:length(Codes_FCD_Type_II)
                
                Codes_FCD_Type_II{i}
                temp_path = Path_FCD_Type_II{i};
                temp_code = Codes_FCD_Type_II{i}; temp_code = temp_code(1:3);
                if(strcmp(temp_code, '080'))
                    temp_code = '080_2';
                end
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), length(Anatomical_Label_structs_idx));
                
                for roi = 1:length(Anatomical_Label_structs_idx)
                    ts_parcel(:,roi) = mean(ts(:,rois==Anatomical_Label_structs_idx(roi)),2);
                end
                ts_parcel_FCD_Type_II(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            ts_parcel_FCD_Type_II_micro = zeros(145, max(microAAL_surf_data_both), num_of_FCD_Type_II);
            for i = 1:length(Codes_FCD_Type_II)
                
                Codes_FCD_Type_II{i}
                temp_path = Path_FCD_Type_II{i};
                temp_code = Codes_FCD_Type_II{i}; temp_code = temp_code(1:3);
                if(strcmp(temp_code, '080'))
                    temp_code = '080_2';
                end
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), max(microAAL_surf_data_both));
                
                for roi = 1:max(microAAL_surf_data_both)
                    ts_parcel(:,roi) = mean(ts(:,rois_micro==roi),2);
                end
                ts_parcel_FCD_Type_II_micro(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            %         ts_parcel_FCD_Type_II_res = zeros(size(ts_parcel_FCD_Type_II));
            %         for i = 1:length(Anatomical_Label_structs_idx)
            %             for j = 1:145
            %
            %                 Y = squeeze(ts_parcel_FCD_Type_II(j, i, :));
            %                 gms = squeeze(mean(mean(ts_parcel_FCD_Type_II, 2)));
            %                 AGE = term(Age_FCD_Type_II);
            %                 SEX = term(Sex_FCD_Type_II);
            %                 GMS = term(gms);
            %                 MODEL = 1 + AGE + SEX + GMS;
            %                 slm = SurfStatLinMod(Y, MODEL);
            %                 ts_parcel_FCD_Type_II_res(j, i, :) = Y -  slm.X*slm.coef;
            %
            %             end
            %         end
            
        end
        
        % HET
        for HET = 1
            
            PREFIX = 'mcd';
            ts_parcel_HET = zeros(145, length(Anatomical_Label_structs_idx), num_of_HET);
            for i = 1:length(Codes_HET)
                
                Codes_HET{i}
                temp_path = Path_HET{i};
                temp_code = Codes_HET{i};
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), length(Anatomical_Label_structs_idx));
                
                for roi = 1:length(Anatomical_Label_structs_idx)
                    ts_parcel(:,roi) = mean(ts(:,rois==Anatomical_Label_structs_idx(roi)),2);
                end
                ts_parcel_HET(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            ts_parcel_HET_micro = zeros(145, max(microAAL_surf_data_both), num_of_HET);
            for i = 1:length(Codes_HET)
                
                Codes_HET{i}
                temp_path = Path_HET{i};
                temp_code = Codes_HET{i};
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), max(microAAL_surf_data_both));
                
                for roi = 1:max(microAAL_surf_data_both)
                    ts_parcel(:,roi) = mean(ts(:,rois_micro==roi),2);
                end
                ts_parcel_HET_micro(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            %         ts_parcel_HET_res = zeros(size(ts_parcel_HET));
            %         for i = 1:length(Anatomical_Label_structs_idx)
            %             for j = 1:145
            %
            %                 Y = squeeze(ts_parcel_HET(j, i, :));
            %                 gms = squeeze(mean(mean(ts_parcel_HET, 2)));
            %                 AGE = term(Age_HET);
            %                 SEX = term(Sex_HET);
            %                 GMS = term(gms);
            %                 MODEL = 1 + AGE + SEX + GMS;
            %                 slm = SurfStatLinMod(Y, MODEL);
            %                 ts_parcel_HET_res(j, i, :) = Y -  slm.X*slm.coef;
            %
            %             end
            %         end
            
        end
        
        % PMG
        for PMG = 1
            
            PREFIX = 'mcd';
            ts_parcel_PMG = zeros(145, length(Anatomical_Label_structs_idx), num_of_PMG);
            for i = 1:length(Codes_PMG)
                
                Codes_PMG{i}
                temp_path = Path_PMG{i};
                temp_code = Codes_PMG{i};
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), length(Anatomical_Label_structs_idx));
                
                for roi = 1:length(Anatomical_Label_structs_idx)
                    ts_parcel(:,roi) = mean(ts(:,rois==Anatomical_Label_structs_idx(roi)),2);
                end
                ts_parcel_PMG(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            ts_parcel_PMG_micro = zeros(145, max(microAAL_surf_data_both), num_of_PMG);
            for i = 1:length(Codes_PMG)
                
                Codes_PMG{i}
                temp_path = Path_PMG{i};
                temp_code = Codes_PMG{i};
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), max(microAAL_surf_data_both));
                
                for roi = 1:max(microAAL_surf_data_both)
                    ts_parcel(:,roi) = mean(ts(:,rois_micro==roi),2);
                end
                ts_parcel_PMG_micro(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
            end
            
            %         ts_parcel_PMG_res = zeros(size(ts_parcel_PMG));
            %         for i = 1:length(Anatomical_Label_structs_idx)
            %             for j = 1:145
            %
            %                 Y = squeeze(ts_parcel_PMG(j, i, :));
            %                 gms = squeeze(mean(mean(ts_parcel_PMG, 2)));
            %                 AGE = term(Age_PMG);
            %                 SEX = term(Sex_PMG);
            %                 GMS = term(gms);
            %                 MODEL = 1 + AGE + SEX + GMS;
            %                 slm = SurfStatLinMod(Y, MODEL);
            %                 ts_parcel_PMG_res(j, i, :) = Y -  slm.X*slm.coef;
            %
            %             end
            %         end
            
        end
        
    end
    
    %% origianl AAL
    for original_AAL = 1
        
        ts_parcel_cont_res          = zeros(size(ts_parcel_cont));
        ts_parcel_FCD_Type_II_res   = zeros(size(ts_parcel_FCD_Type_II));
        ts_parcel_HET_res           = zeros(size(ts_parcel_HET));
        ts_parcel_PMG_res           = zeros(size(ts_parcel_PMG));
        
        ts_parcel_cont_GMS          = squeeze(mean(ts_parcel_cont, 2));
        ts_parcel_FCD_Type_II_GMS   = squeeze(mean(ts_parcel_FCD_Type_II, 2));
        ts_parcel_HET_GMS           = squeeze(mean(ts_parcel_HET, 2));
        ts_parcel_PMG_GMS           = squeeze(mean(ts_parcel_PMG, 2));
        
        gms = [ ts_parcel_cont_GMS ts_parcel_FCD_Type_II_GMS ts_parcel_HET_GMS ts_parcel_PMG_GMS ]';
        age = [ Age_cont; Age_FCD_Type_II; Age_HET; Age_PMG ];
        sex = [ Sex_cont; Sex_FCD_Type_II; Sex_HET; Sex_PMG ];
        AGE = term(age);
        SEX = term(sex);
        
        temp_idx = [ 1 num_of_cont; num_of_cont+1 num_of_cont+num_of_FCD_Type_II; num_of_cont+num_of_FCD_Type_II+1 num_of_cont+num_of_FCD_Type_II+num_of_HET; num_of_cont+num_of_FCD_Type_II+num_of_HET+1 num_of_cont+num_of_FCD_Type_II+num_of_HET+num_of_PMG ]
        
        for i = 1 : 145
            for j = 1 : length(Anatomical_Label_structs_idx)
                
                gms_temp = gms(:, i);
                GMS = term(gms_temp);
                Y = [ squeeze(ts_parcel_cont(i, j, :)); squeeze(ts_parcel_FCD_Type_II(i, j, :)); squeeze(ts_parcel_HET(i, j, :)); squeeze(ts_parcel_PMG(i, j, :)) ];
%                 MODEL = 1 + AGE + SEX + GMS;
                MODEL = 1 + AGE + SEX;
                slm = SurfStatLinMod(Y, MODEL);
                temp =  Y -  slm.X*slm.coef;
                ts_parcel_cont_res(i, j, :)         = temp(temp_idx(1, 1):temp_idx(1, 2));
                ts_parcel_FCD_Type_II_res(i, j, :)  = temp(temp_idx(2, 1):temp_idx(2, 2));
                ts_parcel_HET_res(i, j, :)          = temp(temp_idx(3, 1):temp_idx(3, 2));
                ts_parcel_PMG_res(i, j, :)          = temp(temp_idx(4, 1):temp_idx(4, 2));
                
            end
        end
        
    end
    
    %% micro AAL
    for micro_AAL = 1
        
        parpool(24);
        
        ts_parcel_cont_micro_res          = cell(size(ts_parcel_cont_micro, 1), size(ts_parcel_cont_micro, 2));
        ts_parcel_FCD_Type_II_micro_res   = cell(size(ts_parcel_FCD_Type_II_micro, 1), size(ts_parcel_FCD_Type_II_micro, 2));
        ts_parcel_HET_micro_res           = cell(size(ts_parcel_HET_micro, 1), size(ts_parcel_HET_micro, 2));
        ts_parcel_PMG_micro_res           = cell(size(ts_parcel_PMG_micro, 1), size(ts_parcel_PMG_micro, 2));
        
        ts_parcel_cont_micro_res_temp          = zeros(size(ts_parcel_cont_micro));
        ts_parcel_FCD_Type_II_micro_res_temp   = zeros(size(ts_parcel_FCD_Type_II_micro));
        ts_parcel_HET_micro_res_temp           = zeros(size(ts_parcel_HET_micro));
        ts_parcel_PMG_micro_res_temp           = zeros(size(ts_parcel_PMG_micro));        
        
        ts_parcel_cont_micro_GMS          = squeeze(mean(ts_parcel_cont_micro, 2));
        ts_parcel_FCD_Type_II_micro_GMS   = squeeze(mean(ts_parcel_FCD_Type_II_micro, 2));
        ts_parcel_HET_micro_GMS           = squeeze(mean(ts_parcel_HET_micro, 2));
        ts_parcel_PMG_micro_GMS           = squeeze(mean(ts_parcel_PMG_micro, 2));
        
        gms = [ ts_parcel_cont_micro_GMS ts_parcel_FCD_Type_II_micro_GMS ts_parcel_HET_micro_GMS ts_parcel_PMG_micro_GMS ]';
        age = [ Age_cont; Age_FCD_Type_II; Age_HET; Age_PMG ];
        sex = [ Sex_cont; Sex_FCD_Type_II; Sex_HET; Sex_PMG ];
        AGE = term(age);
        SEX = term(sex);
        
        temp_idx = [ 1 num_of_cont; num_of_cont+1 num_of_cont+num_of_FCD_Type_II; num_of_cont+num_of_FCD_Type_II+1 num_of_cont+num_of_FCD_Type_II+num_of_HET; num_of_cont+num_of_FCD_Type_II+num_of_HET+1 num_of_cont+num_of_FCD_Type_II+num_of_HET+num_of_PMG ]
        
        for i = 1 : 145
            i
            parfor j = 1 : max(microAAL_surf_data_both)
                
                gms_temp = gms(:, i);
                GMS = term(gms_temp);
                Y = [ squeeze(ts_parcel_cont_micro(i, j, :)); squeeze(ts_parcel_FCD_Type_II_micro(i, j, :)); squeeze(ts_parcel_HET_micro(i, j, :)); squeeze(ts_parcel_PMG_micro(i, j, :)) ];
%                 MODEL = 1 + AGE + SEX + GMS;
                MODEL = 1 + AGE + SEX;
                slm = SurfStatLinMod(Y, MODEL);
                temp =  Y -  slm.X*slm.coef;                
                ts_parcel_cont_micro_res{i, j}         = temp(temp_idx(1, 1):temp_idx(1, 2));
                ts_parcel_FCD_Type_II_micro_res{i, j}  = temp(temp_idx(2, 1):temp_idx(2, 2));
                ts_parcel_HET_micro_res{i, j}          = temp(temp_idx(3, 1):temp_idx(3, 2));
                ts_parcel_PMG_micro_res{i, j}          = temp(temp_idx(4, 1):temp_idx(4, 2));
                
            end
                        
            for j = 1 : max(microAAL_surf_data_both)
                ts_parcel_cont_micro_res_temp(i, j, :) = ts_parcel_cont_micro_res{i, j}';
                ts_parcel_FCD_Type_II_micro_res_temp(i, j, :) = ts_parcel_FCD_Type_II_micro_res{i, j}';
                ts_parcel_HET_micro_res_temp(i, j, :) = ts_parcel_HET_micro_res{i, j}';
                ts_parcel_PMG_micro_res_temp(i, j, :) = ts_parcel_PMG_micro_res{i, j}';
            end
            
        end
        
        ts_parcel_cont_micro_res        = ts_parcel_cont_micro_res_temp;
        ts_parcel_FCD_Type_II_micro_res = ts_parcel_FCD_Type_II_micro_res_temp;
        ts_parcel_HET_micro_res         = ts_parcel_HET_micro_res_temp;
        ts_parcel_PMG_micro_res         = ts_parcel_PMG_micro_res_temp;
        
    end
    
end

%% construct association matrices
for construct_asso_mat = 1
    
    load([ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_GSR_Age_Sex1_micro.mat' ]);
    FDR_thres = 0.05;
    
    %% construct association matrices
    for original_AAL = 1
        
        adjacent_matrix_control1              = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_cont);
        adjacent_matrix_FCD_Type_II1          = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
        adjacent_matrix_HET1                  = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_HET);
        adjacent_matrix_PMG1                  = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_PMG);
        zadjacent_matrix_control1             = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_cont);
        zadjacent_matrix_FCD_Type_II1         = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
        zadjacent_matrix_HET1                 = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_HET);
        zadjacent_matrix_PMG1                 = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_PMG);
        
        adjacent_matrix_control2              = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_cont);
        adjacent_matrix_FCD_Type_II2          = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
        adjacent_matrix_HET2                  = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_HET);
        adjacent_matrix_PMG2                  = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_PMG);
        zadjacent_matrix_control2             = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_cont);
        zadjacent_matrix_FCD_Type_II2         = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
        zadjacent_matrix_HET2                 = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_HET);
        zadjacent_matrix_PMG2                 = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_PMG);
        
        adjacent_matrix_control3              = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_cont);
        adjacent_matrix_FCD_Type_II3          = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
        adjacent_matrix_HET3                  = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_HET);
        adjacent_matrix_PMG3                  = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_PMG);
        zadjacent_matrix_control3             = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_cont);
        zadjacent_matrix_FCD_Type_II3         = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
        zadjacent_matrix_HET3                 = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_HET);
        zadjacent_matrix_PMG3                 = zeros(length(Anatomical_Label_structs_idx), length(Anatomical_Label_structs_idx), num_of_PMG);
        
        FDR_thres_final_cont            = [];
        FDR_thres_final_FCD_Type_II     = [];
        FDR_thres_final_HET             = [];
        FDR_thres_final_PMG             = [];
        
        min_thres_r_cont            = [];
        min_thres_r_FCD_Type_II     = [];
        min_thres_r_HET             = [];
        min_thres_r_PMG             = [];
        
        min_thres_z_cont            = [];
        min_thres_z_FCD_Type_II     = [];
        min_thres_z_HET             = [];
        min_thres_z_PMG             = [];
        
        % control
        for i = 1:num_of_cont
            
            i
            ts_temp = ts_parcel_cont_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_cont(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_cont(i) strcmp(Sex_cont(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(length(Anatomical_Label_structs_idx))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_cont = [ FDR_thres_final_cont FDR_thres_final ];
            adjacent_matrix_control1(:, :, i)   = r;
            zadjacent_matrix_control1(:, :, i)  = z;
            adjacent_matrix_control2(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_control2(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_control3(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_control3(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_cont = [ min_thres_r_cont min(temp_r(temp_r~=0)) ];
            min_thres_z_cont = [ min_thres_z_cont min(temp_z(temp_z~=0)) ];
            
        end
        
        % FCD Type-II
        for i = 1:num_of_FCD_Type_II
            
            i
            ts_temp = ts_parcel_FCD_Type_II_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_FCD_Type_II(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_FCD_Type_II(i) strcmp(Sex_FCD_Type_II(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(length(Anatomical_Label_structs_idx))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_FCD_Type_II = [ FDR_thres_final_FCD_Type_II FDR_thres_final ];
            adjacent_matrix_FCD_Type_II1(:, :, i)   = r;
            zadjacent_matrix_FCD_Type_II1(:, :, i)  = z;
            adjacent_matrix_FCD_Type_II2(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_FCD_Type_II2(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_FCD_Type_II3(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_FCD_Type_II3(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_FCD_Type_II = [ min_thres_r_FCD_Type_II min(temp_r(temp_r~=0)) ];
            min_thres_z_FCD_Type_II = [ min_thres_z_FCD_Type_II min(temp_z(temp_z~=0)) ];
            
        end
        
        % HET
        for i = 1:num_of_HET
            
            i
            ts_temp = ts_parcel_HET_res(:, :, i);
            if(i==6)
                ts_temp = ts_temp(1:110, :);
            end
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_HET(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_HET(i) strcmp(Sex_HET(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(length(Anatomical_Label_structs_idx))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_HET = [ FDR_thres_final_HET FDR_thres_final ];
            adjacent_matrix_HET1(:, :, i)   = r;
            zadjacent_matrix_HET1(:, :, i)  = z;
            adjacent_matrix_HET2(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_HET2(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_HET3(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_HET3(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_HET = [ min_thres_r_HET min(temp_r(temp_r~=0)) ];
            min_thres_z_HET = [ min_thres_z_HET min(temp_z(temp_z~=0)) ];
            
        end
        
        % PMG
        for i = 1:num_of_PMG
            
            i
            ts_temp = ts_parcel_PMG_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_PMG(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_PMG(i) strcmp(Sex_PMG(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(length(Anatomical_Label_structs_idx))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_PMG = [ FDR_thres_final_PMG FDR_thres_final ];
            adjacent_matrix_PMG1(:, :, i)   = r;
            zadjacent_matrix_PMG1(:, :, i)  = z;
            adjacent_matrix_PMG2(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_PMG2(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_PMG3(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_PMG3(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_PMG = [ min_thres_r_PMG min(temp_r(temp_r~=0)) ];
            min_thres_z_PMG = [ min_thres_z_PMG min(temp_z(temp_z~=0)) ];
            
        end
        
    end
    
    for micro_AAL = 1
        
        adjacent_matrix_control1_micro              = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_cont);
        adjacent_matrix_FCD_Type_II1_micro          = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_FCD_Type_II);
        adjacent_matrix_HET1_micro                  = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_HET);
        adjacent_matrix_PMG1_micro                  = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_PMG);
        zadjacent_matrix_control1_micro             = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_cont);
        zadjacent_matrix_FCD_Type_II1_micro         = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_FCD_Type_II);
        zadjacent_matrix_HET1_micro                 = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_HET);
        zadjacent_matrix_PMG1_micro                 = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_PMG);
        
        adjacent_matrix_control2_micro              = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_cont);
        adjacent_matrix_FCD_Type_II2_micro          = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_FCD_Type_II);
        adjacent_matrix_HET2_micro                  = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_HET);
        adjacent_matrix_PMG2_micro                  = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_PMG);
        zadjacent_matrix_control2_micro             = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_cont);
        zadjacent_matrix_FCD_Type_II2_micro         = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_FCD_Type_II);
        zadjacent_matrix_HET2_micro                 = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_HET);
        zadjacent_matrix_PMG2_micro                 = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_PMG);
        
        adjacent_matrix_control3_micro              = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_cont);
        adjacent_matrix_FCD_Type_II3_micro          = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_FCD_Type_II);
        adjacent_matrix_HET3_micro                  = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_HET);
        adjacent_matrix_PMG3_micro                  = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_PMG);
        zadjacent_matrix_control3_micro             = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_cont);
        zadjacent_matrix_FCD_Type_II3_micro         = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_FCD_Type_II);
        zadjacent_matrix_HET3_micro                 = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_HET);
        zadjacent_matrix_PMG3_micro                 = zeros(max(microAAL_surf_data_both), max(microAAL_surf_data_both), num_of_PMG);
        
        FDR_thres_final_cont_micro            = [];
        FDR_thres_final_FCD_Type_II_micro     = [];
        FDR_thres_final_HET_micro             = [];
        FDR_thres_final_PMG_micro             = [];
        
        min_thres_r_cont_micro            = [];
        min_thres_r_FCD_Type_II_micro     = [];
        min_thres_r_HET_micro             = [];
        min_thres_r_PMG_micro             = [];
        
        min_thres_z_cont_micro            = [];
        min_thres_z_FCD_Type_II_micro     = [];
        min_thres_z_HET_micro             = [];
        min_thres_z_PMG_micro             = [];
        
        % control
        for i = 1:num_of_cont
            
            i
            ts_temp = ts_parcel_cont_micro_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_cont(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_cont(i) strcmp(Sex_cont(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(max(microAAL_surf_data_both))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_cont_micro = [ FDR_thres_final_cont_micro FDR_thres_final ];
            adjacent_matrix_control1_micro(:, :, i)   = r;
            zadjacent_matrix_control1_micro(:, :, i)  = z;
            adjacent_matrix_control2_micro(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_control2_micro(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_control3_micro(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_control3_micro(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_cont_micro = [ min_thres_r_cont_micro min(temp_r(temp_r~=0)) ];
            min_thres_z_cont_micro = [ min_thres_z_cont_micro min(temp_z(temp_z~=0)) ];
            
        end
        
        % FCD Type-II
        for i = 1:num_of_FCD_Type_II
            
            i
            ts_temp = ts_parcel_FCD_Type_II_micro_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_FCD_Type_II(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_FCD_Type_II(i) strcmp(Sex_FCD_Type_II(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(max(microAAL_surf_data_both))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_FCD_Type_II_micro = [ FDR_thres_final_FCD_Type_II_micro FDR_thres_final ];
            adjacent_matrix_FCD_Type_II1_micro(:, :, i)   = r;
            zadjacent_matrix_FCD_Type_II1_micro(:, :, i)  = z;
            adjacent_matrix_FCD_Type_II2_micro(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_FCD_Type_II2_micro(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_FCD_Type_II3_micro(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_FCD_Type_II3_micro(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_FCD_Type_II_micro = [ min_thres_r_FCD_Type_II_micro min(temp_r(temp_r~=0)) ];
            min_thres_z_FCD_Type_II_micro = [ min_thres_z_FCD_Type_II_micro min(temp_z(temp_z~=0)) ];
            
        end
        
        % HET
        for i = 1:num_of_HET
            
            i
            ts_temp = ts_parcel_HET_micro_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_HET(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_HET(i) strcmp(Sex_HET(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(max(microAAL_surf_data_both))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_HET_micro = [ FDR_thres_final_HET_micro FDR_thres_final ];
            adjacent_matrix_HET1_micro(:, :, i)   = r;
            zadjacent_matrix_HET1_micro(:, :, i)  = z;
            adjacent_matrix_HET2_micro(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_HET2_micro(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_HET3_micro(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_HET3_micro(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_HET_micro = [ min_thres_r_HET_micro min(temp_r(temp_r~=0)) ];
            min_thres_z_HET_micro = [ min_thres_z_HET_micro min(temp_z(temp_z~=0)) ];
            
        end
        
        % PMG
        for i = 1:num_of_PMG
            
            i
            ts_temp = ts_parcel_PMG_micro_res(:, :, i);
            [ r, pval ] = corr(ts_temp);
            %         ts_temp = ts_parcel_PMG(:, :, i);
            %         [ r, pval ] = partialcorr(ts_temp, [ repmat([ Age_PMG(i) strcmp(Sex_PMG(i), 'male') ], size(ts_temp, 1), 1) mean(ts_temp, 2) ]);
            z = 0.5*log((1+r)./(1-r));
            z(logical(eye(max(microAAL_surf_data_both))))=0;
            FDR_thres_final  = FDR(pval, FDR_thres);
            FDR_thres_final_PMG_micro = [ FDR_thres_final_PMG_micro FDR_thres_final ];
            adjacent_matrix_PMG1_micro(:, :, i)   = r;
            zadjacent_matrix_PMG1_micro(:, :, i)  = z;
            adjacent_matrix_PMG2_micro(:, :, i)   = r.*((pval <= FDR_thres_final) .* (r > 0));
            zadjacent_matrix_PMG2_micro(:, :, i)  = z.*((pval <= FDR_thres_final) .* (r > 0));
            adjacent_matrix_PMG3_micro(:, :, i)   = r.*(pval <= FDR_thres_final);
            zadjacent_matrix_PMG3_micro(:, :, i)  = z.*(pval <= FDR_thres_final);
            
            temp_r = r.*((pval <= FDR_thres_final) .* (r > 0));
            temp_z = z.*((pval <= FDR_thres_final) .* (z > 0));
            min_thres_r_PMG_micro = [ min_thres_r_PMG_micro min(temp_r(temp_r~=0)) ];
            min_thres_z_PMG_micro = [ min_thres_z_PMG_micro min(temp_z(temp_z~=0)) ];
            
        end
        
    end
        
    %% check which density has a fully connected matrix: 0.11
    iternation_num = 100;
    interval_thres = 0.01;
    myspar = 0.05:interval_thres:0.4;
    AAL_index_iter = length(Anatomical_Label_structs_idx);
    type_matrix = 1;
    
    for full_conn = 1
        
        node_degree_uw_pos_cont             = zeros(length(myspar), AAL_index_iter);
        node_degree_uw_pos_FCD_Type_II      = zeros(length(myspar), AAL_index_iter);
        node_degree_uw_pos_HET              = zeros(length(myspar), AAL_index_iter);
        node_degree_uw_pos_PMG              = zeros(length(myspar), AAL_index_iter);        
        node_strength_uw_pos_cont           = zeros(length(myspar), AAL_index_iter);
        node_strength_uw_pos_FCD_Type_II    = zeros(length(myspar), AAL_index_iter);
        node_strength_uw_pos_HET            = zeros(length(myspar), AAL_index_iter);
        node_strength_uw_pos_PMG            = zeros(length(myspar), AAL_index_iter);
        
        if(type_matrix == 1)        % raw connectivity matrix
            mat_cont        = zadjacent_matrix_control1;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II1;
            mat_HET         = zadjacent_matrix_HET1;
            mat_PMG         = zadjacent_matrix_PMG1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction            
            mat_cont        = zadjacent_matrix_control2;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II2;
            mat_HET         = zadjacent_matrix_HET2;
            mat_PMG         = zadjacent_matrix_PMG2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            mat_cont        = zadjacent_matrix_control3;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II3;
            mat_HET         = zadjacent_matrix_HET3;
            mat_PMG         = zadjacent_matrix_PMG3;
        end
        
        node_fully_connected                = ones(1, length(myspar));
        
        for i = 1 : length(myspar)
            
            %% controls
            wmatrix                                  = threshold_proportional(mean(mat_cont, 3), myspar(i));
            node_degree_uw_pos_cont(i, :)            = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_cont(i, :)          = strengths_und(double(wmatrix>0));            
                        
            %% FCD Type-II
            wmatrix                                  = threshold_proportional(mean(mat_FCD_Type_II, 3), myspar(i));
            node_degree_uw_pos_FCD_Type_II(i, :)     = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_FCD_Type_II(i, :)   = strengths_und(double(wmatrix>0));
            
            %% HET
            wmatrix                                  = threshold_proportional(mean(mat_HET, 3), myspar(i));
            node_degree_uw_pos_HET(i, :)             = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_HET(i, :)           = strengths_und(double(wmatrix>0));
            
            %% PMG
            wmatrix                                  = threshold_proportional(mean(mat_PMG, 3), myspar(i));
            node_degree_uw_pos_PMG(i, :)             = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_PMG(i, :)           = strengths_und(double(wmatrix>0));
            
            %% check if fully connected
            if((sum(node_degree_uw_pos_cont(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_II(i, :)==0) + sum(node_degree_uw_pos_HET(i, :)==0) + sum(node_degree_uw_pos_PMG(i, :)==0)) > 0)
                node_fully_connected(i) = 0;
            end
            
        end   
                
        node_degree_uw_pos_cont_micro             = zeros(length(myspar), max(microAAL_surf_data_both));
        node_degree_uw_pos_FCD_Type_II_micro      = zeros(length(myspar), max(microAAL_surf_data_both));
        node_degree_uw_pos_HET_micro              = zeros(length(myspar), max(microAAL_surf_data_both));
        node_degree_uw_pos_PMG_micro              = zeros(length(myspar), max(microAAL_surf_data_both));        
        node_strength_uw_pos_cont_micro           = zeros(length(myspar), max(microAAL_surf_data_both));
        node_strength_uw_pos_FCD_Type_II_micro    = zeros(length(myspar), max(microAAL_surf_data_both));
        node_strength_uw_pos_HET_micro            = zeros(length(myspar), max(microAAL_surf_data_both));
        node_strength_uw_pos_PMG_micro            = zeros(length(myspar), max(microAAL_surf_data_both));
        
        if(type_matrix == 1)        % raw connectivity matrix
            mat_cont_micro        = zadjacent_matrix_control1_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II1_micro;
            mat_HET_micro         = zadjacent_matrix_HET1_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG1_micro;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction            
            mat_cont_micro        = zadjacent_matrix_control2_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II2_micro;
            mat_HET_micro         = zadjacent_matrix_HET2_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG2_micro;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            mat_cont_micro        = zadjacent_matrix_control3_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II3_micro;
            mat_HET_micro         = zadjacent_matrix_HET3_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG3_micro;
        end
        
        node_fully_connected_micro                = ones(1, length(myspar));
        
        for i = 1 : length(myspar)
            
            %% controls
            wmatrix                                        = threshold_proportional(mean(mat_cont_micro, 3), myspar(i));
            node_degree_uw_pos_cont_micro(i, :)            = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_cont_micro(i, :)          = strengths_und(double(wmatrix>0));            
                        
            %% FCD Type-II
            wmatrix                                        = threshold_proportional(mean(mat_FCD_Type_II_micro, 3), myspar(i));
            node_degree_uw_pos_FCD_Type_II_micro(i, :)     = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_FCD_Type_II_micro(i, :)   = strengths_und(double(wmatrix>0));
            
            %% HET
            wmatrix                                        = threshold_proportional(mean(mat_HET_micro, 3), myspar(i));
            node_degree_uw_pos_HET_micro(i, :)             = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_HET_micro(i, :)           = strengths_und(double(wmatrix>0));
            
            %% PMG
            wmatrix                                        = threshold_proportional(mean(mat_PMG_micro, 3), myspar(i));
            node_degree_uw_pos_PMG_micro(i, :)             = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_PMG_micro(i, :)           = strengths_und(double(wmatrix>0));
            
            %% check if fully connected
            if((sum(node_degree_uw_pos_cont_micro(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_II_micro(i, :)==0) + sum(node_degree_uw_pos_HET_micro(i, :)==0) + sum(node_degree_uw_pos_PMG_micro(i, :)==0)) > 0)
                node_fully_connected_micro(i) = 0;
            end
            
        end           
        
    end
    
    %% visualize the matrix and save files
    if(printfigs == 1)

        savefigs = 0;
        load('colormap_asso_mat.mat');
        load('colormap_asso_mat2.mat');
        struct_idx = [ 0 13 17 22 29 33 39 52 56 61 68 72 78 ];
        vertical_line_idx_x = [ struct_idx; struct_idx ];
        vertical_line_idx_y = [ zeros(1, length(struct_idx)); ones(1, length(struct_idx))*78 ];
        
        for i = 1 : size(struct_idx, 2)
            
            if(struct_idx(i) == 0)
                micro_struct_idx = 0;
            else
                temp = [];
                aal_idx = Anatomical_Label_structs_idx((1+struct_idx(i-1)):struct_idx(i));
                for j = 1 : size(aal_idx, 2)
                    temp = [ temp length(unique(microAAL_surf_data_both(AAL_surf_data_both == aal_idx(j)))) ];
                end
                micro_struct_idx = [ micro_struct_idx max(micro_struct_idx)+sum(temp) ];
            end
            
        end
        
        micro_vertical_line_idx_x = [ micro_struct_idx; micro_struct_idx ];
        micro_vertical_line_idx_y = [ zeros(1, length(micro_struct_idx)); ones(1, length(micro_struct_idx))*max(microAAL_surf_data_both) ];
        type_matrix = 2;
        if(type_matrix == 1)        % raw connectivity matrix
            mat_cont        = zadjacent_matrix_control1;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II1;
            mat_HET         = zadjacent_matrix_HET1;
            mat_PMG         = zadjacent_matrix_PMG1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction            
            mat_cont        = zadjacent_matrix_control2;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II2;
            mat_HET         = zadjacent_matrix_HET2;
            mat_PMG         = zadjacent_matrix_PMG2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            mat_cont        = zadjacent_matrix_control3;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II3;
            mat_HET         = zadjacent_matrix_HET3;
            mat_PMG         = zadjacent_matrix_PMG3;
        elseif(type_matrix == 4)    % all connectivity matrix after FDR correction
            mat_cont        = zadjacent_matrix_control;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II;
            mat_HET         = zadjacent_matrix_HET;
            mat_PMG         = zadjacent_matrix_PMG;
        end        
        if(type_matrix == 1)        % raw connectivity matrix
            mat_cont_micro        = zadjacent_matrix_control1_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II1_micro;
            mat_HET_micro         = zadjacent_matrix_HET1_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG1_micro;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            mat_cont_micro        = zadjacent_matrix_control2_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II2_micro;
            mat_HET_micro         = zadjacent_matrix_HET2_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG2_micro;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            mat_cont_micro        = zadjacent_matrix_control3_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II3_micro;
            mat_HET_micro         = zadjacent_matrix_HET3_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG3_micro;
        end
       
        for original_AAL = 1
            
            %% 1) control
            % - original AAL
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_cont, 3),{},[-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_control_rs_fMRI.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_cont, 0, 3),{},[0 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
            %% 2) FCD Type-II
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_FCD_Type_II, 3),{},[-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_FCD_TypeII_rs_fMRI.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_FCD_Type_II, 0, 3),{},[0 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
            %% 3) HET
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_HET, 3),{},[-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_HET_rs_fMRI.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_HET, 0, 3),{},[0 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
            %% 4) PMG
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_PMG, 3),{},[-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_PMG2.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_PMG, 0, 3),{},[0 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
        end
        
        for micro_AAL = 1
            
            %% 1) control
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_cont_micro, 3),{},[-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_control_rs_fMRI_micro.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_cont_micro, 0, 3),{},[0 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
            %% 2) FCD Type-II
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_FCD_Type_II_micro, 3),{},[-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_FCD_Type_II_rs_fMRI_micro.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_FCD_Type_II_micro, 0, 3),{},[0 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
            %% 3) HET
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_HET_micro, 3),{},[-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_HET_rs_fMRI_micro.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_HET_micro, 0, 3),{},[0 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
            %% 4) PMG
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(mean(mat_PMG_micro, 3),{},[-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            
            if(savefigs == 1)
                exportfigbo(gcf,[OUTPATH 'adj_mat_PMG_rs_fMRI_micro.png'], 'png', 13); close(gcf);
            end
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(std(mat_PMG_micro, 0, 3),{},[0 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            
        end
        
    end 
       
    %% save variables
    if(exist([ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_GSR_Age_Sex1.mat' ])~=2)
        save([ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_Age_Sex1.mat' ], ...
              'adjacent_matrix_control1',  'adjacent_matrix_FCD_Type_II1', 'adjacent_matrix_HET1', 'adjacent_matrix_PMG1', ...
              'zadjacent_matrix_control1', 'zadjacent_matrix_FCD_Type_II1', 'zadjacent_matrix_HET1', 'zadjacent_matrix_PMG1', ...
              'adjacent_matrix_control2',  'adjacent_matrix_FCD_Type_II2', 'adjacent_matrix_HET2', 'adjacent_matrix_PMG2', ...
              'zadjacent_matrix_control2', 'zadjacent_matrix_FCD_Type_II2', 'zadjacent_matrix_HET2', 'zadjacent_matrix_PMG2', ...
              'adjacent_matrix_control3',  'adjacent_matrix_FCD_Type_II3', 'adjacent_matrix_HET3', 'adjacent_matrix_PMG3', ...
              'zadjacent_matrix_control3', 'zadjacent_matrix_FCD_Type_II3', 'zadjacent_matrix_HET3', 'zadjacent_matrix_PMG3', ...
              'Anatomical_Label_structs', 'Anatomical_Label_structs_idx', ...
              'ts_parcel_cont', 'ts_parcel_FCD_Type_II', 'ts_parcel_HET', 'ts_parcel_PMG');
    end
    if(exist([ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_GSR_Age_Sex1_micro.mat' ])~=2)
        save([ MATDIR '/01_generate_association_matrix_MCD_rsfMRI_Age_Sex1_micro.mat' ], ...
              'adjacent_matrix_control1_micro',  'adjacent_matrix_FCD_Type_II1_micro', 'adjacent_matrix_HET1_micro', 'adjacent_matrix_PMG1_micro', ...
              'zadjacent_matrix_control1_micro', 'zadjacent_matrix_FCD_Type_II1_micro', 'zadjacent_matrix_HET1_micro', 'zadjacent_matrix_PMG1_micro', ...
              'adjacent_matrix_control2_micro',  'adjacent_matrix_FCD_Type_II2_micro', 'adjacent_matrix_HET2_micro', 'adjacent_matrix_PMG2_micro', ...
              'zadjacent_matrix_control2_micro', 'zadjacent_matrix_FCD_Type_II2_micro', 'zadjacent_matrix_HET2_micro', 'zadjacent_matrix_PMG2_micro', ...
              'adjacent_matrix_control3_micro',  'adjacent_matrix_FCD_Type_II3_micro', 'adjacent_matrix_HET3_micro', 'adjacent_matrix_PMG3_micro', ...
              'zadjacent_matrix_control3_micro', 'zadjacent_matrix_FCD_Type_II3_micro', 'zadjacent_matrix_HET3_micro', 'zadjacent_matrix_PMG3_micro', ...
              'Anatomical_Label_structs', 'Anatomical_Label_structs_idx', 'microAAL_surf_data_both', ...
              'ts_parcel_cont_micro', 'ts_parcel_FCD_Type_II_micro', 'ts_parcel_HET_micro', 'ts_parcel_PMG_micro');
    end    
    
end

%% statistical test for each connectivity using permutation and FDR
for statistical_test = 1
    
    for original_AAL_ver1 = 1
        
        parpool(24);
        myspar = [ 0.05:0.01:0.4 ];
        iteration = 1000;
        
        type_matrix = 2;
        
        %% 1) control vs. FCD Type-II
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = zadjacent_matrix_control1;
            group2 = zadjacent_matrix_FCD_Type_II1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control2;
            group2 = zadjacent_matrix_FCD_Type_II2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control3;
            group2 = zadjacent_matrix_FCD_Type_II3;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
       
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_corr_1 = mean(tgroup_Reg(:,:,rp(1:N1)), 3);
            rerand_corr_2 = mean(tgroup_Reg(:,:,rp((N1+1):N_tgroup)), 3);
            rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
            rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
            
            wmatrix1_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            wmatrix2_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            diff_wmatrix_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            for s = 1 : length(myspar)
                
                wmatrix1_temp(:, :, s) = threshold_proportional(rerand_corr_1, myspar(s));
                wmatrix2_temp(:, :, s) = threshold_proportional(rerand_corr_2, myspar(s));
                
                diff_wmatrix_temp(:, :, s) = wmatrix1_temp(:, :, s)-wmatrix2_temp(:, :, s);
                
            end
            wmatrix1_cell{i}        = wmatrix1_temp;
            wmatrix2_cell{i}        = wmatrix2_temp;
            diff_wmatrix_cell{i}    = diff_wmatrix_temp;
            
        end
        wmatrix1        = reshape(cell2mat(wmatrix1_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        wmatrix2        = reshape(cell2mat(wmatrix2_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        for s = 1 : length(myspar)
            
            fprintf('\n s = %d\n',s);
            real_wmatrix1 = threshold_proportional(mean(group1, 3), myspar(s));
            real_wmatrix2 = threshold_proportional(mean(group2, 3), myspar(s));
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 : length(Anatomical_Label_structs_idx)
                for j = 1 : length(Anatomical_Label_structs_idx)
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
        
        diff_real_wmatrix_cont_vs_FCD_TypeII = diff_real_wmatrix;
        pmap_mat_cont_vs_FCD_TypeII = pmap_mat;
        
        %% 2) control vs. HET
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = zadjacent_matrix_control1;
            group2 = zadjacent_matrix_HET1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control2;
            group2 = zadjacent_matrix_HET2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control3;
            group2 = zadjacent_matrix_HET3;
        end
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
       
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_corr_1 = mean(tgroup_Reg(:,:,rp(1:N1)), 3);
            rerand_corr_2 = mean(tgroup_Reg(:,:,rp((N1+1):N_tgroup)), 3);
            rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
            rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
            
            wmatrix1_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            wmatrix2_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            diff_wmatrix_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));            
            for s = 1 : length(myspar)
                
                wmatrix1_temp(:, :, s) = threshold_proportional(rerand_corr_1, myspar(s));
                wmatrix2_temp(:, :, s) = threshold_proportional(rerand_corr_2, myspar(s));
                
                diff_wmatrix_temp(:, :, s) = wmatrix1_temp(:, :, s)-wmatrix2_temp(:, :, s);
                
            end
            wmatrix1_cell{i}        = wmatrix1_temp;
            wmatrix2_cell{i}        = wmatrix2_temp;
            diff_wmatrix_cell{i}    = diff_wmatrix_temp;
            
        end
        wmatrix1        = reshape(cell2mat(wmatrix1_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        wmatrix2        = reshape(cell2mat(wmatrix2_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        for s = 1 : length(myspar)
            
            fprintf('\n s = %d\n',s);
            real_wmatrix1 = threshold_proportional(mean(group1, 3), myspar(s));
            real_wmatrix2 = threshold_proportional(mean(group2, 3), myspar(s));
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 : length(Anatomical_Label_structs_idx)
                for j = 1 : length(Anatomical_Label_structs_idx)
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
        
        diff_real_wmatrix_cont_vs_HET = diff_real_wmatrix;
        pmap_mat_cont_vs_HET = pmap_mat;
        
        %% 3) control vs. PMG
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = zadjacent_matrix_control1;
            group2 = zadjacent_matrix_PMG1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control2;
            group2 = zadjacent_matrix_PMG2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control3;
            group2 = zadjacent_matrix_PMG3;
        end
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
       
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_corr_1 = mean(tgroup_Reg(:,:,rp(1:N1)), 3);
            rerand_corr_2 = mean(tgroup_Reg(:,:,rp((N1+1):N_tgroup)), 3);
            rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
            rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
            
            wmatrix1_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            wmatrix2_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            diff_wmatrix_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));            
            for s = 1 : length(myspar)
                
                wmatrix1_temp(:, :, s) = threshold_proportional(rerand_corr_1, myspar(s));
                wmatrix2_temp(:, :, s) = threshold_proportional(rerand_corr_2, myspar(s));
                
                diff_wmatrix_temp(:, :, s) = wmatrix1_temp(:, :, s)-wmatrix2_temp(:, :, s);
                
            end
            wmatrix1_cell{i}        = wmatrix1_temp;
            wmatrix2_cell{i}        = wmatrix2_temp;
            diff_wmatrix_cell{i}    = diff_wmatrix_temp;
            
        end
        wmatrix1        = reshape(cell2mat(wmatrix1_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        wmatrix2        = reshape(cell2mat(wmatrix2_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        for s = 1 : length(myspar)
            
            fprintf('\n s = %d\n',s);
            real_wmatrix1 = threshold_proportional(mean(group1, 3), myspar(s));
            real_wmatrix2 = threshold_proportional(mean(group2, 3), myspar(s));
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 : length(Anatomical_Label_structs_idx)
                for j = 1 : length(Anatomical_Label_structs_idx)
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
        
        diff_real_wmatrix_cont_vs_PMG = diff_real_wmatrix;
        pmap_mat_cont_vs_PMG = pmap_mat;
        
        %% 4) visualize difference
        if(printfigs == 1)
            
            %% variable setup
            load('colormap_asso_mat.mat');
            load('colormap_asso_mat2.mat');
            load('colormap_asso_mat3.mat');
            struct_idx = [ 0 13 17 22 29 33 39 52 56 61 68 72 78 ];
            vertical_line_idx_x = [ struct_idx; struct_idx ];
            vertical_line_idx_y = [ zeros(1, length(struct_idx)); ones(1, length(struct_idx))*78 ];
            spar = 7; %min(find(node_fully_connected));
            FDRQ = 0.05;
            threshold_p = 0.025;
            
            S_new = S;
            S_new.coord(1, 1:40962)     = S_new.coord(1, 1:40962)-10;
            S_new.coord(1, 40963:81924) = S_new.coord(1, 40963:81924)+10;
            
            %% node position of each parcel in AAL
            centroid = zeros(1, length(Anatomical_Label_structs_idx));
            for i = 1 : length(Anatomical_Label_structs_idx)
                
                vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                vert_coord = S.coord(:, vert_idx);
                [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                centroid(i) = vert_idx(b);
                
            end
            
            nodePos = S_new.coord(:, centroid);
            
            %% 1) control vs. FCD Type-II
            pmap_mat            = pmap_mat_cont_vs_FCD_TypeII(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_FCD_TypeII(:, :, spar);
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_FCD_Type_II_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_FCD_Type_II_fMRI.png'], 'png', 13); close(gcf);
            
            %% 2) control vs. HET
            pmap_mat            = pmap_mat_cont_vs_HET(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_HET(:, :, spar);
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1] = NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
            
            %% 3) control vs. PMG
            pmap_mat            = pmap_mat_cont_vs_PMG(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_PMG(:, :, spar);
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.4, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
            
        end
        
    end
    
    for original_AAL_ver2 = 1
        
        parpool(24);
        myspar = [ 0.05:0.01:0.4 ];
        iteration = 1000;
        
        type_matrix = 2;

        if(type_matrix == 1)        % raw connectivity matrix
            mat_cont        = zadjacent_matrix_control1;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II1;
            mat_HET         = zadjacent_matrix_HET1;
            mat_PMG         = zadjacent_matrix_PMG1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction            
            mat_cont        = zadjacent_matrix_control2;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II2;
            mat_HET         = zadjacent_matrix_HET2;
            mat_PMG         = zadjacent_matrix_PMG2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            mat_cont        = zadjacent_matrix_control3;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II3;
            mat_HET         = zadjacent_matrix_HET3;
            mat_PMG         = zadjacent_matrix_PMG3;
        elseif(type_matrix == 4)    % all connectivity matrix after FDR correction
            mat_cont        = zadjacent_matrix_control;
            mat_FCD_Type_II = zadjacent_matrix_FCD_Type_II;
            mat_HET         = zadjacent_matrix_HET;
            mat_PMG         = zadjacent_matrix_PMG;
        end
        
        mat_thres = 0.5;
        
        %% 0) binarized averaged matrix construction
        averaged_mat_cont        = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
        averaged_mat_FCD_Type_II = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
        averaged_mat_HET         = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
        averaged_mat_PMG         = zeros(AAL_index_iter ,AAL_index_iter, length(myspar));
        
        node_fully_connected = ones(1, length(myspar));
        for i = 1 : length(myspar)
           
            %% control
            mat_cont_temp = mat_cont*0;
            for j = 1 : num_of_cont
                 temp = threshold_proportional(mat_cont(:, :, j), myspar(i));
                 mat_cont_temp(:, :, j) = temp > 0;
                 mat_cont_temp2(:, :, j) = temp;
            end
%             averaged_mat_cont(:, :, i) = (sum(mat_cont_temp, 3)/num_of_cont)>mat_thres;
            averaged_mat_cont(:, :, i) = mean(mat_cont_temp2, 3).*(mean(mat_cont_temp, 3)>mat_thres);
            
            %% FCD Type II
            mat_FCD_Type_II_temp = mat_FCD_Type_II*0;
            for j = 1 : num_of_FCD_Type_II
                temp = threshold_proportional(mat_FCD_Type_II(:, :, j), myspar(i));
                mat_FCD_Type_II_temp(:, :, j) = temp > 0;
                mat_FCD_Type_II_temp2(:, :, j) = temp;
            end
%             averaged_mat_FCD_Type_II(:, :, i) = (sum(mat_FCD_Type_II_temp, 3)/num_of_FCD_Type_II)>mat_thres;
            averaged_mat_FCD_Type_II(:, :, i) = mean(mat_FCD_Type_II_temp2, 3).*(mean(mat_FCD_Type_II_temp, 3)>mat_thres);
            
            %% HET
            mat_HET_temp = mat_HET*0;
            for j = 1 : num_of_HET
                temp = threshold_proportional(mat_HET(:, :, j), myspar(i));
                mat_HET_temp(:, :, j) = temp > 0;
                mat_HET_temp2(:, :, j) = temp;
            end
%             averaged_mat_HET(:, :, i) = (sum(mat_HET_temp, 3)/num_of_HET)>mat_thres;
            averaged_mat_HET(:, :, i) = mean(mat_HET_temp2, 3).*(mean(mat_HET_temp, 3)>mat_thres);
            
            %% PMG
            mat_PMG_temp = mat_PMG*0;
            for j = 1 : num_of_PMG
                temp = threshold_proportional(mat_PMG(:, :, j), myspar(i));
                mat_PMG_temp(:, :, j) = temp > 0;
                mat_PMG_temp2(:, :, j) = temp;
            end
%             averaged_mat_PMG(:, :, i) = (sum(mat_PMG_temp, 3)/num_of_PMG)>mat_thres;
            averaged_mat_PMG(:, :, i) = mean(mat_PMG_temp2, 3).*(mean(mat_PMG_temp, 3)>mat_thres);

            %% Connectedness checking
            if(sum(sum(averaged_mat_cont(:, :, i)) == 0) > 0 | ...
               sum(sum(averaged_mat_FCD_Type_II(:, :, i)) == 0) > 0 | ...
               sum(sum(averaged_mat_HET(:, :, i)) == 0) > 0 | ...
               sum(sum(averaged_mat_PMG(:, :, i)) == 0) > 0)
           
               node_fully_connected(i) = 0;
            end
                 
        end
        
        %% 1) control vs. FCD Type-II
        type_matrix = 1;
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = mat_cont;
            group2 = mat_FCD_Type_II;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = mat_cont;
            group2 = mat_HET;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = mat_cont;
            group2 = mat_PMG;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
       
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_group1 = tgroup_Reg(:,:,rp(1:N1));
            rerand_group2 = tgroup_Reg(:,:,rp((N1+1):N_tgroup));
            
            wmatrix1_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            wmatrix2_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            diff_wmatrix_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));            
            for s = 1 : length(myspar)
                
                mat_temp  = rerand_group1*0;
                mat_temp2 = rerand_group1*0;
                for j = 1 : N1
                    
                    temp = threshold_proportional(rerand_group1(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_1 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                mat_temp  = rerand_group2*0;
                mat_temp2 = rerand_group2*0;
                for j = 1 : N2
                    
                    temp = threshold_proportional(rerand_group2(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_2 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                
                wmatrix1_temp(:, :, s) = threshold_proportional(rerand_corr_1, myspar(s));
                wmatrix2_temp(:, :, s) = threshold_proportional(rerand_corr_2, myspar(s));
                
                diff_wmatrix_temp(:, :, s) = wmatrix1_temp(:, :, s)-wmatrix2_temp(:, :, s);
                
            end
            wmatrix1_cell{i}        = wmatrix1_temp;
            wmatrix2_cell{i}        = wmatrix2_temp;
            diff_wmatrix_cell{i}    = diff_wmatrix_temp;
            
        end
        wmatrix1        = reshape(cell2mat(wmatrix1_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        wmatrix2        = reshape(cell2mat(wmatrix2_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        for s = 1 : length(myspar)
            
            fprintf('\n s = %d\n',s);
            real_wmatrix1 = averaged_mat_cont(:, :, s);
            real_wmatrix2 = averaged_mat_FCD_Type_II(:, :, s);
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 : length(Anatomical_Label_structs_idx)
                for j = 1 : length(Anatomical_Label_structs_idx)
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
        
        diff_real_wmatrix_cont_vs_FCD_TypeII = diff_real_wmatrix;
        pmap_mat_cont_vs_FCD_TypeII = pmap_mat;
        
        %% 2) control vs. HET
        type_matrix = 2;
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = mat_cont;
            group2 = mat_FCD_Type_II;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = mat_cont;
            group2 = mat_HET;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = mat_cont;
            group2 = mat_PMG;
        end
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
       
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_group1 = tgroup_Reg(:,:,rp(1:N1));
            rerand_group2 = tgroup_Reg(:,:,rp((N1+1):N_tgroup));
            
            wmatrix1_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            wmatrix2_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            diff_wmatrix_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));            
            for s = 1 : length(myspar)
                
                mat_temp  = rerand_group1*0;
                mat_temp2 = rerand_group1*0;
                for j = 1 : N1
                    
                    temp = threshold_proportional(rerand_group1(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_1 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                mat_temp  = rerand_group2*0;
                mat_temp2 = rerand_group2*0;
                for j = 1 : N2
                    
                    temp = threshold_proportional(rerand_group2(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_2 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                
                wmatrix1_temp(:, :, s) = threshold_proportional(rerand_corr_1, myspar(s));
                wmatrix2_temp(:, :, s) = threshold_proportional(rerand_corr_2, myspar(s));
                
                diff_wmatrix_temp(:, :, s) = wmatrix1_temp(:, :, s)-wmatrix2_temp(:, :, s);
                
            end
            wmatrix1_cell{i}        = wmatrix1_temp;
            wmatrix2_cell{i}        = wmatrix2_temp;
            diff_wmatrix_cell{i}    = diff_wmatrix_temp;
            
        end
        wmatrix1        = reshape(cell2mat(wmatrix1_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        wmatrix2        = reshape(cell2mat(wmatrix2_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        for s = 1 : length(myspar)
            
            fprintf('\n s = %d\n',s);
            real_wmatrix1 = averaged_mat_cont(:, :, s);
            real_wmatrix2 = averaged_mat_HET(:, :, s);
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 : length(Anatomical_Label_structs_idx)
                for j = 1 : length(Anatomical_Label_structs_idx)
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
        
        diff_real_wmatrix_cont_vs_HET = diff_real_wmatrix;
        pmap_mat_cont_vs_HET = pmap_mat;
        
        %% 3) control vs. PMG
        type_matrix = 3;
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = mat_cont;
            group2 = mat_FCD_Type_II;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = mat_cont;
            group2 = mat_HET;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = mat_cont;
            group2 = mat_PMG;
        end
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
       
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_group1 = tgroup_Reg(:,:,rp(1:N1));
            rerand_group2 = tgroup_Reg(:,:,rp((N1+1):N_tgroup));
            
            wmatrix1_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            wmatrix2_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
            diff_wmatrix_temp = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));            
            for s = 1 : length(myspar)
                
                mat_temp  = rerand_group1*0;
                mat_temp2 = rerand_group1*0;
                for j = 1 : N1
                    
                    temp = threshold_proportional(rerand_group1(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_1 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                mat_temp  = rerand_group2*0;
                mat_temp2 = rerand_group2*0;
                for j = 1 : N2
                    
                    temp = threshold_proportional(rerand_group2(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_2 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                
                wmatrix1_temp(:, :, s) = threshold_proportional(rerand_corr_1, myspar(s));
                wmatrix2_temp(:, :, s) = threshold_proportional(rerand_corr_2, myspar(s));
                
                diff_wmatrix_temp(:, :, s) = wmatrix1_temp(:, :, s)-wmatrix2_temp(:, :, s);
                
            end
            wmatrix1_cell{i}        = wmatrix1_temp;
            wmatrix2_cell{i}        = wmatrix2_temp;
            diff_wmatrix_cell{i}    = diff_wmatrix_temp;
            
        end
        wmatrix1        = reshape(cell2mat(wmatrix1_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        wmatrix2        = reshape(cell2mat(wmatrix2_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
        for s = 1 : length(myspar)
            
            fprintf('\n s = %d\n',s);
            real_wmatrix1 = averaged_mat_cont(:, :, s);
            real_wmatrix2 = averaged_mat_PMG(:, :, s);
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 : length(Anatomical_Label_structs_idx)
                for j = 1 : length(Anatomical_Label_structs_idx)
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, s, :)-wmatrix2(i, j, s, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
        
        diff_real_wmatrix_cont_vs_PMG = diff_real_wmatrix;
        pmap_mat_cont_vs_PMG = pmap_mat;
        
        %% 4) visualize difference
        if(printfigs == 1)
            
            %% variable setup
            load('colormap_asso_mat.mat');
            load('colormap_asso_mat2.mat');
            load('colormap_asso_mat3.mat');
            struct_idx = [ 0 13 17 22 29 33 39 52 56 61 68 72 78 ];
            vertical_line_idx_x = [ struct_idx; struct_idx ];
            vertical_line_idx_y = [ zeros(1, length(struct_idx)); ones(1, length(struct_idx))*78 ];
            spar = 7; %min(find(node_fully_connected));
            FDRQ = 0.05;
            threshold_p = 0.05;
            
            S_new = S;
            S_new.coord(1, 1:40962)     = S_new.coord(1, 1:40962)-3;
            S_new.coord(1, 40963:81924) = S_new.coord(1, 40963:81924)+3;
            
            %% node position of each parcel in AAL
            centroid = zeros(1, length(Anatomical_Label_structs_idx));
            for i = 1 : length(Anatomical_Label_structs_idx)
                
                vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                vert_coord = S.coord(:, vert_idx);
                [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                centroid(i) = vert_idx(b);
                
            end
            
            nodePos = S_new.coord(:, centroid);
            
            %% 1) control vs. FCD Type-II
            pmap_mat            = pmap_mat_cont_vs_FCD_TypeII(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_FCD_TypeII(:, :, spar);
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_FCD_Type_II_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_FCD_Type_II_fMRI.png'], 'png', 13); close(gcf);
            
            %% 2) control vs. HET
            pmap_mat            = pmap_mat_cont_vs_HET(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_HET(:, :, spar);
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1] = NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
            
            %% 3) control vs. PMG
            pmap_mat            = pmap_mat_cont_vs_PMG(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_PMG(:, :, spar);
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
            
        end
        
    end
    
    for micro_AAL_ver2 = 1
        
        parpool(24);
        myspar = [ 0.11:0.01:0.11 ];
        iteration = 1000;
        
        type_matrix = 2;

        if(type_matrix == 1)        % raw connectivity matrix
            mat_cont_micro        = zadjacent_matrix_control1_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II1_micro;
            mat_HET_micro         = zadjacent_matrix_HET1_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG1_micro;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            mat_cont_micro        = zadjacent_matrix_control2_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II2_micro;
            mat_HET_micro         = zadjacent_matrix_HET2_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG2_micro;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            mat_cont_micro        = zadjacent_matrix_control3_micro;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II3_micro;
            mat_HET_micro         = zadjacent_matrix_HET3_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG3_micro;       
        elseif(type_matrix == 4)    % all connectivity matrix after FDR correction
            mat_cont_micro        = zadjacent_matrix_control;
            mat_FCD_Type_II_micro = zadjacent_matrix_FCD_Type_II_micro;
            mat_HET_micro         = zadjacent_matrix_HET_micro;
            mat_PMG_micro         = zadjacent_matrix_PMG_micro;
        end
        
        mat_thres = 0.5;
        micro_AAL_index_iter = 1166;
        
        %% 0) binarized averaged matrix construction
        averaged_mat_cont_micro        = zeros(micro_AAL_index_iter ,micro_AAL_index_iter, length(myspar));
        averaged_mat_FCD_Type_II_micro = zeros(micro_AAL_index_iter ,micro_AAL_index_iter, length(myspar));
        averaged_mat_HET_micro         = zeros(micro_AAL_index_iter ,micro_AAL_index_iter, length(myspar));
        averaged_mat_PMG_micro         = zeros(micro_AAL_index_iter ,micro_AAL_index_iter, length(myspar));
        
        node_fully_connected_micro = ones(1, length(myspar));
        for i = 1 : length(myspar)
           
            i
            
            %% control
            mat_cont_temp = mat_cont_micro*0;
            for j = 1 : num_of_cont
                 temp = threshold_proportional(mat_cont_micro(:, :, j), myspar(i));
                 mat_cont_temp(:, :, j) = temp > 0;
                 mat_cont_temp2(:, :, j) = temp;
            end
%             averaged_mat_cont(:, :, i) = (sum(mat_cont_temp, 3)/num_of_cont)>mat_thres;
            averaged_mat_cont_micro(:, :, i) = mean(mat_cont_temp2, 3).*(mean(mat_cont_temp, 3)>mat_thres);
            
            %% FCD Type II
            mat_FCD_Type_II_temp = mat_FCD_Type_II_micro*0;
            for j = 1 : num_of_FCD_Type_II
                temp = threshold_proportional(mat_FCD_Type_II_micro(:, :, j), myspar(i));
                mat_FCD_Type_II_temp(:, :, j) = temp > 0;
                mat_FCD_Type_II_temp2(:, :, j) = temp;
            end
%             averaged_mat_FCD_Type_II(:, :, i) = (sum(mat_FCD_Type_II_temp, 3)/num_of_FCD_Type_II)>mat_thres;
            averaged_mat_FCD_Type_II_micro(:, :, i) = mean(mat_FCD_Type_II_temp2, 3).*(mean(mat_FCD_Type_II_temp, 3)>mat_thres);
            
            %% HET
            mat_HET_temp = mat_HET_micro*0;
            for j = 1 : num_of_HET
                temp = threshold_proportional(mat_HET_micro(:, :, j), myspar(i));
                mat_HET_temp(:, :, j) = temp > 0;
                mat_HET_temp2(:, :, j) = temp;
            end
%             averaged_mat_HET(:, :, i) = (sum(mat_HET_temp, 3)/num_of_HET)>mat_thres;
            averaged_mat_HET_micro(:, :, i) = mean(mat_HET_temp2, 3).*(mean(mat_HET_temp, 3)>mat_thres);
            
            %% PMG
            mat_PMG_temp = mat_PMG_micro*0;
            for j = 1 : num_of_PMG
                temp = threshold_proportional(mat_PMG_micro(:, :, j), myspar(i));
                mat_PMG_temp(:, :, j) = temp > 0;
                mat_PMG_temp2(:, :, j) = temp;
            end
%             averaged_mat_PMG(:, :, i) = (sum(mat_PMG_temp, 3)/num_of_PMG)>mat_thres;
            averaged_mat_PMG_micro(:, :, i) = mean(mat_PMG_temp2, 3).*(mean(mat_PMG_temp, 3)>mat_thres);

            %% Connectedness checking
            if(sum(sum(averaged_mat_cont_micro(:, :, i)) == 0) > 0 | ...
               sum(sum(averaged_mat_FCD_Type_II_micro(:, :, i)) == 0) > 0 | ...
               sum(sum(averaged_mat_HET_micro(:, :, i)) == 0) > 0 | ...
               sum(sum(averaged_mat_PMG_micro(:, :, i)) == 0) > 0)
           
               node_fully_connected_micro(i) = 0;
            end
                 
        end
        
        %% 1) control vs. FCD Type-II
        type_matrix = 1;
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = mat_cont_micro;
            group2 = mat_FCD_Type_II_micro;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = mat_cont_micro;
            group2 = mat_HET_micro;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = mat_cont_micro;
            group2 = mat_PMG_micro;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
                
        diff_real_wmatrix = zeros(micro_AAL_index_iter,micro_AAL_index_iter,length(myspar));
        pmap_mat = zeros(micro_AAL_index_iter,micro_AAL_index_iter,length(myspar));;
                    
        rp_set = zeros(iteration, N_tgroup);
        for i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp_set(i, :) = randperm(N_tgroup);
            
        end
        
        for s = 1 : length(myspar)
            
            wmatrix1_temp = cell(iteration, 1);
            wmatrix2_temp = cell(iteration, 1);
            diff_wmatrix_temp = cell(iteration, 1);
            
            parfor i = 1 : iteration
                
                fprintf('\n s = %d & i = %d\n', s, i);
                rp = rp_set(i, :);
                
                rerand_group1 = single(tgroup_Reg(:,:,rp(1:N1)));
                rerand_group2 = single(tgroup_Reg(:,:,rp((N1+1):N_tgroup)));
                      
                mat_temp  = rerand_group1*0;
                mat_temp2 = rerand_group1*0;
                for j = 1 : N1
                    
                    temp = threshold_proportional(rerand_group1(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_1 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                mat_temp  = rerand_group2*0;
                mat_temp2 = rerand_group2*0;
                for j = 1 : N2
                    
                    temp = threshold_proportional(rerand_group2(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_2 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                
                wmatrix1_temp{i} = rerand_corr_1;
                wmatrix2_temp{i} = rerand_corr_2;               
                diff_wmatrix_temp{i} = wmatrix1_temp{i}-wmatrix2_temp{i};
                
            end 
            
            wmatrix1     = zeros(1166, 1166, iteration);
            wmatrix2     = zeros(1166, 1166, iteration);
            diff_wmatrix = zeros(1166, 1166, iteration);
            for i = 1 : iteration
                
                wmatrix1(:, :, i) = wmatrix1_temp{i};
                wmatrix2(:, :, i) = wmatrix2_temp{i};
                diff_wmatrix(:, :, i) = diff_wmatrix_temp{i};
                
            end
       
            real_wmatrix1 = averaged_mat_cont_micro(:, :, s);
            real_wmatrix2 = averaged_mat_FCD_Type_II_micro(:, :, s);
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 :micro_AAL_index_iter
                for j = 1 : micro_AAL_index_iter
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, :)-wmatrix2(i, j, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, :)-wmatrix2(i, j, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
       
        diff_real_wmatrix_cont_vs_FCD_TypeII_micro = diff_real_wmatrix;
        pmap_mat_cont_vs_FCD_TypeII_micro = pmap_mat;
        
        %% 2) control vs. HET
        type_matrix = 2;
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = mat_cont_micro;
            group2 = mat_FCD_Type_II_micro;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = mat_cont_micro;
            group2 = mat_HET_micro;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = mat_cont_micro;
            group2 = mat_PMG_micro;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
                
        diff_real_wmatrix = zeros(micro_AAL_index_iter,micro_AAL_index_iter,length(myspar));
        pmap_mat = zeros(micro_AAL_index_iter,micro_AAL_index_iter,length(myspar));;
                    
        rp_set = zeros(iteration, N_tgroup);
        for i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp_set(i, :) = randperm(N_tgroup);
            
        end
        
        for s = 1 : length(myspar)
            
            wmatrix1_temp = cell(iteration, 1);
            wmatrix2_temp = cell(iteration, 1);
            diff_wmatrix_temp = cell(iteration, 1);
            
            parfor i = 1 : iteration
                
                fprintf('\n s = %d & i = %d\n', s, i);
                rp = rp_set(i, :);
                
                rerand_group1 = single(tgroup_Reg(:,:,rp(1:N1)));
                rerand_group2 = single(tgroup_Reg(:,:,rp((N1+1):N_tgroup)));
                      
                mat_temp  = rerand_group1*0;
                mat_temp2 = rerand_group1*0;
                for j = 1 : N1
                    
                    temp = threshold_proportional(rerand_group1(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_1 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                mat_temp  = rerand_group2*0;
                mat_temp2 = rerand_group2*0;
                for j = 1 : N2
                    
                    temp = threshold_proportional(rerand_group2(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_2 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                
                wmatrix1_temp{i} = rerand_corr_1;
                wmatrix2_temp{i} = rerand_corr_2;               
                diff_wmatrix_temp{i} = wmatrix1_temp{i}-wmatrix2_temp{i};
                
            end 
            
            wmatrix1     = zeros(1166, 1166, iteration);
            wmatrix2     = zeros(1166, 1166, iteration);
            diff_wmatrix = zeros(1166, 1166, iteration);
            for i = 1 : iteration
                
                wmatrix1(:, :, i) = wmatrix1_temp{i};
                wmatrix2(:, :, i) = wmatrix2_temp{i};
                diff_wmatrix(:, :, i) = diff_wmatrix_temp{i};
                
            end
       
            real_wmatrix1 = averaged_mat_cont_micro(:, :, s);
            real_wmatrix2 = averaged_mat_HET_micro(:, :, s);
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 :micro_AAL_index_iter
                for j = 1 : micro_AAL_index_iter
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, :)-wmatrix2(i, j, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, :)-wmatrix2(i, j, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
               
        diff_real_wmatrix_cont_vs_HET_micro = diff_real_wmatrix;
        pmap_mat_cont_vs_HET_micro = pmap_mat;
        
        %% 3) control vs. PMG
        type_matrix = 3;
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = mat_cont_micro;
            group2 = mat_FCD_Type_II_micro;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = mat_cont_micro;
            group2 = mat_HET_micro;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = mat_cont_micro;
            group2 = mat_PMG_micro;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
                
        diff_real_wmatrix = zeros(micro_AAL_index_iter,micro_AAL_index_iter,length(myspar));
        pmap_mat = zeros(micro_AAL_index_iter,micro_AAL_index_iter,length(myspar));;
                    
        rp_set = zeros(iteration, N_tgroup);
        for i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp_set(i, :) = randperm(N_tgroup);
            
        end
        
        for s = 1 : length(myspar)
            
            wmatrix1_temp = cell(iteration, 1);
            wmatrix2_temp = cell(iteration, 1);
            diff_wmatrix_temp = cell(iteration, 1);
            
            parfor i = 1 : iteration
                
                fprintf('\n s = %d & i = %d\n', s, i);
                rp = rp_set(i, :);
                
                rerand_group1 = single(tgroup_Reg(:,:,rp(1:N1)));
                rerand_group2 = single(tgroup_Reg(:,:,rp((N1+1):N_tgroup)));
                      
                mat_temp  = rerand_group1*0;
                mat_temp2 = rerand_group1*0;
                for j = 1 : N1
                    
                    temp = threshold_proportional(rerand_group1(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_1 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                mat_temp  = rerand_group2*0;
                mat_temp2 = rerand_group2*0;
                for j = 1 : N2
                    
                    temp = threshold_proportional(rerand_group2(:, :, j), myspar(s));
                    mat_temp(:, :, j) = temp > 0;
                    mat_temp2(:, :, j) = temp;
                    
                end
                rerand_corr_2 = mean(mat_temp2, 3).*(mean(mat_temp, 3)>mat_thres);
                
                rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
                rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
                
                wmatrix1_temp{i} = rerand_corr_1;
                wmatrix2_temp{i} = rerand_corr_2;               
                diff_wmatrix_temp{i} = wmatrix1_temp{i}-wmatrix2_temp{i};
                
            end 
            
            wmatrix1     = zeros(1166, 1166, iteration);
            wmatrix2     = zeros(1166, 1166, iteration);
            diff_wmatrix = zeros(1166, 1166, iteration);
            for i = 1 : iteration
                
                wmatrix1(:, :, i) = wmatrix1_temp{i};
                wmatrix2(:, :, i) = wmatrix2_temp{i};
                diff_wmatrix(:, :, i) = diff_wmatrix_temp{i};
                
            end
       
            real_wmatrix1 = averaged_mat_cont_micro(:, :, s);
            real_wmatrix2 = averaged_mat_PMG_micro(:, :, s);
            
            diff_real_wmatrix(:, :, s) = real_wmatrix1 - real_wmatrix2;
            
            for i = 1 :micro_AAL_index_iter
                for j = 1 : micro_AAL_index_iter
                    
                    if(diff_real_wmatrix(i, j, s)>0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) < squeeze(wmatrix1(i, j, :)-wmatrix2(i, j, :)))/iteration;
                    elseif(diff_real_wmatrix(i, j, s)<0)
                        pmap_mat(i, j, s) = sum(diff_real_wmatrix(i, j, s) > squeeze(wmatrix1(i, j, :)-wmatrix2(i, j, :)))/iteration;
                    else
                        pmap_mat(i, j, s) = 1;
                    end
                    
                end
            end
            
        end
               
        diff_real_wmatrix_cont_vs_PMG_micro = diff_real_wmatrix;
        pmap_mat_cont_vs_PMG_micro = pmap_mat;
        
        %% 4) visualize difference
        if(printfigs == 1)
            
            %% variable setup
            load('colormap_asso_mat.mat');
            load('colormap_asso_mat2.mat');
            load('colormap_asso_mat3.mat');
            struct_idx = [ 0 13 17 22 29 33 39 52 56 61 68 72 78 ];
            vertical_line_idx_x = [ struct_idx; struct_idx ];
            vertical_line_idx_y = [ zeros(1, length(struct_idx)); ones(1, length(struct_idx))*78 ];
            
            for i = 1 : size(struct_idx, 2)
                
                if(struct_idx(i) == 0)
                    micro_struct_idx = 0;
                else
                    temp = [];
                    aal_idx = Anatomical_Label_structs_idx((1+struct_idx(i-1)):struct_idx(i));
                    for j = 1 : size(aal_idx, 2)
                        temp = [ temp length(unique(microAAL_surf_data_both(AAL_surf_data_both == aal_idx(j)))) ];
                    end
                    micro_struct_idx = [ micro_struct_idx max(micro_struct_idx)+sum(temp) ];
                end
                
            end
            micro_vertical_line_idx_x = [ micro_struct_idx; micro_struct_idx ];
            micro_vertical_line_idx_y = [ zeros(1, length(micro_struct_idx)); ones(1, length(micro_struct_idx))*max(microAAL_surf_data_both) ];
                        
            spar = 1 %min(find(micro_node_fully_connected));
            FDRQ = 0.05;
            
            S_new = S;
            S_new.coord(1, 1:40962)     = S_new.coord(1, 1:40962)-3;
            S_new.coord(1, 40963:81924) = S_new.coord(1, 40963:81924)+3;
            
            %% node position of each parcel in AAL
            centroid = zeros(1, max(microAAL_surf_data_both));
            for i = 1 : max(microAAL_surf_data_both)
                
                vert_idx = find(microAAL_surf_data_both == i);
                vert_coord = S.coord(:, vert_idx);
                [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                centroid(i) = vert_idx(b);
                
            end
            
            nodePos = S_new.coord(:, centroid);
            
            %% 1) control vs. FCD Type-II
            pmap_mat            = pmap_mat_cont_vs_FCD_TypeII_micro(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_FCD_TypeII_micro(:, :, spar);
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = diff_real_wmatrix.*(pmap_mat<=0.1);
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_FCD_Type_II_fMRI_micro2.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf_micro(diff_real_wmatrix, nodePos, (1:1166), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'node_radius', 1, 'links_visible', 0.01, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_FCD_Type_II_fMRI_micro.png'], 'png', 13); close(gcf);
            
            %% 2) control vs. HET
            pmap_mat            = pmap_mat_cont_vs_HET_micro(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_HET_micro(:, :, spar);
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));
            temp_mat = diff_real_wmatrix.*(pmap_mat<=0.1);
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_HET_fMRI_micro2.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1] =NOELNetVisuAssociateMatrixOnSurf_micro(diff_real_wmatrix, nodePos, (1:1166), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'node_radius', 1, 'links_visible', 0.01, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_HET_fMRI_micro.png'], 'png', 13); close(gcf);
            
            %% 3) control vs. PMG
            pmap_mat            = pmap_mat_cont_vs_PMG_micro(:, :, spar);
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_PMG_micro(:, :, spar);
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(temp(1:round(size(pmap_mat, 1)*size(pmap_mat, 2)*myspar(spar))), FDRQ));
%             temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ));            
            temp_mat = diff_real_wmatrix.*(pmap_mat<=0.1);
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(micro_vertical_line_idx_x, micro_vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(micro_vertical_line_idx_y, micro_vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', micro_struct_idx );
            set(gca, 'XTick', micro_struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_PMG_fMRI_micro3.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1] =NOELNetVisuAssociateMatrixOnSurf_micro(diff_real_wmatrix, nodePos, (1:1166), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'node_radius', 1, 'links_visible', 0.012, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_PMG_fMRI_micro.png'], 'png', 13); close(gcf);
            
        end
        
    end
        
    for original_AAL_new = 1
        
        parpool(24);
        myspar = [ 0.05:0.01:0.4 ];
        iteration = 1000;
        
        type_matrix = 2;
        
        %% 1) control vs. FCD Type-II
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = zadjacent_matrix_control1;
            group2 = zadjacent_matrix_FCD_Type_II1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control2;
            group2 = zadjacent_matrix_FCD_Type_II2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control3;
            group2 = zadjacent_matrix_FCD_Type_II3;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
        
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_corr_1 = mean(tgroup_Reg(:,:,rp(1:N1)), 3);
            rerand_corr_2 = mean(tgroup_Reg(:,:,rp((N1+1):N_tgroup)), 3);
            rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
            rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
            rerand_corr_1 = rerand_corr_1.*(rerand_corr_1>0);
            rerand_corr_2 = rerand_corr_2.*(rerand_corr_2>0);
        
            diff_wmatrix_cell{i}    = rerand_corr_1 - rerand_corr_2;
            
        end        
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx));
                
        diff_real_wmatrix = mean(group1, 3) - mean(group2, 3);
        
        for i = 1 : length(Anatomical_Label_structs_idx)
            for j = 1 : length(Anatomical_Label_structs_idx)
                
                if(diff_real_wmatrix(i, j)>0)
                    pmap_mat(i, j) = sum(diff_real_wmatrix(i, j) < diff_wmatrix(i, j, :))/iteration;
                elseif(diff_real_wmatrix(i, j)<0)
                    pmap_mat(i, j) = sum(diff_real_wmatrix(i, j) > diff_wmatrix(i, j, :))/iteration;
                else
                    pmap_mat(i, j) = 1;
                end
                
            end
        end
         
        diff_real_wmatrix_cont_vs_FCD_TypeII = diff_real_wmatrix;
        pmap_mat_cont_vs_FCD_TypeII = pmap_mat;
        
        %% 2) control vs. HET
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = zadjacent_matrix_control1;
            group2 = zadjacent_matrix_HET1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control2;
            group2 = zadjacent_matrix_HET2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control3;
            group2 = zadjacent_matrix_HET3;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
        
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_corr_1 = mean(tgroup_Reg(:,:,rp(1:N1)), 3);
            rerand_corr_2 = mean(tgroup_Reg(:,:,rp((N1+1):N_tgroup)), 3);
            rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
            rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
            rerand_corr_1 = rerand_corr_1.*(rerand_corr_1>0);
            rerand_corr_2 = rerand_corr_2.*(rerand_corr_2>0);
        
            diff_wmatrix_cell{i}    = rerand_corr_1 - rerand_corr_2;
            
        end        
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx));
                
        diff_real_wmatrix = mean(group1, 3) - mean(group2, 3);
        
        for i = 1 : length(Anatomical_Label_structs_idx)
            for j = 1 : length(Anatomical_Label_structs_idx)
                
                if(diff_real_wmatrix(i, j)>0)
                    pmap_mat(i, j) = sum(diff_real_wmatrix(i, j) < diff_wmatrix(i, j, :))/iteration;
                elseif(diff_real_wmatrix(i, j)<0)
                    pmap_mat(i, j) = sum(diff_real_wmatrix(i, j) > diff_wmatrix(i, j, :))/iteration;
                else
                    pmap_mat(i, j) = 1;
                end
                
            end
        end
        
        diff_real_wmatrix_cont_vs_HET = diff_real_wmatrix;
        pmap_mat_cont_vs_HET = pmap_mat;
        
        %% 3) control vs. PMG
        if(type_matrix == 1)        % raw connectivity matrix
            group1 = zadjacent_matrix_control1;
            group2 = zadjacent_matrix_PMG1;
        elseif(type_matrix == 2)    % positive connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control2;
            group2 = zadjacent_matrix_PMG2;
        elseif(type_matrix == 3)    % all connectivity matrix after FDR correction
            group1 = zadjacent_matrix_control3;
            group2 = zadjacent_matrix_PMG3;
        end
        
        group1(isinf(group1)) = 0;
        group2(isinf(group2)) = 0;
        
        N1 = size(group1,3);
        N2 = size(group2,3);
        N_tgroup = N1 + N2;
        tgroup_Reg = cat(3, group1,group2);
        wmatrix1_cell       = cell(iteration, 1);
        wmatrix2_cell       = cell(iteration, 1);
        diff_wmatrix_cell   = cell(iteration, 1);
        
        diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
        rp_set = zeros(iteration, N_tgroup);
        
        parfor i = 1 : iteration
            
            fprintf('\n i = %d\n',i);
            rp = randperm(N_tgroup);
            
            rerand_corr_1 = mean(tgroup_Reg(:,:,rp(1:N1)), 3);
            rerand_corr_2 = mean(tgroup_Reg(:,:,rp((N1+1):N_tgroup)), 3);
            rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
            rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
            rerand_corr_1 = rerand_corr_1.*(rerand_corr_1>0);
            rerand_corr_2 = rerand_corr_2.*(rerand_corr_2>0);
        
            diff_wmatrix_cell{i}    = rerand_corr_1 - rerand_corr_2;
            
        end        
        diff_wmatrix    = reshape(cell2mat(diff_wmatrix_cell), length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),iteration);
        
        pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx));
        diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx));
                
        diff_real_wmatrix = mean(group1, 3) - mean(group2, 3);
        
        for i = 1 : length(Anatomical_Label_structs_idx)
            for j = 1 : length(Anatomical_Label_structs_idx)
                
                if(diff_real_wmatrix(i, j)>0)
                    pmap_mat(i, j) = sum(diff_real_wmatrix(i, j) < diff_wmatrix(i, j, :))/iteration;
                elseif(diff_real_wmatrix(i, j)<0)
                    pmap_mat(i, j) = sum(diff_real_wmatrix(i, j) > diff_wmatrix(i, j, :))/iteration;
                else
                    pmap_mat(i, j) = 1;
                end
                
            end
        end
        
        diff_real_wmatrix_cont_vs_PMG = diff_real_wmatrix;
        pmap_mat_cont_vs_PMG = pmap_mat;
        
        %% 4) visualize difference
        if(printfigs == 1)
            
            %% variable setup
            load('colormap_asso_mat.mat');
            load('colormap_asso_mat2.mat');
            load('colormap_asso_mat3.mat');
            struct_idx = [ 0 13 17 22 29 33 39 52 56 61 68 72 78 ];
            vertical_line_idx_x = [ struct_idx; struct_idx ];
            vertical_line_idx_y = [ zeros(1, length(struct_idx)); ones(1, length(struct_idx))*78 ];
            spar = 7; %min(find(node_fully_connected));
            FDRQ = 0.05;
            cVID = 1;
            threshold_p = 0.025;
            
            S_new = S;
            S_new.coord(1, 1:40962)     = S_new.coord(1, 1:40962)-10;
            S_new.coord(1, 40963:81924) = S_new.coord(1, 40963:81924)+10;
            
            %% node position of each parcel in AAL
            centroid = zeros(1, length(Anatomical_Label_structs_idx));
            for i = 1 : length(Anatomical_Label_structs_idx)
                
                vert_idx = find(AAL_surf_data_both == Anatomical_Label_structs_idx(i));
                vert_coord = S.coord(:, vert_idx);
                [a b] = min(sum((vert_coord - repmat(mean(vert_coord, 2), 1, length(vert_idx))).^2, 1));
                centroid(i) = vert_idx(b);
                
            end
            
            nodePos = S_new.coord(:, centroid);
            
            %% 1) control vs. FCD Type-II
            pmap_mat            = pmap_mat_cont_vs_FCD_TypeII;
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_FCD_TypeII;            
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ, cVID));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-1 1]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_FCD_Type_II_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_FCD_Type_II_fMRI.png'], 'png', 13); close(gcf);
            
            %% 2) control vs. HET
            pmap_mat            = pmap_mat_cont_vs_HET;
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_HET;         
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ, cVID));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-0.3 0.3]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1] = NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.3, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_HET_fMRI.png'], 'png', 13); close(gcf);
            
            %% 3) control vs. PMG
            pmap_mat            = pmap_mat_cont_vs_PMG;
            temp = sort(pmap_mat(:));
            diff_real_wmatrix   = diff_real_wmatrix_cont_vs_PMG;            
            temp_mat = diff_real_wmatrix.*(pmap_mat<=FDR(pmap_mat(:), FDRQ, cVID));
            temp_mat = (triu(temp_mat)'+tril(temp_mat)+tril(temp_mat)'+triu(temp_mat))/2;
            
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(temp_mat, '', [-0.3 0.3]);
            hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1]);
            hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1]);
            set(gca, 'TickLength', [ -0.03 -0.03 ] );
            set(gca, 'YTick', struct_idx );
            set(gca, 'XTick', struct_idx );
            colorbar('off'); axes(ax1); colorbar;
            colormap(colormap_asso_mat3);
            text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
            text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
            text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
            exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
            
            figure;
            [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.4, 'colormap', colormap_asso_mat3);
            exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_PMG_fMRI.png'], 'png', 13); close(gcf);
            
        end
        
    end
    
end
