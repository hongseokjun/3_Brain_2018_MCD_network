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
CIVET_ver1 = 2; % 1 -> new data processed by 1.2.1 | 2 -> old data processed by 1.2.0 (control, FCD Type-II)
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
    
    fid = fopen('DemographicNGrouping_for_controls_Data_3T.csv');
    C = textscan(fid, '%s%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);    
    Path_cont_thickness             = C{1}(:, 2);    
    fclose(fid);
    
    fid = fopen('DemographicNGrouping_for_controls_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%s%f%s%s', 'Delimiter', ',', 'CollectOutput', 1);
    Codes_cont                      = C{1}(:, 1);    
    Path_cont_rsfMRI                = C{1}(:, 2);
    Path_cont_thickness             = strcat(Path_cont_thickness(1, :), repmat(['ver' num2str(CIVET_ver1) '/' ], length(Codes_cont), 1));
    Age_cont                        = C{2};                 % mean/sd:     29.2/7.1
    Sex_cont                        = C{3}(:, 1);           % male/female:   18/15    
    fclose(fid);
    
end

%% read FCD Type-II demographics
for read_FCD_Type_II = 1
    
    fid = fopen('DemographicNGrouping_for_FCD_Type-II_Data_3T.csv');
    C = textscan(fid, '%s%f%s%d%d%d%d%s%d%d%d%d','Delimiter',',','CollectOutput', 1);
    Path_FCD_Type_II_thickness      = C{5};
    fclose(fid);
    
    fid = fopen('DemographicNGrouping_for_FCD_Type-II_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%f%s%d%d%d%d%s%d%d%d%d','Delimiter',',','CollectOutput', 1);
    Codes_FCD_Type_II               = C{1};
    Age_FCD_Type_II                 = C{2};                 % mean/sd: 27.3/8.6
    Sex_FCD_Type_II                 = C{3};                 % male/female: 14/13        
    Left_FCD_Type_II                = C{4}(:, 1);           % 13
    Right_FCD_Type_II               = C{4}(:, 2);           % 14
    Pre_Central_FCD_Type_II         = C{4}(:, 3);           % 17
    Post_Central_FCD_Type_II        = C{4}(:, 4);           % 10    
    Path_FCD_Type_II_rsfMRI         = C{5};
    Path_FCD_Type_II_thickness       = strcat(Path_FCD_Type_II_thickness(1, :), repmat(['ver' num2str(CIVET_ver1) '/' ], length(Codes_FCD_Type_II), 1));
    Duration_FCD_Type_II            = double(C{6}(:, 1));   % mean/sd: 20.9/11.9
    Operation_FCD_Type_II           = double(C{6}(:, 2));   % 27
    Scanner_FCD_Type_II             = double(C{6}(:, 3));
    Lesion_Volume                   = double(C{6}(:, 4));    
    
    LabelPath_FCD = '/host/gypsy/local_raid/seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label_3T/';
    fclose(fid); 
    
end

%% read HET demographics
for read_HET = 1
    
    fid = fopen('DemographicNGrouping_for_HET_Data_3T.csv');
    C = textscan(fid, '%s%f%s%d%s%s%d','Delimiter',',','CollectOutput', 1);  
    Path_HET_thickness = C{5}(:, 1);
    
    fid = fopen('DemographicNGrouping_for_HET_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%f%s%d%s%s%d','Delimiter',',','CollectOutput', 1);
    Inclusion                       = C{4}(:, 1);
    Codes_HET                       = C{1}(Inclusion>0);
    Age_HET                         = C{2}(Inclusion>0);    % mean/sd: 31.4/10.6
    Sex_HET                         = C{3}(Inclusion>0);    % male/female: 8/6
    Left_HET                        = C{4}(Inclusion>0, 1) == 1; 
    Right_HET                       = C{4}(Inclusion>0, 1) == 2;
    Bilateral_HET                   = C{4}(Inclusion>0, 1) == 3;  % Left/Right/Bilateral: 4/4/6  
    Path_HET_rsfMRI                 = C{5}(Inclusion>0, 1);
    Path_HET_thickness              = strcat(Path_HET_thickness(1, :), repmat(['ver' num2str(CIVET_ver1) '/' ], length(Codes_HET), 1));
    Subtype_HET                     = C{5}(Inclusion>0, 2); % PVNH/SUBC: 13/1
    res_type_HET                    = C{6}(Inclusion>0);    % 1: low res (4x4x4), 2: high res (2x2x4)
    
    LabelPath_HET = '/host/gypsy/local_raid/seokjun/01_project/02_Morphometry_MRI_pos_neg/02_Org_data/03_lesion_label/';
    fclose(fid);
    
end

%% read PMG demographics
for read_PMG = 1
    
    fid = fopen('DemographicNGrouping_for_PMG_Data_3T.csv');
    C = textscan(fid, '%s%f%s%d%s%s','Delimiter',',','CollectOutput', 1); 
    Path_PMG_thickness = C{5}(:, 1);
    
    fid = fopen('DemographicNGrouping_for_PMG_Data_3T_fMRI.csv');
    C = textscan(fid, '%s%f%s%d%s%s%d','Delimiter',',','CollectOutput', 1);
    Inclusion                       = C{4}(:, 1);
    Codes_PMG                       = C{1}(Inclusion>0);
    Age_PMG                         = C{2}(Inclusion>0);     % mean/sd: 29.5/11.4
    Sex_PMG                         = C{3}(Inclusion>0);     % male/female: 6/5      
    Left_PMG                        = C{4}(Inclusion>0, 1) == 1; 
    Right_PMG                       = C{4}(Inclusion>0, 1) == 2;
    Bilateral_PMG                   = C{4}(Inclusion>0, 1) == 3;  % Left/Right/Bilateral: 2/1/8  
    Path_PMG_rsfMRI                 = C{5}(Inclusion>0, 1);
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
    % Control      32      29.2/7.1      18/15         N/A               N/A            
    % FCD Type-II  27      27.3/8.6      14/13       15/15/0          20.9/11.9   
    % HET          14      31.4/10.6      8/6        6/10/11          Not yet   
    % PMG          11      29.5/11.4      6/5        6/5/10           Not yet   
    
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

%% read files
for read_measures = 1
    
    rois = AAL_surf_data_both;
        
    % control
    for control = 1
        
        PREFIX = 'TLE';
        
        for thickness = 1
            
            names_ct_left_cont = strcat('TLE_', Codes_cont, '_native_rms_rsl_tlink_20mm_left.txt');
            names_ct_right_cont = strcat('TLE_', Codes_cont, '_native_rms_rsl_tlink_20mm_right.txt');
            
            T_control_left = zeros(size(Codes_cont, 1), 81924/2); cont_nofile = [];
            for i = 1 : size(Codes_cont, 1)
                if(isempty(Path_cont_thickness{i, 1})) continue; end
                if(exist([Path_cont_thickness{i, 1}  names_ct_left_cont{i, :}], 'file'))
                    disp([Path_cont_thickness{i, 1}  names_ct_left_cont{i, :}]);
                else
                    cont_nofile = [cont_nofile, i]; continue;
                end
                T_control_left(i, :) = SurfStatReadData( [{[Path_cont_thickness{i, 1}  names_ct_left_cont{i, :}]}] );
            end
            
            T_control_right = zeros(size(Codes_cont, 1), 81924/2); cont_nofile = [];
            for i = 1 : size(Codes_cont, 1)
                if(isempty(Path_cont_thickness{i, 1})) continue; end
                if(exist([Path_cont_thickness{i, 1}  names_ct_right_cont{i, :}], 'file'))
                    disp([Path_cont_thickness{i, 1}  names_ct_right_cont{i, :}]);
                else
                    cont_nofile = [cont_nofile, i]; continue;
                end
                T_control_right(i, :) = SurfStatReadData( [{[Path_cont_thickness{i, 1}  names_ct_right_cont{i, :}]}] );
            end
            
            T_control = zeros(size(Codes_cont, 1), 81924); cont_nofile = [];
            for i = 1 : size(Codes_cont, 1)
                if(isempty(Path_cont_thickness{i, 1})) continue; end
                if(exist([Path_cont_thickness{i, 1}  names_ct_left_cont{i, :}], 'file') && exist([Path_cont_thickness{i, 1}  names_ct_right_cont{i, :}], 'file'))
                    disp([Path_cont_thickness{i, 1}  names_ct_left_cont{i, :} ' + ' Path_cont_thickness{i, 1}  names_ct_right_cont{i, :}]);
                else
                    cont_nofile = [cont_nofile, i]; continue;
                end
                T_control(i, :) = SurfStatReadData( [{[Path_cont_thickness{i, 1}  names_ct_left_cont{i, :}]}, {[Path_cont_thickness{i, 1}  names_ct_right_cont{i, :}]}] );
            end
            
        end
        
        for rs_fMRI = 1
            ts_parcel_cont = zeros(145, length(Anatomical_Label_structs_idx), num_of_cont);
            for i = 1:length(Codes_cont)
                
                Codes_cont{i}
                temp_path = Path_cont_rsfMRI{i};
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
        
    end
    
    % FCD_Type_II
    for FCD_Type_II = 1
        
        PREFIX = 'mcd';
        
        for thickness = 1
           
            Codes_FCD_Type_II_org = Codes_FCD_Type_II;
            for i = 1 : length(Codes_FCD_Type_II_org)
                Codes_FCD_Type_II{i} = Codes_FCD_Type_II{i}(1:3);
                if(strcmp(Codes_FCD_Type_II{i}, '080'))
                    Codes_FCD_Type_II{i} = '080_1';
                end
            end
            names_ct_left_mcd = strcat('mcd_', Codes_FCD_Type_II, '_native_rms_rsl_tlink_20mm_left.txt');
            names_ct_right_mcd = strcat('mcd_', Codes_FCD_Type_II, '_native_rms_rsl_tlink_20mm_right.txt');
            
            T_FCD_Type_II_whole_left = zeros(size(Codes_FCD_Type_II, 1), 81924/2); FCD_Type_II_nofile = [];
            for i = 1 : size(Codes_FCD_Type_II, 1)
                if(isempty(Path_FCD_Type_II_thickness{i, 1})) continue; end
                if(exist([Path_FCD_Type_II_thickness{i, 1}  names_ct_left_mcd{i, :}], 'file'))
                    disp([Path_FCD_Type_II_thickness{i, 1}  names_ct_left_mcd{i, :}]);
                else
                    FCD_Type_II_nofile = [FCD_Type_II_nofile, i]; continue;
                end
                T_FCD_Type_II_whole_left(i, :) = SurfStatReadData( {[Path_FCD_Type_II_thickness{i, 1}  names_ct_left_mcd{i, :}]} );
            end
            
            T_FCD_Type_II_whole_right = zeros(size(Codes_FCD_Type_II, 1), 81924/2); FCD_Type_II_nofile = [];
            for i = 1 : size(Codes_FCD_Type_II, 1)
                if(isempty(Path_FCD_Type_II_thickness{i, 1})) continue; end
                if(exist([Path_FCD_Type_II_thickness{i, 1}  names_ct_right_mcd{i, :}], 'file'))
                    disp([Path_FCD_Type_II_thickness{i, 1}  names_ct_right_mcd{i, :}]);
                else
                    FCD_Type_II_nofile = [FCD_Type_II_nofile, i]; continue;
                end
                T_FCD_Type_II_whole_right(i, :) = SurfStatReadData( {[Path_FCD_Type_II_thickness{i, 1}  names_ct_right_mcd{i, :}]} );
            end
            
            T_FCD_Type_II_whole = zeros(size(Codes_FCD_Type_II, 1), 81924); FCD_Type_II_nofile = [];
            for i = 1 : size(Codes_FCD_Type_II, 1)
                if(isempty(Path_FCD_Type_II_thickness{i, 1})) continue; end
                if(exist([Path_FCD_Type_II_thickness{i, 1}  names_ct_left_mcd{i, :}], 'file') && exist([Path_FCD_Type_II_thickness{i, 1}  names_ct_right_mcd{i, :}], 'file'))
                    disp([Path_FCD_Type_II_thickness{i, 1}  names_ct_left_mcd{i, :} ' + ' Path_FCD_Type_II_thickness{i, 1}  names_ct_right_mcd{i, :}]);
                else
                    FCD_Type_II_nofile = [FCD_Type_II_nofile, i]; continue;
                end
                T_FCD_Type_II_whole(i, :) = SurfStatReadData( [{[Path_FCD_Type_II_thickness{i, 1}  names_ct_left_mcd{i, :}]}, {[Path_FCD_Type_II_thickness{i, 1}  names_ct_right_mcd{i, :}]}] );
            end
            Codes_FCD_Type_II = Codes_FCD_Type_II_org;
            
        end
        
        for rs_fMRI = 1
            
            ts_parcel_FCD_Type_II = zeros(145, length(Anatomical_Label_structs_idx), num_of_FCD_Type_II);
            for i = 1:length(Codes_FCD_Type_II)
                
                Codes_FCD_Type_II{i}
                temp_path = Path_FCD_Type_II_rsfMRI{i};
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
        
    end
    
    % HET
    for HET = 1
        
        PREFIX = 'mcd';
        
        for thickness = 1
            
            names_ct_left_pg  = strcat('mcd_', Codes_HET, '_native_rms_rsl_tlink_20mm_left.txt');
            names_ct_right_pg = strcat('mcd_', Codes_HET, '_native_rms_rsl_tlink_20mm_right.txt');
            
            T_HET_whole_left = zeros(size(Codes_HET, 1), 81924/2); HET_nofile = [];
            for i = 1 : size(Codes_HET, 1)
                if(isempty(Path_HET_thickness{i, 1})) continue; end
                if(exist([Path_HET_thickness{i, 1}  names_ct_left_pg{i, :}], 'file'))
                    disp([Path_HET_thickness{i, 1}  names_ct_left_pg{i, :}]);
                else
                    HET_nofile = [HET_nofile, i]; continue;
                end
                T_HET_whole_left(i, :) = SurfStatReadData( {[Path_HET_thickness{i, 1}  names_ct_left_pg{i, :}]} );
            end
            
            T_HET_whole_right = zeros(size(Codes_HET, 1), 81924/2); HET_nofile = [];
            for i = 1 : size(Codes_HET, 1)
                if(isempty(Path_HET_thickness{i, 1})) continue; end
                if(exist([Path_HET_thickness{i, 1}  names_ct_right_pg{i, :}], 'file'))
                    disp([Path_HET_thickness{i, 1}  names_ct_right_pg{i, :}]);
                else
                    HET_nofile = [HET_nofile, i]; continue;
                end
                T_HET_whole_right(i, :) = SurfStatReadData( {[Path_HET_thickness{i, 1}  names_ct_right_pg{i, :}]} );
            end
            
            T_HET_whole = zeros(size(Codes_HET, 1), 81924); HET_nofile = [];
            for i = 1 : size(Codes_HET, 1)
                if(isempty(Path_HET_thickness{i, 1})) continue; end
                if(exist([Path_HET_thickness{i, 1}  names_ct_left_pg{i, :}], 'file') && exist([Path_HET_thickness{i, 1}  names_ct_right_pg{i, :}], 'file'))
                    disp([Path_HET_thickness{i, 1}  names_ct_left_pg{i, :} ' + ' Path_HET_thickness{i, 1}  names_ct_right_pg{i, :}]);
                else
                    HET_nofile = [HET_nofile, i]; continue;
                end
                T_HET_whole(i, :) = SurfStatReadData( [{[Path_HET_thickness{i, 1}  names_ct_left_pg{i, :}]}, {[Path_HET_thickness{i, 1}  names_ct_right_pg{i, :}]}] );
            end
            
        end
        
        for rs_fMRI = 1
            
            ts_parcel_HET = zeros(145, length(Anatomical_Label_structs_idx), num_of_HET);
            for i = 1:length(Codes_HET)
                
                Codes_HET{i}
                temp_path = Path_HET_rsfMRI{i};
                temp_code = Codes_HET{i};
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), length(Anatomical_Label_structs_idx));
                
                for roi = 1:length(Anatomical_Label_structs_idx)
                    ts_parcel(:,roi) = mean(ts(:,rois==Anatomical_Label_structs_idx(roi)),2);
                end
                ts_parcel_HET(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
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
        
    end
    
    % PMG
    for PMG = 1
        
        PREFIX = 'mcd';
        
        for thickness = 1
            
            names_ct_left_PMG  = strcat('mcd_', Codes_PMG, '_native_rms_rsl_tlink_20mm_left.txt');
            names_ct_right_PMG = strcat('mcd_', Codes_PMG, '_native_rms_rsl_tlink_20mm_right.txt');
            
            T_PMG_whole_left = zeros(size(Codes_PMG, 1), 81924/2); PMG_nofile = [];
            for i = 1 : size(Codes_PMG, 1)
                if(isempty(Path_PMG_thickness{i, 1})) continue; end
                if(exist([Path_PMG_thickness{i, 1}  names_ct_left_PMG{i, :}], 'file'))
                    disp([Path_PMG_thickness{i, 1}  names_ct_left_PMG{i, :}]);
                else
                    PMG_nofile = [PMG_nofile, i]; continue;
                end
                T_PMG_whole_left(i, :) = SurfStatReadData( {[Path_PMG_thickness{i, 1}  names_ct_left_PMG{i, :}]} );
            end
            
            T_PMG_whole_right = zeros(size(Codes_PMG, 1), 81924/2); PMG_nofile = [];
            for i = 1 : size(Codes_PMG, 1)
                if(isempty(Path_PMG_thickness{i, 1})) continue; end
                if(exist([Path_PMG_thickness{i, 1}  names_ct_right_PMG{i, :}], 'file'))
                    disp([Path_PMG_thickness{i, 1}  names_ct_right_PMG{i, :}]);
                else
                    PMG_nofile = [PMG_nofile, i]; continue;
                end
                T_PMG_whole_right(i, :) = SurfStatReadData( {[Path_PMG_thickness{i, 1}  names_ct_right_PMG{i, :}]} );
            end
            
            T_PMG_whole = zeros(size(Codes_PMG, 1), 81924); PMG_nofile = [];
            for i = 1 : size(Codes_PMG, 1)
                if(isempty(Path_PMG_thickness{i, 1})) continue; end
                if(exist([Path_PMG_thickness{i, 1}  names_ct_left_PMG{i, :}], 'file') && exist([Path_PMG_thickness{i, 1}  names_ct_right_PMG{i, :}], 'file'))
                    disp([Path_PMG_thickness{i, 1}  names_ct_left_PMG{i, :} ' + ' Path_PMG_thickness{i, 1}  names_ct_right_PMG{i, :}]);
                else
                    PMG_nofile = [PMG_nofile, i]; continue;
                end
                T_PMG_whole(i, :) = SurfStatReadData( [{[Path_PMG_thickness{i, 1}  names_ct_left_PMG{i, :}]}, {[Path_PMG_thickness{i, 1}  names_ct_right_PMG{i, :}]}] );
            end
            
        end
        
        for rs_fMRI = 1
            
            ts_parcel_PMG = zeros(145, length(Anatomical_Label_structs_idx), num_of_PMG);
            for i = 1:length(Codes_PMG)
                
                Codes_PMG{i}
                temp_path = Path_PMG_rsfMRI{i};
                temp_code = Codes_PMG{i};
                
                ts_name     = [temp_path '/' PREFIX '_' temp_code '/surf_ts/' PREFIX '_' temp_code '_4D_sm_5_rsl.mgh'];
                ts          = SurfStatReadData_mgh1(ts_name);
                ts_parcel   = zeros(size(ts,1), length(Anatomical_Label_structs_idx));
                
                for roi = 1:length(Anatomical_Label_structs_idx)
                    ts_parcel(:,roi) = mean(ts(:,rois==Anatomical_Label_structs_idx(roi)),2);
                end
                ts_parcel_PMG(1:size(ts_parcel, 1), :, i) = ts_parcel;
                
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
            MODEL = 1 + AGE + SEX + GMS;
            slm = SurfStatLinMod(Y, MODEL);
            temp =  Y -  slm.X*slm.coef;
            ts_parcel_cont_res(i, j, :)         = temp(temp_idx(1, 1):temp_idx(1, 2));
            ts_parcel_FCD_Type_II_res(i, j, :)  = temp(temp_idx(2, 1):temp_idx(2, 2));
            ts_parcel_HET_res(i, j, :)          = temp(temp_idx(3, 1):temp_idx(3, 2));
            ts_parcel_PMG_res(i, j, :)          = temp(temp_idx(4, 1):temp_idx(4, 2));
            
        end
    end
    
end

%% construct association matrices
for construct_asso_mat = 1
    
    for thickness = 1
        
        %% variable setup
        ControlCoding                = cell(size(Codes_cont, 1), 1);                ControlCoding(:, 1)              = {'Control'};
        FCD_Type_II_PatientCoding    = cell(size(Codes_FCD_Type_II, 1), 1);         FCD_Type_II_PatientCoding(:, 1)  = {'FCD_Type_II'};
        HETPatientCoding             = cell(size(Codes_HET, 1), 1);                 HETPatientCoding(:, 1)           = {'HET'};
        PMGPatientCoding             = cell(size(Codes_PMG, 1), 1);                 PMGPatientCoding(:, 1)           = {'PMG'};
        GroupCoding                  = [ ControlCoding; FCD_Type_II_PatientCoding; HETPatientCoding; PMGPatientCoding ];
        
        excluded_index = [ 37, 38, 41, 42, 71, 72, 73, 74, 75, 76, 77, 78 ];
        
        num_of_cont        = length(ControlCoding);
        num_of_FCD_Type_II = length(FCD_Type_II_PatientCoding);
        num_of_HET         = length(HETPatientCoding);
        num_of_PMG         = length(PMGPatientCoding);
        
        asso_mat_cont_val            = zeros(num_of_cont, length(Anatomical_Label_structs_idx));
        asso_mat_FCD_Type_II_val     = zeros(num_of_FCD_Type_II, length(Anatomical_Label_structs_idx));
        asso_mat_HET_val             = zeros(num_of_HET, length(Anatomical_Label_structs_idx));
        asso_mat_PMG_val             = zeros(num_of_PMG, length(Anatomical_Label_structs_idx));
        
        micro_asso_mat_cont_val            = zeros(num_of_cont,         max(microAAL_surf_data_both));
        micro_asso_mat_FCD_Type_II_val     = zeros(num_of_FCD_Type_II,  max(microAAL_surf_data_both));
        micro_asso_mat_HET_val             = zeros(num_of_HET,          max(microAAL_surf_data_both));
        micro_asso_mat_PMG_val             = zeros(num_of_PMG,          max(microAAL_surf_data_both));
        
        count = 1;
        count2 = 1;
        myspar = [ 0.05:0.01:0.4 ];
        
        %% aal-based mean cortical thickness
        for i = Anatomical_Label_structs_idx
            
            i
            
            if(sum(ismember(excluded_index, i)) )
                continue;
            end
            
            asso_mat_cont_val(:, count)          = mean(T_control(:, AAL_surf_data_both == i), 2);
            asso_mat_FCD_Type_II_val(:, count)   = mean(T_FCD_Type_II_whole(:, AAL_surf_data_both == i), 2);
            asso_mat_HET_val(:, count)           = mean(T_HET_whole(:, AAL_surf_data_both == i), 2);
            asso_mat_PMG_val(:, count)           = mean(T_PMG_whole(:, AAL_surf_data_both == i), 2);
            
            microAAL_idx = unique(microAAL_surf_data_both(AAL_surf_data_both == i));
            
            for j = 1 : size(microAAL_idx, 2)
                
                micro_asso_mat_cont_val(:, count2)               = mean(T_control(:, microAAL_surf_data_both == microAAL_idx(j)), 2);
                micro_asso_mat_FCD_Type_II_val(:, count2)        = mean(T_FCD_Type_II_whole(:, microAAL_surf_data_both == microAAL_idx(j)), 2);
                micro_asso_mat_HET_val(:, count2)                = mean(T_HET_whole(:, microAAL_surf_data_both == microAAL_idx(j)), 2);
                micro_asso_mat_PMG_val(:, count2)                = mean(T_PMG_whole(:, microAAL_surf_data_both == microAAL_idx(j)), 2);
                
                count2 = count2 + 1;
                
            end
            count = count + 1;
            
        end
        
        %% regress out the effects of age, gender and mean cortical thickness
        for regress_out_nuisance = 1
            
            Sex_group = [ Sex_cont; Sex_FCD_Type_II; Sex_HET; Sex_PMG ];
            Age_group = [ Age_cont; Age_FCD_Type_II; Age_HET; Age_PMG ];
            T         = [ T_control; T_FCD_Type_II_whole; T_HET_whole; T_PMG_whole ];
            mean_T = zeros(size(T, 1), 1);
            
            for i = 1 : size(T, 1)
                mean_T(i) = mean(T(i, mask==1));
            end
            
            GENDER = term(Sex_group);
            AGE = term(Age_group);
            model = 1 + AGE + GENDER + term(mean_T);
            
            asso_mat_cont_res_val         = zeros(num_of_cont,        length(Anatomical_Label_structs_idx));
            asso_mat_FCD_Type_II_res_val  = zeros(num_of_FCD_Type_II, length(Anatomical_Label_structs_idx));
            asso_mat_HET_res_val          = zeros(num_of_HET,         length(Anatomical_Label_structs_idx));
            asso_mat_PMG_res_val          = zeros(num_of_PMG,         length(Anatomical_Label_structs_idx));
            
            micro_asso_mat_cont_res_val            = zeros(num_of_cont,         max(microAAL_surf_data_both));
            micro_asso_mat_FCD_Type_II_res_val     = zeros(num_of_FCD_Type_II,  max(microAAL_surf_data_both));
            micro_asso_mat_HET_res_val             = zeros(num_of_HET,          max(microAAL_surf_data_both));
            micro_asso_mat_PMG_res_val             = zeros(num_of_PMG,          max(microAAL_surf_data_both));
            
            count  = 1;
            count2 = 1;
            for i = Anatomical_Label_structs_idx
                
                if(sum(ismember(excluded_index, i)) )
                    continue;
                end
                
                ROI_thick = [ asso_mat_cont_val(:, count); asso_mat_FCD_Type_II_val(:, count); asso_mat_HET_val(:, count); asso_mat_PMG_val(:, count) ];
                slm = SurfStatLinMod(ROI_thick, model);
                ROI_res = ROI_thick - slm.X*slm.coef;
                
                asso_mat_cont_res_val(:, count)           = ROI_res(strcmp(GroupCoding, 'Control'));
                asso_mat_FCD_Type_II_res_val(:, count)    = ROI_res(strcmp(GroupCoding, 'FCD_Type_II'));
                asso_mat_HET_res_val(:, count)            = ROI_res(strcmp(GroupCoding, 'HET'));
                asso_mat_PMG_res_val(:, count)            = ROI_res(strcmp(GroupCoding, 'PMG'));
                
                microAAL_idx = unique(microAAL_surf_data_both(AAL_surf_data_both == i));
                
                for j = 1 : size(microAAL_idx, 2)
                    
                    ROI_thick = [ micro_asso_mat_cont_val(:, count2); micro_asso_mat_FCD_Type_II_val(:, count2); micro_asso_mat_HET_val(:, count2); micro_asso_mat_PMG_val(:, count2) ];
                    slm = SurfStatLinMod(ROI_thick, model);
                    ROI_res = ROI_thick - slm.X*slm.coef;
                    
                    micro_asso_mat_cont_res_val(:, count2)               = ROI_res(strcmp(GroupCoding, 'Control'));
                    micro_asso_mat_FCD_Type_II_res_val(:, count2)        = ROI_res(strcmp(GroupCoding, 'FCD_Type_II'));
                    micro_asso_mat_HET_res_val(:, count2)                = ROI_res(strcmp(GroupCoding, 'HET'));
                    micro_asso_mat_PMG_res_val(:, count2)                = ROI_res(strcmp(GroupCoding, 'PMG'));
                    
                    count2 = count2 + 1;
                    
                end
                count = count + 1;
                
            end
            
        end
        
        %% construct association matrices
        load( [ MATDIR '/07_generate_association_matrix_MCD2_3T_HET_with_NODULE.mat' ] )
        
        adjacent_matrix_control_thickness             = corr(asso_mat_cont_res_val);
        adjacent_matrix_FCD_Type_II_thickness         = corr(asso_mat_FCD_Type_II_res_val);
        adjacent_matrix_HET_thickness                 = corr(asso_mat_HET_res_val);
        adjacent_matrix_PMG_thickness                 = corr(asso_mat_PMG_res_val);
        
        zadjacent_matrix_control_thickness            = 0.5*(log((1+adjacent_matrix_control_thickness)./(1-adjacent_matrix_control_thickness)));
        zadjacent_matrix_FCD_Type_II_thickness        = 0.5*(log((1+adjacent_matrix_FCD_Type_II_thickness)./(1-adjacent_matrix_FCD_Type_II_thickness)));
        zadjacent_matrix_HET_thickness                = 0.5*(log((1+adjacent_matrix_HET_thickness)./(1-adjacent_matrix_HET_thickness)));
        zadjacent_matrix_PMG_thickness                = 0.5*(log((1+adjacent_matrix_PMG_thickness)./(1-adjacent_matrix_PMG_thickness)));
        
        micro_adjacent_matrix_control_thickness       = corr(micro_asso_mat_cont_res_val);
        micro_adjacent_matrix_FCD_Type_II_thickness   = corr(micro_asso_mat_FCD_Type_II_res_val);
        micro_adjacent_matrix_HET_thickness           = corr(micro_asso_mat_HET_res_val);
        micro_adjacent_matrix_PMG_thickness           = corr(micro_asso_mat_PMG_res_val);
        
        micro_zadjacent_matrix_control_thickness      = 0.5*(log((1+micro_adjacent_matrix_control_thickness)./(1-micro_adjacent_matrix_control_thickness)));
        micro_zadjacent_matrix_FCD_Type_II_thickness  = 0.5*(log((1+micro_adjacent_matrix_FCD_Type_II_thickness)./(1-micro_adjacent_matrix_FCD_Type_II_thickness)));
        micro_zadjacent_matrix_HET_thickness          = 0.5*(log((1+micro_adjacent_matrix_HET_thickness)./(1-micro_adjacent_matrix_HET_thickness)));
        micro_zadjacent_matrix_PMG_thickness          = 0.5*(log((1+micro_adjacent_matrix_PMG_thickness)./(1-micro_adjacent_matrix_PMG_thickness)));
        
        %% check which density has a fully connected matrix: 0.11
        iternation_num = 100;
        interval_thres = 0.01;
        myspar = 0.05:interval_thres:0.4;
        AAL_index_iter = length(Anatomical_Label_structs_idx);
        for full_conn = 1
            
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
                wmatrix                                  = threshold_proportional(zadjacent_matrix_control_thickness, myspar(i));
                node_degree_uw_pos_cont(i, :)            = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_cont(i, :)          = strengths_und(double(wmatrix>0));
                
                %% FCD Type-II
                wmatrix                                  = threshold_proportional(zadjacent_matrix_FCD_Type_II_thickness, myspar(i));
                node_degree_uw_pos_FCD_Type_II(i, :)     = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_FCD_Type_II(i, :)   = strengths_und(double(wmatrix>0));
                
                %% HET
                wmatrix                                  = threshold_proportional(zadjacent_matrix_HET_thickness, myspar(i));
                node_degree_uw_pos_HET(i, :)             = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_HET(i, :)           = strengths_und(double(wmatrix>0));
                
                %% PMG
                wmatrix                                  = threshold_proportional(zadjacent_matrix_PMG_thickness, myspar(i));
                node_degree_uw_pos_PMG(i, :)             = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_PMG(i, :)           = strengths_und(double(wmatrix>0));
                
                %% check if fully connected
                if((sum(node_degree_uw_pos_cont(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_II(i, :)==0) + sum(node_degree_uw_pos_HET(i, :)==0) + sum(node_degree_uw_pos_PMG(i, :)==0)) > 0)
                    node_fully_connected(i) = 0;
                end
                
            end
            
            node_degree_uw_pos_cont             = zeros(length(myspar), max(microAAL_surf_data_both));
            node_degree_uw_pos_FCD_Type_II      = zeros(length(myspar), max(microAAL_surf_data_both));
            node_degree_uw_pos_HET              = zeros(length(myspar), max(microAAL_surf_data_both));
            node_degree_uw_pos_PMG              = zeros(length(myspar), max(microAAL_surf_data_both));
            
            node_strength_uw_pos_cont             = zeros(length(myspar), max(microAAL_surf_data_both));
            node_strength_uw_pos_FCD_Type_II      = zeros(length(myspar), max(microAAL_surf_data_both));
            node_strength_uw_pos_HET            = zeros(length(myspar), max(microAAL_surf_data_both));
            node_strength_uw_pos_PMG            = zeros(length(myspar), max(microAAL_surf_data_both));
            
            micro_node_fully_connected                = ones(1, length(myspar));
            
            for i = 1 : length(myspar)
                
                %% controls
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_control_thickness, myspar(i));
                node_degree_uw_pos_cont(i, :)            = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_cont(i, :)          = strengths_und(double(wmatrix>0));
                
                %% FCD Type-II
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II_thickness, myspar(i));
                node_degree_uw_pos_FCD_Type_II(i, :)     = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_FCD_Type_II(i, :)   = strengths_und(double(wmatrix>0));
                
                %% HET
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_HET_thickness, myspar(i));
                node_degree_uw_pos_HET(i, :)             = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_HET(i, :)           = strengths_und(double(wmatrix>0));
                
                %% PMG
                wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_PMG_thickness, myspar(i));
                node_degree_uw_pos_PMG(i, :)             = degrees_und(double(wmatrix>0));
                node_strength_uw_pos_PMG(i, :)           = strengths_und(double(wmatrix>0));
                
                %% check if fully connected
                if((sum(node_degree_uw_pos_cont(i, :)==0) + sum(node_degree_uw_pos_FCD_Type_II(i, :)==0) + sum(node_degree_uw_pos_HET(i, :)==0) + sum(node_degree_uw_pos_PMG(i, :)==0)) > 0)
                    micro_node_fully_connected(i) = 0;
                end
                
            end
            
        end
        
        %% visualize the matrix and save files
        if(printfigs == 1)
            
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
            
            %% 1) control
            % - original AAL
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(zadjacent_matrix_control_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_control_3T.png'], 'png', 13); close(gcf);
            
            % - micro
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(micro_zadjacent_matrix_control_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_control_micro_3T.png'], 'png', 13); close(gcf);
            
            %% 2) FCD Type-II
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(zadjacent_matrix_FCD_Type_II_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_FCD_TypeII_3T.png'], 'png', 13); close(gcf);
            
            % - micro
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(micro_zadjacent_matrix_FCD_Type_II_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_FCD_Type_II_micro_3T.png'], 'png', 13); close(gcf);
            
            %% 3) HET
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(adjacent_matrix_HET_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_HET_3T.png'], 'png', 13); close(gcf);
            
            f = figure;
            for i = 1 : length(myspar)
                
                wmatrix = threshold_proportional(zadjacent_matrix_HET_thickness, myspar(i));
                
                ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
                ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
                plotCorrelationMatrix(double(wmatrix>0),{},[0 1]);
                hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [0.8 0.8 0.8]);
                hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [0.8 0.8 0.8]);
                set(gca, 'TickLength', [ -0.03 -0.03 ] );
                set(gca, 'YTick', struct_idx );
                set(gca, 'XTick', struct_idx );
                colorbar('off'); axes(ax1); colorbar;
                colormap(colormap_asso_mat2);
                text(0.09, 0.91, 'F'); text(0.09, 0.805, 'C'); text(0.09, 0.745, 'P'); text(0.09, 0.67, 'O'); text(0.09, 0.6, 'L'); text(0.09, 0.54, 'T');
                text(0.09, 0.91-0.5, 'F'); text(0.09, 0.805-0.5, 'C'); text(0.09, 0.745-0.5, 'P'); text(0.09, 0.67-0.5, 'O'); text(0.09, 0.6-0.5, 'L'); text(0.09, 0.54-0.5, 'T');
                text(0.04, 0.75, 'L'); text(0.04, 0.25, 'R');
                exportfigbo(gcf,[OUTPATH 'adj_mat_HET_' num2str(myspar(i)) '.png'], 'png', 13);
                
            end
            close(gcf);
            
            % - micro
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(micro_zadjacent_matrix_HET_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_HET_micro_3T.png'], 'png', 13); close(gcf);
            
            %% 4) PMG
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(zadjacent_matrix_PMG_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_PMG_3T_4.png'], 'png', 13); close(gcf);
                        
            % - micro
            f = figure;
            ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            plotCorrelationMatrix(micro_zadjacent_matrix_PMG_thickness,{},[-1 1]);
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_PMG_micro_3T.png'], 'png', 13); close(gcf);
            
        end
        
    end
    
    for rs_fMRI = 1
        
        FDR_thres = 0.05;
        
        %% construct association matrices
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
        
        %% check which density has a fully connected matrix: 0.11
        iternation_num = 100;
        interval_thres = 0.01;
        myspar = 0.05:interval_thres:0.4;
        AAL_index_iter = length(Anatomical_Label_structs_idx);
        type_matrix = 2;
        
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
            
        end
        
        %% visualize the matrix and save files
        if(printfigs == 1)
            
            savefigs = 1;
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
            type_matrix = 3;
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
            
            %         f = figure;
            %         ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            %         ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            %         plotCorrelationMatrix(std(mat_cont, 0, 3),{},[0 1]);
            %         hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            %         hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            %         set(gca, 'TickLength', [ -0.03 -0.03 ] );
            %         set(gca, 'YTick', struct_idx );
            %         set(gca, 'XTick', struct_idx );
            %         colorbar('off'); axes(ax1); colorbar;
            
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
            
            %         f = figure;
            %         ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            %         ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            %         plotCorrelationMatrix(std(mat_FCD_Type_II, 0, 3),{},[0 1]);
            %         hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            %         hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            %         set(gca, 'TickLength', [ -0.03 -0.03 ] );
            %         set(gca, 'YTick', struct_idx );
            %         set(gca, 'XTick', struct_idx );
            %         colorbar('off'); axes(ax1); colorbar;
            
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
            
            %         f = figure;
            %         ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            %         ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            %         plotCorrelationMatrix(std(mat_HET, 0, 3),{},[0 1]);
            %         hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            %         hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            %         set(gca, 'TickLength', [ -0.03 -0.03 ] );
            %         set(gca, 'YTick', struct_idx );
            %         set(gca, 'XTick', struct_idx );
            %         colorbar('off'); axes(ax1); colorbar;
            
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
            
            %         f = figure;
            %         ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off');
            %         ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
            %         plotCorrelationMatrix(std(mat_PMG, 0, 3),{},[0 1]);
            %         hold on; plot(vertical_line_idx_x, vertical_line_idx_y, 'Color', [1 1 1], 'LineWidth', 2);
            %         hold on; plot(vertical_line_idx_y, vertical_line_idx_x, 'Color', [1 1 1], 'LineWidth', 2);
            %         set(gca, 'TickLength', [ -0.03 -0.03 ] );
            %         set(gca, 'YTick', struct_idx );
            %         set(gca, 'XTick', struct_idx );
            %         colorbar('off'); axes(ax1); colorbar;
            
        end
        
    end
    
end

%% compute coupling coefficient of structural and functional network
for coupling_s_f_network = 1
    
    OUTPATH = '/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/02_result/95_FCD_PMG_HET/07_s_f_coupling/';
    temp_mat = ones(78, 78);
    iteration_num = 1000;
    visualization = 0;
    
    lobe_wise_str = { 'frontal', 'parietal', 'occipital', 'temporal', 'limbic', 'insula' };
    lobe_wise_L = { [ 1:13 15 ], [ 16:22 ], [ 23:29 ], [ 30:33 ], [ 34:39 ], [ 14 ] };
    lobe_wise_R = { 39 + [ 1:13 15 ], 39 + [ 16:22 ], 39 + [ 23:29 ], 39 + [ 30:33 ], 39 + [ 34:39 ], 39 + [ 14 ] };
    
    for control = 1
        
        zadjacent_matrix_control_rsfMRI = mean(zadjacent_matrix_control1, 3);
        pos_idx = zadjacent_matrix_control_thickness>0 & zadjacent_matrix_control_rsfMRI>0 & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_control_thickness(pos_idx), zadjacent_matrix_control_rsfMRI(pos_idx));
        
        pos_idx = ~isinf(zadjacent_matrix_control_thickness) & ~isinf(zadjacent_matrix_control_rsfMRI) & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_control_thickness(pos_idx), zadjacent_matrix_control_rsfMRI(pos_idx))
        figure; scatter(zadjacent_matrix_control_thickness(pos_idx), zadjacent_matrix_control_rsfMRI(pos_idx), 'o', 'filled', 'MarkerFaceColor', [0.118, 0.565, 1.000], 'MarkerEdgeColor', [0.878, 1.000, 1.000]); xlim([-1.5 1.5]); ylim([-1 2]);
        
        figure; plotCorrelationMatrix(zadjacent_matrix_control_thickness, '', [-1 1]);
        figure; plotCorrelationMatrix(zadjacent_matrix_control_rsfMRI, '', [-1 1]);
        
    end
    
    for FCD_Type_II = 1
        
        zadjacent_matrix_FCD_Type_II_rsfMRI = mean(zadjacent_matrix_FCD_Type_II1, 3);
        pos_idx = zadjacent_matrix_FCD_Type_II_thickness>0 & zadjacent_matrix_FCD_Type_II_rsfMRI>0 & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_FCD_Type_II_thickness(pos_idx), zadjacent_matrix_FCD_Type_II_rsfMRI(pos_idx));
        
        pos_idx = ~isinf(zadjacent_matrix_FCD_Type_II_thickness) & ~isinf(zadjacent_matrix_FCD_Type_II_rsfMRI) & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_FCD_Type_II_thickness(pos_idx), zadjacent_matrix_FCD_Type_II_rsfMRI(pos_idx))
        figure; scatter(zadjacent_matrix_FCD_Type_II_thickness(pos_idx), zadjacent_matrix_FCD_Type_II_rsfMRI(pos_idx), 'o', 'filled', 'MarkerFaceColor', [1.000, 0.647, 0.000], 'MarkerEdgeColor', [1.000, 0.973, 0.863]); xlim([-1.5 1.5]); ylim([-1 2]);
        
        figure; plotCorrelationMatrix(zadjacent_matrix_FCD_Type_II_thickness, '', [-1 1]);
        figure; plotCorrelationMatrix(zadjacent_matrix_FCD_Type_II_rsfMRI, '', [-1 1]);
        
    end
    
    for HET = 1
        
        zadjacent_matrix_HET_rsfMRI = mean(zadjacent_matrix_HET1, 3);
        pos_idx = zadjacent_matrix_HET_thickness>0 & zadjacent_matrix_HET_rsfMRI>0 & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_HET_thickness(pos_idx), zadjacent_matrix_HET_rsfMRI(pos_idx));
        
        pos_idx = ~isinf(zadjacent_matrix_HET_thickness) & ~isinf(zadjacent_matrix_HET_rsfMRI) & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_HET_thickness(pos_idx), zadjacent_matrix_HET_rsfMRI(pos_idx))
        figure; scatter(zadjacent_matrix_HET_thickness(pos_idx), zadjacent_matrix_HET_rsfMRI(pos_idx), 'o', 'filled', 'MarkerFaceColor', [0.294, 0.000, 0.510], 'MarkerEdgeColor', [0.902, 0.902, 0.980]); xlim([-1.5 1.5]); ylim([-1 2]);
        
        figure; plotCorrelationMatrix(zadjacent_matrix_HET_thickness, '', [-1 1]);
        figure; plotCorrelationMatrix(zadjacent_matrix_HET_rsfMRI, '', [-1 1]);
        
    end
    
    for PMG = 1
        
        zadjacent_matrix_PMG_rsfMRI = mean(zadjacent_matrix_PMG1, 3);
        pos_idx = zadjacent_matrix_PMG_thickness>0 & zadjacent_matrix_PMG_rsfMRI>0 & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_PMG_thickness(pos_idx), zadjacent_matrix_PMG_rsfMRI(pos_idx));
        
        pos_idx = ~isinf(zadjacent_matrix_PMG_thickness) & ~isinf(zadjacent_matrix_PMG_rsfMRI) & ~eye(78, 78) & tril(temp_mat);
        [r p] = corr(zadjacent_matrix_PMG_thickness(pos_idx), zadjacent_matrix_PMG_rsfMRI(pos_idx))
        figure; scatter(zadjacent_matrix_PMG_thickness(pos_idx), zadjacent_matrix_PMG_rsfMRI(pos_idx), 'o', 'filled', 'MarkerFaceColor', [0.000, 0.502, 0.000], 'MarkerEdgeColor', [0.596, 0.984, 0.596]); xlim([-1.5 1.5]); ylim([-1 2]);
        
        figure; plotCorrelationMatrix(zadjacent_matrix_PMG_thickness, '', [-1 1]);
        figure; plotCorrelationMatrix(zadjacent_matrix_PMG_rsfMRI, '', [-1 1]);
        
    end
        
    for control_vs_FCD_Type_II = 1
        
        feature1_thickness = zadjacent_matrix_control_thickness;
        feature1_rsfMRI    = mean(zadjacent_matrix_control1, 3);
        feature2_thickness = zadjacent_matrix_FCD_Type_II_thickness;
        feature2_rsfMRI    = mean(zadjacent_matrix_FCD_Type_II1, 3);
        
        pos_idx1 = find(~isinf(feature1_thickness) & ~isinf(feature1_rsfMRI) & ~eye(78, 78) & tril(temp_mat));
        pos_idx2 = find(~isinf(feature2_thickness) & ~isinf(feature2_rsfMRI) & ~eye(78, 78) & tril(temp_mat));
        pos_idx_length = min([ length(pos_idx1) length(pos_idx2) ]);
        
        feature1_thickness_temp = feature1_thickness(pos_idx1);
        feature1_rsfMRI_temp    = feature1_rsfMRI(pos_idx1);
        feature2_thickness_temp = feature2_thickness(pos_idx2);
        feature2_rsfMRI_temp    = feature2_rsfMRI(pos_idx2);
        
        [h,p,ksstat,cv] = kstest(feature1_thickness_temp)
        [h,p,ksstat,cv] = kstest(feature1_rsfMRI_temp)
        [h,p,ksstat,cv] = kstest(feature2_thickness_temp)
        [h,p,ksstat,cv] = kstest(feature2_rsfMRI_temp)
        
        [r1 p] = corr(feature1_thickness_temp, feature1_rsfMRI_temp, 'type', 'spearman')
        [r2 p] = corr(feature2_thickness_temp, feature2_rsfMRI_temp, 'type', 'spearman')
        
        THICK1 = term(feature1_thickness_temp);
        MODEL1 = 1 + THICK1;
        slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
        slm1 = SurfStatT(slm1, feature1_thickness_temp); slm1.t       
        
        THICK2 = term(feature2_thickness_temp);
        MODEL2 = 1 + THICK2;
        slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
        slm2 = SurfStatT(slm2, feature2_thickness_temp); slm2.t
        
        lobe_wise_corr        = cell(length(lobe_wise_L), 1);
        lobe_wise_overallcorr = cell(length(lobe_wise_L), 1);
        for i = 1 : length(lobe_wise_L)
            
            lobe_wise_idx = [ lobe_wise_L{i} lobe_wise_R{i} ];
            
            curr_lobe_r_p1 = [];
            curr_lobe_r_p2 = [];
            for j = 1 : length(lobe_wise_idx)
                
                feature1_thickness_temp2 = feature1_thickness(lobe_wise_idx(j), :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
                feature1_rsfMRI_temp2    = feature1_rsfMRI(lobe_wise_idx(j), :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
                feature2_thickness_temp2 = feature2_thickness(lobe_wise_idx(j), :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
                feature2_rsfMRI_temp2    = feature2_rsfMRI(lobe_wise_idx(j), :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
                
                [r1_lobe p1_lobe] = corr(feature1_thickness_temp2', feature1_rsfMRI_temp2');
                [r2_lobe p2_lobe] = corr(feature2_thickness_temp2', feature2_rsfMRI_temp2');
                
                curr_lobe_r_p1 = [ curr_lobe_r_p1; r1_lobe p1_lobe ];
                curr_lobe_r_p2 = [ curr_lobe_r_p2; r2_lobe p2_lobe ];
                
            end
            
            lobe_wise_corr{i}.g1 = curr_lobe_r_p1(:, 1)';
            lobe_wise_corr{i}.g2 = curr_lobe_r_p2(:, 1)';
            
            [ mean(lobe_wise_corr{i}.g1) mean(lobe_wise_corr{i}.g2) std(lobe_wise_corr{i}.g1) std(lobe_wise_corr{i}.g2) ]
            
            feature1_thickness_temp2 = feature1_thickness(lobe_wise_idx, :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
            feature1_rsfMRI_temp2    = feature1_rsfMRI(lobe_wise_idx, :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
            feature2_thickness_temp2 = feature2_thickness(lobe_wise_idx, :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
            feature2_rsfMRI_temp2    = feature2_rsfMRI(lobe_wise_idx, :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
            
            [r1_olobe p1_olobe] = corr(feature1_thickness_temp2(:), feature1_rsfMRI_temp2(:));
            [r2_olobe p2_olobe] = corr(feature2_thickness_temp2(:), feature2_rsfMRI_temp2(:));
                        
            z = (fisherz(r1_olobe)-fisherz(r2_olobe))/sqrt((1/(length(feature1_thickness_temp2(:))-3))+(1/(length(feature1_thickness_temp2(:)))))
            onet = 1-normcdf(z,0,1);
            twot = 2*onet;
            
            lobe_wise_overallcorr{i} = [ r1_olobe r2_olobe length(feature1_thickness_temp2(:)) length(feature1_thickness_temp2(:)) z onet twot ];
            
        end
        lobe_wise_overallcorr_mat = cell2mat(lobe_wise_overallcorr);
        
        slm1.coef(2) - slm2.coef(2) % 0.0340
        real_r_diff = r1 - r2                     % 0.0255
        z = (fisherz(r1)-fisherz(r2))/sqrt((1/(pos_idx_length-3))+(1/(pos_idx_length-3))) % 1.1248
        onet = 1-normcdf(z,0,1)                                                           % 0.13
        twot = 2*onet                                                                     % 0.26
        
        if(visualization)
                        
            figure; xlim([-1.5 2]); ylim([-0.7 1.5]); hold on;
            pos_idx_length = min([ length(pos_idx1) length(pos_idx2) ]);
            for i = 1 : 100 : pos_idx_length
                i
                range_idx = i : i + 100;
                if(sum(range_idx > pos_idx_length) > 0)
                    range_idx = i : pos_idx_length;
                end
                scatter_patches(feature1_thickness(pos_idx1(range_idx)), feature1_rsfMRI(pos_idx1(range_idx)), 10, 'FaceColor', [0.118, 0.565, 1.000], 'EdgeColor', [0.878, 1.000, 1.000], 'FaceAlpha', 0.5);
                scatter_patches(feature2_thickness(pos_idx2(range_idx)), feature2_rsfMRI(pos_idx2(range_idx)), 10, 'FaceColor', [1.000, 0.647, 0.000], 'EdgeColor', [1.000, 0.973, 0.863], 'FaceAlpha', 0.5);
            end
            
            THICK1 = term(feature1_thickness_temp);
            MODEL1 = 1 + THICK1;
            slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
            slm1 = SurfStatT(slm1, feature1_thickness_temp); slm1.t
            
            x_val = -1.3:0.1:1.7;
            y_val = slm1.coef(2)*x_val + slm1.coef(1);
            
            plot(x_val, y_val, 'Color', [0.118, 0.565, 1.000], 'LineWidth', 2)
            
            THICK2 = term(feature2_thickness_temp);
            MODEL2 = 1 + THICK2;
            slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
            
            x_val = -1.3:0.1:1.7;
            y_val = slm2.coef(2)*x_val + slm2.coef(1);
            plot(x_val, y_val, 'Color', [1.000, 0.647, 0.000], 'LineWidth', 2)
            
            exportfigbo(gcf,[OUTPATH 's_f_coupling_control_FCD_Type_II.png'], 'png', 13); close(gcf);
            
            order_idx = [ 3 2 1 6 5 4 ];
            spider(lobe_wise_overallcorr_mat(order_idx , 1:2), 'structural-functional coupling', repmat([0 0.5],  length(lobe_wise_L), 1), lobe_wise_str(order_idx));
            exportfigbo(gcf,[OUTPATH 's_f_coupling_control_FCD_Type_II_lobewise.png'], 'png', 13); close(gcf);
            
        end
        
        for statistics = 1
            
            group1_thickness   = asso_mat_cont_res_val;
            group1_rsfMRI      = zadjacent_matrix_control1;
            group2_thickness   = asso_mat_FCD_Type_II_res_val;
            group2_rsfMRI      = zadjacent_matrix_FCD_Type_II1;
            
            num_group1 = size(group1_rsfMRI, 3);
            num_group2 = size(group2_rsfMRI, 3);
            N_tgroup   = num_group1 + num_group2;
            num_group_min = min([ num_group1 num_group2 ]);
            
            group_thickness = [ group1_thickness; group2_thickness ];
            group_rsfMRI    = cat(3, group1_rsfMRI, group2_rsfMRI);
            
            rerand_diff_rho  = zeros(1, iteration_num);
            rerand_diff_beta = zeros(1, iteration_num);
            rp_set = zeros(iteration_num, N_tgroup);
                        
            rerand_lobe_wise_overallcorr = cell(length(lobe_wise_L), iteration_num);
            
            for i = 1 : iteration_num
                
                rp = randperm(N_tgroup);
                rp_set(i, :) = rp;
                
                % thickness
                rerand_thickness_group1 = group_thickness(rp(1:num_group_min), :);
                rerand_thickness_group2 = group_thickness(rp(end-num_group_min+1:end), :);                
                rerand_scn_group1 = corr(rerand_thickness_group1);
                rerand_scn_group2 = corr(rerand_thickness_group2);
                
                % functional connectivity
                rerand_rsfMRI_group1 = group_rsfMRI(:, :, rp(1:num_group_min));
                rerand_rsfMRI_group2 = group_rsfMRI(:, :, rp(end-num_group_min+1:end));
                rerand_fn_group1 = mean(rerand_rsfMRI_group1, 3); 
                rerand_fn_group2 = mean(rerand_rsfMRI_group2, 3);
                
                rp1_temp = randperm(size(rerand_fn_group1, 1));
                rp2_temp = randperm(size(rerand_fn_group2, 1));
                rerand_fn_group1 = rerand_fn_group1(rp1_temp, rp1_temp);
                rerand_fn_group2 = rerand_fn_group2(rp2_temp, rp2_temp);
                
                pos_idx1 = find(~isinf(rerand_scn_group1) & ~isinf(rerand_fn_group1) & ~eye(78, 78) & tril(temp_mat));
                pos_idx2 = find(~isinf(rerand_scn_group2) & ~isinf(rerand_fn_group2) & ~eye(78, 78) & tril(temp_mat));
                
                feature1_thickness_temp = rerand_scn_group1(pos_idx1);
                feature1_rsfMRI_temp    = rerand_fn_group1(pos_idx1);
                feature2_thickness_temp = rerand_scn_group2(pos_idx2);
                feature2_rsfMRI_temp    = rerand_fn_group2(pos_idx2);
               
                THICK1 = term(feature1_thickness_temp);
                MODEL1 = 1 + THICK1;
                slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
                
                THICK2 = term(feature2_thickness_temp);
                MODEL2 = 1 + THICK2;
                slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
                
                rerand_diff_beta(i) = slm1.coef(2) - slm2.coef(2);
                [r1 p] = corr(feature1_thickness_temp, feature1_rsfMRI_temp, 'type', 'spearman');
                [r2 p] = corr(feature2_thickness_temp, feature2_rsfMRI_temp, 'type', 'spearman');
                rerand_diff_rho(i) = r1 - r2;
                
                for k = 1 : length(lobe_wise_L)
                    
                    lobe_wise_idx = [ lobe_wise_L{k} lobe_wise_R{k} ];
                                        
                    feature1_thickness_temp2 = rerand_scn_group1(lobe_wise_idx, :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
                    feature1_rsfMRI_temp2    = rerand_fn_group1(lobe_wise_idx, :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
                    feature2_thickness_temp2 = rerand_scn_group2(lobe_wise_idx, :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
                    feature2_rsfMRI_temp2    = rerand_fn_group1(lobe_wise_idx, :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
                    
                    [r1_olobe p1_olobe] = corr(feature1_thickness_temp2(:), feature1_rsfMRI_temp2(:));
                    [r2_olobe p2_olobe] = corr(feature2_thickness_temp2(:), feature2_rsfMRI_temp2(:));
                    
                    z = (fisherz(r1_olobe)-fisherz(r2_olobe))/sqrt((1/(length(feature1_thickness_temp2(:))-3))+(1/(length(feature1_thickness_temp2(:)))));
                    onet = 1-normcdf(z,0,1);
                    twot = 2*onet;
                    
                    rerand_lobe_wise_overallcorr{k, i} = [ r1_olobe r2_olobe length(feature1_thickness_temp2(:)) length(feature1_thickness_temp2(:)) z onet twot ];
                    
                end
                
            end
            
            rerand_lobe_wise_overallcorr_temp = cell(length(lobe_wise_L), 1);
            p_FCD_Type_II_temp = [];
            for k = 1 : length(lobe_wise_L)
                
                rerand_lobe_wise_overallcorr_temp{k} = cell2mat(rerand_lobe_wise_overallcorr(k, :)');
                diff_real = lobe_wise_overallcorr_mat(k, 1) - lobe_wise_overallcorr_mat(k, 2);
                diff_rand = rerand_lobe_wise_overallcorr_temp{k}(:, 1) - rerand_lobe_wise_overallcorr_temp{k}(:, 2);
                if(diff_real>0)
                    p_FCD_Type_II_temp =  [p_FCD_Type_II_temp; sum(diff_rand>diff_real)/iteration_num; ];
                elseif(diff_real<0)
                    p_FCD_Type_II_temp =  [p_FCD_Type_II_temp; sum(diff_rand<diff_real)/iteration_num; ];
                end
                
            end
            lobe_wise_overallcorr_mat = [ lobe_wise_overallcorr_mat p_FCD_Type_II_temp ];
            lobe_wise_overallcorr_mat(:, [6 8])            
            
            if(real_r_diff>0)
                p_FCD_Type_II = sum(rerand_diff_rho>real_r_diff)/iteration_num
            elseif(real_r_diff<0)
                p_FCD_Type_II = sum(rerand_diff_rho<real_r_diff)/iteration_num
            end
            
        end
        
    end
    
    for control_vs_HET = 1
        
        feature1_thickness = zadjacent_matrix_control_thickness;
        feature1_rsfMRI    = mean(zadjacent_matrix_control1, 3);
        feature2_thickness = zadjacent_matrix_HET_thickness;
        feature2_rsfMRI    = mean(zadjacent_matrix_HET1, 3);
        
        pos_idx1 = find(~isinf(feature1_thickness) & ~isinf(feature1_rsfMRI) & ~eye(78, 78) & tril(temp_mat));
        pos_idx2 = find(~isinf(feature2_thickness) & ~isinf(feature2_rsfMRI) & ~eye(78, 78) & tril(temp_mat));
        pos_idx_length = min([ length(pos_idx1) length(pos_idx2) ]);
        
        feature1_thickness_temp = feature1_thickness(pos_idx1);
        feature1_rsfMRI_temp    = feature1_rsfMRI(pos_idx1);
        feature2_thickness_temp = feature2_thickness(pos_idx2);
        feature2_rsfMRI_temp    = feature2_rsfMRI(pos_idx2);
        
        [h,p,ksstat,cv] = kstest(feature1_thickness_temp)
        [h,p,ksstat,cv] = kstest(feature1_rsfMRI_temp)
        [h,p,ksstat,cv] = kstest(feature2_thickness_temp)
        [h,p,ksstat,cv] = kstest(feature2_rsfMRI_temp)
        
        [r1 p] = corr(feature1_thickness_temp, feature1_rsfMRI_temp, 'type', 'spearman')
        [r2 p] = corr(feature2_thickness_temp, feature2_rsfMRI_temp, 'type', 'spearman')
        
        THICK1 = term(feature1_thickness_temp);
        MODEL1 = 1 + THICK1;
        slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
        slm1 = SurfStatT(slm1, feature1_thickness_temp); slm1.t       
        
        THICK2 = term(feature2_thickness_temp);
        MODEL2 = 1 + THICK2;
        slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
        slm2 = SurfStatT(slm2, feature2_thickness_temp); slm2.t
        
        lobe_wise_corr        = cell(length(lobe_wise_L), 1);
        lobe_wise_overallcorr = cell(length(lobe_wise_L), 1);
        for i = 1 : length(lobe_wise_L)
            
            lobe_wise_idx = [ lobe_wise_L{i} lobe_wise_R{i} ];
            
            curr_lobe_r_p1 = [];
            curr_lobe_r_p2 = [];
            for j = 1 : length(lobe_wise_idx)
                
                feature1_thickness_temp2 = feature1_thickness(lobe_wise_idx(j), :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
                feature1_rsfMRI_temp2    = feature1_rsfMRI(lobe_wise_idx(j), :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
                feature2_thickness_temp2 = feature2_thickness(lobe_wise_idx(j), :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
                feature2_rsfMRI_temp2    = feature2_rsfMRI(lobe_wise_idx(j), :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
                
                [r1_lobe p1_lobe] = corr(feature1_thickness_temp2', feature1_rsfMRI_temp2');
                [r2_lobe p2_lobe] = corr(feature2_thickness_temp2', feature2_rsfMRI_temp2');
                
                curr_lobe_r_p1 = [ curr_lobe_r_p1; r1_lobe p1_lobe ];
                curr_lobe_r_p2 = [ curr_lobe_r_p2; r2_lobe p2_lobe ];
                
            end
            
            lobe_wise_corr{i}.g1 = curr_lobe_r_p1(:, 1)';
            lobe_wise_corr{i}.g2 = curr_lobe_r_p2(:, 1)';
            
            [ mean(lobe_wise_corr{i}.g1) mean(lobe_wise_corr{i}.g2) std(lobe_wise_corr{i}.g1) std(lobe_wise_corr{i}.g2) ]
            
            feature1_thickness_temp2 = feature1_thickness(lobe_wise_idx, :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
            feature1_rsfMRI_temp2    = feature1_rsfMRI(lobe_wise_idx, :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
            feature2_thickness_temp2 = feature2_thickness(lobe_wise_idx, :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
            feature2_rsfMRI_temp2    = feature2_rsfMRI(lobe_wise_idx, :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
            
            [r1_olobe p1_olobe] = corr(feature1_thickness_temp2(:), feature1_rsfMRI_temp2(:));
            [r2_olobe p2_olobe] = corr(feature2_thickness_temp2(:), feature2_rsfMRI_temp2(:));
            
            z = (fisherz(r1_olobe)-fisherz(r2_olobe))/sqrt((1/(length(feature1_thickness_temp2(:))-3))+(1/(length(feature1_thickness_temp2(:)))))
            onet = 1-normcdf(z,0,1);
            twot = 2*onet;
            
            lobe_wise_overallcorr{i} = [ r1_olobe r2_olobe length(feature1_thickness_temp2(:)) length(feature1_thickness_temp2(:)) z onet twot ];
            
        end
        lobe_wise_overallcorr_mat = cell2mat(lobe_wise_overallcorr);
        
        slm1.coef(2) - slm2.coef(2) % 0.2057
        real_r_diff = r1 - r2                     % 0.0587
        z = (fisherz(r1)-fisherz(r2))/sqrt((1/(pos_idx_length-3))+(1/(pos_idx_length-3))) % 2.5570
        onet = 1-normcdf(z,0,1)                                                           % 0.0053
        twot = 2*onet                                                                     % 0.0106
        
        if(visualization)
                        
            figure; xlim([-1.5 2]); ylim([-0.7 1.5]); hold on;
            pos_idx_length = min([ length(pos_idx1) length(pos_idx2) ]);
            for i = 1 : 100 : pos_idx_length
                i
                range_idx = i : i + 100;
                if(sum(range_idx > pos_idx_length) > 0)
                    range_idx = i : pos_idx_length;
                end
                scatter_patches(feature1_thickness(pos_idx1(range_idx)), feature1_rsfMRI(pos_idx1(range_idx)), 10, 'FaceColor', [0.118, 0.565, 1.000], 'EdgeColor', [0.878, 1.000, 1.000], 'FaceAlpha', 0.5);
                scatter_patches(feature2_thickness(pos_idx2(range_idx)), feature2_rsfMRI(pos_idx2(range_idx)), 8,  'FaceColor', [0.294, 0.000, 0.510], 'EdgeColor', [0.902, 0.902, 0.980], 'FaceAlpha', 0.5);
            end
            
            THICK1 = term(feature1_thickness_temp);
            MODEL1 = 1 + THICK1;
            slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
            slm1 = SurfStatT(slm1, feature1_thickness_temp); slm1.t
            
            x_val = -1.3:0.1:1.7;
            y_val = slm1.coef(2)*x_val + slm1.coef(1);
            
            plot(x_val, y_val, 'Color', [0.118, 0.565, 1.000], 'LineWidth', 2)
            
            THICK2 = term(feature2_thickness_temp);
            MODEL2 = 1 + THICK2;
            slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
            
            x_val = -1.3:0.1:1.7;
            y_val = slm2.coef(2)*x_val + slm2.coef(1);
            plot(x_val, y_val, 'Color', [0.294, 0.000, 0.510], 'LineWidth', 2)
            
            exportfigbo(gcf,[OUTPATH 's_f_coupling_control_HET.png'], 'png', 13); close(gcf);
            
            order_idx = [ 3 2 1 6 5 4 ];
            spider(lobe_wise_overallcorr_mat(order_idx , 1:2), 'structural-functional coupling', repmat([0 0.5],  length(lobe_wise_L), 1), lobe_wise_str(order_idx));
            exportfigbo(gcf,[OUTPATH 's_f_coupling_control_HET_lobewise.png'], 'png', 13); close(gcf);
            
        end
        
        for statistics = 1
            
            group1_thickness   = asso_mat_cont_res_val;
            group1_rsfMRI      = zadjacent_matrix_control1;
            group2_thickness   = asso_mat_HET_res_val;
            group2_rsfMRI      = zadjacent_matrix_HET1;
            
            num_group1 = size(group1_rsfMRI, 3);
            num_group2 = size(group2_rsfMRI, 3);
            N_tgroup   = num_group1 + num_group2;
            num_group_min = min([ num_group1 num_group2 ]);
            
            group_thickness = [ group1_thickness; group2_thickness ];
            group_rsfMRI    = cat(3, group1_rsfMRI, group2_rsfMRI);
            
            rerand_diff_rho  = zeros(1, iteration_num);
            rerand_diff_beta = zeros(1, iteration_num);
            rp_set = zeros(iteration_num, N_tgroup);
                        
            rerand_lobe_wise_overallcorr = cell(length(lobe_wise_L), iteration_num);
            
            for i = 1 : iteration_num
                
                rp = randperm(N_tgroup);
                rp_set(i, :) = rp;
                
                % thickness
                rerand_thickness_group1 = group_thickness(rp(1:num_group_min), :);
                rerand_thickness_group2 = group_thickness(rp(end-num_group_min+1:end), :);                
                rerand_scn_group1 = corr(rerand_thickness_group1);
                rerand_scn_group2 = corr(rerand_thickness_group2);
                
                % functional connectivity
                rerand_rsfMRI_group1 = group_rsfMRI(:, :, rp(1:num_group_min));
                rerand_rsfMRI_group2 = group_rsfMRI(:, :, rp(end-num_group_min+1:end));
                rerand_fn_group1 = mean(rerand_rsfMRI_group1, 3); 
                rerand_fn_group2 = mean(rerand_rsfMRI_group2, 3);
                
                rp1_temp = randperm(size(rerand_fn_group1, 1));
                rp2_temp = randperm(size(rerand_fn_group2, 1));
                rerand_fn_group1 = rerand_fn_group1(rp1_temp, rp1_temp);
                rerand_fn_group2 = rerand_fn_group2(rp2_temp, rp2_temp);
                
                pos_idx1 = find(~isinf(rerand_scn_group1) & ~isinf(rerand_fn_group1) & ~eye(78, 78) & tril(temp_mat));
                pos_idx2 = find(~isinf(rerand_scn_group2) & ~isinf(rerand_fn_group2) & ~eye(78, 78) & tril(temp_mat));
                
                feature1_thickness_temp = rerand_scn_group1(pos_idx1);
                feature1_rsfMRI_temp    = rerand_fn_group1(pos_idx1);
                feature2_thickness_temp = rerand_scn_group2(pos_idx2);
                feature2_rsfMRI_temp    = rerand_fn_group2(pos_idx2);
               
                THICK1 = term(feature1_thickness_temp);
                MODEL1 = 1 + THICK1;
                slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
                
                THICK2 = term(feature2_thickness_temp);
                MODEL2 = 1 + THICK2;
                slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
                
                rerand_diff_beta(i) = slm1.coef(2) - slm2.coef(2);
                [r1 p] = corr(feature1_thickness_temp, feature1_rsfMRI_temp, 'type', 'spearman');
                [r2 p] = corr(feature2_thickness_temp, feature2_rsfMRI_temp, 'type', 'spearman');
                rerand_diff_rho(i) = r1 - r2;
                
                for k = 1 : length(lobe_wise_L)
                    
                    lobe_wise_idx = [ lobe_wise_L{k} lobe_wise_R{k} ];
                                        
                    feature1_thickness_temp2 = rerand_scn_group1(lobe_wise_idx, :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
                    feature1_rsfMRI_temp2    = rerand_fn_group1(lobe_wise_idx, :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
                    feature2_thickness_temp2 = rerand_scn_group2(lobe_wise_idx, :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
                    feature2_rsfMRI_temp2    = rerand_fn_group1(lobe_wise_idx, :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
                    
                    [r1_olobe p1_olobe] = corr(feature1_thickness_temp2(:), feature1_rsfMRI_temp2(:));
                    [r2_olobe p2_olobe] = corr(feature2_thickness_temp2(:), feature2_rsfMRI_temp2(:));
                    
                    z = (fisherz(r1_olobe)-fisherz(r2_olobe))/sqrt((1/(length(feature1_thickness_temp2(:))-3))+(1/(length(feature1_thickness_temp2(:)))));
                    onet = 1-normcdf(z,0,1);
                    twot = 2*onet;
                    
                    rerand_lobe_wise_overallcorr{k, i} = [ r1_olobe r2_olobe length(feature1_thickness_temp2(:)) length(feature1_thickness_temp2(:)) z onet twot ];
                    
                end
                
            end
            
            rerand_lobe_wise_overallcorr_temp = cell(length(lobe_wise_L), 1);
            p_HET_temp = [];
            for k = 1 : length(lobe_wise_L)
                
                rerand_lobe_wise_overallcorr_temp{k} = cell2mat(rerand_lobe_wise_overallcorr(k, :)');
                diff_real = lobe_wise_overallcorr_mat(k, 1) - lobe_wise_overallcorr_mat(k, 2);
                diff_rand = rerand_lobe_wise_overallcorr_temp{k}(:, 1) - rerand_lobe_wise_overallcorr_temp{k}(:, 2);
                if(diff_real>0)
                    p_HET_temp =  [p_HET_temp; sum(diff_rand>diff_real)/iteration_num; ];
                elseif(diff_real<0)
                    p_HET_temp =  [p_HET_temp; sum(diff_rand<diff_real)/iteration_num; ];
                end
                
            end
            lobe_wise_overallcorr_mat = [ lobe_wise_overallcorr_mat p_HET_temp ];
            lobe_wise_overallcorr_mat(:, [6 8])  
            
            if(real_r_diff>0)
                p_HET = sum(rerand_diff_rho>real_r_diff)/iteration_num
            elseif(real_r_diff<0)
                p_HET = sum(rerand_diff_rho<real_r_diff)/iteration_num
            end
            
        end
        
    end
    
    for control_vs_PMG = 1
        
        feature1_thickness = zadjacent_matrix_control_thickness;
        feature1_rsfMRI    = mean(zadjacent_matrix_control1, 3);
        feature2_thickness = zadjacent_matrix_PMG_thickness;
        feature2_rsfMRI    = mean(zadjacent_matrix_PMG1, 3);
        
        pos_idx1 = find(~isinf(feature1_thickness) & ~isinf(feature1_rsfMRI) & ~eye(78, 78) & tril(temp_mat));
        pos_idx2 = find(~isinf(feature2_thickness) & ~isinf(feature2_rsfMRI) & ~eye(78, 78) & tril(temp_mat));
        pos_idx_length = min([ length(pos_idx1) length(pos_idx2) ]);
        
        feature1_thickness_temp = feature1_thickness(pos_idx1);
        feature1_rsfMRI_temp    = feature1_rsfMRI(pos_idx1);
        feature2_thickness_temp = feature2_thickness(pos_idx2);
        feature2_rsfMRI_temp    = feature2_rsfMRI(pos_idx2);
        
        [h,p,ksstat,cv] = kstest(feature1_thickness_temp)
        [h,p,ksstat,cv] = kstest(feature1_rsfMRI_temp)
        [h,p,ksstat,cv] = kstest(feature2_thickness_temp)
        [h,p,ksstat,cv] = kstest(feature2_rsfMRI_temp)
        
        [r1 p] = corr(feature1_thickness_temp, feature1_rsfMRI_temp, 'type', 'spearman')
        [r2 p] = corr(feature2_thickness_temp, feature2_rsfMRI_temp, 'type', 'spearman')
        
        THICK1 = term(feature1_thickness_temp);
        MODEL1 = 1 + THICK1;
        slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
        slm1 = SurfStatT(slm1, feature1_thickness_temp); slm1.t       
        
        THICK2 = term(feature2_thickness_temp);
        MODEL2 = 1 + THICK2;
        slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
        slm2 = SurfStatT(slm2, feature2_thickness_temp); slm2.t
        
        lobe_wise_corr        = cell(length(lobe_wise_L), 1);
        lobe_wise_overallcorr = cell(length(lobe_wise_L), 1);
        for i = 1 : length(lobe_wise_L)
            
            lobe_wise_idx = [ lobe_wise_L{i} lobe_wise_R{i} ];
            
            curr_lobe_r_p1 = [];
            curr_lobe_r_p2 = [];
            for j = 1 : length(lobe_wise_idx)
                
                feature1_thickness_temp2 = feature1_thickness(lobe_wise_idx(j), :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
                feature1_rsfMRI_temp2    = feature1_rsfMRI(lobe_wise_idx(j), :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
                feature2_thickness_temp2 = feature2_thickness(lobe_wise_idx(j), :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
                feature2_rsfMRI_temp2    = feature2_rsfMRI(lobe_wise_idx(j), :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
                
                [r1_lobe p1_lobe] = corr(feature1_thickness_temp2', feature1_rsfMRI_temp2');
                [r2_lobe p2_lobe] = corr(feature2_thickness_temp2', feature2_rsfMRI_temp2');
                
                curr_lobe_r_p1 = [ curr_lobe_r_p1; r1_lobe p1_lobe ];
                curr_lobe_r_p2 = [ curr_lobe_r_p2; r2_lobe p2_lobe ];
                
            end
            
            lobe_wise_corr{i}.g1 = curr_lobe_r_p1(:, 1)';
            lobe_wise_corr{i}.g2 = curr_lobe_r_p2(:, 1)';
            
            [ mean(lobe_wise_corr{i}.g1) mean(lobe_wise_corr{i}.g2) std(lobe_wise_corr{i}.g1) std(lobe_wise_corr{i}.g2) ]
            
            feature1_thickness_temp2 = feature1_thickness(lobe_wise_idx, :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
            feature1_rsfMRI_temp2    = feature1_rsfMRI(lobe_wise_idx, :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
            feature2_thickness_temp2 = feature2_thickness(lobe_wise_idx, :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
            feature2_rsfMRI_temp2    = feature2_rsfMRI(lobe_wise_idx, :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
            
            [r1_olobe p1_olobe] = corr(feature1_thickness_temp2(:), feature1_rsfMRI_temp2(:));
            [r2_olobe p2_olobe] = corr(feature2_thickness_temp2(:), feature2_rsfMRI_temp2(:));
            
            z = (fisherz(r1_olobe)-fisherz(r2_olobe))/sqrt((1/(length(feature1_thickness_temp2(:))-3))+(1/(length(feature1_thickness_temp2(:)))))
            onet = 1-normcdf(z,0,1);
            twot = 2*onet;
            
            lobe_wise_overallcorr{i} = [ r1_olobe r2_olobe length(feature1_thickness_temp2(:)) length(feature1_thickness_temp2(:)) z onet twot ];
            
        end
        lobe_wise_overallcorr_mat = cell2mat(lobe_wise_overallcorr);
        
        slm1.coef(2) - slm2.coef(2) % 0.2903
        real_r_diff = r1 - r2                     % 0.1019
        z = (fisherz(r1)-fisherz(r2))/sqrt((1/(pos_idx_length-3))+(1/(pos_idx_length-3))) % 4.3721
        onet = 1-normcdf(z,0,1)                                                           % 0.000001
        twot = 2*onet                                                                     % 0.00001
        
        if(visualization)
                        
            figure; xlim([-1.5 2]); ylim([-0.7 1.5]); hold on;
            pos_idx_length = min([ length(pos_idx1) length(pos_idx2) ]);
            for i = 1 : 100 : pos_idx_length
                i
                range_idx = i : i + 100;
                if(sum(range_idx > pos_idx_length) > 0)
                    range_idx = i : pos_idx_length;
                end
                scatter_patches(feature1_thickness(pos_idx1(range_idx)), feature1_rsfMRI(pos_idx1(range_idx)), 8, 'FaceColor', [0.118, 0.565, 1.000], 'EdgeColor', [0.878, 1.000, 1.000], 'FaceAlpha', 0.5);
                scatter_patches(feature2_thickness(pos_idx2(range_idx)), feature2_rsfMRI(pos_idx2(range_idx)), 8, 'FaceColor', [0.000, 0.502, 0.000], 'EdgeColor', [0.596, 0.984, 0.596], 'FaceAlpha', 0.5);
            end
            
            THICK1 = term(feature1_thickness_temp);
            MODEL1 = 1 + THICK1;
            slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
            slm1 = SurfStatT(slm1, feature1_thickness_temp); slm1.t
            
            x_val = -1.3:0.1:1.7;
            y_val = slm1.coef(2)*x_val + slm1.coef(1);
            
            plot(x_val, y_val, 'Color', [0.118, 0.565, 1.000], 'LineWidth', 2)
            
            THICK2 = term(feature2_thickness_temp);
            MODEL2 = 1 + THICK2;
            slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
            
            x_val = -1.3:0.1:1.7;
            y_val = slm2.coef(2)*x_val + slm2.coef(1);
            plot(x_val, y_val, 'Color', [0.000, 0.502, 0.000], 'LineWidth', 2)
            
            exportfigbo(gcf,[OUTPATH 's_f_coupling_control_PMG.png'], 'png', 13); close(gcf);
            
            order_idx = [ 3 2 1 6 5 4 ];
            spider(lobe_wise_overallcorr_mat(order_idx , 1:2), 'structural-functional coupling', repmat([0 0.5],  length(lobe_wise_L), 1), lobe_wise_str(order_idx));
            exportfigbo(gcf,[OUTPATH 's_f_coupling_control_PMG_lobewise.png'], 'png', 13); close(gcf);
            
        end
        
        for statistics = 1
            
            group1_thickness   = asso_mat_cont_res_val;
            group1_rsfMRI      = zadjacent_matrix_control1;
            group2_thickness   = asso_mat_PMG_res_val;
            group2_rsfMRI      = zadjacent_matrix_PMG1;
            
            num_group1 = size(group1_rsfMRI, 3);
            num_group2 = size(group2_rsfMRI, 3);
            N_tgroup   = num_group1 + num_group2;
            num_group_min = min([ num_group1 num_group2 ]);
            
            group_thickness = [ group1_thickness; group2_thickness ];
            group_rsfMRI    = cat(3, group1_rsfMRI, group2_rsfMRI);
            
            rerand_diff_rho  = zeros(1, iteration_num);
            rerand_diff_beta = zeros(1, iteration_num);
            rp_set = zeros(iteration_num, N_tgroup);
                        
            rerand_lobe_wise_overallcorr = cell(length(lobe_wise_L), iteration_num);
            
            for i = 1 : iteration_num
                
                rp = randperm(N_tgroup);
                rp_set(i, :) = rp;
                
                % thickness
                rerand_thickness_group1 = group_thickness(rp(1:num_group_min), :);
                rerand_thickness_group2 = group_thickness(rp(end-num_group_min+1:end), :);                
                rerand_scn_group1 = corr(rerand_thickness_group1);
                rerand_scn_group2 = corr(rerand_thickness_group2);
                
                % functional connectivity
                rerand_rsfMRI_group1 = group_rsfMRI(:, :, rp(1:num_group_min));
                rerand_rsfMRI_group2 = group_rsfMRI(:, :, rp(end-num_group_min+1:end));
                rerand_fn_group1 = mean(rerand_rsfMRI_group1, 3); 
                rerand_fn_group2 = mean(rerand_rsfMRI_group2, 3);
                
                rp1_temp = randperm(size(rerand_fn_group1, 1));
                rp2_temp = randperm(size(rerand_fn_group2, 1));
                rerand_fn_group1 = rerand_fn_group1(rp1_temp, rp1_temp);
                rerand_fn_group2 = rerand_fn_group2(rp2_temp, rp2_temp);
                
                pos_idx1 = find(~isinf(rerand_scn_group1) & ~isinf(rerand_fn_group1) & ~eye(78, 78) & tril(temp_mat));
                pos_idx2 = find(~isinf(rerand_scn_group2) & ~isinf(rerand_fn_group2) & ~eye(78, 78) & tril(temp_mat));
                
                feature1_thickness_temp = rerand_scn_group1(pos_idx1);
                feature1_rsfMRI_temp    = rerand_fn_group1(pos_idx1);
                feature2_thickness_temp = rerand_scn_group2(pos_idx2);
                feature2_rsfMRI_temp    = rerand_fn_group2(pos_idx2);
               
                THICK1 = term(feature1_thickness_temp);
                MODEL1 = 1 + THICK1;
                slm1 = SurfStatLinMod(feature1_rsfMRI_temp, MODEL1);
                
                THICK2 = term(feature2_thickness_temp);
                MODEL2 = 1 + THICK2;
                slm2 = SurfStatLinMod(feature2_rsfMRI_temp, MODEL2);
                
                rerand_diff_beta(i) = slm1.coef(2) - slm2.coef(2);
                [r1 p] = corr(feature1_thickness_temp, feature1_rsfMRI_temp, 'type', 'spearman');
                [r2 p] = corr(feature2_thickness_temp, feature2_rsfMRI_temp, 'type', 'spearman');
                rerand_diff_rho(i) = r1 - r2;
                
                for k = 1 : length(lobe_wise_L)
                    
                    lobe_wise_idx = [ lobe_wise_L{k} lobe_wise_R{k} ];
                                        
                    feature1_thickness_temp2 = rerand_scn_group1(lobe_wise_idx, :); feature1_thickness_temp2(isinf(feature1_thickness_temp2)) = 0;
                    feature1_rsfMRI_temp2    = rerand_fn_group1(lobe_wise_idx, :);    feature1_rsfMRI_temp2(isinf(feature1_rsfMRI_temp2)) = 0;
                    feature2_thickness_temp2 = rerand_scn_group2(lobe_wise_idx, :); feature2_thickness_temp2(isinf(feature2_thickness_temp2)) = 0;
                    feature2_rsfMRI_temp2    = rerand_fn_group1(lobe_wise_idx, :);    feature2_rsfMRI_temp2(isinf(feature2_rsfMRI_temp2)) = 0;
                    
                    [r1_olobe p1_olobe] = corr(feature1_thickness_temp2(:), feature1_rsfMRI_temp2(:));
                    [r2_olobe p2_olobe] = corr(feature2_thickness_temp2(:), feature2_rsfMRI_temp2(:));
                    
                    z = (fisherz(r1_olobe)-fisherz(r2_olobe))/sqrt((1/(length(feature1_thickness_temp2(:))-3))+(1/(length(feature1_thickness_temp2(:)))));
                    onet = 1-normcdf(z,0,1);
                    twot = 2*onet;
                    
                    rerand_lobe_wise_overallcorr{k, i} = [ r1_olobe r2_olobe length(feature1_thickness_temp2(:)) length(feature1_thickness_temp2(:)) z onet twot ];
                    
                end
                
            end
            
            rerand_lobe_wise_overallcorr_temp = cell(length(lobe_wise_L), 1);
            p_PMG_temp = [];
            for k = 1 : length(lobe_wise_L)
                
                rerand_lobe_wise_overallcorr_temp{k} = cell2mat(rerand_lobe_wise_overallcorr(k, :)');
                diff_real = lobe_wise_overallcorr_mat(k, 1) - lobe_wise_overallcorr_mat(k, 2);
                diff_rand = rerand_lobe_wise_overallcorr_temp{k}(:, 1) - rerand_lobe_wise_overallcorr_temp{k}(:, 2);
                if(diff_real>0)
                    p_PMG_temp =  [p_PMG_temp; sum(diff_rand>diff_real)/iteration_num; ];
                elseif(diff_real<0)
                    p_PMG_temp =  [p_PMG_temp; sum(diff_rand<diff_real)/iteration_num; ];
                end
                
            end
            lobe_wise_overallcorr_mat = [ lobe_wise_overallcorr_mat p_PMG_temp ];
            lobe_wise_overallcorr_mat(:, [6 8])  
            
            if(real_r_diff>0)
                p_PMG = sum(rerand_diff_rho>real_r_diff)/iteration_num
            elseif(real_r_diff<0)
                p_PMG = sum(rerand_diff_rho<real_r_diff)/iteration_num
            end
            
        end
        
    end
    
    for aal_lobe_wise_brain = 1
       
        lobe_wise_L_AAL = { [ frontal_lobe_lateral_surface_L_index, frontal_lobe_medial_surface_L_index, frontal_lobe_orbital_surace_L_index, central_region_L_index(1) ], ...
                            [ parietal_lobe_lateral_surface_L_index, parietal_lobe_medial_surface_L_index central_region_L_index(2:3) ], ...
                            [ occipital_lobe_lateral_surface_L_index, occipital_lobe_medial_inferior_surface_L_index ], ...
                            [ temporal_lobe_lateral_surface_L_index ], ...
                            [ limbic_lobe_L_index ], ...
                            [ insula_L_index ] };
        AAL_surf_data_lobewise_left = AAL_surf_data_left * 0;
        
        for i = 1 : length(lobe_wise_L_AAL)
            
            for j = 1 : length(lobe_wise_L_AAL{i})
                lobe_wise_L_AAL{i}(j)
                AAL_surf_data_lobewise_left(AAL_surf_data_left==lobe_wise_L_AAL{i}(j)) = i;
            end
            
        end
        S_left_inf = SurfStatInflate(S_left, 0.2);
        figure; SurfStatView(AAL_surf_data_lobewise_left, S_left_inf); SurfStatColLim([0 6]);
        exportfigbo(gcf,[OUTPATH 'AAL_lobe_wise_on_inf_surf.png'], 'png', 13); close(gcf);
        
    end
    
end
