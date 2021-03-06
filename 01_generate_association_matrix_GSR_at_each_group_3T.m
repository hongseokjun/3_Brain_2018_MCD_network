+clear all
close all

%% variable setup
addpath('/data/noel/noel2/local/matlab/surfstat_chicago');
addpath(genpath('/host/gypsy/local_raid/seokjun/01_project/03_NetworkAnalysis_FCD_NLES/01_analysis/00_toolbox/'));
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

%% read cortical thickness
for read_CT = 1
    
    % control
    for control = 1
        
        names_ct_left_cont = strcat('TLE_', Codes_cont, '_native_rms_rsl_tlink_20mm_left.txt');
        names_ct_right_cont = strcat('TLE_', Codes_cont, '_native_rms_rsl_tlink_20mm_right.txt');
        
        T_control_left = zeros(size(Codes_cont, 1), 81924/2); cont_nofile = [];
        for i = 1 : size(Codes_cont, 1)
            if(isempty(Path_cont{i, 1})) continue; end
            if(exist([Path_cont{i, 1}  names_ct_left_cont{i, :}], 'file'))
                disp([Path_cont{i, 1}  names_ct_left_cont{i, :}]);
            else
                cont_nofile = [cont_nofile, i]; continue;
            end
            T_control_left(i, :) = SurfStatReadData( [{[Path_cont{i, 1}  names_ct_left_cont{i, :}]}] );
        end
        
        T_control_right = zeros(size(Codes_cont, 1), 81924/2); cont_nofile = [];
        for i = 1 : size(Codes_cont, 1)
            if(isempty(Path_cont{i, 1})) continue; end
            if(exist([Path_cont{i, 1}  names_ct_right_cont{i, :}], 'file'))
                disp([Path_cont{i, 1}  names_ct_right_cont{i, :}]);
            else
                cont_nofile = [cont_nofile, i]; continue;
            end
            T_control_right(i, :) = SurfStatReadData( [{[Path_cont{i, 1}  names_ct_right_cont{i, :}]}] );
        end
        
        T_control = zeros(size(Codes_cont, 1), 81924); cont_nofile = [];
        for i = 1 : size(Codes_cont, 1)
            if(isempty(Path_cont{i, 1})) continue; end
            if(exist([Path_cont{i, 1}  names_ct_left_cont{i, :}], 'file') && exist([Path_cont{i, 1}  names_ct_right_cont{i, :}], 'file'))
                disp([Path_cont{i, 1}  names_ct_left_cont{i, :} ' + ' Path_cont{i, 1}  names_ct_right_cont{i, :}]);
            else
                cont_nofile = [cont_nofile, i]; continue;
            end
            T_control(i, :) = SurfStatReadData( [{[Path_cont{i, 1}  names_ct_left_cont{i, :}]}, {[Path_cont{i, 1}  names_ct_right_cont{i, :}]}] );
        end
        
    end
   
    % FCD_Type_II
    for FCD_Type_II = 1
        
        names_ct_left_mcd = strcat('mcd_', Codes_FCD_Type_II, '_native_rms_rsl_tlink_20mm_left.txt');
        names_ct_right_mcd = strcat('mcd_', Codes_FCD_Type_II, '_native_rms_rsl_tlink_20mm_right.txt');
        
        T_FCD_Type_II_whole_left = zeros(size(Codes_FCD_Type_II, 1), 81924/2); FCD_Type_II_nofile = [];
        for i = 1 : size(Codes_FCD_Type_II, 1)
            if(isempty(Path_FCD_Type_II{i, 1})) continue; end
            if(exist([Path_FCD_Type_II{i, 1}  names_ct_left_mcd{i, :}], 'file'))
                disp([Path_FCD_Type_II{i, 1}  names_ct_left_mcd{i, :}]);
            else
                FCD_Type_II_nofile = [FCD_Type_II_nofile, i]; continue;
            end
            T_FCD_Type_II_whole_left(i, :) = SurfStatReadData( {[Path_FCD_Type_II{i, 1}  names_ct_left_mcd{i, :}]} );
        end
        
        T_FCD_Type_II_whole_right = zeros(size(Codes_FCD_Type_II, 1), 81924/2); FCD_Type_II_nofile = [];
        for i = 1 : size(Codes_FCD_Type_II, 1)
            if(isempty(Path_FCD_Type_II{i, 1})) continue; end
            if(exist([Path_FCD_Type_II{i, 1}  names_ct_right_mcd{i, :}], 'file'))
                disp([Path_FCD_Type_II{i, 1}  names_ct_right_mcd{i, :}]);
            else
                FCD_Type_II_nofile = [FCD_Type_II_nofile, i]; continue;
            end
            T_FCD_Type_II_whole_right(i, :) = SurfStatReadData( {[Path_FCD_Type_II{i, 1}  names_ct_right_mcd{i, :}]} );
        end
        
        T_FCD_Type_II_whole = zeros(size(Codes_FCD_Type_II, 1), 81924); FCD_Type_II_nofile = [];
        for i = 1 : size(Codes_FCD_Type_II, 1)
            if(isempty(Path_FCD_Type_II{i, 1})) continue; end
            if(exist([Path_FCD_Type_II{i, 1}  names_ct_left_mcd{i, :}], 'file') && exist([Path_FCD_Type_II{i, 1}  names_ct_right_mcd{i, :}], 'file'))
                disp([Path_FCD_Type_II{i, 1}  names_ct_left_mcd{i, :} ' + ' Path_FCD_Type_II{i, 1}  names_ct_right_mcd{i, :}]);
            else
                FCD_Type_II_nofile = [FCD_Type_II_nofile, i]; continue;
            end
            T_FCD_Type_II_whole(i, :) = SurfStatReadData( [{[Path_FCD_Type_II{i, 1}  names_ct_left_mcd{i, :}]}, {[Path_FCD_Type_II{i, 1}  names_ct_right_mcd{i, :}]}] );
        end
        
    end
    
    % HET
    for HET = 1
        
        names_ct_left_pg  = strcat('mcd_', Codes_HET, '_native_rms_rsl_tlink_20mm_left.txt');
        names_ct_right_pg = strcat('mcd_', Codes_HET, '_native_rms_rsl_tlink_20mm_right.txt');
        
        T_HET_whole_left = zeros(size(Codes_HET, 1), 81924/2); HET_nofile = [];
        for i = 1 : size(Codes_HET, 1)
            if(isempty(Path_HET{i, 1})) continue; end
            if(exist([Path_HET{i, 1}  names_ct_left_pg{i, :}], 'file'))
                disp([Path_HET{i, 1}  names_ct_left_pg{i, :}]);
            else
                HET_nofile = [HET_nofile, i]; continue;
            end
            T_HET_whole_left(i, :) = SurfStatReadData( {[Path_HET{i, 1}  names_ct_left_pg{i, :}]} );
        end
        
        T_HET_whole_right = zeros(size(Codes_HET, 1), 81924/2); HET_nofile = [];
        for i = 1 : size(Codes_HET, 1)
            if(isempty(Path_HET{i, 1})) continue; end
            if(exist([Path_HET{i, 1}  names_ct_right_pg{i, :}], 'file'))
                disp([Path_HET{i, 1}  names_ct_right_pg{i, :}]);
            else
                HET_nofile = [HET_nofile, i]; continue;
            end
            T_HET_whole_right(i, :) = SurfStatReadData( {[Path_HET{i, 1}  names_ct_right_pg{i, :}]} );
        end
        
        T_HET_whole = zeros(size(Codes_HET, 1), 81924); HET_nofile = [];
        for i = 1 : size(Codes_HET, 1)
            if(isempty(Path_HET{i, 1})) continue; end
            if(exist([Path_HET{i, 1}  names_ct_left_pg{i, :}], 'file') && exist([Path_HET{i, 1}  names_ct_right_pg{i, :}], 'file'))
                disp([Path_HET{i, 1}  names_ct_left_pg{i, :} ' + ' Path_HET{i, 1}  names_ct_right_pg{i, :}]);
            else
                HET_nofile = [HET_nofile, i]; continue;
            end
            T_HET_whole(i, :) = SurfStatReadData( [{[Path_HET{i, 1}  names_ct_left_pg{i, :}]}, {[Path_HET{i, 1}  names_ct_right_pg{i, :}]}] );
        end
        
    end
    
    % PMG
    for PMG = 1
        
        names_ct_left_PMG  = strcat('mcd_', Codes_PMG, '_native_rms_rsl_tlink_20mm_left.txt');
        names_ct_right_PMG = strcat('mcd_', Codes_PMG, '_native_rms_rsl_tlink_20mm_right.txt');
        
        T_PMG_whole_left = zeros(size(Codes_PMG, 1), 81924/2); PMG_nofile = [];
        for i = 1 : size(Codes_PMG, 1)
            if(isempty(Path_PMG{i, 1})) continue; end
            if(exist([Path_PMG{i, 1}  names_ct_left_PMG{i, :}], 'file'))
                disp([Path_PMG{i, 1}  names_ct_left_PMG{i, :}]);
            else
                PMG_nofile = [PMG_nofile, i]; continue;
            end
            T_PMG_whole_left(i, :) = SurfStatReadData( {[Path_PMG{i, 1}  names_ct_left_PMG{i, :}]} );
        end
        
        T_PMG_whole_right = zeros(size(Codes_PMG, 1), 81924/2); PMG_nofile = [];
        for i = 1 : size(Codes_PMG, 1)
            if(isempty(Path_PMG{i, 1})) continue; end
            if(exist([Path_PMG{i, 1}  names_ct_right_PMG{i, :}], 'file'))
                disp([Path_PMG{i, 1}  names_ct_right_PMG{i, :}]);
            else
                PMG_nofile = [PMG_nofile, i]; continue;
            end
            T_PMG_whole_right(i, :) = SurfStatReadData( {[Path_PMG{i, 1}  names_ct_right_PMG{i, :}]} );
        end
        
        T_PMG_whole = zeros(size(Codes_PMG, 1), 81924); PMG_nofile = [];
        for i = 1 : size(Codes_PMG, 1)
            if(isempty(Path_PMG{i, 1})) continue; end
            if(exist([Path_PMG{i, 1}  names_ct_left_PMG{i, :}], 'file') && exist([Path_PMG{i, 1}  names_ct_right_PMG{i, :}], 'file'))
                disp([Path_PMG{i, 1}  names_ct_left_PMG{i, :} ' + ' Path_PMG{i, 1}  names_ct_right_PMG{i, :}]);
            else
                PMG_nofile = [PMG_nofile, i]; continue;
            end
            T_PMG_whole(i, :) = SurfStatReadData( [{[Path_PMG{i, 1}  names_ct_left_PMG{i, :}]}, {[Path_PMG{i, 1}  names_ct_right_PMG{i, :}]}] );
        end
        
    end
    
    num_of_cont         = size(T_control, 1);
    num_of_FCD_Type_II  = size(T_FCD_Type_II_whole, 1);
    num_of_HET          = size(T_HET_whole, 1);
    num_of_PMG          = size(T_PMG_whole, 1);
        
    if printfigs == 1
                
        figure; SurfStatView(mean(T_control, 1), S);            SurfStatColLim([1 4.5]);
        exportfigbo(gcf,[OUTPATH 'mean_CT_control_3T.png'], 'png', 13); close(gcf);
        
        figure; SurfStatView(mean(T_FCD_Type_II_whole, 1), S);  SurfStatColLim([1 4.5]);
        exportfigbo(gcf,[OUTPATH 'mean_CT_FCD_Type_II_3T.png'], 'png', 13); close(gcf);
     
        figure; SurfStatView(mean(T_HET_whole, 1), S);          SurfStatColLim([1 4.5]);
        exportfigbo(gcf,[OUTPATH 'mean_CT_HET_3T.png'], 'png', 13); close(gcf);
        
        figure; SurfStatView(mean(T_PMG_whole, 1), S);          SurfStatColLim([1 4.5]);
        exportfigbo(gcf,[OUTPATH 'mean_CT_PMG_3T.png'], 'png', 13); close(gcf);
       
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

%% construct association matrices, correcting for age, gender and mean thickness
for construct_asso_mat = 1
     
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
        
        %% control
        Sex_group = [ Sex_cont; ];
        Age_group = [ Age_cont; ];
        T         = [ T_control; ];
        mean_T = zeros(size(T, 1), 1);
        
        for i = 1 : size(T, 1)
            mean_T(i) = mean(T(i, mask==1));
        end
        
        GENDER = term(Sex_group);
        AGE = term(Age_group);
        model = 1 + AGE + GENDER + term(mean_T);
        
        asso_mat_cont_res_val         = zeros(num_of_cont, length(Anatomical_Label_structs_idx));
        micro_asso_mat_cont_res_val   = zeros(num_of_cont, max(microAAL_surf_data_both));
        
        count = 1;
        count2 = 1;
        for i = Anatomical_Label_structs_idx
            
            if(sum(ismember(excluded_index, i)) )
                continue;
            end
            
            ROI_thick = [ asso_mat_cont_val(:, count); ];
            slm = SurfStatLinMod(ROI_thick, model);
            ROI_res = ROI_thick - slm.X*slm.coef;
            
            asso_mat_cont_res_val(:, count)           = ROI_res;     
            
            microAAL_idx = unique(microAAL_surf_data_both(AAL_surf_data_both == i));
            
            for j = 1 : size(microAAL_idx, 2)
                
                ROI_thick = [ micro_asso_mat_cont_val(:, count2); ];
                slm = SurfStatLinMod(ROI_thick, model);
                ROI_res = ROI_thick - slm.X*slm.coef;
                
                micro_asso_mat_cont_res_val(:, count2) = ROI_res;
                
                count2 = count2 + 1;
                
            end
            count = count + 1;
            
        end
        
        %% FCD Type II
        Sex_group = [ Sex_FCD_Type_II; ];
        Age_group = [ Age_FCD_Type_II; ];
        T         = [ T_FCD_Type_II_whole; ];
        mean_T = zeros(size(T, 1), 1);
        
        for i = 1 : size(T, 1)
            mean_T(i) = mean(T(i, mask==1));
        end
        
        GENDER = term(Sex_group);
        AGE = term(Age_group);
        model = 1 + AGE + GENDER + term(mean_T);        
        
        asso_mat_FCD_Type_II_res_val         = zeros(num_of_FCD_Type_II,  length(Anatomical_Label_structs_idx));
        micro_asso_mat_FCD_Type_II_res_val   = zeros(num_of_FCD_Type_II,   max(microAAL_surf_data_both));
        
        count = 1;
        count2 = 1;
        for i = Anatomical_Label_structs_idx
            
            if(sum(ismember(excluded_index, i)) )
                continue;
            end
            
            ROI_thick = [ asso_mat_FCD_Type_II_val(:, count); ];
            slm = SurfStatLinMod(ROI_thick, model);
            ROI_res = ROI_thick - slm.X*slm.coef;            
            
            asso_mat_FCD_Type_II_res_val(:, count)     = ROI_res;
             microAAL_idx = unique(microAAL_surf_data_both(AAL_surf_data_both == i));
            
            for j = 1 : size(microAAL_idx, 2)
                
                ROI_thick = [ micro_asso_mat_FCD_Type_II_val(:, count2); ];
                slm = SurfStatLinMod(ROI_thick, model);
                ROI_res = ROI_thick - slm.X*slm.coef;
                
                micro_asso_mat_FCD_Type_II_res_val(:, count2) = ROI_res;
                
                count2 = count2 + 1;
                
            end
            count = count + 1;
            
        end
        
        %% HET
        Sex_group = [ Sex_HET; ];
        Age_group = [ Age_HET; ];
        T         = [ T_HET_whole; ];
        mean_T = zeros(size(T, 1), 1);
        
        for i = 1 : size(T, 1)
            mean_T(i) = mean(T(i, mask==1));
        end
        
        GENDER = term(Sex_group);
        AGE = term(Age_group);
        model = 1 + AGE + GENDER + term(mean_T);        
        
        asso_mat_HET_res_val         = zeros(num_of_HET,  length(Anatomical_Label_structs_idx));
        micro_asso_mat_HET_res_val   = zeros(num_of_HET,   max(microAAL_surf_data_both));
        
        count = 1;
        count2 = 1;
        for i = Anatomical_Label_structs_idx
            
            if(sum(ismember(excluded_index, i)) )
                continue;
            end
            
            ROI_thick = [ asso_mat_HET_val(:, count); ];
            slm = SurfStatLinMod(ROI_thick, model);
            ROI_res = ROI_thick - slm.X*slm.coef;            
            
            asso_mat_HET_res_val(:, count)     = ROI_res;
             microAAL_idx = unique(microAAL_surf_data_both(AAL_surf_data_both == i));
            
            for j = 1 : size(microAAL_idx, 2)
                
                ROI_thick = [ micro_asso_mat_HET_val(:, count2); ];
                slm = SurfStatLinMod(ROI_thick, model);
                ROI_res = ROI_thick - slm.X*slm.coef;
                
                micro_asso_mat_HET_res_val(:, count2) = ROI_res;
                
                count2 = count2 + 1;
                
            end
            count = count + 1;
            
        end
        
        %% PMG
        Sex_group = [ Sex_PMG; ];
        Age_group = [ Age_PMG; ];
        T         = [ T_PMG_whole; ];
        mean_T = zeros(size(T, 1), 1);
        
        for i = 1 : size(T, 1)
            mean_T(i) = mean(T(i, mask==1));
        end
        
        GENDER = term(Sex_group);
        AGE = term(Age_group);
        model = 1 + AGE + GENDER + term(mean_T);        
        
        asso_mat_PMG_res_val         = zeros(num_of_PMG,  length(Anatomical_Label_structs_idx));
        micro_asso_mat_PMG_res_val   = zeros(num_of_PMG,   max(microAAL_surf_data_both));
        
        count = 1;
        count2 = 1;
        for i = Anatomical_Label_structs_idx
            
            if(sum(ismember(excluded_index, i)) )
                continue;
            end
            
            ROI_thick = [ asso_mat_PMG_val(:, count); ];
            slm = SurfStatLinMod(ROI_thick, model);
            ROI_res = ROI_thick - slm.X*slm.coef;            
            
            asso_mat_PMG_res_val(:, count)     = ROI_res;
             microAAL_idx = unique(microAAL_surf_data_both(AAL_surf_data_both == i));
            
            for j = 1 : size(microAAL_idx, 2)
                
                ROI_thick = [ micro_asso_mat_PMG_val(:, count2); ];
                slm = SurfStatLinMod(ROI_thick, model);
                ROI_res = ROI_thick - slm.X*slm.coef;
                
                micro_asso_mat_PMG_res_val(:, count2) = ROI_res;
                
                count2 = count2 + 1;
                
            end
            count = count + 1;
            
        end
        
    end
        
    %% construct association matrices
    adjacent_matrix_control             = corr(asso_mat_cont_res_val);    
    adjacent_matrix_FCD_Type_II         = corr(asso_mat_FCD_Type_II_res_val);
    adjacent_matrix_HET                 = corr(asso_mat_HET_res_val);
    adjacent_matrix_PMG                 = corr(asso_mat_PMG_res_val);
    
    zadjacent_matrix_control            = 0.5*(log((1+adjacent_matrix_control)./(1-adjacent_matrix_control)));    
    zadjacent_matrix_FCD_Type_II        = 0.5*(log((1+adjacent_matrix_FCD_Type_II)./(1-adjacent_matrix_FCD_Type_II)));
    zadjacent_matrix_HET                = 0.5*(log((1+adjacent_matrix_HET)./(1-adjacent_matrix_HET)));
    zadjacent_matrix_PMG                = 0.5*(log((1+adjacent_matrix_PMG)./(1-adjacent_matrix_PMG)));
    
    micro_adjacent_matrix_control       = corr(micro_asso_mat_cont_res_val);    
    micro_adjacent_matrix_FCD_Type_II   = corr(micro_asso_mat_FCD_Type_II_res_val);
    micro_adjacent_matrix_HET           = corr(micro_asso_mat_HET_res_val);
    micro_adjacent_matrix_PMG           = corr(micro_asso_mat_PMG_res_val);
    
    micro_zadjacent_matrix_control      = 0.5*(log((1+micro_adjacent_matrix_control)./(1-micro_adjacent_matrix_control)));    
    micro_zadjacent_matrix_FCD_Type_II  = 0.5*(log((1+micro_adjacent_matrix_FCD_Type_II)./(1-micro_adjacent_matrix_FCD_Type_II)));
    micro_zadjacent_matrix_HET          = 0.5*(log((1+micro_adjacent_matrix_HET)./(1-micro_adjacent_matrix_HET)));
    micro_zadjacent_matrix_PMG          = 0.5*(log((1+micro_adjacent_matrix_PMG)./(1-micro_adjacent_matrix_PMG)));    
        
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
            wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_control, myspar(i));
            node_degree_uw_pos_cont(i, :)            = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_cont(i, :)          = strengths_und(double(wmatrix>0));
            
            %% FCD Type-II
            wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_FCD_Type_II, myspar(i));
            node_degree_uw_pos_FCD_Type_II(i, :)     = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_FCD_Type_II(i, :)   = strengths_und(double(wmatrix>0));
            
            %% HET
            wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_HET, myspar(i));
            node_degree_uw_pos_HET(i, :)             = degrees_und(double(wmatrix>0));
            node_strength_uw_pos_HET(i, :)           = strengths_und(double(wmatrix>0));
            
            %% PMG
            wmatrix                                  = threshold_proportional(micro_zadjacent_matrix_PMG, myspar(i));
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
        plotCorrelationMatrix(zadjacent_matrix_control,{},[-1 1]); 
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
        
        f = figure;
        for i = 1 : length(myspar)
            
            wmatrix = threshold_proportional(zadjacent_matrix_control, myspar(i));
            
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_control_' num2str(myspar(i)) '.png'], 'png', 13); 
            
        end
        close(gcf);
        
        % - micro
        f = figure;
        ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off'); 
        ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
        plotCorrelationMatrix(micro_zadjacent_matrix_control,{},[-1 1]); 
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
        plotCorrelationMatrix(zadjacent_matrix_FCD_Type_II,{},[-1 1]); 
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
        
        f = figure;
        for i = 1 : length(myspar)
            
            wmatrix = threshold_proportional(zadjacent_matrix_FCD_Type_II, myspar(i));
            
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_FCD_TypeII_' num2str(myspar(i)) '.png'], 'png', 13); 
            
        end
        close(gcf);
        
        % - micro
        f = figure;
        ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off'); 
        ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
        plotCorrelationMatrix(micro_zadjacent_matrix_FCD_Type_II,{},[-1 1]); 
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
        plotCorrelationMatrix(zadjacent_matrix_HET,{},[-1 1]); 
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
            
            wmatrix = threshold_proportional(zadjacent_matrix_HET, myspar(i));
            
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
        plotCorrelationMatrix(micro_zadjacent_matrix_HET,{},[-1 1]); 
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
        plotCorrelationMatrix(zadjacent_matrix_PMG,{},[-1 1]); 
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
        
        f = figure;
        for i = 1 : length(myspar)
            
            wmatrix = threshold_proportional(zadjacent_matrix_PMG, myspar(i));
            
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
            exportfigbo(gcf,[OUTPATH 'adj_mat_PMG_' num2str(myspar(i)) '.png'], 'png', 13); 
            
        end
        close(gcf);     
        
        % - micro
        f = figure;
        ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off'); 
        ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
        plotCorrelationMatrix(micro_zadjacent_matrix_PMG,{},[-1 1]); 
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
       
    %% save variables
    save_file = 0;
    if(save_file)
        if(exist([ MATDIR '/07_generate_association_matrix_MCD2_3T.mat' ])~=2)
            save([ MATDIR '/07_generate_association_matrix_MCD2_3T_HET_with_NODULE.mat' ], ...
                'adjacent_matrix_control', 'adjacent_matrix_FCD_Type_II', 'adjacent_matrix_HET', 'adjacent_matrix_PMG', ...
                'zadjacent_matrix_control', 'zadjacent_matrix_FCD_Type_II', 'zadjacent_matrix_HET', 'zadjacent_matrix_PMG', ...
                'Anatomical_Label_structs', 'Anatomical_Label_structs_idx', ...
                'asso_mat_cont_res_val', 'asso_mat_FCD_Type_II_res_val', 'asso_mat_HET_res_val', 'asso_mat_PMG_res_val', ...
                'asso_mat_cont_val', 'asso_mat_FCD_Type_II_val', 'asso_mat_HET_val', 'asso_mat_PMG_val');
        end
        
        if(exist([ MATDIR '/01_generate_micro_association_matrix_MCD2_3T.mat' ])~=2)
            save([ MATDIR '/01_generate_micro_association_matrix_MCD2_3T.mat' ], ...
                'micro_adjacent_matrix_control', 'micro_adjacent_matrix_FCD_Type_I', 'micro_adjacent_matrix_FCD_Type_II', 'micro_adjacent_matrix_HET', 'micro_adjacent_matrix_PMG', ...
                'micro_zadjacent_matrix_control', 'micro_zadjacent_matrix_FCD_Type_I', 'micro_zadjacent_matrix_FCD_Type_II', 'micro_zadjacent_matrix_HET', 'micro_zadjacent_matrix_PMG', ...
                'Anatomical_Label_structs', 'Anatomical_Label_structs_idx', ...
                'micro_asso_mat_cont_res_val', 'micro_asso_mat_FCD_Type_I_res_val', 'micro_asso_mat_FCD_Type_II_res_val', 'micro_asso_mat_HET_res_val', 'micro_asso_mat_PMG_res_val', ...
                'micro_asso_mat_cont_val', 'micro_asso_mat_FCD_Type_I_val', 'micro_asso_mat_FCD_Type_II_val', 'micro_asso_mat_HET_val', 'micro_asso_mat_PMG_val');
        end
    end
end

%% statistical test for each connectivity using permutation and FDR
for statistical_test = 1
    
    myspar = [ 0.05:0.01:0.4 ];
    iteration = 1000;
    
    %% 1) control vs. FCD Type-II
    group1 = zadjacent_matrix_control; 
    group2 = zadjacent_matrix_FCD_Type_II;
    group1(isinf(group1)) = 0;
    group2(isinf(group2)) = 0;
    
    N1 = size(group1,1);
    N2 = size(group2,1);
    N_tgroup = N1 + N2;
    tgroup_Reg = [group1;group2];
    wmatrix1 = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    wmatrix2 = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    
    for i = 1 : iteration
        
        fprintf('\n i = %d\n',i);
        rp = randperm(N_tgroup);
        
        rp_set(i, :) = rp;
        
        rerand_corr_1 = corrcoef(tgroup_Reg(rp(1:N1),:));
        rerand_corr_2 = corrcoef(tgroup_Reg(rp((N1+1):N_tgroup),:));
        rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
        rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
        
        for s = 1 : length(myspar)
            
            wmatrix1(:, :, s, i) = threshold_proportional(rerand_corr_1, myspar(s)); 
            wmatrix2(:, :, s, i) = threshold_proportional(rerand_corr_2, myspar(s));
            
            diff_wmatrix(:, :, s, i) = wmatrix1(:, :, s, i)-wmatrix2(:, :, s, i);
            
        end        
        
    end
    
    pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
    diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
    for s = 1 : length(myspar)
        
        fprintf('\n s = %d\n',s);
        real_wmatrix1 = threshold_proportional(group1, myspar(s));
        real_wmatrix2 = threshold_proportional(group2, myspar(s));
        
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
    group1 = zadjacent_matrix_control; 
    group2 = zadjacent_matrix_HET;
    group1(isinf(group1)) = 0;
    group2(isinf(group2)) = 0;
    
    N1 = size(group1,1);
    N2 = size(group2,1);
    N_tgroup = N1 + N2;
    tgroup_Reg = [group1;group2];
    wmatrix1 = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    wmatrix2 = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    
    for i = 1 : iteration
        
        fprintf('\n i = %d\n',i);
        rp = randperm(N_tgroup);
        
        rp_set(i, :) = rp;
        
        rerand_corr_1 = corrcoef(tgroup_Reg(rp(1:N1),:));
        rerand_corr_2 = corrcoef(tgroup_Reg(rp((N1+1):N_tgroup),:));
        rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
        rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
        
        for s = 1 : length(myspar)
            
            wmatrix1(:, :, s, i) = threshold_proportional(rerand_corr_1, myspar(s)); 
            wmatrix2(:, :, s, i) = threshold_proportional(rerand_corr_2, myspar(s));
            
            diff_wmatrix(:, :, s, i) = wmatrix1(:, :, s, i)-wmatrix2(:, :, s, i);
            
        end        
        
    end
    
    pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
    diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
    for s = 1 : length(myspar)
        
        fprintf('\n s = %d\n',s);
        real_wmatrix1 = threshold_proportional(group1, myspar(s));
        real_wmatrix2 = threshold_proportional(group2, myspar(s));
        
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
    group1 = zadjacent_matrix_control; 
    group2 = zadjacent_matrix_PMG;
    group1(isinf(group1)) = 0;
    group2(isinf(group2)) = 0;
    
    N1 = size(group1,1);
    N2 = size(group2,1);
    N_tgroup = N1 + N2;
    tgroup_Reg = [group1;group2];
    wmatrix1 = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    wmatrix2 = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    diff_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar),iteration);
    
    for i = 1 : iteration
        
        fprintf('\n i = %d\n',i);
        rp = randperm(N_tgroup);
        
        rp_set(i, :) = rp;
        
        rerand_corr_1 = corrcoef(tgroup_Reg(rp(1:N1),:));
        rerand_corr_2 = corrcoef(tgroup_Reg(rp((N1+1):N_tgroup),:));
        rerand_corr_1 = rerand_corr_1 - diag(diag(rerand_corr_1));
        rerand_corr_2 = rerand_corr_2 - diag(diag(rerand_corr_2));
        
        for s = 1 : length(myspar)
            
            wmatrix1(:, :, s, i) = threshold_proportional(rerand_corr_1, myspar(s)); 
            wmatrix2(:, :, s, i) = threshold_proportional(rerand_corr_2, myspar(s));
            
            diff_wmatrix(:, :, s, i) = wmatrix1(:, :, s, i)-wmatrix2(:, :, s, i);
            
        end        
        
    end
    
    pmap_mat = ones(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
    diff_real_wmatrix = zeros(length(Anatomical_Label_structs_idx),length(Anatomical_Label_structs_idx),length(myspar));
    for s = 1 : length(myspar)
        
        fprintf('\n s = %d\n',s);
        real_wmatrix1 = threshold_proportional(group1, myspar(s));
        real_wmatrix2 = threshold_proportional(group2, myspar(s));
        
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
    
    %% 5) visualize difference
    if(printfigs == 1)

        %% variable setup
        load('colormap_asso_mat.mat');
        load('colormap_asso_mat2.mat');
        load('colormap_asso_mat3.mat');
        struct_idx = [ 0 13 17 22 29 33 39 52 56 61 68 72 78 ];
        vertical_line_idx_x = [ struct_idx; struct_idx ];
        vertical_line_idx_y = [ zeros(1, length(struct_idx)); ones(1, length(struct_idx))*78 ];
        spar =  7
        FDRQ = 0.05;
        
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
        diff_real_wmatrix   = diff_real_wmatrix_cont_vs_FCD_TypeII(:, :, spar);

        f = figure;
        ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off'); 
        ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
        plotCorrelationMatrix(diff_real_wmatrix.*(pmap_mat<FDR(pmap_mat(:), FDRQ)), '', [-1 1]);
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
        exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_FCD_Type_II_3T_2.png'], 'png', 13); close(gcf);
        
        figure;
        [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.1, 'colormap', colormap_asso_mat3);
        exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_FCD_Type_II.png'], 'png', 13); close(gcf);     
        
        %% 2) control vs. HET
        pmap_mat            = pmap_mat_cont_vs_HET(:, :, spar);
        diff_real_wmatrix   = diff_real_wmatrix_cont_vs_HET(:, :, spar);

        f = figure;
        ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off'); 
        ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
        plotCorrelationMatrix(diff_real_wmatrix.*(pmap_mat<FDR(pmap_mat(:), FDRQ)), '', [-1 1]);
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
        exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_HET_3T_2.png'], 'png', 13); close(gcf);
        
        figure;
        [aa1, cb1] = NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.9, 'colormap', colormap_asso_mat3);
        exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_HET.png'], 'png', 13); close(gcf);     
        
        %% 3) control vs. PMG
        pmap_mat            = pmap_mat_cont_vs_PMG(:, :, spar);
        diff_real_wmatrix   = diff_real_wmatrix_cont_vs_PMG(:, :, spar);

        f = figure;
        ax1 = axes('Position',[0 0.1100 0.92 0.8150],'Visible','off'); 
        ax2 = axes('Position',[0.1129 0.1100 0.6731 0.8150]); axes(ax2);
        plotCorrelationMatrix(diff_real_wmatrix.*(pmap_mat<FDR(pmap_mat(:), FDRQ)), '', [-1 1]);
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
        exportfigbo(gcf,[OUTPATH 'adj_diff_mat_cont_vs_PMG_3T_2.png'], 'png', 13); close(gcf);  
        
        figure;
        [aa1, cb1]=NOELNetVisuAssociateMatrixOnSurf(diff_real_wmatrix, nodePos, (1:78), S_new, 'white', 'link_width', 1, 'link_width_weight', 'no', 'links_visible', 0.1, 'colormap', colormap_asso_mat3);
        exportfigbo(gcf,[OUTPATH 'brain_adj_diff_mat_cont_vs_PMG.png'], 'png', 13); close(gcf);     
        
    end 
    
end
