
%% 1.5T
PREFIX = 'PMG';
aal_data  = SurfStatReadData('aal_both_rsl_final.txt');

ID='234_1';
system(['gunzip -vf /local_raid/noel4_backup/CT' PREFIX '-1.1.5_2NUC/' ID '/surfaces/mcd_' ID '_mid_surface_*_81920.obj.gz']);
surf_data = SurfStatReadSurf({ [ '/local_raid/noel4_backup/CT' PREFIX '-1.1.5_2NUC/' ID '/surfaces/mcd_' ID '_mid_surface_left_81920.obj' ], [ '/local_raid/noel4_backup/CT' PREFIX '-1.1.5_2NUC/' ID '/surfaces/mcd_' ID '_mid_surface_right_81920.obj' ]});
system(['gzip -vf /local_raid/noel4_backup/CT' PREFIX '-1.1.5_2NUC/' ID '/surfaces/mcd_' ID '_mid_surface_*_81920.obj']);

figure; SurfStatView(aal_data, surf_data);

%% 3T
PREFIX = 'PMG';
aal_data  = SurfStatReadData('aal_both_rsl_final.txt');
BASEDIR = [ '/data/noel/noel7/CT' PREFIX '-3T_CIVET-1.2.1/NoelCIVET/' ];
ID='337'; % HET: 198, 201; PMG: 301_3T, 306_3T (x), 319_3T, 335, 337 (HET?), 342

[r, s] = system(['echo $$']); s = s(1:end-1);
TEMPDIR = [ '/tmp/mcd_' ID '_' num2str(s) ];
mkdir(TEMPDIR);
[r, s] = system([ 'sphere_resample_obj ' [ BASEDIR '/' ID '/surfaces/mcd_' ID '_mid_surface_left_81920.obj '  ], [ BASEDIR '/' ID '/transforms/surfreg/mcd_' ID '_left_surfmap.sm '  ], [ TEMPDIR '/mcd_' ID '_mid_surface_left_81920_rsl.obj'  ] ]);
[r, s] = system([ 'sphere_resample_obj ' [ BASEDIR '/' ID '/surfaces/mcd_' ID '_mid_surface_right_81920.obj ' ], [ BASEDIR '/' ID '/transforms/surfreg/mcd_' ID '_right_surfmap.sm ' ], [ TEMPDIR '/mcd_' ID '_mid_surface_right_81920_rsl.obj' ] ]);

surf_data = SurfStatReadSurf({ [ TEMPDIR '/mcd_' ID '_mid_surface_left_81920_rsl.obj'  ], [ TEMPDIR '/mcd_' ID '_mid_surface_right_81920_rsl.obj' ]});

figure; SurfStatView(aal_data, surf_data);
rmdir(TEMPDIR, 's');


