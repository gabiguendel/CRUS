clear; close all;

pathname = 'E:\Exp_data\CRUS\19-5-2016';

experiment = '23';
filenumber = '1';
dead_time = '125';

full_path_FB = [pathname, '\', 'data_FB_', experiment, '_', filenumber,...
    '_', 'dt', dead_time, '.mat'];
full_path_RT = [pathname, '\', 'data_RT_', experiment, '_', filenumber,...
    '_', 'dt', dead_time, '.mat'];
load(full_path_FB);
load(full_path_RT);


%% Switches

period_iteration = false;
convert_fb_to_rt = true;
create_rt_tot = true;
create_fb_tot = true;
save_tot_files = true;

%%

