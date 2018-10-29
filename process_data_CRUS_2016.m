% clear all
pathname = 'E:\Exp_data\CRUS\19-5-2016';

% available experiments:

% 23    20/05/2016, Bb, critical coupling, CRUS measurement high power
% 24    20/05/2016, Bb, critical coupling, CRUS measurement low power

experiment = 23;
input = false;
background = false;

separate = true;
file_num = 1;

deadtime = 125; % in ns

%% Switches

convert_spcx_to_fb = true;
save_fb_file = true;
plot_fb_counts_vs_sweep = true;%false;
save_counts_vs_sweeps = false;
plot_fb_pulse_folding = true;%false;
period_iteration = false;
convert_fb_to_rt = true;
create_rt_tot = true;
create_fb_tot = true;
save_tot_files = true;

%%
switch experiment
% pulse direction: 1=forward, 2=backward
    % backward  -->  orange curve
    % forward   -->  blue curve
        
    case 23
        % CRUS experiment high power
        filename = {
            '20-05-16 09-29-33'
            '20-05-16 10-33-55'
            '20-05-16 11-32-00'
            '20-05-16 11-46-39'
            '20-05-16 12-47-58'
            '20-05-16 13-50-32'
            };
        inputfilename = {
            '20-05-16 09-27-32'
            '20-05-16 10-27-50'
%             '20-05-16 11-21-01'
            '20-05-16 12-43-14'
            };      % off resonance no MOT
        
        width_temp = 128;
        period =  4859.9875 * width_temp/64;			% 4859.98 is for 64 FPGA ticks (in TDC units)
                                                        % for NI PXI-7853R FPGA card
                                                        
%         period_lim = [4859.9859 4859.9892] * width_temp/64;
        period_lim = [4859.985 4859.989] * width_temp/64;
        dperiod = 0.00005;

        delay = 3400;
        
        pulse_delay = delay + [390 1470 2820 3900];
        pulse_hi = [360  760  360  760];
        pulse_dir = [1  1  1  1];

        delper = 0.93;
        
    case 24
        % CRUS experiment low power
        filename = {
            '20-05-16 14-27-50'
            '20-05-16 15-06-51'
            '20-05-16 15-54-06'
            '20-05-16 16-25-46'
            '20-05-16 17-03-57'
            '20-05-16 17-42-44'
            };
        inputfilename = {
            '20-05-16 14-22-07'
            '20-05-16 15-45-06'
            '20-05-16 16-59-41'
            };      % off resonance no MOT
        
        width_temp = 64;
%         width_temp = 128; % double the period for 1-3, 1-5 and 1-7 pulse # postselection methods
        
        % 4859.98 is for 64 FPGA ticks (in TDC units)
%         period =  4859.988 * width_temp/64;
        period =  4859.9877 * width_temp/64;

%         period_lim = [4859.9859 4859.9892] * width_temp/64;
%         period_lim = [4859.985 4859.99] * width_temp/64;
%         period_lim = [4859.986 4859.99] * width_temp/64;
        period_lim = [4859.985 4859.989] * width_temp/64;
        
        dperiod = 0.00005;

        delay = 3400;
        
        pulse_delay = delay + [390 1470 2820 3900];
        pulse_hi = [360  820  360  820];
        pulse_dir = [1  1  1  1];
        
%         pulse_delay = delay + [390 1720 2820 4120 390+period/2 1720+period/2 2840+period/2 4120+period/2];
%         pulse_hi = [360  360  360  360  360  360  360  360];
%         pulse_dir = [1  1  1  1  1  1  1  1];

        delper = 0.93;
        
    otherwise
        error('Unknown experiment number');
end

        fit = [0 1 0 1];
        threshL = repmat(0.01, 1, length(pulse_dir)); % Setting threshold for noise L
        threshR = repmat(0.01, 1, length(pulse_dir)); % Setting threshold for noise R

if input % Switch for traeting the input
    rt_file_tot = fullfile(pathname, sprintf('inputdata_RT_%1.0f_dt%1.0f.mat', experiment, deadtime));
    fb_file_tot = fullfile(pathname, sprintf('inputdata_FB_%1.0f_dt%1.0f.mat', experiment, deadtime));
    filename = inputfilename;
    
    if separate
        rt_file_tot = fullfile(pathname, sprintf('inputdata_RT_%1.0f_%1.0f_dt%1.0f.mat', experiment, file_num, deadtime));
        fb_file_tot = fullfile(pathname, sprintf('inputdata_FB_%1.0f_%1.0f_dt%1.0f.mat', experiment, file_num, deadtime));
    end
    
elseif background % Switch for traeting the background
    rt_file_tot = fullfile(pathname, sprintf('bgdata_RT_%1.0f_dt%1.0f.mat', experiment, deadtime));
    fb_file_tot = fullfile(pathname, sprintf('bgdata_FB_%1.0f_dt%1.0f.mat', experiment, deadtime));
    
    if separate
        rt_file_tot = fullfile(pathname, sprintf('bgdata_RT_%1.0f_%1.0f_dt%1.0f.mat', experiment, file_num, deadtime));
        fb_file_tot = fullfile(pathname, sprintf('bgdata_FB_%1.0f_%1.0f_dt%1.0f.mat', experiment, file_num, deadtime));
    end
    
elseif separate
    rt_file_tot = fullfile(pathname, sprintf('data_RT_%1.0f_%1.0f_dt%1.0f.mat', experiment, file_num, deadtime));
    fb_file_tot = fullfile(pathname, sprintf('data_FB_%1.0f_%1.0f_dt%1.0f.mat', experiment, file_num, deadtime));
    
else
    rt_file_tot = fullfile(pathname, sprintf('data_RT_%1.0f_dt%1.0f.mat', experiment, deadtime));
    fb_file_tot = fullfile(pathname, sprintf('data_FB_%1.0f_dt%1.0f.mat', experiment, deadtime));
    
end

%% Go over all files

num_pulses = length(pulse_dir);

Rtot = {};
Ttot = {};

D_Rtot = {};
D_Ttot = {};

Btot = {};
Ftot = {};

D_Btot = {};
D_Ftot = {};

if separate
    filelist = file_num;
else
    filelist = 1:length(filename);
end

for i = filelist
    
    spcx_file = fullfile(pathname, [filename{i} '.spcx']);
    
    if input
        rt_file = fullfile(pathname, [filename{i} sprintf('_RT_inputdt%1.0f.mat',deadtime)]);
        fb_file = fullfile(pathname, [filename{i} sprintf('_FB_inputdt%1.0f.mat',deadtime)]);
        
    elseif background
        rt_file = fullfile(pathname, [filename{i} sprintf('_RT_bgdt%1.0f.mat',deadtime)]);
        fb_file = fullfile(pathname, [filename{i} sprintf('_FB_bgdt%1.0f.mat',deadtime)]);
        
    else
        rt_file = fullfile(pathname, [filename{i} sprintf('_RT_dt%1.0f.mat',deadtime)]);
        fb_file = fullfile(pathname, [filename{i} sprintf('_FB_dt%1.0f.mat',deadtime)]);
        
    end

    %% Convert the spcx files to FB files
    if convert_spcx_to_fb
        
        [data, tpb] = read_spcx(spcx_file, deadtime, 10);
        
        % total acquisition time in us
        acq_time = 30e3;

        tpb_ns = tpb*1e-6;

%         if input
%             num_sweeps = 1;
%         if background
%             num_sweeps = 10;
%         else
            num_sweeps = size(data,1);
%         end

        F = cell(num_sweeps,1);
        B = cell(num_sweeps,1);
        D_F = cell(num_sweeps,1);
        D_B = cell(num_sweeps,1);

        newF = cell(num_sweeps,1);
        newB = cell(num_sweeps,1);
        newD_F = cell(num_sweeps,1);
        newD_B = cell(num_sweeps,1);

        for j = 1:num_sweeps

            FF{j} = sortrows([data{j,1} ones(length(data{j,1}),1); data{j,2}  repmat(2,length(data{j,2}),1); data{j,3} repmat(3,length(data{j,3}),1); data{j,4} repmat(4,length(data{j,4}),1); data{j,5} repmat(5,length(data{j,5}),1);...
                data{j,6} repmat(6,length(data{j,6}),1); data{j,7} repmat(7,length(data{j,7}),1); data{j,8} repmat(8,length(data{j,8}),1); data{j,9} repmat(9,length(data{j,9}),1); data{j,10} repmat(10,length(data{j,10}),1)]);
            
            %             F{j} = sortrows([data{j,1} ones(length(data{j,1}),1); data{j,2}  repmat(2,length(data{j,2}),1); data{j,3} repmat(3,length(data{j,3}),1); data{j,4} repmat(4,length(data{j,4}),1); data{j,5} repmat(5,length(data{j,5}),1)]);
            %             B{j} = sortrows([data{j,6} repmat(6,length(data{j,6}),1); data{j,7} repmat(7,length(data{j,7}),1); data{j,8} repmat(8,length(data{j,8}),1); data{j,9} repmat(9,length(data{j,9}),1); data{j,10} repmat(10,length(data{j,10}),1)]);
                        
            D_F{j} = FF{j}(:,2);
%             D_B{j} = B{j}(:,2);
            D_B{j} = [];
            
            F{j} = FF{j}(:,1);
%             B{j} = B{j}(:,1);
            
        end
        
        e64 = uint64([]);   % empty array of type uint64

%         if ~background
            
            Fdel = F;
            
            switch filename{i}
                case '20-05-16 13-50-32'
                    F{180} = e64;
                case '20-05-16 17-42-44'
                    F{101} = e64;
                    F{139} = e64;
                    F{165} = e64;
                    F{655} = e64;
                    F{892} = e64;
                    F{894} = e64;
                    F{978} = e64;
                    F{938} = e64;
                    F{986} = e64;
                    F{1007} = e64;
                    F{1080} = e64;
                    F{1217} = e64;
                    F{1245} = e64;
                    F{1411} = e64;
                    F{1559} = e64;
                    F{1566} = e64;
                    F{1586} = e64;
                    F{1597} = e64;
                    F{1670} = e64;
                    F{1731} = e64;
                    F{1859} = e64;
                otherwise
                     for jjj = 1:num_sweeps
                        if length(Fdel{jjj}) < mean(cellfun(@length, Fdel(max(jjj-100,1):min(jjj+100,num_sweeps))))*delper
                            Fdel{jjj} = e64;
                        elseif length(Fdel{jjj}) > mean(cellfun(@length, Fdel(max(jjj-100,1):min(jjj+100,num_sweeps))))*(1+1-delper)
                            Fdel{jjj} = e64;
                        end
                    end
            end
%         end

        % backward --> green curve
        % forward --> blue curve

%         num_sweeps = length(F);
        
        if save_fb_file
                        
            Ftemp = F;
%             Btemp = B;

            F = Fdel;
%             B = Bdel;
            
            fprintf('Writing %s\n', fb_file);
            save(fb_file, 'F', 'B', 'D_F', 'D_B', 'acq_time', 'deadtime', 'tpb_ns','-v7.3');
            
            pause(1)
            
            F = Ftemp;
%             B = Btemp;
            
        end
        
    elseif plot_fb_counts_vs_sweep||plot_fb_pulse_folding||convert_fb_to_rt||create_fb_tot

        fprintf('Reading %s... ', fb_file);
        load(fb_file);
        num_sweeps = length(F);
        fprintf('%d sweeps\n', num_sweeps);

    end

    %% Show total counts vs. sweeps on FB files
    if plot_fb_counts_vs_sweep
        
        figure(i)
        clf
        set(gcf, 'Name', fb_file);
        
        if convert_spcx_to_fb
            subplot(2,1,2)
            plot(1:num_sweeps, [cellfun(@length,Fdel)]);

            subplot(2,1,1)
        end
        plot(1:num_sweeps, [cellfun(@length,F)]);
        
        xlabel 'Sweep'
        ylabel 'Total Counts'
        legend('F', 'B')
        title(filename{i});
        set(gcf, 'position', [10 420 1661 420])
        
        if save_counts_vs_sweeps
            savefig(gcf, strrep(rt_file, '.mat', 'Counts.fig'),'compact')
            export_fig(gcf, strrep(rt_file, '.mat', 'Counts.png'), '-r300');
        end

        if convert_spcx_to_fb
            F = Fdel;
%             B = Bdel;
        end
        
    elseif convert_spcx_to_fb
        F = Fdel;
%         B = Bdel;
    end
    
    %% Display pulse folding for each FB file
    if plot_fb_pulse_folding

        psum = zeros(num_pulses,1);			% area of each pulse
        psum2 = zeros(num_pulses,1);        % background of each pulse (area in the opposite detector)
        ER = zeros(num_pulses,1);
        rise_time = zeros(num_pulses,1);
        rise_time_2 = zeros(num_pulses,1);
        fwhm = zeros(num_pulses,1);
        g2_SG2 = cell(num_pulses,1);
        
        pulse_gate_refined = cell(1, num_pulses);
        pulse_delay_refined = zeros(num_pulses,1);
        pulse_hi_refined = zeros(num_pulses,1);
        
        g = zeros(num_pulses,3);
        
        ws = warning('off','all');  % Turn off fit warnings

        file_div = length(F); % dividing run file to sweep chunks
        Fnew = F(1:find(~cellfun(@isempty, F),1,'last'));
        Ftemp = F;
        F = Fnew;
        num_sweeps = length(F);

        if num_sweeps >= file_div
            F_reshaped = reshape(F(1:file_div*floor(num_sweeps/file_div)), [file_div, floor(num_sweeps/file_div)]);
%             B_reshaped = reshape(B(1:file_div*floor(num_sweeps/file_div)), [file_div, floor(num_sweeps/file_div)]);
            if floor(num_sweeps/file_div) ~= num_sweeps/file_div
                F_reshaped(1:length(F)-(file_div*floor(num_sweeps/file_div)),floor(num_sweeps/file_div)+1) = F(file_div*floor(num_sweeps/file_div)+1:end);
%                 B_reshaped(1:length(B)-(file_div*floor(num_sweeps/file_div)),floor(num_sweeps/file_div)+1) = B(file_div*floor(num_sweeps/file_div)+1:end);
            end
        else
            F_reshaped = F;
%             B_reshaped = B;
        end

        pulse_fit = (1:num_pulses).*fit;
        pulse_fit = pulse_fit(~pulse_fit == 0);
        
        all_periods = zeros(size(F_reshaped,2),1);
        all_rise_time_data = zeros(size(F_reshaped,2),1);
        all_pulse_gate = cell(num_pulses,size(F_reshaped,2));
        
        for o = 1:size(F_reshaped,2)
            
            if period_iteration
                period_arr = period_lim(1):dperiod:period_lim(2);
            elseif exist(rt_file, 'file') == 2
                period_load = load(rt_file_tot, 'period');
                period_arr = period_load.period;
                fprintf('Loading factor from RT file %s.\n', rt_file);
            else
                period_arr = period;
            end
            
                rise_time_data = zeros(num_pulses, length(period_arr));
            
            for pp = 1:length(period_arr)
                
                period = period_arr(pp);
                
                %% Define pulse window functions
                pulse_gate = cell(num_pulses,1);
                
                for j = 1:num_pulses
                    pulse_gate{j} = @(t) find_pulse_train_level(t, pulse_hi(j), period-pulse_hi(j), pulse_delay(j));
                end
                
                % define pulse_forward and pulse_backward functions
                if isequal(pulse_dir, [1 1 1 1])
                    % for Ff, spectrum measurements
                    pulse_forward = @(x) (pulse_gate{1}(x) | pulse_gate{2}(x) | pulse_gate{3}(x) | pulse_gate{4}(x));
                    pulse_backward = @(x) [];
                    % background taking only information from pulse 2 and 4
                    pulse_forward_meas = @(x) (pulse_gate{2}(x) | pulse_gate{4}(x));
                    pulse_backward_meas = @(x) 0;
                    
                elseif isequal(pulse_dir, [1 1 1 1 1 1 1 1])
                    % for Ff, spectrum measurements
                    pulse_forward = @(x) (pulse_gate{1}(x) | pulse_gate{2}(x) | pulse_gate{3}(x) | pulse_gate{4}(x) |...
                        pulse_gate{5}(x) | pulse_gate{6}(x) | pulse_gate{7}(x) | pulse_gate{8}(x));
                    pulse_backward = @(x) [];
                    % background taking only information from pulses 2, 4, 6, 8
                    pulse_forward_meas = @(x) (pulse_gate{2}(x) | pulse_gate{4}(x) | pulse_gate{6}(x) | pulse_gate{8}(x));
                    pulse_backward_meas = @(x) 0;

                else
                    error('Unknown pulse sequence\n');
                end

%%
                bin = 10;
                tmax = 33000*period; 
                tmin = 0*period;
                
                Ffold = cell2mat(F_reshaped(~cellfun(@isempty, F_reshaped(:,o)),o));
                Bfold = cell2mat(B);
                
                Ffold = mod(Ffold(Ffold>tmin & Ffold<tmax) - delay, period);
                Bfold = mod(Bfold(Bfold>tmin & Bfold<tmax) - delay, period);
                
                [Hfold, tfold] = my_histogram2({Ffold, Bfold}, bin, 0, period);
%                 [Hfold, tfold] = my_histogram2(Ffold, bin, 0, period);
                
                pfoldd = zeros(length(tfold), 1);	% logical level of pulses, folded
                pfold_refined = zeros(length(tfold), 1);	% logical level of pulses, folded
                
                Hfold = Hfold / num_sweeps / ((tmax-tmin)/period);		% normalize to per sweep per pulse. 
                
                u = max(Hfold(:));
                
                N2 = 24;
                SG2 = @(tt0, sig, t) exp(-(t-tt0).^N2/(2*sig^N2));
                SG2_fun = @(a,t)(a(1)*SG2(a(2),a(3),t) + a(4));
                
                for j = 1:num_pulses
                    
                    pfold = pulse_gate{j}(tfold + delay);
                    
                    if ~fit(j) && input
                        ps = find(pfold,1);
                        pe = find(pfold,1,'last');
                        pfold = [false(ps-round((pe-ps)/1.5),1) ; true(round((pe-ps)*2.2),1) ; false(length(pfold)-ps+round((pe-ps)/1.5)-round((pe-ps)*2.2),1)];
                    end
                    
                    % to refine the positions/widths of the detection windows:
                    [maxi, position_maxi] = max(Hfold(:, pulse_dir(j)).*pfold);	% maximum+location of each pulse
                    
                    % finding refined pulse windows
                    noise_thresholdL = threshL(j) * maxi; % Left threshold
                    noise_thresholdR = threshR(j) * maxi; % Right threshold
                    
                    [thresL, pos_thresL] = min(abs(Hfold(1:position_maxi, pulse_dir(j)).*pfold(1:position_maxi,1)-noise_thresholdL));	% to define the frame position to the left\
                    if pos_thresL == 1
                        findthresL = abs(Hfold(1:position_maxi, pulse_dir(j)).*pfold(1:position_maxi,1)-noise_thresholdL);
                        findthresL1 = findthresL(findthresL~=0);
                        pos_thresL = find(findthresL == findthresL1(1));
                        pos_thresL = pos_thresL(end);
                    end
                    [thresR, pos_thresR] = min(abs(Hfold(position_maxi:length(tfold), pulse_dir(j)).*pfold(position_maxi:length(tfold),1)-noise_thresholdR)); %to define the frame position to the right
                    pos_thresR = pos_thresR + position_maxi;
                    
                    % new refined windows:
                    pulse_delay_refined(j) = delay + tfold(pos_thresL);
                    pulse_hi_refined(j) = tfold(pos_thresR) - tfold(pos_thresL);
                    pulse_gate_refined{j} = @(t) find_pulse_train_level(t, pulse_hi_refined(j), period-pulse_hi_refined(j), pulse_delay_refined(j));
                    pfold_refined = pfold_refined + pulse_gate_refined{j}(tfold + delay);
                    
                    psum(j) = sum(Hfold(pulse_gate_refined{j}(tfold + delay), pulse_dir(j)));
                    psum2(j) = sum(Hfold(pulse_gate_refined{j}(tfold + delay), 3 - pulse_dir(j)));
                    
                    if fit(j)
                        
                        ttfold0 = tfold(pulse_gate_refined{j}(tfold+delay));

                        [y50,pos_y50] = min(abs(0.5*maxi-Hfold(1:position_maxi,pulse_dir(j)).*pfold(1:position_maxi,1)));
                        [y50_2,pos_y50_2] = min(abs(0.5*maxi-Hfold(position_maxi:length(tfold),pulse_dir(j)).*pfold(position_maxi:length(tfold),1)));
                        pos_y50_2 = pos_y50_2 + position_maxi;
                        fwhm(j) = (tfold(pos_y50_2) - tfold(pos_y50))*tpb_ns;
                        
                        g0_SG2 = [maxi;(ttfold0(end)+ttfold0(1))/2;pulse_hi_refined(j)/2;thresL];
                        g2_SG2{j} = nlinfit(ttfold0, Hfold(pulse_gate_refined{j}(tfold+delay),pulse_dir(j)), SG2_fun, g0_SG2);
                        
                        t = tfold(pulse_gate_refined{j}(tfold+delay));
%                         tt = tfold(pulse_gate{j}(tfold+delay));
                        
                        eta = 20; %to be sure to average the  max value on the plateau
                        plateau = mean(Hfold(pos_thresL+eta:pos_thresR-eta));
                        
                        %rise time
                        [y10_data, pos_y10_data] = min(abs(0.1*plateau - Hfold(1:position_maxi,pulse_dir(j)).*pfold(1:position_maxi,1)));
                        [y90_data, pos_y90_data] = min(abs(0.9*plateau - Hfold(1:position_maxi,pulse_dir(j)).*pfold(1:position_maxi,1)));
                        rise_time_data(j, pp) = (tfold(pos_y90_data) - tfold(pos_y10_data))*tpb_ns;
                        
                        %extinction ratio
                        ER(j) = Hfold(pos_thresL)/max(SG2_fun(g2_SG2{j}, t));
                        
                    else
                        pfold = pulse_gate{j}(tfold + delay);
                        
                    end
                    pfoldd = pfoldd + pfold;
                    
                end
                
                pfold = pfoldd;
                plot_bool = false;
                
                if pp > 1
                    if sum(rise_time_data(:,pp)) < all_rise_time_data(o)
                        all_periods(o) = period_arr(pp);
                        all_rise_time_data(o) = sum(rise_time_data(:,pp));
                        all_pulse_gate(:,o) = pulse_gate;
                        plot_bool = true;
                        rise_time = sum(rise_time_data(:,pp))/2;
                    end
                else
                    all_periods(o) = period_arr(pp);
                    all_rise_time_data(o) = sum(rise_time_data(:,pp));
                    all_pulse_gate(:,o) = pulse_gate;
                    plot_bool = true;
                    rise_time = sum(rise_time_data(:,pp))/2;
                end
            
                if plot_bool
                    % plotting
                    figure(i+100*o)
                    clf
                    set(gcf, 'Name', fb_file);
                    axes('Position',[0.0436473429951689 0.11952380952381 0.936425120772947 0.815000000000001]);
                    
                    hold all
                    plot(tfold*tpb_ns, [Hfold u*pfold]);
                    plot(tfold*tpb_ns, u*pfold_refined, 'k--');
                    xlabel 't (ns)'
                    xlim([0 period*tpb_ns])
                    title(sprintf('Folded Pulses %s', filename{i}));
                    
                    legend('F', 'B');
                    
                    for j = 1:num_pulses
                        if fit(j)
                            t = tfold(pulse_gate_refined{j}(tfold+delay));
                            plot(t*tpb_ns, SG2_fun(g2_SG2{j}, t), 'm--')
                            txxt = text((pulse_delay(j)-delay+pulse_hi(j)/2)*tpb_ns, u*5/6, ...
                                sprintf('%.2g \n fwhm %.1f ns \n rise time DATA %.1f ns \n ER %.1f%% (%.1f dB)', ...
                                psum(j), fwhm(j), rise_time_data(j, pp), ER(j)*100, 10*log10(ER(j))), 'HorizontalAlignment', 'center');
                            
                        else
                            txxt = text((pulse_delay(j)-delay+pulse_hi(j)/2)*tpb_ns, u*5/6, ...
                                sprintf('%.2g', psum(j)), 'HorizontalAlignment', 'center');
                        end
                        
                        txxt.FontSize = 9;
                        txxt.FontWeight = 'bold';
                    end
                    
                    set(gcf, 'position', [10 420 1661 420])
                    savefig(gcf, strrep(rt_file, '.mat', sprintf('Pulses_%g.fig', o)),'compact')
                    export_fig(gcf, strrep(rt_file, '.mat', sprintf('Pulses_%g.png', o)), '-r300');
                    drawnow()
                end
                
            end

        end
        
        warning(ws)  % Turn fit warnings back on
        
        pulse_gate = all_pulse_gate;
        period = all_periods;
        
        % redefine pulse_forward and pulse_backward functions
        if isequal(pulse_dir, [1 1 1 1])
            % for Ff, spectrum measurements
            pulse_forward = @(x) (pulse_gate{1}(x) | pulse_gate{2}(x) | pulse_gate{3}(x) | pulse_gate{4}(x));
            pulse_backward = @(x) [];
            % background taking only information from pulse 2 and 4
            pulse_forward_meas = @(x) (pulse_gate{2}(x) | pulse_gate{4}(x));
            pulse_backward_meas = @(x) 0;
        else
            error('Unknown pulse sequence\n');
        end
                
    end

    %% Convert FB file to RT file
    if convert_fb_to_rt
        
        if ~plot_fb_pulse_folding%% removed ~
            load(rt_file, 'period', 'pulse_gate');
            
            if isequal(pulse_dir, [1 1 1 1])
                % for Ff, spectrum measurements
                pulse_forward = @(x) (pulse_gate{1}(x) | pulse_gate{2}(x) | pulse_gate{3}(x) | pulse_gate{4}(x));
                pulse_backward = @(x) [];
                % background taking only information from pulse 2 and 4
                pulse_forward_meas = @(x) (pulse_gate{2}(x) | pulse_gate{4}(x));
                pulse_backward_meas = @(x) 0;
            else
                error('Unknown pulse sequence\n');
            end
            
        end
        
    	R = cell(num_sweeps,1);
        T = cell(num_sweeps,1);
        
    	D_R = cell(num_sweeps,1);
        D_T = cell(num_sweeps,1);

        for j = 1:num_sweeps

            % transmitted forward pulses: photons detected on forward detector when there was forward pulse
            TF = F{j}(pulse_forward(F{j}));
            % numbers of detectors triggered
            D_TF = D_F{j}(pulse_forward(F{j}));
            
            % reflected backward pulses: photons detected on forward detector when there was backward pulse
            RB = F{j}(pulse_backward(F{j}));
            % numbers of detectors triggered
            D_RB = D_F{j}(pulse_backward(F{j}));
            
            % transmitted backward pulses: photons detected on backward detector when there was backward pulse
            TB = B{j}(pulse_backward(B{j}));
            % numbers of detectors triggered
            D_TB = D_B{j}(pulse_backward(B{j}));
            
            % reflected forward pulses: photons detected on backward detector when there was forward pulse
            RF = B{j}(pulse_forward(B{j}));
            % numbers of detectors triggered
            D_RF = D_B{j}(pulse_forward(B{j}));
        
            % combine reflected/transmitted from both F and B
            R{j} = sortrows([RF D_RF; RB D_RB]);
            T{j} = sortrows([TF D_TF; TB D_TB]);

            if ~isempty(R{j})
                D_R{j} = R{j}(:,2);
                R{j} = R{j}(:,1);

            else
                D_R{j} = [];
                R{j} = [];

            end
            
            if ~isempty(T{j})
                D_T{j} = T{j}(:,2);
                T{j} = T{j}(:,1);

            else
                D_T{j} = [];
                T{j} = [];

            end
            
            % for each event, compute in which period it occured
            Rn = floor((double(R{j}) - delay)/period);
            Tn = floor((double(T{j}) - delay)/period);

            Bn = floor((double(B{j}) - delay)/period);
            Fn = floor((double(F{j}) - delay)/period);
            
            tpb_us = tpb_ns*1e-3;
            
            if ~input && ~background
%                 retain data only in periods where transmission occured
                R{j} = R{j}(ismember(Rn, Tn));
                D_R{j} = D_R{j}(ismember(Rn, Tn));
                
                F{j} = F{j}(ismember(Fn, Tn));
                B{j} = B{j}(ismember(Bn, Tn));
                
                D_F{j} = D_F{j}(ismember(Fn, Tn));
                D_B{j} = D_B{j}(ismember(Bn, Tn));
                
            elseif background
%                 % retain only a small part of the data
                R{j} = R{j}(R{j} >= 26e3/tpb_us);
                T{j} = T{j}(T{j} >= 26e3/tpb_us);

                D_R{j} = D_R{j}(R{j} >= 26e3/tpb_us);
                D_T{j} = D_T{j}(T{j} >= 26e3/tpb_us);

                F{j} = F{j}(F{j} >= 26e3/tpb_us);
                B{j} = B{j}(B{j} >= 26e3/tpb_us);

                D_F{j} = D_F{j}(F{j} >= 26e3/tpb_us);
                D_B{j} = D_B{j}(B{j} >= 26e3/tpb_us);
            end

        end
        
        fprintf('Writing %s.\n', rt_file);
        save(rt_file, 'R', 'T', 'D_R', 'D_T', 'tpb_ns', 'deadtime', 'acq_time', 'period', 'delay', 'pulse_hi', ...
            'pulse_delay', 'pulse_dir', 'pulse_gate', 'rise_time', '-v7.3');
    
    elseif create_rt_tot
            fprintf('Reading %s... ', rt_file);
            load(rt_file);
            num_sweeps = length(R);
            fprintf('%d sweeps\n', num_sweeps);
    end

    if create_rt_tot
        Rtot = [Rtot; R];
        Ttot = [Ttot; T];
        
        D_Rtot = [D_Rtot; D_R];
        D_Ttot = [D_Ttot; D_T];        
    end
       
    if create_fb_tot
        Btot = [Btot; B];
        Ftot = [Ftot; F];
        
        D_Btot = [D_Btot; D_B];
        D_Ftot = [D_Ftot; D_F];
    end
    
end