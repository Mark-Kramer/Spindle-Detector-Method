function [duration_sec, energy, peak_val, peak_loc, edge1_store, edge2_store, peak_freq, mean_sigma] = fun_Calculate_Spindle_Features(C_dur, Fs)

if size(C_dur,1)==0
    duration_sec = NaN;
    energy       = NaN;
    edge1_store  = NaN;
    edge2_store  = NaN;
    peak_val     = NaN;
    peak_loc     = NaN;
    peak_freq    = NaN;
    mean_sigma   = Nan;
    return
end

% C = 3 sec segments of EEG spindles, e.g. 360spindles * 1201 samples
for x=1:size(C_dur,1) % complete these operations for every segment in the data
    
    % (ignore 1st spindle since there might not be 2 s before its detection
    spindle=C_dur(x,:); % extract x segment as "spindle" and perform all operations
    
    if or(sum(spindle)==0,isnan(sum(spindle)))
        
        duration_sec(x) = NaN;
        energy(x)       = NaN;
        edge1_store(x)  = NaN;
        edge2_store(x)  = NaN;
        peak_val(x)     = NaN;
        peak_loc(x)     = NaN;
        peak_freq(x)    = NaN;
        mean_sigma(x)   = NaN;
    else
        
        %% Step 2: Obtain the complex morlet wavelet coefficients
        Num_scales = 60;
        scales     = (1:Num_scales); % define vector of scales to use for wavelet tranform
        coefs      = cwt(spindle,scales,'cmor1-1.5'); % Transform using complex Morlet 1-1.5
        
        % Baseline correction procedure below: Initial 1 second of data is taken as baseline.
        % For each scale, the mean value of the baseline period is subtracted from the data.
        baseline        = coefs(:,1:Fs); % Define the length of the baseline period as 1st second
        meanbase        = mean(baseline,2); % Take the baseline period average for each of 40 scales
        baselined_coefs = zeros(Num_scales,length(spindle)); % preallocate array for the baseline corrected coefficients
        for y=1:Num_scales
            baselined_coefs(y,:)=coefs(y,:)-meanbase(y,1);
        end
        
        %% Step 3: Restrict dataset to a broad spindle range
        sss = scal2frq(scales,'cmor1-1.5',1/Fs); % check out frequencies that correspond to selected scales (7.5 to 300 Hz)
        [~,low_fedge]  = min(abs(sss-12));
        [~,high_fedge] = min(abs(sss-15));
        coefs         = coefs(high_fedge:low_fedge,:); % define untransformed coefficients (coefs from 10 to 15.79 Hz)
        b_coefs       = baselined_coefs(high_fedge:low_fedge, Fs+1:end-Fs); % define transformed coefficents
        
        %% Find peak Frequency
        Y = fft(mean(real(b_coefs)));
        L = length(b_coefs);
        f = Fs*(0:round(L/2))/L;
        Y = Y(1:round(L/2)+1);
        [~,PeakFreqIdx] = max(abs(Y));
        
        peak_freq(x)  = f(PeakFreqIdx);
        
        %% Find sigma power
        Y = (1/(Fs*L)) * abs(Y).^2; % Calculate uV^2/Hz
        [~,low_fedge]  = min(abs(f-12));
        [~,high_fedge] = min(abs(f-15));
        mean_sigma(x)  = mean(Y(low_fedge:high_fedge));
        
        %% Plot FFT of coefficients
% % %         for sc=1:1:20
% % %             plot(real(b_coefs(sc,:)));
% % %             y=real(b_coefs(sc,:));
% % %             Y=fft(y);
% % %             L=length(y);
% % %             f = Fs*(0:(L/2))/L;
% % %             [~,idx] = max(abs(Y(1:L/2+1)));
% % %             mx(sc) = f(idx);
% % %             figure; plot(f,abs(Y(1:L/2+1)))
% % %             pause
% % %         end
        
        %% Step 4: Extract mean (real) wavelet coefficient for each timepoint & smooth
%         mean_coef=mean(abs(real(coefs(:,Fs+1:end)))); % compute the mean of the abolute value real numbers in the wavelet coefficients for spindle frequencies
        mean_coef = mean(abs(real(coefs)));
        window = ones(round(Fs/5),1)/round(Fs/5);  % smooth the data with sliding window of Fs/5
        smoothed_mean = filtfilt(window,1,mean_coef); 
        smoothed_mean = smoothed_mean(Fs+1:end-Fs);
        
        %% Step 5: Find the local maxima closest to center of smoothed mean data
        % PICK THE MAXIMAL VALUE AS THE PEAK
        peak = max(smoothed_mean); % value of peak
        indexx = find(smoothed_mean==peak); % location of peak
        
        % % OPTIONAL: Plot location of selected peak
% % %                 figure;
% % %                 sp_plot = eegfilt_mine(spindle, Fs, 12, 15);
% % %                 val=(ones(1,1));
% % % %                 plot(smoothed_mean, 'DisplayName', 'Smoothed Mean (whole)','Color',[0.8 0.5 0.8]);
% % %                 hold on
% % %                 plot([zeros(1,Fs) smoothed_mean], 'DisplayName', 'Smoothed Mean', 'Color',[0 0.8 1]);
% % %                 scatter(indexx+Fs, val.*peak, 'DisplayName', 'Peak','sizedata', 20,'MarkerFaceColor','flat','CData',[0.5 0.8 0]);
% % %                 plot(sp_plot,'DisplayName', 'Spindle filtered at [9 16]','Color',[1 0.8 0]);
% % %                 scatter(Fs*2, val.*smoothed_mean(1), 'DisplayName', 'Detection point', 'sizedata', 50, 'CData',[1 0 0]);

        %% Step 6: Calculate DURATION
        
        dropmax=.50; % define where edge of peak will be (e.g. .5 of max)
        
        val_peak=peak;
        normed=smoothed_mean/val_peak; % normalize data so max=1 and min=0
        %figure; plot(normed);
        
        centerindex=indexx; % begin search for edges at predetermined peak
        
        i=centerindex;  % find ending edge
        while normed(i)>dropmax
            if i==(length(normed))
                break
            else i=i+1;
            end
        end
        edge2=i;
        
        j=centerindex;  % find beginning edge
        while normed(j)>dropmax
            if j==1
                break
            else
                j=j-1;
            end
        end
        edge1=j;
        
        %             % more plotting
        %             plot([zeros(1,Fs+edge1) sp_plot(Fs+edge1+1:Fs+edge2-1)],'DisplayName','Detected spindle','Color','k')
        %             stem(edge1+Fs, val.*smoothed_mean(edge1), 'DisplayName', 'Spindle start', 'Color',[0 0 0],'MarkerFaceColor','k','LineStyle',':');
        %             stem(edge2+Fs, val.*smoothed_mean(edge2), 'DisplayName', 'Spindle end', 'Color',[0 0 0],'MarkerFaceColor','k','LineStyle',':');
        %             title('Spindle duration detection - Original version');
        
        duration(x)     = edge2-edge1;
        duration_sec(x) = duration(x)/Fs;
        edge2_store(x)  = edge2;
        edge1_store(x)  = edge1;
        peak_val(x)     = peak;
        peak_loc(x)     = indexx;
        energy(x)       = 0;
        
        %% Step 7: Find where this same peak drops to 50% height on either side in
        % frequency space, and find the area under this curve
        coefs2     = b_coefs(:,(edge1:edge2));  %Extract the previsouly-defined spindle from wavelet coefficients
        imagecoefs = wscalogram('', coefs2); % obtain wscalogram coefficients (corresponds to plot)
        % figure;
        % imagecoefs2=wscalogram('image', coefs2);
        % figure; imagesc(abs(coefs2).^2);
        
        dropmax2     =.50; %How frequency edge of spindle will be defined (e.g. 50% of max value at each time slice)
        areas        = zeros(1,length(coefs2));
        edge1_vector = zeros(1,length(coefs2));
        edge2_vector = zeros(1,length(coefs2));
        
        for k=1:(size(coefs2,2)) % For all timepoints in the spindle
            
            col_mx            = max(coefs2(:,k)); % Find maximum value at that timepoint (column)
            [index_x,index_y] = find(coefs2==col_mx); % Return indices of maxvalue
            
            normed2     = (imagecoefs(:,k)); % normalize data so max=1 and min=0
            col_mx2     = max(normed2); % get max value of image data
            normed2     = (imagecoefs(:,k))./col_mx2;
            centerindex = index_x; % begin search for frequency-domain edges at predetermined peak
            
            m = centerindex;  % find ending edge
            while normed2(m)>dropmax2
                if m==(length(normed2))
                    break
                else
                    m=m+1;
                end
            end
            f_edge2=m;
            edge2_vector(1,k)=m; %save edge data
            
            n=centerindex;  % find beginning edge
            while normed2(n)>dropmax2
                if n==1
                    break
                else
                    n=n-1;
                end
            end
            f_edge1=n;
            edge1_vector(1,k)=n; %save edge data
            
            areas(1,k)=trapz(imagecoefs((f_edge1:f_edge2),k)); % Calculate area under the curve
            
        end
        
        % %Optional: Plot graphic of energy results
        % figure
        % plot(edge1_vector, 'DisplayName', 'edge1_vector', 'YDataSource', 'edge1_vector', 'color', [1 0 1], 'linewidth', 5); hold on;
        % plot(edge2_vector, 'DisplayName', 'edge2_vector', 'YDataSource', 'edge2_vector','color', [1 0 1], 'linewidth', 5);
        % title 'Result of Area Boundaries';
        % set(gcf, 'color', [1 1 1]);
        % box off;
        % hold off;
        
        %% Step 9: Sum the areas to obtain the spindle power
        energy(x) = sum(areas);
        
    end % End the case there signal is not truncated
end % End the loop for all the detection points