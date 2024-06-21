%% Mouth-Based IMU Feature Engineering for Future Classification Models
% Utilises Gabler 2020 (https://    )

% This scripts stands as a theoretical pipeline through which to define
% features from mouth-based IMU sensors. The specific classification these
% features are applied to is dependent upon the specific aims of the
% project, and will require appropriately labeled files to determine these
% classifications.

% Utilises four separate feature types:

% 1. Statistical Summaries (Mean, min, max, SD, variance, ratio of lin/angresultant)
% 2. Pulse Size Summaries (Width and AUC of resultant traces)
% 3. PSD Summaries
% 4. Kinematic Based Summaries

% Data are then summarised, additional normalization applied and exports

clc
clearvars 
close all

%% Data Import
% Assumes Data is formatted in individual .csv files and saved within file
% directory. Data stored as cell array, with the assume column structure
% of column 1 = time; columns 2-5 rotational acceleration
% (x,y,z,resultant); columns 6-9 rotational velocity (x,y,z,resultant);
% columns 10-13 linear acceleration (x,y,z,resultant). 

% DATA MUST BE RAW

files = dir(fullfile('*.csv'));
fs = ; %DEFINE SAMPLING FREQ

for i = 1:length(files)
    
    % Import 
    data = readmatrix(fullfile(files(i).name));

    % Extract specific axes and calculate ratio of linear/rotational
    iMG{i,1} = data(:,2:13) % extracts x,y,z,res for rot_acc, rot_vel and linear
    iMG{i,1}(:,14) = iMG{i,1}(:,13) / iMG{i,1}(:5); %ratio of lin_res/rotacc_res

    % time (assumes each datasheet has a time column in column 1)
    time = data(:,1);

    clear data

end

%% Statistical Summaries
% Defines trial mean, min, max, SD, variance and range of all variables,
% and ratio of linear to angular resultant

% Statistics
for i = 1:14
    for ii = 1:n
        trialmean(ii,i) = mean(iMG{ii,i}); % Trial mean
        trialmin(ii,i) = min(iMG{ii,i}); % Trial min
        [trialmax(ii,i), trialmaxloc(ii,i)] = max(iMG{ii,i}); % Trial max
        trialSD(ii,i) = std(iMG{ii,i}); % Trial stdev
        trialvar(ii,i) = var(iMG{ii,i}); % Trial variance
        trialrange(ii,i) = range(iMG{ii,i}); % Trial range
    end
end


%% Pulse Size Summaries
% Defines pulse duration uses threshold of trial mean as start and end
% point. Also calculate AUC for linear and rotational acceleration. 

for i = 1:length(files)
 
    % Linear Pulse Width
    impactstartlin(i,1) = find(iMG{i,4}>trialmean(i,4),1); % LinRes start time
    % Determine end of impact - finds widest possible
    potential_end_index = find(iMG{i,4}(impactstartlin(i,1):end,1) < trialmean(i,4), 1) + impactstartlin(i,1);

    % Check if the signal goes back below the mean
    if isempty(potential_end_index)
        % Set to the minimal value after the impact start
        [~, min_index] = min(iMG{i,4}(impactstartlin(i,1):end,1));
        impactendlin(i,1) = min_index + impactstartlin(i,1) - 1;
    else
        % Set to the calculated end index
        impactendlin(i,1) = potential_end_index;
    end

    linearpulsewidth(i,1) = impactendlin(i,1) - impactstartlin(i,1);
    clear potential_end_index min_index

    % Linear Pulse Width
    impactstartang(i,1) = find(iMG{i,12}>trialmean(i,12),1); % LinRes start time
    % Determine end of impact - finds widest possible
    potential_end_index = find(iMG{i,12}(impactstartang(i,1):end,1) < trialmean(i,12), 1) + impactstartang(i,1);

    % Check if the signal goes back below the mean
    if isempty(potential_end_index)
        % Set to the minimal value after the impact start
        [~, min_index] = min(iMG{i,12}(impactstartang(i,1):end,1));
        impactendang(i,1) = min_index + impactstartang(i,1) - 1;
    else
        % Set to the calculated end index
        impactendang(i,1) = potential_end_index;
    end

    angularpulsewidth(i,1) = impactendang(i,1) - impactstartang(i,1);
    clear potential_end_index min_index

    % Linear AUC
    linAUC(i,1) = trapz(time(impactstartlin(i,1):impactendlin(i,1)),iMG{i,4}(impactstartlin(i,1):impactendlin(i,1))) ;  
    angAUC(i,1) = trapz(time(impactstartang(i,1):impactendang(i,1)),iMG{i,4}(impactstartang(i,1):impactendang(i,1)));
end
 
%% PSD Summaries
% Calculate power spectral density summaries. This includes calculating
% the values at which 95% power within the signal is achieve; then splits the
% PSD into 20Hz 'bins' up to 200Hz+ to determine signal frequency content;
% and calculates AUC of these bins. 

for i = 1:n
    
    % Linear PSD
    NFFT = 500;
    FFT = fft(iMG{i,4}-mean(iMG{i,4}),NFFT);
    Pyy = FFT.*conj(FFT);
    f = fs/2*linspace(0,1,NFFT/2);
    PSDlin{i,1} = 2*abs(Pyy(1:NFFT/2));
    PSDlin{i,1} = PSDlin{i,1}./sum(PSDlin{i,1});
    clear NFFT FFT Pyy

    % Angular PSD
    NFFT = 500;
    FFT = fft(iMG{i,12}-mean(iMG{i,12}),NFFT);
    Pyy = FFT.*conj(FFT);
    PSDang{i,1} = 2*abs(Pyy(1:NFFT/2));
    PSDang{i,1} = PSDang{i,1}./sum(PSDang{i,1});
    clear NFFT FFT Pyy

    % PSD features Linear
    [PSDlinpeak(i,1), PSDlinpeakloc] = max(PSDlin{i,1}); % Max PSD absolute value
    PSDlinpeakfreq(i,1) = f(PSDlinpeakloc); % Dominant frequency lin
    cslin(:,i) = cumsum(PSDlin{i,1}); % Cumulative sum of linear PSD
    normcslin = cumsum(PSDlin{i,1})./sum(PSDlin{i,1}); % Normalised CS
    ind1 = find(normcslin <=0.95,1,'last'); %finds frequnecy at which 95% power if achieved
    linpsdninefive(i,1) = round(f(ind1));
    linpsdbins(i,1) = sum(PSDlin{i,1}(1:11,1)); % each column represents a 20Hz bin 1 to 25Hz
    linpsdAUCbins(i,1) = trapz(f(1,1:11),PSDlin{i,1}(1:11,1)); % each column represents a 20Hz bin 1 to 20
    for ii = 2:10
        linpsdbins(i,ii) = sum(PSDlin{i,1}((ii*10-8):(ii*10+1),1)); % each column represents a 20Hz bin 1 to 20Hz
        linpsdAUCbins(i,ii) = trapz(f(1,(ii*10-8):(ii*10+1)),PSDlin{i,1}((ii*10-8):(ii*10+1),1));
    end
    linpsdbins(i,11) = sum(PSDlin{i,1}(102:end,1)); % 200Hz+ bin
    linpsdAUCbins(i,11) = trapz(f(1,102:end),PSDlin{i,1}(102:end,1)); % 200Hz+ AUC
    linpsdskew(i,1) = skewness(PSDlin{i,1}); % Skewness of PSD curve 

    % PSD features angular
    [PSDangpeak(i,1), PSDangpeakloc] = max(PSDang{i,1}); % Max PSD absolute value
    PSDangpeakfreq(i,1) = f(PSDangpeakloc); % Dominant frequency ang
    csang(:,i) = cumsum(PSDang{i,1}); % Cumulative sum of angear PSD
    normcsang = cumsum(PSDang{i,1})./sum(PSDang{i,1}); % Normalised CS
    ind1 = find(normcsang <=0.95,1,'last'); %finds frequnecy at which 95% power if achieved
    angpsdninefive(i,1) = round(f(ind1));
    angpsdbins(i,1) = sum(PSDang{i,1}(1:11,1)); % each column represents a 20Hz bin 1 to 25Hz
    angpsdAUCbins(i,1) = trapz(f(1,1:11),PSDang{i,1}(1:11,1)); % each column represents a 20Hz bin 1 to 20
    for ii = 2:10
        angpsdbins(i,ii) = sum(PSDang{i,1}((ii*10-8):(ii*10+1),1)); % each column represents a 20Hz bin 1 to 20Hz
        angpsdAUCbins(i,ii) = trapz(f(1,(ii*10-8):(ii*10+1)),PSDang{i,1}((ii*10-8):(ii*10+1),1));
    end
    angpsdbins(i,11) = sum(PSDang{i,1}(102:end,1)); % 200Hz+ bin
    angpsdAUCbins(i,11) = trapz(f(1,102:end),PSDang{i,1}(102:end,1)); % 200Hz+ AUC
    angpsdskew(i,1) = skewness(PSDang{i,1});

end


%% Output Data Table

% External file import from .csv of all column headings
headings = readcell("Headings.xlsx");
headings = headings(:,3)';

% constructs output dataframe
out = [impactloc,trialmean,trialmin,trialmax,trialmaxloc,trialSD,trialvar,trialrange,linearpulsewidth,angularpulsewidth,linAUC,angAUC,PSDlinpeak,PSDlinpeakfreq,linpsdninefive,linpsdbins, linpsdskew,PSDangpeak,PSDangpeakfreq,angpsdninefive,angpsdbins,angpsdskew,linpsdAUCbins, angpsdAUCbins];

% Normalizes all values for additional features
outz = zeros(n,width(out));
outz(:,1) = impactloc;
for i = 2:width(out)
    outz(:, i) = (out(:, i) - mean(out(:, i), 'omitnan')) / std(out(:, i), 'omitnan');
end

% combines non-normlized and normalized data
out = [out,outz(:,2:end)];
out = num2cell(out(:,:));
output = [headings(1,1:148);out(:,1:148)];

% Output Data
writecell(output,"RawBig.csv")


