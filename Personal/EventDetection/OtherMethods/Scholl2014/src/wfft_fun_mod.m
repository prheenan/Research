function [good noises]=wfft_fun(proteinType,recordingType,~,...
                                noiseLevel,minPeakDistance,...
                                minPeakHeight,bestLength)
%
% wfft - windowed fast fourier transform
% Syntax = [good noises] =
%  wfft(proteinType,recordingType,toplot,noiseLevel)
% 
% Input:
% proteinType = vector designating which force regimes 
%   I27 = 3, SNase = 2, Ank = 1.  I.e. Ank+I27 would be [1 3]
% recordingType = name of the file to open. I.e. 'Ni10C_3.afm'
% toplot = will plot if set to 1
% noiseLevel = level of noise to add
%
% Output: Classification of force-rupture peaks
% good = true if recording is selected with heuristics
% noises = measured level of noise
%
% Zack Scholl, January, 2013
% Revised 04/26/2014 (stepWindow based off window size)
% Revised 02/21/2017 by Patrick Heenan to take minPeakDistance, 
%  minPeakHeight,bestLength,file as parameters and remove plotting

goodPeaks = [];

% load a recording
file=sprintf('%s',recordingType);
gooddata = importdata(file);
%% Add noise to the recording
noise=noiseLevel*sin(rand(length(gooddata(:,2)),1))+noiseLevel*rand(length(gooddata(:,2)),1);
noise = noise-mean(noise);
noises=sqrt(mean(gooddata(4900:5000,2).^2));
gooddata(:,2)=gooddata(:,2)+noise;
% Compare noise levels
noises=sqrt(mean(gooddata(4900:5000,2).^2))-noises;

% Start timing
tic;

% Define vectors
forceData=gooddata(:,2)';
extensionData=gooddata(:,1);
x=size(extensionData);
X=forceData;


%% Determine the window length, i.e. find datapoints/nm *roughly*
conversion=sum(gooddata(:,1)>10)/(max(gooddata(:,1))-10);
windowLength = round(bestLength*conversion);
stepWindow=  round(conversion*0.5);
% Get length of data vector
dataLength = length(X);
windowedDataLength = ceil(dataLength/stepWindow);
if (windowLength/2 == round(windowLength/2))
	fftLength = (windowLength/2) + 1;
else
	fftLength = ceil(windowLength/2);
end

% Initialize vector for holding windowed FT
wfftSignal = zeros(windowedDataLength,fftLength);
X = [zeros(1,windowLength/2) X zeros(1,windowLength/2)]; % padding the signal with zeros
ii = 0;
% Loop through each window, taking FT
for i = 1:stepWindow:dataLength
	ii = ii + 1;
	dataFrame = X(i:(i + windowLength - 1));  
	fftValues = abs(fft(dataFrame, windowLength));
	wfftSignal(ii,:) = fftValues(1:fftLength);
end


% Sum up the odd coefficients from the FT
coefficientSum=zeros(length(wfftSignal(:,1)),1);
for i=1:2:size(wfftSignal,2)
    coefficientSum=coefficientSum+wfftSignal(:,i);
end
coefficientSum = coefficientSum/windowLength;

% Two cases: Signal ends in rupture or cyclic 
% For rupture, simply baseline FT sum by taking average of final points
% (optimal)
% For cyclic, baseline FT by taking average of minimum of points
foo = sort(coefficientSum);
if (std(coefficientSum(end-50:end)) < 20)
    coefficientSum = coefficientSum-mean(coefficientSum(end-50:end));
else
    coefficientSum = coefficientSum-mean(foo(1:30));
end

% Cutoff first 10nm as they contain nonspecific events
forget=floor(sum(gooddata(:,1)<10)/stepWindow);
PeakSig = (coefficientSum(forget:end));

% Find peaks using the protein database parameters
[pks,locs] = findpeaks(PeakSig,'SORTSTR','descend','minpeakdistance',minPeakDistance,'minpeakheight',minPeakHeight);
% Limit the number of peaks
if (length(pks)> 20)
   pks = pks(1:20);
   locs = locs(1:20);
end

% Correct for the elimination of the nonspecific events
newlocs = [];
xx=gooddata(1:stepWindow:end,1);
for j=1:length(locs)
   newlocs = [newlocs; xx(locs(j)+forget) pks(j)]; 
   newlocs(j,1)=mean(newlocs(j,1)-1:newlocs(j,1));
end
goods = 0;
if (size(newlocs,1)>0)
newlocs=sortrows(newlocs);
end

extrapolatedData = [];
goodPeaks = zeros(size(newlocs,1),1);
previous = 0;   
for j=1:size(newlocs,1)
ii=(gooddata(:,1)>newlocs(j,1)-newlocs(j,1)*.01 & gooddata(:,1)<newlocs(j,1)+newlocs(j,1)*.011);
if (sum(ii)>0)
    extrapolatedForce = mean(gooddata(ii,2));
    extrapolatedData = [extrapolatedData; newlocs(j,1) extrapolatedForce];
    if (size(extrapolatedData,1)>1)
        jj=size(extrapolatedData,1);
        contour = extrapolatedData(jj,1)-extrapolatedData(jj-1,1);
        noIdea = 1;
        
        for iii=1:length(proteinType)
        if (contour >= proteinDatabase{proteinType(iii)}.contourLength(1) ...
            && contour <= proteinDatabase{proteinType(iii)}.contourLength(2) ...
            && extrapolatedData(jj-1,2) >= proteinDatabase{proteinType(iii)}.force(1) ...
            && extrapolatedData(jj-1,2) <= proteinDatabase{proteinType(iii)}.force(2))
            
                type = proteinDatabase{proteinType(iii)}.name;
                proteinDatabase{proteinType(iii)}.num = proteinDatabase{proteinType(iii)}.num +1;
                  goodPeaks(jj-1)=1;
              end 
              noIdea = 0;
              previous=iii;
        elseif (iii==previous ...
            && contour >= proteinDatabase{proteinType(iii)}.contourLength(1) ...
            && extrapolatedData(jj-1,2) >= proteinDatabase{proteinType(iii)}.force(1) ...
            && extrapolatedData(jj-1,2) <= proteinDatabase{proteinType(iii)}.force(2))
            
                type = proteinDatabase{proteinType(iii)}.name;
                proteinDatabase{proteinType(iii)}.num = proteinDatabase{proteinType(iii)}.num +1;
              noIdea = 0;
              previous=iii;
            
        end
        end
        
    end

end
end
      
end
