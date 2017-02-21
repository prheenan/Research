def wfft_fun(proteinType,recordingType,toplot,noiseLevel):
    #
    # wfft - windowed fast fourier transform
    # Syntax = [good noises] =
    #  wfft(proteinType,recordingType,toplot,noiseLevel)
    # 
    # Input:
    # proteinType = vector designating which force regimes 
    #   I27 = 3, SNase = 2, Ank = 1.  I.e. Ank+I27 would be [1 3]
    # recordingType = name of the file to open. I.e. 'Ni10C_3.afm'
    # toplot = will plot if set to 1
    # noiseLevel = level of noise to add
    #
    # Output: Classification of force-rupture peaks
    # good = true if recording is selected with heuristics
    # noises = measured level of noise
    #
    # Zack Scholl, January, 2013
    # Revised 04/26/2014 (stepWindow based off window size)

    warning off
    ## Protein database in order of increasing force
    proteinDatabase{3}.name='I27'
    proteinDatabase{3}.contourLength=[15 35]
    proteinDatabase{3}.force = [95 300]
    proteinDatabase{3}.num = 0
    proteinDatabase{2}.name='SNase'
    proteinDatabase{2}.contourLength=[23 70]
    proteinDatabase{2}.force = [10 85]
    proteinDatabase{2}.num = 0
    proteinDatabase{1}.name='Ank'
    proteinDatabase{1}.contourLength=[4 20]
    proteinDatabase{1}.force = [6.5 50]
    proteinDatabase{1}.num = 0
    goodPeaks = []

    # Set the window length
    bestLength = 5 #nm
    # Use the protein database to deterimine the peak finding parameters
    minPeakDistance = 0
    minPeakHeight = 0

    # Define vectors
    forceData=gooddata(:,2)'
    extensionData=gooddata(:,1)
    x=size(extensionData)
    X=forceData


    ## Determine the window length, i.e. find datapoints/nm *roughly*
    conversion=sum(gooddata(:,1)>10)/(max(gooddata(:,1))-10)
    windowLength = round(bestLength*conversion)
    stepWindow=  round(conversion*0.5)
    # Get length of data vector
    dataLength = length(X)
    windowedDataLength = ceil(dataLength/stepWindow)
    if (windowLength/2 == round(windowLength/2))
            fftLength = (windowLength/2) + 1
    else
            fftLength = ceil(windowLength/2)
    end

    # Initialize vector for holding windowed FT
    wfftSignal = zeros(windowedDataLength,fftLength)
    X = [zeros(1,windowLength/2) X zeros(1,windowLength/2)] # padding the signal with zeros
    ii = 0
    # Loop through each window, taking FT
    for i = 1:stepWindow:dataLength
            ii = ii + 1
            dataFrame = X(i:(i + windowLength - 1))  
            fftValues = abs(fft(dataFrame, windowLength))
            wfftSignal(ii,:) = fftValues(1:fftLength)
    end


    # Sum up the odd coefficients from the FT
    coefficientSum=zeros(length(wfftSignal(:,1)),1)
    for i=1:2:size(wfftSignal,2)
        coefficientSum=coefficientSum+wfftSignal(:,i)
    end
    coefficientSum = coefficientSum/windowLength

    # Two cases: Signal ends in rupture or cyclic 
    # For rupture, simply baseline FT sum by taking average of final points
    # (optimal)
    # For cyclic, baseline FT by taking average of minimum of points
    foo = sort(coefficientSum)
    if (std(coefficientSum(end-50:end)) < 20)
        coefficientSum = coefficientSum-mean(coefficientSum(end-50:end))
    else
        coefficientSum = coefficientSum-mean(foo(1:30))
    end

    # Cutoff first 10nm as they contain nonspecific events
    forget=floor(sum(gooddata(:,1)<10)/stepWindow)
    PeakSig = (coefficientSum(forget:end))

    # Find peaks using the protein database parameters
    [pks,locs] = findpeaks(PeakSig,'SORTSTR','descend','minpeakdistance',minPeakDistance,'minpeakheight',minPeakHeight)
    # Limit the number of peaks
    if (length(pks)> 20)
       pks = pks(1:20)
       locs = locs(1:20)
    end

    # Correct for the elimination of the nonspecific events
    newlocs = []
    xx=gooddata(1:stepWindow:end,1)
    for j=1:length(locs)
       newlocs = [newlocs xx(locs(j)+forget) pks(j)] 
       newlocs(j,1)=mean(newlocs(j,1)-1:newlocs(j,1))
    end
    goods = 0
    if (size(newlocs,1)>0)
    newlocs=sortrows(newlocs)
    end


    timed=toc

    ## Plot coefficientSum if toplot==2
    xxx = (1:stepWindow:dataLength)
    xxx = (xxx-500)*.0935
    if (toplot==2)
        figure(4)
        subplot(2,1,1)
        plot(xxx,coefficientSum(1:end),'b','LineWidth',1) 
        hold on
        for j=1:length(locs)
        plot(xxx(locs(j)+forget),pks(j)+15,'kv','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
        end
        xlabel('extension [nm]')
        ylabel('sum_{n,odd} X_n(k)')
        axis([-5 360 -5 450])

    end


    ## Determine which peaks are events matching up with protein database
    ## parameters
    if (toplot==1)
        hold on
        plot(gooddata(:,1),gooddata(:,2),'b','LineWidth',0.5)
        plot(-1000,-1000,'kv','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
        plot(-1000,-1000,'kv','markerfacecolor',[30/255 144/255 1],'markeredgecolor',[30/255 144/255 1])
        plot(-1000,-1000,'kv','markerfacecolor',[1 1 0],'markeredgecolor',[1 1 0]) 
        plot(-1000,-1000,'kv','markerfacecolor',[124/255 252/255 0],'markeredgecolor',[124/255 252/255 0])           
        legend('Data','Nonspecific','I27','SNase','Ank')
    end
    extrapolatedData = []
    goodPeaks = zeros(size(newlocs,1),1)
    previous = 0   
    for j=1:size(newlocs,1)
    ii=(gooddata(:,1)>newlocs(j,1)-newlocs(j,1)*.01 & gooddata(:,1)<newlocs(j,1)+newlocs(j,1)*.011)
    if (sum(ii)>0)
        extrapolatedForce = mean(gooddata(ii,2))
        extrapolatedData = [extrapolatedData newlocs(j,1) extrapolatedForce]
        if (size(extrapolatedData,1)>1)
            jj=size(extrapolatedData,1)
            contour = extrapolatedData(jj,1)-extrapolatedData(jj-1,1)
            noIdea = 1

            for iii=1:length(proteinType)
            if (contour >= proteinDatabase{proteinType(iii)}.contourLength(1) ...
                && contour <= proteinDatabase{proteinType(iii)}.contourLength(2) ...
                && extrapolatedData(jj-1,2) >= proteinDatabase{proteinType(iii)}.force(1) ...
                && extrapolatedData(jj-1,2) <= proteinDatabase{proteinType(iii)}.force(2))

                    type = proteinDatabase{proteinType(iii)}.name
                    proteinDatabase{proteinType(iii)}.num = proteinDatabase{proteinType(iii)}.num +1
                  if (toplot==1) 
                        if (proteinType(iii)==3)
                      plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[30/255 144/255 1],'markeredgecolor',[30/255 144/255 1])
                        end
                        if (proteinType(iii)==2)
                      plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[1 1 0],'markeredgecolor',[1 1 0]) 
                        end
                        if (proteinType(iii)==1)
                      plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[124/255 252/255 0],'markeredgecolor',[124/255 252/255 0])                        
                        end
                      goodPeaks(jj-1)=1
                  end 
                  noIdea = 0
                  previous=iii
            elseif (iii==previous ...
                && contour >= proteinDatabase{proteinType(iii)}.contourLength(1) ...
                && extrapolatedData(jj-1,2) >= proteinDatabase{proteinType(iii)}.force(1) ...
                && extrapolatedData(jj-1,2) <= proteinDatabase{proteinType(iii)}.force(2))

                    type = proteinDatabase{proteinType(iii)}.name
                    proteinDatabase{proteinType(iii)}.num = proteinDatabase{proteinType(iii)}.num +1
                  if (toplot==1) 
                        if (proteinType(iii)==3)
                      plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[30/255 144/255 1],'markeredgecolor',[30/255 144/255 1])
                        end
                        if (proteinType(iii)==2)
                      plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[1 1 0],'markeredgecolor',[1 1 0]) 
                        end
                        if (proteinType(iii)==1)
                      plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[124/255 252/255 0],'markeredgecolor',[124/255 252/255 0])                        
                        end
                      goodPeaks(jj-1)=1
                  end 
                  noIdea = 0
                  previous=iii

            end
            end
                if (extrapolatedForce > 10 && contour>5 && noIdea == 1)
                   if (toplot==1)
                       plot(newlocs(jj-1,1),extrapolatedData(jj-1,2)+20,'kv','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])            
                   end
                end  

        end

    end
    end

    # # Nice plotting
    if (toplot ==1 )
        hXLabel=xlabel('extension [nm]')
         hYLabel=ylabel('force [pN]') 
            axis([-5  max(extrapolatedData(:,1))*1.5 -25 max(gooddata(:,2))*1.1])
        #axis([-5  350 -30 350])

    end

    bestRun=0
    longestRun = 0
    for j=1:length(goodPeaks)
        if (goodPeaks(j)>0)
           longestRun = longestRun+1
        else
            if (longestRun > bestRun)
                bestRun = longestRun
            end
            longestRun=0
        end
    end
    good = 0
    if ((sum(proteinType==1)>0 && proteinDatabase{1}.num >= 4) ...
            || (sum(proteinType==2)>0 && sum(proteinType==3)>0 && proteinDatabase{2}.num+proteinDatabase{3}.num >= 3) ...
            || (sum(proteinType==3)>0 && proteinDatabase{3}.num >= 3))
        good = 1
    end

    end  
