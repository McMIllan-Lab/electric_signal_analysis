#!/usr/bin/R

# NOTE: SOME CALCULATIONS WILL STILL THROW WARNINGS - RUN  options(warn=-1) to stop seeing them
# Example usage (20/11/2019):
# in_recording<-readMP3("SomeRecording.mp3")
# in_recording<-mono(in_recording,"left")
# peaksTable<-signalSummaryTable(in_recording,100)

library(tuneR)

relufy<-function(invect,lowestNumber=0) { # Converts a vector into a rectified unit (i.e. converts all values below lowestNumber to lowestNumber - a classical rectifier changes all negative numbers to zero, for example).
  outvect<-sapply(invect,max,lowestNumber)
  return(outvect)
}

monofy<-function(inSignal) {
  outSound<-mono(inSignal,which="left")
  return(outSound)
}

samplesToSeconds<-function(numSamples,samplingRate) {
  seconds<-numSamples/samplingRate
  return(seconds)
}

secondsToSamples<-function(numSeconds,samplingRate) {
  samples<-round(numSeconds*samplingRate)
  return(samples)
}

clipIncompleteWaves<-function(inSignal) {
  signalSize<-length(inSignal@left)

  maxPeaks<-which(inSignal@left > (max(inSignal@left) * 0.8))
  minPeaks<-which(inSignal@left < (min(inSignal@left) * 0.8))

  diffMaxPeaks<-diff(maxPeaks)
  avgStepSize<-round(mean(diffMaxPeaks[diffMaxPeaks > mean(diffMaxPeaks)]))
  halfStep<-round(avgStepSize/2)
  sizeMinusHalfStep<-signalSize-halfStep

# A warning will pop up if the signal cut in the edge has only one of the peaks. Maybe fix one day.
  if (any(maxPeaks > sizeMinusHalfStep) || any(minPeaks > sizeMinusHalfStep)) {
    upperSignalBound<-min(min(maxPeaks[maxPeaks > sizeMinusHalfStep]), min(minPeaks[minPeaks > sizeMinusHalfStep]))
    upperCutoff<-upperSignalBound - halfStep
    inSignal<-inSignal[1:upperCutoff]
  }

  if (any(maxPeaks < halfStep) || any(minPeaks < halfStep)) {
    lowerSignalBound<-max(max(maxPeaks[maxPeaks < halfStep]), max(minPeaks[minPeaks < halfStep]))
    lowerCutoff<-lowerSignalBound + halfStep
    inSignal<-inSignal[lowerCutoff:length(inSignal@left)]
  }

  return(inSignal)
}

burnTails<-function(inSignal,burnSeconds=1) {
#  print(paste("Burning",burnSeconds,"seconds from beginning and end of recording."))
  if (inSignal@stereo) {
   print("Input in stereo converted to mono")
   inSignal<-monofy(inSignal)
  }

  samplesToBurn<-secondsToSamples(burnSeconds,samplingRate=inSignal@samp.rate)
  startSample<-samplesToBurn
  endSample<-length(inSignal)-samplesToBurn
  inSignal<-inSignal[startSample:endSample]
  inSignal<-clipIncompleteWaves(inSignal)
  return(inSignal)
}

frequencyEstimate<-function(inSignal) {
#  peakRegionsIndices<-which(inSignal@left > quantile(inSignal@left,c(.99))) #Not usable across the board
  peakRegionsIndices<-which(inSignal@left > (max(inSignal@left) * 0.8))
  diffPRI<-diff(peakRegionsIndices)
  peakCount<-sum(diffPRI > 1)
  outValue<-round(peakCount/samplesToSeconds(length(inSignal),inSignal@samp.rate))
  return(outValue)
}

fourierTransfTable<-function(inSignal) {
  options(warn=-1)
  estimatedFreq<-frequencyEstimate(inSignal)
  minFreqRange<-round(estimatedFreq*0.5)
  maxFreqRange<-round(estimatedFreq*1.5)

  if (inSignal@stereo == TRUE) {
    print("Input signal is stereo. Converting to mono.")
    inSignal<-monofy(inSignal)
  }
  perGram<-periodogram(inSignal,frqRange=c(minFreqRange,maxFreqRange))
  mainDF<-data.frame(perGram@freq,perGram@spec)
  colnames(mainDF)<-c("freq","spec")
  mainDF<-mainDF[order(-mainDF$spec),]
  return(mainDF)
}

identifyPeaks<-function(inSignal) {
  if (inSignal@stereo) {
    print("Input signal is stereo. Converting to mono.")
    inSignal<-monofy(inSignal)
  }

  peakRegionsIndicesHigh<-which(inSignal@left > (max(inSignal@left) * 0.8))
  diffMaxPeaks<-diff(peakRegionsIndicesHigh)
  avgStepSize<-round(mean(diffMaxPeaks[diffMaxPeaks > mean(diffMaxPeaks)]))
  halfStep<-round(avgStepSize/2)

  highestLowest<-max(diffMaxPeaks[diffMaxPeaks < mean(diffMaxPeaks)])
  stepVect<-c(peakRegionsIndicesHigh[1])
  resultsDF<-data.frame(High.Peak.Index=NA,High.Peak.Value=NA,Low.Peak.Index=NA,Low.Peak.Value=NA)
  inRowNum<-1

  for (i in 1:length(diffMaxPeaks)) {
    if (diffMaxPeaks[i] <= highestLowest) {
      stepVect<-c(stepVect,peakRegionsIndicesHigh[i+1])
    } else {
      stepMaxIndex<-stepVect[which.max(inSignal@left[stepVect])]
      stepMaxValue<-max(inSignal@left[stepVect])

      stepMinIndex<-(stepMaxIndex-1) + (which.min(inSignal@left[stepMaxIndex:(stepMaxIndex + halfStep)]))
      stepMinValue<-inSignal@left[stepMinIndex]

      resultsDF[inRowNum,]<-c(stepMaxIndex,stepMaxValue,stepMinIndex,stepMinValue)
      inRowNum<-inRowNum+1
      stepVect<-c(peakRegionsIndicesHigh[i+1])
    }
  }
  return(resultsDF)
}


extractSingleWaveform<-function(inSignal,waveformNumber=1) {
  inSignal<-monofy(inSignal)
  fullPeaksTable<-identifyPeaks(inSignal)
  chosenSignal<-fullPeaksTable[waveformNumber,]

  waveformHighPeak<-chosenSignal$High.Peak.Index
  waveformLowPeak<-chosenSignal$Low.Peak.Index

  meanStepSize<-round(mean(diff(fullPeaksTable$High.Peak.Index)))

  diffHighLows<-abs(waveformHighPeak-waveformLowPeak)
  padding<-round((meanStepSize-diffHighLows)/2)

  if (padding > min(waveformHighPeak,waveformLowPeak)) {
    padding<-min(waveformHighPeak,waveformLowPeak)
  }
  lowerValue<-min(waveformHighPeak,waveformLowPeak)-padding
  higherValue<-max(waveformHighPeak,waveformLowPeak)+padding

  return(inSignal[lowerValue:higherValue])
}


identifyLowerBound<-function(waveform,windowSize=10,stepSize=10) {
  lowestPeak<-min(which.max(waveform@left),which.min(waveform@left))
  lowerBoundWindowStart<-lowestPeak-windowSize
  lowerBoundWindowEnd<-lowestPeak
#print(paste("Window Start:",lowerBoundWindowStart,"window ends:",lowerBoundWindowEnd))
  testAgain=TRUE
  while (testAgain) {
#print("Entered while loop")
#    windowVals<-tryCatch(waveform@left[lowerBoundWindowStart:lowerBoundWindowEnd], error=function(e) { return(NA) })
    windowVals<-waveform@left[lowerBoundWindowStart:lowerBoundWindowEnd]
    minVal<-min(windowVals)
    maxVal<-max(windowVals)
#print(paste("Window Bounds:",lowerBoundWindowStart,lowerBoundWindowEnd))
#print(windowVals)
    timeVect<-c(1:length(windowVals))
    correlation<-cor.test(windowVals,timeVect,method="spearman")	
    if (correlation$p.value > 0.05) {
#print("Entered SECOND if statement")
      lowerBound<-round(abs(lowerBoundWindowStart+lowerBoundWindowEnd)/2)
      testAgain=FALSE
    }
    # Update bounds values and stepNum value
    lowerBoundWindowStart=lowerBoundWindowStart-stepSize
    lowerBoundWindowEnd=lowerBoundWindowEnd-stepSize
#print(paste(lowerBoundWindowStart,lowerBoundWindowEnd))
  }
  return(lowerBound)
}


identifyUpperBound<-function(waveform,windowSize=10,stepSize=10) {
  upperPeak<-max(which.max(waveform@left),which.min(waveform@left))
  upperBoundWindowStart<-upperPeak
  upperBoundWindowEnd<-upperPeak+windowSize

  testAgain=TRUE
  while (testAgain) {
#print("Entered while loop")
    windowVals<-waveform@left[upperBoundWindowStart:upperBoundWindowEnd]
    minVal<-min(windowVals)
    maxVal<-max(windowVals)
#print(paste("Window Bounds:",lowerBoundWindowStart,lowerBoundWindowEnd))
#print(windowVals)
    timeVect<-c(1:length(windowVals))
    correlation<-cor.test(windowVals,timeVect,method="spearman")
    if (correlation$p.value > 0.05) {
#print("Entered SECOND if statement")
      upperBound<-round(abs(upperBoundWindowStart+upperBoundWindowEnd)/2)
      testAgain=FALSE
    }
    # Update bounds values and stepNum value
    upperBoundWindowStart=upperBoundWindowStart+stepSize
    upperBoundWindowEnd=upperBoundWindowEnd+stepSize
#print(paste(lowerBoundWindowStart,lowerBoundWindowEnd))
  }
  return(upperBound)
}


findBounds<-function(waveform,title="No Title") { # NOTE: must contain a single waveform
  maxValIndex<-which.max(waveform@left)
  minValIndex<-which.min(waveform@left)
  minMaxInterval<-abs(maxValIndex-minValIndex)

  #pdf(paste(title,".pdf",sep=""))
  plot(waveform,main=title,xlab="Time (s)", ylab="Amplitude")
  abline(h=0,col="red")

  # Identify lower bound
  lowerBoundsVector<-c()
  for (i in (round(minMaxInterval/2)):minMaxInterval) {
    newEstimate<-identifyLowerBound(waveform,i,1)
#print(paste("Lower Bound Estimate Num ",i,": ",newEstimate))
    lowerBoundsVector<-c(lowerBoundsVector,newEstimate)
  }
  lowerBound<-round(median(lowerBoundsVector))
  lowestPeak<-min(which.max(waveform@left),which.min(waveform@left))
  lowPeakVal<-waveform@left[lowestPeak]
  distFromLowerPeak<-lowestPeak-lowerBound

  abline(v=samplesToSeconds(lowerBound,waveform@samp.rate),col="green")
  points(samplesToSeconds(lowestPeak,waveform@samp.rate),lowPeakVal,pch=23,col="blue")

  #Identify upper bound
  upperBoundsVector<-c()
  for (i in (round(minMaxInterval/2)):round(minMaxInterval*1.5)) { # Larger window sizes give better reads for the smaller slope
    newEstimate<-identifyUpperBound(waveform,i,1)
#print(paste("Upper Bound Estimate Num ",i,": ",newEstimate))
    upperBoundsVector<-c(upperBoundsVector,newEstimate)
  }
  upperBound<-round(median(upperBoundsVector))
  upperPeak<-max(which.max(waveform@left),which.min(waveform@left))
  upperPeakVal<-waveform@left[upperPeak]
  distFromUpperPeak<-upperBound-upperPeak

  abline(v=samplesToSeconds(upperBound,waveform@samp.rate),col="green")
  points(samplesToSeconds(upperPeak,waveform@samp.rate),upperPeakVal,pch=23,col="blue")
#  dev.off()
#  return(c(lowerBound,upperBound))
  return(c(distFromLowerPeak,distFromUpperPeak))
}


signalSummaryTable<-function(inSignal,howMany="ALL",location="Undetermined",individual="Undetermined",gain="?") {
  if (inSignal@stereo == TRUE) {
    print("Input signal is stereo. Converting to mono.")
    inSignal<-monofy(inSignal)
  }
  peaksTable<-identifyPeaks(inSignal)
  totalWaveforms<-nrow(peaksTable)

  if (howMany == "ALL") {
    print(paste("Sampling all",totalWaveforms,"waveforms in recorded signal"))
    peakNumbers<-c(1:totalWaveforms)

  } else if (is.numeric(howMany)) {
    howMany<-min(round(howMany),totalWaveforms)
    print(paste("Sampling",howMany,"random waveforms."))
    randVector<-sample(1:totalWaveforms,howMany)
    randVector<-randVector[order(randVector)]

    peakNumbers<-randVector
  } else {
    print(paste("Your input for 'howMany',",howMany,"was not understood"))
    quit()
  }

  resultsDF<-data.frame(Peak.Number=peakNumbers,Lower.Bound.Index=NA,P1.Index=peaksTable$High.Peak.Index[peakNumbers],Intercept.Index=NA,P2.Index=peaksTable$Low.Peak.Index[peakNumbers],Upper.Bound.Index=NA,P1.Value=peaksTable$High.Peak.Value[peakNumbers],P2.Value=peaksTable$Low.Peak.Value[peakNumbers],P1.P2.Ratio=NA,Time.To.P1=NA,Time.To.Switch=NA,Time.To.Repolarize=NA,Decay.Lambda=NA)
  resultsDF$P1.P2.Ratio<-abs(resultsDF$P1.Value/resultsDF$P2.Value)

  interceptsVect<-c()
  lowerBoundsVect<-c()
  upperBoundsVect<-c()
  lambdasVect<-c()
  inputIndex<-1

  pdffilename<-paste(location,"_",individual,"_G",gain,".pdf",sep="")
  pdf(pdffilename)

  for (i in 1:nrow(resultsDF)) { # CAUTION: This section will produce incorrect results if signal is inverted. Avoid this by always using the channel that records the first peak as positive.
    waveformNum<-resultsDF$Peak.Number[i]
    identifier<-paste(location,"_",individual,"_","gain-",gain,"_","pulse-",waveformNum,"_",i,"-",nrow(resultsDF),sep="")
    print(paste(i,": Summarizing peak",waveformNum))

    waveform<-extractSingleWaveform(inSignal,waveformNum)
    highPeak<-resultsDF$P1.Index[i]
    lowPeak<-resultsDF$P2.Index[i]

    indexVector<-c(highPeak:lowPeak)
    interceptIndex<-indexVector[which.min(abs(inSignal@left[highPeak:lowPeak]))]
    interceptsVect[i]<-interceptIndex

    bounds<-findBounds(waveform,identifier)

#print(paste("Bounds were found at:",bounds))
    lowerBound<-highPeak-bounds[1]
    upperBound<-lowPeak+bounds[2]
    lowerBoundsVect[i]<-lowerBound
    upperBoundsVect[i]<-upperBound

    repolarizationSlice<-inSignal@left[lowPeak:upperBound] # Caution: this assumes the negative peak goes after the positive peak
    repoSliceList<-list(sample=0:(length(repolarizationSlice)-1),value=repolarizationSlice)
    fitModel<-tryCatch(nls(value ~ SSasymp(sample, Asym, R0, lrc), data=repoSliceList),error=function(e) { return(NA) })

    if ( is.na(fitModel)) {
      lambdasVect[i]<-NA
    } else {
      lambdasVect[i]<-exp(coef(fitModel)[["lrc"]])
    }
  }

  dev.off()

  resultsDF$Intercept.Index<-interceptsVect
  resultsDF$Lower.Bound.Index<-lowerBoundsVect
  resultsDF$Upper.Bound.Index<-upperBoundsVect
  resultsDF$Time.To.P1<-samplesToSeconds(resultsDF$P1.Index-resultsDF$Lower.Bound.Index,inSignal@samp.rate)
  resultsDF$Time.To.Switch<-samplesToSeconds(resultsDF$P2.Index-resultsDF$P1.Index,inSignal@samp.rate)
  resultsDF$Time.To.Repolarize<-samplesToSeconds(resultsDF$Upper.Bound.Index-resultsDF$P2.Index,inSignal@samp.rate)
  resultsDF$Decay.Lambda<-lambdasVect

  return(resultsDF)
}


# The function that takes MP3 files with their collection data, already organized in a table (see Table_Brachys_singleGain.csv), and gives out the average values of the signalSummaryTable for each individual)
totalStudyAnalysis<-function(inTable,pulsesPerIndiv=10) {
  outDFMeans<-data.frame(ID=NA,Location=NA,P1.Value=NA, P2.Value=NA,P1.P2.Ratio=NA,Time.To.P1=NA,Time.To.Switch=NA,Time.To.Repolarize=NA,Decay.Lambda=NA,Num.Pulses=NA)
  rawData<-data.frame()

  for (i in 1:length(rownames(inTable))) {
    fileName<-as.character(inTable[i,"Recording"])
    gain<-inTable[i,"Gain"]
    idNum<-as.character(inTable[i,"ID"])
    location<-as.character(inTable[i,"Location"])
    print(paste("Analyzing",fileName))

    recording<-readMP3(fileName)
    burntInSound<-burnTails(recording,1)
    burntInSound@left<-burntInSound@left/gain
    indivTable<-signalSummaryTable(burntInSound,pulsesPerIndiv,location,idNum,gain)
    indivTable$Location<-location
    indivTable$ID<-idNum
    rawData<-rbind(rawData,indivTable)
    outDFMeans[i,]<-c(idNum,location,colMeans(indivTable[,7:13]),pulsesPerIndiv)
#    outDFSds[i,]<-c(idNum,location,colMeans(indivTable[,7:13]))
  }

  return(list(outDFMeans,rawData))
}









