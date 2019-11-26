#!/usr/bin/R

# NOTE: SOME CALCULATIONS WILL STILL THROW WARNINGS - RUN  options(warn=-1) to stop seeing them
# Example usage (20/11/2019):
# in_recording<-readMP3("SomeRecording.mp3")
# in_recording<-mono(in_recording,"left")
# peaksTable<-signalSummaryTable(in_recording,100)

library(tuneR)

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

burnTails<-function(inSignal,burnSeconds=1) {
  print(paste("Burning",burnSeconds,"seconds from beginning and end of recording."))
  if (inSignal@stereo) {
   print("Input in stereo converted to mono")
   inSignal<-monofy(inSignal)
  }
  samplesToBurn<-secondsToSamples(burnSeconds,samplingRate=inSignal@samp.rate)
  startSample<-samplesToBurn
  endSample<-length(inSignal)-samplesToBurn
  return(inSignal[startSample:endSample])
}

frequencyEstimate<-function(inSignal) {
  peakRegionsIndices<-which(inSignal@left > quantile(inSignal@left,c(.99))) #Not usable across the board
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

identifyPeaks<-function(inSignal,type="high") {
  if (inSignal@stereo) {
    print("Input signal is stereo. Converting to mono.")
    inSignal<-monofy(inSignal)
  }

  if (type == "high") {
    peakRegionsIndices<-which(inSignal@left > quantile(inSignal@left,c(.99)))
  } else if (type == "low") {
    peakRegionsIndices<-which(inSignal@left < quantile(inSignal@left,c(.01)))
  } else {
    paste("Type of peak",type,"unrecognized. Please input  \"high\" or \"low\"")
    return()
  }

  diffPeakRegionsIndices<-diff(peakRegionsIndices)
  stepVect<-c(peakRegionsIndices[1])
  resultsDF<-data.frame(peakIndex=NA,value=NA)
  inRowNum<-1

  for (i in 1:length(diffPeakRegionsIndices)) {
    if (diffPeakRegionsIndices[i] == 1) {
      stepVect<-c(stepVect,peakRegionsIndices[i+1])
    } else {
      stepMaxIndex<-stepVect[which.max(inSignal@left[stepVect])]
      stepMaxValue<-max(inSignal@left[stepVect])
      stepMinIndex<-stepVect[which.min(inSignal@left[stepVect])]
      stepMinValue<-min(inSignal@left[stepVect])

      if (type == "high") {
        resultsDF[inRowNum,]<-c(stepMaxIndex,stepMaxValue)
      } else if (type == "low") {
        resultsDF[inRowNum,]<-c(stepMinIndex,stepMinValue)
      }
      inRowNum<-inRowNum+1
      stepVect<-c(peakRegionsIndices[i+1])
    }
  }
  return(resultsDF)
}

peakStepsInSeconds<-function(inSignal) {
  samplingRate<-inSignal@samp.rate
  highPeaksDF<-identifyPeaks(inSignal,type="high")
  outValue<-samplesToSeconds(diff(highPeaksDF$peakIndex),samplingRate)
  return(outValue)
}

extractSingleWaveform<-function(inSignal,waveformNumber=1) {
  highPeaksDF<-identifyPeaks(inSignal,type="high")
  lowPeaksDF<-identifyPeaks(inSignal,type="low")
  waveformHighPeak<-highPeaksDF[waveformNumber,]
  waveformLowPeak<-lowPeaksDF[waveformNumber,]
  meanStepSize<-round(mean(diff(highPeaksDF$peakIndex)))
  diffHighLows<-abs(waveformHighPeak$peakIndex-waveformLowPeak$peakIndex)
  padding<-round((meanStepSize-diffHighLows)/2)
  if (padding > min(waveformHighPeak$peakIndex,waveformLowPeak$peakIndex)) {
    padding<-min(waveformHighPeak$peakIndex,waveformLowPeak$peakIndex)
  }
  lowerValue<-min(waveformHighPeak$peakIndex,waveformLowPeak$peakIndex)-padding
  higherValue<-max(waveformHighPeak$peakIndex,waveformLowPeak$peakIndex)+padding

  return(inSignal[lowerValue:higherValue])
}


signalSummaryTable<-function(inSignal,howMany="ALL") {
  if (inSignal@stereo == TRUE) {
    print("Input signal is stereo. Converting to mono.")
    inSignal<-monofy(inSignal)
  }

  highPeaks<-identifyPeaks(inSignal,"high")
  lowPeaks<-identifyPeaks(inSignal,"low")
  if (length(highPeaks) != length(lowPeaks)) {
    print("Number of high vs low peaks does not match, this may bring troubles later. Please align sample better.")
    return
  }

  if (howMany == "ALL") {
    print("Sampling all waveforms in recorded signal")
    peakNumbers<-c(1:length(highPeaks$peakIndex))
  } else if (is.numeric(howMany)) {
    howMany<-min(round(howMany),length(highPeaks$peakIndex))
    print(paste("Sampling",howMany,"random waveforms."))
    randVector<-sample.int(length(highPeaks$peakIndex),howMany)
    randVector<-randVector[order(randVector)]
    highPeaks<-highPeaks[randVector,]
    lowPeaks<-lowPeaks[randVector,]
    peakNumbers<-randVector
  } else {
    print(paste("Your input for 'howMany',",howMany,"was not understood"))
  }

  resultsDF<-data.frame(Peak.Number=peakNumbers,Lower.Bound.Index=NA,P1.Index=highPeaks$peakIndex,Intercept.Index=NA,P2.Index=lowPeaks$peakIndex,Upper.Bound.Index=NA,P1.Value=highPeaks$value,P2.Value=lowPeaks$value,P1.P2.Ratio=NA,Time.To.P1=NA,Time.To.Switch=NA,Time.To.Repolarize=NA)
  resultsDF$P1.P2.Ratio<-abs(resultsDF$P1.Value/resultsDF$P2.Value)

  interceptsVect<-c()
  lowerBoundsVect<-c()
  upperBoundsVect<-c()
  inputIndex<-1
  for (i in 1:length(rownames(resultsDF))) { # CAUTION: This section will produce incorrect results if signal is inverted. Avoid this by always using the channel that records the first peak as positive.
print(paste("Summarizing peak",i))
    waveformNum<-i
    waveform<-extractSingleWaveform(inSignal,i)
    highPeak<-resultsDF$P1.Index[i]
    lowPeak<-resultsDF$P2.Index[i]
    indexVector<-c(highPeak:lowPeak)
    interceptIndex<-indexVector[which.min(abs(inSignal@left[highPeak:lowPeak]))]
    interceptsVect[i]<-interceptIndex

    bounds<-findBounds(waveform)
#print(bounds)
    lowerBound<-highPeak-bounds[1]
    upperBound<-lowPeak+bounds[2]
    lowerBoundsVect[i]<-lowerBound
    upperBoundsVect[i]<-upperBound
  }

  resultsDF$Intercept.Index<-interceptsVect
  resultsDF$Lower.Bound.Index<-lowerBoundsVect
  resultsDF$Upper.Bound.Index<-upperBoundsVect
  resultsDF$Time.To.P1<-samplesToSeconds(resultsDF$P1.Index-resultsDF$Lower.Bound.Index,inSignal@samp.rate)
  resultsDF$Time.To.Switch<-samplesToSeconds(resultsDF$P2.Index-resultsDF$P1.Index,inSignal@samp.rate)
  resultsDF$Time.To.Repolarize<-samplesToSeconds(resultsDF$Upper.Bound.Index-resultsDF$P2.Index,inSignal@samp.rate)

#  resultsDF[(length(rownames(resultsDF))+1),]<-colMeans(resultsDF)

  return(resultsDF)
}


identifyLowerBound<-function(waveform,windowSize=10,stepSize=10) {
  lowestPeak<-min(which.max(waveform@left),which.min(waveform@left))
  lowerBoundWindowStart<-lowestPeak-windowSize
  lowerBoundWindowEnd<-lowestPeak

  testAgain=TRUE
  while (testAgain) {
#print("Entered while loop")
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


findBounds<-function(waveform) { # NOTE: must contain a single waveform
  maxValIndex<-which.max(waveform@left)
  minValIndex<-which.min(waveform@left)
  minMaxInterval<-abs(maxValIndex-minValIndex)

  # Identify lower bound
  lowerBoundsVector<-c()
  for (i in (round(minMaxInterval/2)):minMaxInterval) {
    newEstimate<-identifyLowerBound(waveform,i,1)
    lowerBoundsVector<-c(lowerBoundsVector,newEstimate)
  }
  lowerBound<-round(median(lowerBoundsVector))
  lowestPeak<-min(which.max(waveform@left),which.min(waveform@left))
  distFromLowerPeak<-lowestPeak-lowerBound

  #Identify upper bound
  upperBoundsVector<-c()
  for (i in (round(minMaxInterval/2)):round(minMaxInterval*1.5)) { # Larger window sizes give better reads for the smaller slope
    newEstimate<-identifyUpperBound(waveform,i,1)
    upperBoundsVector<-c(upperBoundsVector,newEstimate)
  }
  upperBound<-round(median(upperBoundsVector))
  upperPeak<-max(which.max(waveform@left),which.min(waveform@left))
  distFromUpperPeak<-upperBound-upperPeak

#  return(c(lowerBound,upperBound))
  return(c(distFromLowerPeak,distFromUpperPeak))
}






















if (FALSE) { ############### START DELETED SECTION

signalSlidingWindow<-function(inSignal,windowSize=10,stepSize=10) {
  if (windowSize < length(inSignal) && stepSize < length(inSignal)) {
    outputDF<-data.frame(Window.Start=NA,Window.Stop=NA,Max.Signal=NA,Min.Signal=NA,Mean.Signal=NA,Slope=NA,Is.Nonrandom=NA)
    rowNumber<-1
    windowStart<-1
    windowStop<-windowSize
    while (windowStart < length(inSignal)-windowSize+1) {
#      print(paste("Doing step",windowStart,windowStop))
      window<-inSignal@left[windowStart:windowStop]
      maxSignal<-max(window)
      minSignal<-min(window)
      meanSignal<-mean(window)
      sdSignal<-sd(window)
      numSamples<-length(window)
      timeVector<-c(1:numSamples)
      randomPoints<-rnorm(numSamples,meanSignal,sdSignal)
      tTest<-t.test(window,randomPoints)
      isNonrandom<-tTest$p.value <= 0.05
      print(isNonrandom)
      linearModel<-lm(window ~ timeVector)
      lmSlope<-linearModel$coefficients[2]
      outputDF[rowNumber,]<-c(windowStart,windowStop,maxSignal,minSignal,meanSignal,lmSlope,isNonrandom)
      windowStart<-windowStart+stepSize
      windowStop<-min(windowStop+stepSize,length(inSignal))
      rowNumber<-rowNumber+1
    }
  } else {
    paste("Aborting: Window size",windowSize,"and/or step size",stepSize,"larger than size of dataset",length(inSignal))
    return
  }

  return(outputDF)
}



signalSlidingWindow<-function(inSignal,windowSize=10,stepSize=10) {
  if (windowSize < length(inSignal) && stepSize < length(inSignal)) {
    outputDF<-data.frame(Window.Start=NA,Window.Stop=NA,Max.Signal=NA,Min.Signal=NA,Mean.Signal=NA,Slope=NA)
    rowNumber<-1
    windowStart<-1
    windowStop<-windowSize
    while (windowStart < length(inSignal)-windowSize+1) {
#      print(paste("Doing step",windowStart,windowStop))
      maxSignal<-max(inSignal@left[windowStart:windowStop])
      minSignal<-min(inSignal@left[windowStart:windowStop])
      meanSignal<-mean(inSignal@left[windowStart:windowStop])
      timeVector<-c(1:length(inSignal@left[windowStart:windowStop]))
      amplitudeVector<-inSignal@left[windowStart:windowStop]
      linearModel<-lm(amplitudeVector ~ timeVector)
      lmSlope<-linearModel$coefficients[2]
      outputDF[rowNumber,]<-c(windowStart,windowStop,maxSignal,minSignal,meanSignal,lmSlope)
      windowStart<-windowStart+stepSize
      windowStop<-min(windowStop+stepSize,length(inSignal))
      rowNumber<-rowNumber+1
    }
  } else {
    paste("Aborting: Window size",windowSize,"and/or step size",stepSize,"larger than size of dataset",length(inSignal))
    return
  }

  return(outputDF)
}


extractWaveformBracket<-function(inSignal,centerPulse=1,bracketSizeInPulses=10) {
  highPeaksDF<-identifyPeaks(inSignal,type="high")
  lowPeaksDF<-identifyPeaks(inSignal,type="low")
  meanPulseLength<-round(mean(diff(highPeaksDF$peakIndex)))
  meanDiffHighLows<-round(mean(abs(highPeaksDF$peakIndex - lowPeaksDF$peakIndex)))
  padding<-round(meanPulseLength/2)
  centerOfCenterPulse<-highPeaksDF$peakIndex[centerPulse] + round(meanDiffHighLows/2)

  if (bracketSizeInPulses %% 2 == 0) {
      pulsesBefore<-(bracketSizeInPulses/2)-1
      pulsesAfter<-(bracketSizeInPulses/2)
  } else {
      pulsesBefore<-(bracketSizeInPulses/-1)/2
      pulsesAfter<-(bracketSizeInPulses/-1)/2
  }

  bracketLowerLimit_pulses<-centerPulse-pulsesBefore
  bracketUpperLimit_pulses<-centerPulse+pulsesAfter

  totalPulses<-min(length(highPeaksDF$peakIndex),length(lowPeaksDF$peakIndex))
  totalSamples<-length(inSignal@left)

  if (bracketLowerLimit_pulses < 1) {
    shift<-abs(bracketLowerLimit_pulses)
    bracketLowerLimit_pulses<-1
    bracketUpperLimit_pulses<-bracketUpperLimit_pulses+shift
    paste("Caution: Bracket had to be shifted",shift,"pulses up.")
  }

  if (bracketUpperLimit_pulses > totalPulses) {
    shift<-abs(bracketUpperLimit_pulses-totalPulses)
    bracketUpperLimit_pulses<-totalPulses
    bracketLowerLimit_pulses<-bracketLowerLimit_pulses-shift
    paste("Caution: Bracket had to be shifted",shift,"pulses up.")
  }

  bracketLowerLimit_samples<-max(highPeaksDF$peakIndex[bracketLowerLimit_pulses]+round(meanDiffHighLows/2)-padding,0)
  bracketUpperLimit_samples<-min(highPeaksDF$peakIndex[bracketLowerLimit_pulses]+round(meanDiffHighLows/2)+padding,totalPulses)

  return(inSignal[bracketLowerLimit_samples:bracketUpperLimit_samples])

}

chooseBest<-function(inSignal,what="waveform",bracketSizeInPulses=10) {
  allHighPeaks<-identifyPeaks(inSignal,type="high")
  highestPeakIndex<-which.max(allHighPeaks$value)
  highestPeak<-allHighPeaks[highestPeakIndex,]

  if (what == "waveform") {
    output<-extractSingleWaveform(inSignal,highestPeakIndex)
  } else if (what == "bracket") {
    lowNumBracket<-highestPeakIndex-pulsesBefore
    highNumBracket<-highestPeakIndex+pulsesAfter
    print(lowNumBracket)
    print(highNumBracket)

    if (lowNumBracket < 1) {
      shift<-abs(lowNumBracket)
      lowNumBracket<-1
      highNumBracket<-highNumBracket+shift
      paste("Caution: Bracket had to be shifted",shift,"pulses up.")
    }

    if (highNumBracket > length(allHighPeaks$value)) {
      shift<-highNumBracket-length(allHighPeaks$value)
      highNumBracket<-length(allHighPeaks$value)
      lowNumBracket<-lowNumBracket-shift
      paste("Caution: Bracket had to be shifted",shift,"pulses down.")
    }
    print(lowNumBracket)
    print(highNumBracket)

    output<-inSignal[lowNumBracket:highNumBracket]

  } else {
    paste("Argument 'what':",what,"not understood.\nPleas write \"waveform\" or \"bracket\"")
    return()
  }
  return(output)
}


if (FALSE) {
  peakRegionsIndices_high<-which(inSignal@left > quantile(inSignal@left,c(.99)))
  peakRegionsIndices_low<-which(inSignal@left < quantile(inSignal@left,c(.01)))

  diffPeakRegionsIndices_high<-diff(peakRegionsIndices_high)
  diffPeakRegionsIndices_low<-diff(peakRegionsIndices_low)

  stepVect_high<-c(peakRegionsIndices_high[1])
  stepVect_low<-c(peakRegionsIndices_low[1])

  inRowNum<-1

  for (i in 1:length(diffPeakRegionsIndices_high)) {
    if (diffPeakRegionsIndices_high[i] == 1) {
      stepVect<-c(stepVect_high,peakRegionsIndices_high[i+1])
    } else {
      stepMaxIndex_high<-stepVect_high[which.max(inSignal@left[stepVect_high])]
      stepMaxValue_high<-max(inSignal@left[stepVect_high])
      resultsDF[inRowNum,c(1,2)]<-c(stepMaxIndex_high,stepMaxValue_high)
      inRowNum<-inRowNum+1
      stepVect_high<-c(peakRegionsIndices_high[i+1])
    }
  }

  inRowNum<-1
  for (i in 1:length(diffPeakRegionsIndices_low)) {
    if (diffPeakRegionsIndices_low[i] == 1) {
      stepVect<-c(stepVect_low,peakRegionsIndices_low[i+1])
    } else {
      stepMinIndex_low<-stepVect_low[which.min(inSignal@left[stepVect_low])]
      stepMinValue_low<-min(inSignal@left[stepVect_low])
      resultsDF[inRowNum,c(3,4)]<-c(stepMinIndex_low,stepMinValue_low)
      inRowNum<-inRowNum+1
      stepVect_low<-c(peakRegionsIndices_low[i+1])
    }
  }
}


} ####################################### END DELETED SECTION

