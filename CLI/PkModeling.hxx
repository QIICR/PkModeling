#ifndef __PkModeling_h
#define __PkModeling_h

#include "Configuration.h"

#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMultiThreader.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkPluginUtilities.h"

#include "itkSignalIntensityToConcentrationImageFilter.h"
#include "itkConcentrationToQuantitativeImageFilter.h"

#include "BolusArrivalTimeEstimator.h"
#include "BolusArrivalTimeEstimatorConstant.h"
#include "BolusArrivalTimeEstimatorPeakGradient.h"

#include <sstream>
#include <fstream>
#include <memory>


class PkModeling {
private:
  const Configuration cfg;

public:
  static const unsigned int VectorVolumeDimension = 3;
  typedef itk::VectorImage<float, VectorVolumeDimension>     VectorVolumeType;
  typedef itk::VectorImage<float, VectorVolumeDimension>     FloatVectorVolumeType;
  typedef VectorVolumeType::RegionType              VectorVolumeRegionType;
  typedef itk::ImageFileReader<VectorVolumeType>             VectorVolumeReaderType;
  typedef itk::ImageFileWriter<FloatVectorVolumeType>        VectorVolumeWriterType;

  static const unsigned int MaskVolumeDimension = 3;
  typedef itk::Image<unsigned short, MaskVolumeDimension>      MaskVolumeType;
  typedef itk::ImageFileReader<MaskVolumeType>                 MaskVolumeReaderType;

  typedef itk::Image<float, VectorVolumeDimension> OutputVolumeType;
  typedef itk::ImageFileWriter< OutputVolumeType> OutputVolumeWriterType;

  typedef itk::ResampleImageFilter<MaskVolumeType, MaskVolumeType> ResamplerType;
  typedef itk::NearestNeighborInterpolateImageFunction<MaskVolumeType> InterpolatorType;

  typedef itk::SignalIntensityToConcentrationImageFilter<VectorVolumeType, MaskVolumeType, FloatVectorVolumeType> ConvertFilterType;
  typedef itk::ConcentrationToQuantitativeImageFilter<FloatVectorVolumeType, MaskVolumeType, OutputVolumeType> QuantifierType;


  PkModeling(Configuration config) : cfg(config) {}
  virtual ~PkModeling() {}

#define SimpleAttributeGetMethodMacro(name, key, type)     \
type Get##name(itk::MetaDataDictionary& dictionary)           \
  {\
  type value = type(); \
  if (dictionary.HasKey(key))\
        {\
    /* attributes stored as strings */ \
    std::string valueString; \
    itk::ExposeMetaData(dictionary, key, valueString);  \
    std::stringstream valueStream(valueString); \
    valueStream >> value; \
        }\
        else\
      {\
    itkGenericExceptionMacro("Missing attribute '" key "'.");\
      }\
  return value;\
  }

  //SimpleAttributeGetMethodMacro(EchoTime, "MultiVolume.DICOM.EchoTime", float);
  SimpleAttributeGetMethodMacro(RepetitionTime, "MultiVolume.DICOM.RepetitionTime", float);
  SimpleAttributeGetMethodMacro(FlipAngle, "MultiVolume.DICOM.FlipAngle", float);

  std::vector<float> GetTiming(itk::MetaDataDictionary& dictionary)
  {
    std::vector<float> triggerTimes;

    if (dictionary.HasKey("MultiVolume.FrameIdentifyingDICOMTagName"))
    {
      std::string tag;
      itk::ExposeMetaData(dictionary, "MultiVolume.FrameIdentifyingDICOMTagName", tag);
      if (dictionary.HasKey("MultiVolume.FrameLabels"))
      {
        // Acquisition parameters stored as text, FrameLabels are comma separated
        std::string frameLabelsString;
        itk::ExposeMetaData(dictionary, "MultiVolume.FrameLabels", frameLabelsString);
        std::stringstream frameLabelsStream(frameLabelsString);
        if (tag == "TriggerTime" || tag == "AcquisitionTime" || tag == "SeriesTime" || tag == "ContentTime")
        {
          float t;
          float t0 = 0.0;
          bool first = true;
          while (frameLabelsStream >> t)
          {
            t /= 1000.0;  // convert to seconds
            if (first)
            {
              t0 = t;
              first = false;
            }
            t = t - t0;

            triggerTimes.push_back(t);
            frameLabelsStream.ignore(1); // skip the comma
          }
        }
        else
        {
          itkGenericExceptionMacro("Unrecognized frame identifying DICOM tag name " << tag);
        }
        // what other frame identification methods are there?
      }
      else
      {
        itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameLabels'.")
      }
    }
    else
    {
      itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameIdentifyingDICOMTagName'.");
    }

    return triggerTimes;
  }


  // Read an AIF from a CSV style file.  Columns are timing and concentration.
  bool GetPrescribedAIF(const std::string& fileName, std::vector<float>& timing, std::vector<float>& aif)
  {
    timing.clear();
    aif.clear();

    std::string line;
    std::ifstream csv;
    csv.open(fileName.c_str());
    if (csv.fail())
    {
      std::cout << "Cannot open file " << fileName << std::endl;
      return false;
    }

    while (!csv.eof())
    {
      getline(csv, line);

      if (line[0] == '#')
      {
        continue;
      }

      std::vector<std::string> svalues;
      splitString(line, ",", svalues);  /// from PkModelingCLP.h

      if (svalues.size() < 2)
      {
        // not enough values on the line
        continue;
      }

      // only keep the time and concentration value
      std::stringstream tstream;
      float time = -1.0, value = -1.0;

      tstream << svalues[0];
      tstream >> time;
      if (tstream.fail())
      {
        // not a float, probably the column labels, skip the row
        continue;
      }
      tstream.str("");
      tstream.clear(); // need to clear the flags since at eof of the stream

      tstream << svalues[1];
      tstream >> value;
      if (tstream.fail())
      {
        // not a float, could be column labels, skip the row
        continue;
      }
      tstream.str("");
      tstream.clear(); // need to clear the flags since at eof of the stream

      timing.push_back(time);
      aif.push_back(value);
    }

    csv.close();

    if (timing.size() > 0)
    {
      return true;
    }

    return false;
  }

  MaskVolumeType::Pointer getMaskVolumeOrNull(const std::string& maskFileName)
  {
    MaskVolumeReaderType::Pointer maskVolumeReader = MaskVolumeReaderType::New();
    MaskVolumeType::Pointer maskVolume = NULL;
    if (maskFileName != "")
    {
      maskVolumeReader->SetFileName(maskFileName.c_str());
      maskVolumeReader->Update();
      maskVolume = maskVolumeReader->GetOutput();
    }
    return maskVolume;
  }

  MaskVolumeType::Pointer getResampledMaskVolumeOrNull(const std::string& maskFileName, const VectorVolumeType::Pointer referenceVolume)
  {
    MaskVolumeType::Pointer maskVolume = getMaskVolumeOrNull(maskFileName);
    if (maskVolume.IsNotNull()) {
      ResamplerType::Pointer resampler = ResamplerType::New();
      InterpolatorType::Pointer interpolator = InterpolatorType::New();

      resampler->SetOutputDirection(referenceVolume->GetDirection());
      resampler->SetOutputSpacing(referenceVolume->GetSpacing());
      resampler->SetOutputStartIndex(referenceVolume->GetBufferedRegion().GetIndex());
      resampler->SetSize(referenceVolume->GetBufferedRegion().GetSize());
      resampler->SetOutputOrigin(referenceVolume->GetOrigin());
      resampler->SetInput(maskVolume);
      resampler->SetInterpolator(interpolator);
      resampler->Update();

      maskVolume = resampler->GetOutput();
    }
    return maskVolume;
  }

  VectorVolumeType::Pointer getVectorVolume(const std::string& volumeFileName)
  {
    VectorVolumeReaderType::Pointer multiVolumeReader = VectorVolumeReaderType::New();
    multiVolumeReader->SetFileName(volumeFileName.c_str());
    multiVolumeReader->Update();
    return multiVolumeReader->GetOutput();
  }

  std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> getBatEstimator(const std::string& mode, int defaultConstantBAT)
  {
    std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> batEstimator(new BolusArrivalTime::BolusArrivalTimeEstimatorConstant(defaultConstantBAT));
    if (mode == "PeakGradient") {
      batEstimator.reset(new BolusArrivalTime::BolusArrivalTimeEstimatorPeakGradient());
    }
    return batEstimator;
  }

  void writeVectorVolumeIfFileNameValid(std::string fileName, const FloatVectorVolumeType::Pointer outVolume, const VectorVolumeType::Pointer referenceVolume)
  {
    if (!fileName.empty())
    {
      // this line is needed to make Slicer recognize this as a VectorVolume and not a MultiVolume
      outVolume->SetMetaDataDictionary(referenceVolume->GetMetaDataDictionary());

      VectorVolumeWriterType::Pointer multiVolumeWriter = VectorVolumeWriterType::New();
      multiVolumeWriter->SetFileName(fileName.c_str());
      multiVolumeWriter->SetInput(outVolume);
      multiVolumeWriter->SetUseCompression(1);
      multiVolumeWriter->Update();
    }
  }

  void writeVolumeIfFileNameValid(std::string fileName, const OutputVolumeType::Pointer outVolume)
  {
    if (!fileName.empty())
    {
      OutputVolumeWriterType::Pointer volumeWriter = OutputVolumeWriterType::New();
      volumeWriter->SetInput(outVolume);
      volumeWriter->SetFileName(fileName.c_str());
      volumeWriter->SetUseCompression(1);
      volumeWriter->Update();
    }
  }

  void saveRequestedOutputs(const ConvertFilterType::Pointer signalToConcentrationsConverter, 
                            const QuantifierType::Pointer concentrationsToQuantitativeImageFilter,
                            const VectorVolumeType::Pointer inputVectorVolume)
  {
    writeVectorVolumeIfFileNameValid(cfg.OutputConcentrationsImageFileName, signalToConcentrationsConverter->GetOutput(), inputVectorVolume);
    writeVectorVolumeIfFileNameValid(cfg.OutputFittedDataImageFileName, concentrationsToQuantitativeImageFilter->GetFittedDataOutput(), inputVectorVolume);

    writeVolumeIfFileNameValid(cfg.OutputKtransFileName, concentrationsToQuantitativeImageFilter->GetKTransOutput());
    writeVolumeIfFileNameValid(cfg.OutputVeFileName, concentrationsToQuantitativeImageFilter->GetVEOutput());
    writeVolumeIfFileNameValid(cfg.OutputMaxSlopeFileName, concentrationsToQuantitativeImageFilter->GetMaxSlopeOutput());
    writeVolumeIfFileNameValid(cfg.OutputAUCFileName, concentrationsToQuantitativeImageFilter->GetAUCOutput());
    writeVolumeIfFileNameValid(cfg.OutputRSquaredFileName, concentrationsToQuantitativeImageFilter->GetRSquaredOutput());
    writeVolumeIfFileNameValid(cfg.OutputBolusArrivalTimeImageFileName, concentrationsToQuantitativeImageFilter->GetBATOutput());
    writeVolumeIfFileNameValid(cfg.OutputOptimizerDiagnosticsImageFileName, concentrationsToQuantitativeImageFilter->GetOptimizerDiagnosticsOutput());

    if (cfg.ComputeFpv) {
      writeVolumeIfFileNameValid(cfg.OutputFpvFileName, concentrationsToQuantitativeImageFilter->GetFPVOutput());
    }
  }


  int execute()
  {
    VectorVolumeType::Pointer inputVectorVolume = getVectorVolume(cfg.InputFourDImageFileName);
    std::unique_ptr<BolusArrivalTime::BolusArrivalTimeEstimator> batEstimator = getBatEstimator(cfg.BATCalculationMode, cfg.ConstantBAT);

    std::vector<float> Timing;
    float FAValue = 0.0;
    float TRValue = 0.0;
    try
    {
      Timing = GetTiming(inputVectorVolume->GetMetaDataDictionary());
      FAValue = GetFlipAngle(inputVectorVolume->GetMetaDataDictionary());
      TRValue = GetRepetitionTime(inputVectorVolume->GetMetaDataDictionary());
    }
    catch (itk::ExceptionObject &exc)
    {
      itkGenericExceptionMacro(<< exc.GetDescription()
        << " Image " << cfg.InputFourDImageFileName.c_str()
        << " does not contain sufficient attributes to support algorithms.");
    }

    //Read masks
    MaskVolumeType::Pointer aifMaskVolume = getMaskVolumeOrNull(cfg.AIFMaskFileName);
    MaskVolumeType::Pointer T1MapVolume = getMaskVolumeOrNull(cfg.T1MapFileName);
    MaskVolumeType::Pointer roiMaskVolume = getResampledMaskVolumeOrNull(cfg.ROIMaskFileName, inputVectorVolume);

    //Read prescribed aif
    bool usingPrescribedAIF = false;
    std::vector<float> prescribedAIFTiming;
    std::vector<float> prescribedAIF;
    if (cfg.PrescribedAIFFileName != "")
    {
      usingPrescribedAIF = GetPrescribedAIF(cfg.PrescribedAIFFileName, prescribedAIFTiming, prescribedAIF);
    }

    if (cfg.AIFMaskFileName == "" && !usingPrescribedAIF && !cfg.UsePopulationAIF)
    {
      itkGenericExceptionMacro(<< "Either a mask localizing the region over which to "
        << "calculate the arterial input function or a prescribed "
        << "arterial input function must be specified.");
    }

    /////////////////////////////// PROCESSING /////////////////////

    //Convert to concentration values
    ConvertFilterType::Pointer signalToConcentrationsConverter = ConvertFilterType::New();
    signalToConcentrationsConverter->SetInput(inputVectorVolume);

    if (!usingPrescribedAIF && !cfg.UsePopulationAIF)
    {
      signalToConcentrationsConverter->SetAIFMask(aifMaskVolume);
    }

    if (cfg.ROIMaskFileName != "")
    {
      signalToConcentrationsConverter->SetROIMask(roiMaskVolume);
    }

    signalToConcentrationsConverter->SetT1PreBlood(cfg.T1PreBloodValue);
    signalToConcentrationsConverter->SetT1PreTissue(cfg.T1PreTissueValue);
    signalToConcentrationsConverter->SetTR(TRValue);
    signalToConcentrationsConverter->SetFA(FAValue);
    signalToConcentrationsConverter->SetBatEstimator(batEstimator.get());
    signalToConcentrationsConverter->SetRGD_relaxivity(cfg.RelaxivityValue);
    signalToConcentrationsConverter->SetS0GradThresh(cfg.S0GradValue);

    if (cfg.T1MapFileName != "")
    {
      signalToConcentrationsConverter->SetT1Map(T1MapVolume);
    }

    itk::PluginFilterWatcher watchConverter(signalToConcentrationsConverter, "Concentrations", cfg.CLPProcessInformation, 1.0 / 20.0, 0.0);
    signalToConcentrationsConverter->Update();

    //Calculate parameters
    QuantifierType::Pointer concentrationsToQuantitativeImageFilter = QuantifierType::New();
    concentrationsToQuantitativeImageFilter->SetInput(signalToConcentrationsConverter->GetOutput());
    if (usingPrescribedAIF)
    {
      concentrationsToQuantitativeImageFilter->SetPrescribedAIF(prescribedAIFTiming, prescribedAIF);
      concentrationsToQuantitativeImageFilter->UsePrescribedAIFOn();
    }
    else if (cfg.UsePopulationAIF)
    {
      concentrationsToQuantitativeImageFilter->UsePopulationAIFOn();
      concentrationsToQuantitativeImageFilter->SetAIFMask(aifMaskVolume);
    }
    else
    {
      concentrationsToQuantitativeImageFilter->SetAIFMask(aifMaskVolume);
    }

    concentrationsToQuantitativeImageFilter->SetAUCTimeInterval(cfg.AUCTimeInterval);
    concentrationsToQuantitativeImageFilter->SetTiming(Timing);
    concentrationsToQuantitativeImageFilter->SetfTol(cfg.FTolerance);
    concentrationsToQuantitativeImageFilter->SetgTol(cfg.GTolerance);
    concentrationsToQuantitativeImageFilter->SetxTol(cfg.XTolerance);
    concentrationsToQuantitativeImageFilter->Setepsilon(cfg.Epsilon);
    concentrationsToQuantitativeImageFilter->SetmaxIter(cfg.MaxIter);
    concentrationsToQuantitativeImageFilter->Sethematocrit(cfg.Hematocrit);
    concentrationsToQuantitativeImageFilter->SetBatEstimator(batEstimator.get());
    if (cfg.ROIMaskFileName != "")
    {
      concentrationsToQuantitativeImageFilter->SetROIMask(roiMaskVolume);
    }

    if (cfg.ComputeFpv)
    {
      concentrationsToQuantitativeImageFilter->SetModelType(itk::LMCostFunction::TOFTS_3_PARAMETER);
    }
    else
    {
      concentrationsToQuantitativeImageFilter->SetModelType(itk::LMCostFunction::TOFTS_2_PARAMETER);
    }

    itk::PluginFilterWatcher watchQuantifier(concentrationsToQuantitativeImageFilter, "Quantifying", cfg.CLPProcessInformation, 19.0 / 20.0, 1.0 / 20.0);
    concentrationsToQuantitativeImageFilter->Update();

    ///////////////////////////////////// OUTPUT ////////////////////

    saveRequestedOutputs(signalToConcentrationsConverter, concentrationsToQuantitativeImageFilter, inputVectorVolume);

    return EXIT_SUCCESS;
  }

};


#endif