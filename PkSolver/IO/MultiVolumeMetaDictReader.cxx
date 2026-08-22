#include "MultiVolumeMetaDictReader.h"

#include "itkMetaDataObject.h"
#include "Exceptions.h"


MultiVolumeMetaDictReader::MultiVolumeMetaDictReader(const itk::MetaDataDictionary& dictionary) 
  : m_dictionary(dictionary) 
{}

float MultiVolumeMetaDictReader::get(const std::string& key)
{
  std::string valueString = "";
  bool readSuccess = itk::ExposeMetaData(m_dictionary, key, valueString);
  if (readSuccess) {
    try 
    {
      return std::stof(valueString);
    }
    catch (const std::invalid_argument& ia)
    {}
    throw FailedDictionaryLookup(key);
  }
}

std::vector<float> MultiVolumeMetaDictReader::getTiming()
{
  std::vector<float> triggerTimes;

  if (m_dictionary.HasKey("MultiVolume.FrameIdentifyingDICOMTagName"))
  {
    std::string tag;
    itk::ExposeMetaData(m_dictionary, "MultiVolume.FrameIdentifyingDICOMTagName", tag);
    if (m_dictionary.HasKey("MultiVolume.FrameLabels"))
    {
      // Acquisition parameters stored as text, FrameLabels are comma separated
      std::string frameLabelsString;
      itk::ExposeMetaData(m_dictionary, "MultiVolume.FrameLabels", frameLabelsString);
      std::stringstream frameLabelsStream(frameLabelsString);
      if (tag == "TriggerTime" || tag == "AcquisitionTime" || tag == "SeriesTime" || tag == "ContentTime" || tag == "Time")
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
    }
    else
    {
      itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameLabels'.");
    }
  }
  else
  {
    itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameIdentifyingDICOMTagName'.");
  }
  return triggerTimes;
}

