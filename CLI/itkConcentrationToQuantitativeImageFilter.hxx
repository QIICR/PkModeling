#ifndef _itkConcentrationToQuantitativeImageFilter_hxx
#define _itkConcentrationToQuantitativeImageFilter_hxx
#endif

#include "itkImageRegionConstIterator.h"
#include <vnl/algo/vnl_levenberg_marquardt.h>
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "vnl/vnl_math.h"

// work around compile error on Windows
#define M_PI 3.1415926535897932384626433832795

#include "itkConcentrationToQuantitativeImageFilter.h"

#include "AIF/ArterialInputFunction.h"

namespace itk
{

  template <class TInputImage, class TMaskImage, class TOutputImage>
  ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>::ConcentrationToQuantitativeImageFilter()
  {
    m_T1Pre = 0.0f;
    m_TR = 0.0f;
    m_FA = 0.0f;
    m_RGD_relaxivity = 4.9E-3f;
    m_S0GradThresh = 15.0f;
    m_fTol = 1e-4f;
    m_gTol = 1e-4f;
    m_xTol = 1e-5f;
    m_epsilon = 1e-9f;
    m_maxIter = 200;
    m_hematocrit = 0.4f;
    m_aifAUC = 0.0f;
    m_AIFBATIndex = 0;
    m_ModelType = itk::LMCostFunction::TOFTS_2_PARAMETER;
    this->Superclass::SetNumberOfRequiredInputs(1);
    this->Superclass::SetNthOutput(1, static_cast<TOutputImage*>(this->MakeOutput(1).GetPointer()));  // Ve
    this->Superclass::SetNthOutput(2, static_cast<TOutputImage*>(this->MakeOutput(2).GetPointer()));  // FPV
    this->Superclass::SetNthOutput(3, static_cast<TOutputImage*>(this->MakeOutput(3).GetPointer()));  // Max slope
    this->Superclass::SetNthOutput(4, static_cast<TOutputImage*>(this->MakeOutput(4).GetPointer()));  // AUC
    this->Superclass::SetNthOutput(5, static_cast<TOutputImage*>(this->MakeOutput(5).GetPointer()));  // R^2
    this->Superclass::SetNthOutput(6, static_cast<TOutputImage*>(this->MakeOutput(6).GetPointer()));  // BAT
    this->Superclass::SetNthOutput(7, static_cast<VectorVolumeType*>(this->MakeOutput(7).GetPointer())); // fitted
    this->Superclass::SetNthOutput(8, static_cast<TOutputImage*>(this->MakeOutput(8).GetPointer())); // diagnostics
    this->DynamicMultiThreadingOff();
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  typename ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>::DataObjectPointer
    ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
    ::MakeOutput(DataObjectPointerArraySizeType idx)
  {
    if (idx == 7)
    {
      return VectorVolumeType::New().GetPointer();
    }
    else
    {
      return TOutputImage::New().GetPointer();
    }
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  void
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::SetAIF(const ArterialInputFunction* aif)
  {
    if (aif->getSignalSize() < 2)
    {
      itkExceptionMacro(<< "AIF must contain at least two time points");
    }
    m_aif = aif;
  }

  // Set 3D ROI mask as third input
  template< class TInputImage, class TMaskImage, class TOutputImage >
  void
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::SetROIMask(const TMaskImage* volume)
  {
    this->SetNthInput(2, const_cast<TMaskImage*>(volume));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  const TMaskImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetROIMask() const
  {
    return dynamic_cast<const TMaskImage *>(this->ProcessObject::GetInput(2));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetKTransOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(0));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetVEOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(1));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetFPVOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(2));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetMaxSlopeOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(3));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetAUCOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(4));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetRSquaredOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(5));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetBATOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(6));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TInputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetFittedDataOutput()
  {
    return dynamic_cast<TInputImage *>(this->ProcessObject::GetOutput(7));
  }

  template< class TInputImage, class TMaskImage, class TOutputImage >
  TOutputImage*
    ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
    ::GetOptimizerDiagnosticsOutput()
  {
    return dynamic_cast<TOutputImage *>(this->ProcessObject::GetOutput(8));
  }

  template <class TInputImage, class TMaskImage, class TOutputImage>
  void
    ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
    ::BeforeThreadedGenerateData()
  {
    const VectorVolumeType* inputVectorVolume = this->GetInput();

    std::cout << "Model type: " << m_ModelType << std::endl;

    int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

    // Some of the outputs are optional and may not be calculated.
    // Let's initialize those to all zeros
    OutputVolumeType *fpv = this->GetFPVOutput();
    fpv->FillBuffer(0.0);

    // get AIF signal
    m_AIF = m_aif->getSignalValues();
    
    // Compute the bolus arrival time
    m_AIFBATIndex = m_batEstimator->getBATIndex(m_AIF.size(), &m_AIF[0]);

    // Compute the area under the curve for the AIF
    m_aifAUC = area_under_curve(timeSize, &m_Timing[0], &m_AIF[0], m_AIFBATIndex, m_AUCTimeInterval);
  }

  template <class TInputImage, class TMaskImage, class TOutputImage>
  void
    ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
#if ITK_VERSION_MAJOR < 4
    ::ThreadedGenerateData( const OutputVolumeRegionType & outputRegionForThread, int threadId )
#else
    ::ThreadedGenerateData(const OutputVolumeRegionType& outputRegionForThread, ThreadIdType threadId)
#endif
  {
    VectorVoxelType vectorVoxel, fittedVectorVoxel;

    float tempFpv = 0.0f;
    float tempKtrans = 0.0f;
    float tempVe = 0.0f;
    float tempMaxSlope = 0.0f;
    float tempAUC = 0.0f;
    int   BATIndex = 0;

    const VectorVolumeType* inputVectorVolume = this->GetInput();

    VectorVolumeConstIterType inputVectorVolumeIter(inputVectorVolume, outputRegionForThread);
    OutputVolumeIterType ktransVolumeIter(this->GetKTransOutput(), outputRegionForThread);
    OutputVolumeIterType veVolumeIter(this->GetVEOutput(), outputRegionForThread);
    typename VectorVolumeType::Pointer fitted = this->GetFittedDataOutput();
    VectorVolumeIterType fittedVolumeIter(fitted, outputRegionForThread);
    OutputVolumeIterType diagVolumeIter(this->GetOptimizerDiagnosticsOutput(), outputRegionForThread);

    MaskVolumeConstIterType roiMaskVolumeIter;
    if (this->GetROIMask())
    {
      roiMaskVolumeIter = MaskVolumeConstIterType(this->GetROIMask(), outputRegionForThread);
    }

    OutputVolumeIterType fpvVolumeIter;
    if (m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
    {
      fpvVolumeIter = OutputVolumeIterType(this->GetFPVOutput(), outputRegionForThread);
    }
    OutputVolumeIterType maxSlopeVolumeIter(this->GetMaxSlopeOutput(), outputRegionForThread);
    OutputVolumeIterType aucVolumeIter(this->GetAUCOutput(), outputRegionForThread);
    OutputVolumeIterType rsqVolumeIter(this->GetRSquaredOutput(), outputRegionForThread);
    OutputVolumeIterType batVolumeIter(this->GetBATOutput(), outputRegionForThread);

    //set up optimizer and cost function
    itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New();
    LMCostFunction::Pointer                   costFunction = LMCostFunction::New();
    int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

    std::vector<float> timeMinute;
    timeMinute = m_Timing;
    for (unsigned int i = 0; i < timeMinute.size(); i++)
    {
      timeMinute[i] = m_Timing[i] / 60.0;
    }

    ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

    VectorVoxelType shiftedVectorVoxel(timeSize);
    int shift;
    unsigned int shiftStart = 0, shiftEnd = 0;
    bool success = true;
    while (!ktransVolumeIter.IsAtEnd())
    {
      success = true;
      float optimizerErrorCode = -1;
      tempKtrans = tempVe = tempFpv = tempMaxSlope = tempAUC = 0.0;
      BATIndex = 0;

      if (!this->GetROIMask() || (this->GetROIMask() && roiMaskVolumeIter.Get()))
      {
        vectorVoxel = inputVectorVolumeIter.Get();
        fittedVectorVoxel = inputVectorVolumeIter.Get();
        // dump a specific voxel
        // std::cout << "VectorVoxel = " << vectorVoxel;

        // Compute the bolus arrival time and the max slope parameter
        if (success)
        {
          try {
            BATIndex = m_batEstimator->getBATIndex(timeSize, &vectorVoxel[0], &tempMaxSlope);
          }
          catch (...)
          {
            success = false;
            optimizerErrorCode = BAT_DETECTION_FAILED;
          }
        }


        // Shift the current time course to align with the BAT of the AIF
        // (note the sense of the shift)
        if (success)
        {
          batVolumeIter.Set(BATIndex);
          shift = m_AIFBATIndex - BATIndex;
          shiftedVectorVoxel.Fill(0.0);
          if (shift <= 0)
          {
            // AIF BAT before current BAT, should always be the case
            shiftStart = 0;
            shiftEnd = vectorVoxel.Size() + shift;
          }
          else
          {
            success = false;
            optimizerErrorCode = BAT_BEFORE_AIF_BAT;
          }
        }
        if (success)
        {
          for (unsigned int i = shiftStart; i < shiftEnd; ++i)
          {
            shiftedVectorVoxel[i] = vectorVoxel[i - shift];
          }
        }

        // Calculate parameter ktrans, ve, and fpv
        double rSquared = 0.0;
        if (success)
        {
          optimizerErrorCode = pk_solver(timeSize, &timeMinute[0],
            const_cast<float *>(shiftedVectorVoxel.GetDataPointer()),
            &m_AIF[0],
            tempKtrans, tempVe, tempFpv,
            m_fTol, m_gTol, m_xTol,
            m_epsilon, m_maxIter, m_hematocrit,
            optimizer, costFunction, m_ModelType, m_batEstimator);

          itk::LMCostFunction::ParametersType param(3);
          param[0] = tempKtrans; param[1] = tempVe;
          if (m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
          {
            param[2] = tempFpv;
          }
          itk::LMCostFunction::MeasureType measure =
            costFunction->GetFittedFunction(param);
          for (size_t i = 0; i < fittedVectorVoxel.GetSize(); i++)
          {
            fittedVectorVoxel[i] = measure[i];
          }

          // Shift the current time course to align with the BAT of the AIF
          // (note the sense of the shift)
          shiftedVectorVoxel.Fill(0.0);
          if (shift <= 0)
          {
            // AIF BAT before current BAT, should always be the case
            shiftStart = shift*-1.;
            shiftEnd = vectorVoxel.Size();
            for (unsigned int i = shiftStart; i < shiftEnd; ++i)
            {
              shiftedVectorVoxel[i] = fittedVectorVoxel[i + shift];
            }
          }

          fittedVolumeIter.Set(shiftedVectorVoxel);

          // Only keep the estimated values if the optimization produced a good answer
          // Check R-squared:
          //   R2 = 1 - SSerr / SStot
          // where
          //   SSerr = \sum (y_i - f_i)^2
          //   SStot = \sum (y_i - \bar{y})^2
          //
          // Note: R-squared is not a good metric for nonlinear function
          // fitting. R-squared values are not bound between [0,1] when
          // fitting nonlinear functions.

          // SSerr we can get easily from the optimizer
          double rms = optimizer->GetOptimizer()->get_end_error();
          double SSerr = rms*rms*shiftedVectorVoxel.GetSize();

          // if we couldn't get rms from the optimizer, we would calculate SSerr ourselves
          // LMCostFunction::MeasureType residuals = costFunction->GetValue(optimizer->GetCurrentPosition());
          // double SSerr = 0.0;
          // for (unsigned int i=0; i < residuals.size(); ++i)
          //   {
          //   SSerr += (residuals[i]*residuals[i]);
          //   }

          // SStot we need to calculate
          double sumSquared = 0.0;
          double sum = 0.0;
          for (unsigned int i = 0; i < shiftedVectorVoxel.GetSize(); ++i)
          {
            sum += shiftedVectorVoxel[i];
            sumSquared += (shiftedVectorVoxel[i] * shiftedVectorVoxel[i]);
          }
          double SStot = sumSquared - sum*sum / (double)shiftedVectorVoxel.GetSize();

          rSquared = 1.0 - (SSerr / SStot);

          /*
          double rSquaredThreshold = 0.15;
          if (rSquared < rSquaredThreshold)
          {
          success = false;
          }
          */
        }
        // Calculate parameter AUC, normalized by AIF AUC
        if (success)
        {
          tempAUC =
            (area_under_curve(timeSize, &m_Timing[0], const_cast<float *>(shiftedVectorVoxel.GetDataPointer()), BATIndex, m_AUCTimeInterval)) / m_aifAUC;
        }

        // If we were successful, save the estimated values, otherwise
        // default to zero
        if (success)
        {
          ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(tempKtrans));
          veVolumeIter.Set(static_cast<OutputVolumePixelType>(tempVe));
          maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(tempMaxSlope));
          aucVolumeIter.Set(static_cast<OutputVolumePixelType>(tempAUC));
          if (m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
          {
            fpvVolumeIter.Set(static_cast<OutputVolumePixelType>(tempFpv));
          }
        }
        else
        {
          ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(0));
          veVolumeIter.Set(static_cast<OutputVolumePixelType>(0));
          maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(0));
          aucVolumeIter.Set(static_cast<OutputVolumePixelType>(0));

          batVolumeIter.Set(-1);
          rsqVolumeIter.Set(0.0);
          shiftedVectorVoxel.Fill(0.0);
          fittedVolumeIter.Set(shiftedVectorVoxel);

        }

        // RSquared output volume is always written
        rsqVolumeIter.Set(rSquared);
      }
      else
      {
        ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(0));
        veVolumeIter.Set(static_cast<OutputVolumePixelType>(0));
        maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(0));
        aucVolumeIter.Set(static_cast<OutputVolumePixelType>(0));

        batVolumeIter.Set(-1);
        rsqVolumeIter.Set(0.0);
        shiftedVectorVoxel.Fill(0.0);
        fittedVolumeIter.Set(shiftedVectorVoxel);
      }

      ++ktransVolumeIter;
      ++veVolumeIter;
      ++maxSlopeVolumeIter;
      ++aucVolumeIter;
      ++rsqVolumeIter;
      ++batVolumeIter;
      ++inputVectorVolumeIter;
      ++fittedVolumeIter;

      if (this->GetROIMask())
      {
        ++roiMaskVolumeIter;
      }

      if (m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
      {
        ++fpvVolumeIter;
      }

      diagVolumeIter.Set(static_cast<OutputVolumePixelType>(optimizerErrorCode));
      ++diagVolumeIter;

      progress.CompletedPixel();
    }
  }

  template <class TInputImage, class TMaskImage, class TOutputImage>
  void ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
    ::SetTiming(const std::vector<float>& inputTiming)
  {
    m_Timing = inputTiming;
  }

  template <class TInputImage, class TMaskImage, class TOutputImage>
  const std::vector<float>& ConcentrationToQuantitativeImageFilter < TInputImage, TMaskImage, TOutputImage >
    ::GetTiming()
  {
    return m_Timing;
  }

  template <class TInputImage, class TMaskImage, class TOutputImage>
  void ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
    ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "Function tolerance: " << m_fTol << std::endl;
    os << indent << "Gradient tolerance: " << m_gTol << std::endl;
    os << indent << "Parameter tolerance: " << m_xTol << std::endl;
    os << indent << "Epsilon: " << m_epsilon << std::endl;
    os << indent << "Maximum number of iterations: " << m_maxIter << std::endl;
    os << indent << "Hematocrit: " << m_hematocrit << std::endl;
  }

} // end namespace itk
