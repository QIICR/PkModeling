/*=auto=========================================================================

  Portions (c) Copyright 2009 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: pk_solver.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.13 $

  =========================================================================auto=*/

#ifndef PkSolver_h_
#define PkSolver_h_

#include "itkLevenbergMarquardtOptimizer.h"
#include <math.h>
#include "itkArray.h"
#include <string>
#include <exception>
#include "BAT/BolusArrivalTimeEstimator.h"


// work around compile error on Win
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

class NoSignalException : public std::exception
{
  virtual const char* what() const throw()
  {
    return "The provided signal vector/array was emtpy, no computation can be performed.";
  }
};


// codes defined in ITKv4 vnl_levenberg_marquardt.cxx:386
enum OptimizerDiagnosticCodes
{
  ERROR_FAILURE = 0,
  ERROR_DODGY_INPUT = 1,
  CONVERGED_FTOL = 2,
  CONVERGED_XTOL = 3,
  CONVERGED_XFTOL = 4,
  CONVERGED_GTOL = 5,
  TOO_MANY_ITERATIONS = 6,
  FAILED_FTOL_TOO_SMALL = 7,
  FAILED_XTOL_TOO_SMALL = 8,
  FAILED_GTOL_TOO_SMALL = 9,
  FAILED_UNKNOWN = 10,
  // next are the masks that are specific to the PK modeling process
  FAILED_NOMATCH = 11, // optimizer failed, but diagnostics string was not recognized
  KTRANS_CLAMPED = 0x10, // = 16 Ktrans was clamped to [0..5]
  VE_CLAMPED = 0x20, // = 32 Ve was clamped to [0..1]
  BAT_DETECTION_FAILED = 0x30, // = 48 BAT detection procedure failed
  BAT_BEFORE_AIF_BAT = 0x40 // = 64 BAT at the voxel was before AIF BAT
};

const std::string OptimizerDiagnosticStrings[] =
{
  "failure in leastsquares function",
  "lmdif dodgy input",
  "converged to ftol",
  "converged to xtol",
  "converged nicely",
  "converged via gtol",
  "too many iterations",
  "ftol is too small",
  "xtol is too small",
  "gtol is too small",
  "unkown info code"
};

const unsigned NumOptimizerDiagnosticCodes = 11;

namespace itk
{

  class LMCostFunction : public itk::MultipleValuedCostFunction
  {
  public:
    typedef LMCostFunction                    Self;
    typedef itk::MultipleValuedCostFunction   Superclass;
    typedef itk::SmartPointer<Self>           Pointer;
    typedef itk::SmartPointer<const Self>     ConstPointer;
    itkNewMacro(Self);

    enum { SpaceDimension = 2 };
    unsigned int RangeDimension;

    enum ModelType { TOFTS_2_PARAMETER = 1, TOFTS_3_PARAMETER };

    typedef Superclass::ParametersType              ParametersType;
    typedef Superclass::DerivativeType              DerivativeType;
    typedef Superclass::MeasureType                 MeasureType, ArrayType;
    typedef Superclass::ParametersValueType         ValueType;


    float m_Hematocrit;

    int m_ModelType;

    LMCostFunction()
    {
    }

    void SetHematocrit(float hematocrit)
    {
      m_Hematocrit = hematocrit;
    }

    void SetModelType(int model)
    {
      m_ModelType = model;
    }

    void SetNumberOfValues(unsigned int NumberOfValues)
    {
      RangeDimension = NumberOfValues;
    }

    void SetCb(const float* cb, int sz) //BloodConcentrationCurve.
    {
      Cb.set_size(sz);
      for (int i = 0; i < sz; ++i)
        Cb[i] = cb[i];
      //std::cout << "Cb: " << Cb << std::endl;
    }


    void SetCv(const float* cv, int sz) //Self signal Y
    {
      Cv.set_size(sz);
      for (int i = 0; i < sz; ++i)
        Cv[i] = cv[i];
      //std::cout << "Cv: " << Cv << std::endl;
    }

    void SetTime(const float* cx, int sz) //Self signal X
    {
      Time.set_size(sz);
      for (int i = 0; i < sz; ++i)
        Time[i] = cx[i];
      //std::cout << "Time: " << Time << std::endl;
    }

    MeasureType GetValue(const ParametersType & parameters) const
    {
      MeasureType measure(RangeDimension);

      ValueType Ktrans = parameters[0];
      ValueType Ve = parameters[1];

      ArrayType VeTerm;
      VeTerm = -Ktrans / Ve*Time;
      ValueType deltaT = Time(1) - Time(0);

      if (m_ModelType == TOFTS_3_PARAMETER)
      {
        ValueType f_pv = parameters[2];
        measure = Cv - (1 / (1.0 - m_Hematocrit)*(Ktrans*deltaT*Convolution(Cb, Exponential(VeTerm)) + f_pv*Cb));
      }
      else if (m_ModelType == TOFTS_2_PARAMETER)
      {
        measure = Cv - (1 / (1.0 - m_Hematocrit)*(Ktrans*deltaT*Convolution(Cb, Exponential(VeTerm))));
      }

      return measure;
    }

    MeasureType GetFittedFunction(const ParametersType & parameters) const
    {
      MeasureType measure(RangeDimension);

      ValueType Ktrans = parameters[0];
      ValueType Ve = parameters[1];

      ArrayType VeTerm;
      VeTerm = -Ktrans / Ve*Time;
      ValueType deltaT = Time(1) - Time(0);

      if (m_ModelType == TOFTS_3_PARAMETER)
      {
        ValueType f_pv = parameters[2];
        measure = 1 / (1.0 - m_Hematocrit)*(Ktrans*deltaT*Convolution(Cb, Exponential(VeTerm)) + f_pv*Cb);
      }
      else if (m_ModelType == TOFTS_2_PARAMETER)
      {
        measure = 1 / (1.0 - m_Hematocrit)*(Ktrans*deltaT*Convolution(Cb, Exponential(VeTerm)));
      }

      return measure;
    }

    //Not going to be used
    void GetDerivative(const ParametersType & /* parameters*/,
      DerivativeType  & /*derivative*/) const
    {
    }

    unsigned int GetNumberOfParameters(void) const
    {
      if (m_ModelType == TOFTS_2_PARAMETER)
      {
        return 2;
      }
      else // if(m_ModelType == TOFTS_3_PARAMETER)
      {
        return 3;
      }
    }

    unsigned int GetNumberOfValues(void) const
    {
      return RangeDimension;
    }

  protected:
    virtual ~LMCostFunction(){}
  private:

    ArrayType Cv, Cb, Time;

    ArrayType Convolution(ArrayType X, ArrayType Y) const
    {
      const unsigned int n = X.size();
      ArrayType Z(n, 0.0);
      for (unsigned int i = 0; i < n; ++i)
        {
        ValueType acc = 0;
        const unsigned int jmin = (i + 1 > Y.size()) ? (i + 1 - Y.size()) : 0u;
        for (unsigned int j = jmin; j <= i; ++j)
          {
          acc += X[j] * Y[i - j];
          }
        Z[i] = acc;
        }
      return Z;
    };

    ArrayType Exponential(ArrayType X) const
    {
      ArrayType Z;
      Z.set_size(X.size());
      for (unsigned int i = 0; i < X.size(); i++)
      {
        Z[i] = exp(X(i));
      }
      return Z;
    };

    int constraintFunc(ValueType x) const
    {
      if (x < 0 || x>1)
        return 1;
      else
        return 0;
    };


  };

  class CommandIterationUpdateLevenbergMarquardt : public itk::Command
  {
  public:
    typedef  CommandIterationUpdateLevenbergMarquardt   Self;
    typedef  itk::Command                               Superclass;
    typedef itk::SmartPointer<Self>                     Pointer;
    itkNewMacro(Self);
  protected:
    CommandIterationUpdateLevenbergMarquardt()
    {
      m_IterationNumber = 0;
    }
    virtual ~CommandIterationUpdateLevenbergMarquardt(){}
  public:
    typedef itk::LevenbergMarquardtOptimizer   OptimizerType;
    typedef   const OptimizerType   *          OptimizerPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute((const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      //std::cout << "Observer::Execute() " << std::endl;
      OptimizerPointer optimizer =
        dynamic_cast<OptimizerPointer>(object);
      if (m_FunctionEvent.CheckEvent(&event))
      {
        // std::cout << m_IterationNumber++ << "   ";
        // std::cout << optimizer->GetCachedValue() << "   ";
        // std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
      }
      else if (m_GradientEvent.CheckEvent(&event))
      {
        std::cout << "Gradient " << optimizer->GetCachedDerivative() << "   ";
      }

    }
  private:
    unsigned long m_IterationNumber;

    itk::FunctionEvaluationIterationEvent m_FunctionEvent;
    itk::GradientEvaluationIterationEvent m_GradientEvent;
  };

  // returns diagnostic error code from the VNL optimizer,
  //  as defined by OptimizerDiagnosticCodes, and masked to indicate
  //  wheather Ktrans or Ve were clamped.
  unsigned pk_solver(int signalSize, const float* timeAxis,
    const float* PixelConcentrationCurve, const float* BloodConcentrationCurve,
    float& Ktrans, float& Ve, float& Fpv,
    float fTol, float gTol, float xTol,
    float epsilon, int maxIter, float hematocrit,
    itk::LevenbergMarquardtOptimizer* optimizer,
    LMCostFunction* costFunction,
    int modelType = itk::LMCostFunction::TOFTS_2_PARAMETER,
    const BolusArrivalTime::BolusArrivalTimeEstimator* batEstimator = NULL);

  void pk_report();
  void pk_clear();

  bool convert_signal_to_concentration(unsigned int signalSize,
    const float* SignalIntensityCurve,
    float T1, float TR, float FA,
    float* concentration,
    float relaxivity = 4.9E-3f,
    float s0 = -1.0f,
    float S0GradThresh = 15.0f);

  float area_under_curve(int signalSize, const float* timeAxis, const float* concentration, int BATIndex, float aucTimeInterval);

  float intergrate(float* yValues, float * xValues, int size);

  void compute_derivative(int signalSize, const float* SingnalY, float* YDeriv);

  void compute_derivative_forward(int signalSize, const float* SignalY, float* YDeriv);

  void compute_derivative_backward(int signalSize, const float* SignalY, float* YDeriv);

  void compute_gradient(int signalSize, const float* SignalY, float* SignalGradient);

  void compute_gradient_forward(int signalSize, const float* SignalY, float* SignalGradient);

  void compute_gradient_backward(int signalSize, const float* SignalY, float* SignalGradient);

  float compute_s0_using_sumsignal_properties(int signalSize, const float* SignalY,
    const short* lowGradIndex, int FirstPeak);

  float compute_s0_individual_curve(int signalSize, const float* SignalY, float S0GradThresh, const BolusArrivalTime::BolusArrivalTimeEstimator* batEstimator);

};

#endif
