//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// ClassificationKNN.h
//
// Code generation for function 'ClassificationKNN'
//

#ifndef CLASSIFICATIONKNN_H
#define CLASSIFICATIONKNN_H

// Include files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
namespace coder {
class categorical;

}

// Type Definitions
namespace coder {
namespace classreg {
namespace learning {
namespace coderutils {
enum class Transform : int
{
  Logit = 0, // Default value
  Doublelogit,
  Invlogit,
  Ismax,
  Sign,
  Symmetric,
  Symmetricismax,
  Symmetriclogit,
  Identity
};

}
} // namespace learning
} // namespace classreg
class ClassificationKNN {
public:
  void init();
  void predict(const double Xin[32], categorical *labels) const;
  double OutputPrecision;
  int ClassNamesLength[6];
  classreg::learning::coderutils::Transform ScoreTransform;
  double Prior[6];
  bool ClassLogicalIndices[6];
  double Cost[36];
  double X[7232];
  unsigned int YIdx[226];
  double W[226];
  double Mu[32];
  double Sigma[32];
};

} // namespace coder

#endif
// End of code generation (ClassificationKNN.h)
