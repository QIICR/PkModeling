#include "itkTestMain.h"

// STD includes
#include <iostream>

#if defined(_WIN32) && !defined(MODULE_STATIC)
#define MODULE_IMPORT __declspec(dllimport)
#else
#define MODULE_IMPORT
#endif

extern "C" MODULE_IMPORT int ModuleEntryPoint(int, char *[]);

int ModuleEntryPointExpectFail(int argc, char * argv[])
{
  int errorCode = ModuleEntryPoint(argc, argv);
  if (errorCode != 0) {
    std::cout << "Module execution failed, TEST PASSED!" << std::endl;
    return 0;
  }
  else {
    std::cout << "Module execution passed, TEST FAILED!" << std::endl;
    return 1;
  }
}

int DoNothingAndPass(int argc, char * argv[])
{
  return 0;
}

void RegisterTests()
{
  StringToTestFunctionMap["ModuleEntryPoint"] = ModuleEntryPoint;
  StringToTestFunctionMap["ModuleEntryPointExpectFail"] = ModuleEntryPointExpectFail;
  StringToTestFunctionMap["DoNothingAndPass"] = DoNothingAndPass;
}
