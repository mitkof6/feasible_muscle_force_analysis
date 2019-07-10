#ifdef WIN32
#   ifdef FeasibleMuscleForceAnalysis_EXPORTS
#       define FeasibleMuscleForceAnalysis_API __declspec(dllexport)
#   else
#       define FeasibleMuscleForceAnalysis_API  __declspec(dllimport)
#   endif
#else
#   define FeasibleMuscleForceAnalysis_API
#endif // WIN32