#include "femparameters.h"

// Forward declarations
namespace libMesh {
class UnstructuredMesh;
}
class FEMParameters;

void build_domain (UnstructuredMesh &mesh, FEMParameters &param);
