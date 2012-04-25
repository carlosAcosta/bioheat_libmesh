// C++ include files that we need
#include <parameters.h>
#include <point.h>
#include <vector_value.h>

void read_initial_parameters();
void finish_initialization();

Number initial_value(const Point& p,
                     const Parameters&,
                     const std::string&,
                     const std::string&);

Gradient initial_grad(const Point& p,
                      const Parameters&,
                      const std::string&,
                      const std::string&);
