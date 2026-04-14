#include <pybind11/pybind11.h>

namespace py = pybind11;

int add(int i, int j) {
    return i + j;
}

PYBIND11_MODULE (htsgen, m, py::mod_gil_not_used()) {
  m.doc() = "Python bindings for htsgen, a synthetic data generation library for high-throughput sequencing";

  m.def("add", &add, "A function that adds two numbers");
}
