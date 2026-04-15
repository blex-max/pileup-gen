#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "generate-pileup.hpp"

namespace py = pybind11;

// mod_gil_not_used() is omitted intentionally
// despite being the suggested default in the pybind11 docs.
// PileupReadSet::qpos_cb calls back into Python
// - callbacks into python require the GIL to be held.
PYBIND11_MODULE (htsgen, m) {
  m.doc() = "Python bindings for htsgen, a synthetic data generation library for high-throughput sequencing";

  py::enum_<BaseEvents> (m, "BaseEvents")
    .value("ref", BaseEvents::ref)
    .value("A", BaseEvents::A)
    .value("C", BaseEvents::C)
    .value("G", BaseEvents::G)
    .value("T", BaseEvents::T)
    .value("N", BaseEvents::N)
    .value("deleted", BaseEvents::del)
    .export_values();

  py::class_<EventSpec> (m, "EventSpec")
    .def (
      py::init<BaseEvents, int, std::string>(),
      py::arg("base"), py::arg("indel") = 0, py::arg("ins") = std::string{}
    )
    .def_readwrite("base", &EventSpec::base)
    .def_readwrite("indel", &EventSpec::indel)
    .def_readwrite("ins", &EventSpec::ins);

  // py::init supports multiple overloads via repeated .def(py::init<...>()),
  // so a second constructor taking a PMF (probability mass)
  // array (e.g. std::vector<double>) can be added later
  // without disturbing the python-side callback approach.
  py::class_<PileupReadSet> (m, "PileupReadSet")
    .def (
      py::init<EventSpec, std::function<uint16_t()>>(),
      py::arg("event"), py::arg("qpos_cb")
    )
    .def_readwrite("event", &PileupReadSet::event)
    .def_readwrite("qpos_cb", &PileupReadSet::qpos_cb);

  py::class_<PileupCoordinates> (m, "PileupCoordinates")
    .def (
      py::init<hts_pos_t, hts_pos_t, hts_pos_t, int32_t>(),
      py::arg("gstart"), py::arg("gend"), py::arg("gpos"), py::arg("tid")
    )
    .def_readwrite("gstart", &PileupCoordinates::gstart)
    .def_readwrite("gend", &PileupCoordinates::gend)
    .def_readwrite("gpos", &PileupCoordinates::gpos)
    .def_readwrite("tid", &PileupCoordinates::tid);

  py::class_<PileupParams> (m, "PileupParams")
    .def (
      py::init<PileupCoordinates, std::string_view, uint16_t>(),
      py::arg("coordinates"), py::arg("refseq"), py::arg("readlen")
    )
    .def_readwrite("coordinates", &PileupParams::coord)
    .def_readwrite("refseq", &PileupParams::refseq)
    .def_readwrite("readlen", &PileupParams::readlen);

}
