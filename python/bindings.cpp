#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "generate-pileup.hpp"
#include "accessors.hpp"

namespace py = pybind11;

// mod_gil_not_used() is omitted intentionally
// despite being the suggested default in the pybind11 docs.
// PileupReadSet::qpos_cb calls back into Python
// - my understanding is that callbacks into python
// require the GIL to be held.
PYBIND11_MODULE (htsgen, m) {
  m.doc() = "Python bindings for htsgen, a synthetic data generation library for high-throughput sequencing";

  py::enum_<BaseEvents> (m, "BaseEvents")
    .value("ref", BaseEvents::ref)
    .value("A", BaseEvents::A)
    .value("C", BaseEvents::C)
    .value("G", BaseEvents::G)
    .value("T", BaseEvents::T)
    .value("N", BaseEvents::N)
    .value("deleted", BaseEvents::del);

  py::class_<EventSpec, py::smart_holder> (m, "EventSpec")
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
  py::class_<PileupReadSet, py::smart_holder> (m, "PileupReadSet")
    .def (
      py::init<EventSpec, std::function<uint16_t()>>(),
      py::arg("event"), py::arg("qpos_cb")
    )
    .def_readwrite("event", &PileupReadSet::event)
    .def_readwrite("qpos_cb", &PileupReadSet::qpos_cb);

  py::class_<PileupCoordinates, py::smart_holder> (m, "PileupCoordinates")
    .def (
      py::init<hts_pos_t, hts_pos_t, hts_pos_t, int32_t>(),
      py::arg("gstart"), py::arg("gend"), py::arg("gpos"), py::arg("tid")
    )
    .def_readwrite("gstart", &PileupCoordinates::gstart)
    .def_readwrite("gend", &PileupCoordinates::gend)
    .def_readwrite("gpos", &PileupCoordinates::gpos)
    .def_readwrite("tid", &PileupCoordinates::tid);

  py::class_<PileupParams, py::smart_holder> (m, "PileupParams")
    .def (
      py::init<PileupCoordinates, std::string_view, uint16_t>(),
      py::arg("coordinates"), py::arg("refseq"), py::arg("readlen")
    )
    .def_readwrite("coordinates", &PileupParams::coord)
    .def_readwrite("refseq", &PileupParams::refseq)
    .def_readwrite("readlen", &PileupParams::readlen);

  py::class_<bam_pileup1_t, py::smart_holder> (m, "PileupEntry")
    .def_property_readonly("qpos",       [](const bam_pileup1_t* p) { return htsacc::qpos(p); })
    .def_property_readonly("gstart",     [](const bam_pileup1_t* p) { return htsacc::gstart(p); })
    .def_property_readonly("seq",        [](const bam_pileup1_t* p) { return htsacc::seq(p); })
    .def_property_readonly("base",       [](const bam_pileup1_t* p) { return htsacc::pileup_base(p); })
    .def_property_readonly("base_qual",  [](const bam_pileup1_t* p) { return htsacc::base_qual(p); })
    .def_property_readonly("flag",       [](const bam_pileup1_t* p) { return htsacc::flag(p); })
    .def_property_readonly("indel",      [](const bam_pileup1_t* p) { return p->indel; })
    .def_property_readonly("is_del",     [](const bam_pileup1_t* p) { return static_cast<bool>(p->is_del); })
    .def_property_readonly("is_head",    [](const bam_pileup1_t* p) { return static_cast<bool>(p->is_head); })
    .def_property_readonly("is_tail",    [](const bam_pileup1_t* p) { return static_cast<bool>(p->is_tail); })
    .def_property_readonly("is_refskip", [](const bam_pileup1_t* p) { return static_cast<bool>(p->is_refskip); });

  py::class_<PileupData, py::smart_holder> (m, "PileupData")
    .def_readonly("nread", &PileupData::nread)
    .def("__len__", [](const PileupData& pd) { return pd.nread; })
    .def("__iter__", [](PileupData& pd) {
      return py::make_iterator(pd.p1arr.get(), pd.p1arr.get() + pd.nread);
    }, py::keep_alive<0, 1>())
    .def("__getitem__", [](PileupData& pd, size_t i) -> bam_pileup1_t& {
      if (i >= pd.nread) throw py::index_error();
      return pd.p1arr[i];
    }, py::keep_alive<0, 1>());

  m.def("generate_pileup",
    [](const PileupParams& pars, std::vector<std::pair<size_t, PileupReadSet>> sets) {
      return generate_pileup(pars, sets);
    },
    "generate a synthetic pileup containing sets of reads");

}
