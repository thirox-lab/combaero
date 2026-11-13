#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "acoustics.h"
#include "area_change.h"
#include "turbomachinery.h"
#include "can_annular_solvers.h"
#include "combustion.h"
#include "combustor_loss.h"
#include "common_names.h"
#include "compressible.h"
#include "cooling_correlations.h"
#include "impingement_correlation.h"
#include "pin_fin_correlation.h"
#include "rib_correlation.h"
#include "rib_ratio_correlation.h"
#include "crossflow_cd.h"
#include "side_stream_momentum.h"
#include "fanno_mach.h"
#include "correlation_status.h"
#include "ejector.h"
#include "equilibrium.h"
#include "friction.h"
#include "geometry.h"
#include "heat_transfer.h"
#include "hole_cd.h"
#include "humidair.h"
#include "input_validation.h"
#include "incompressible.h"
#include "ideal_junction.h"
#include "mpce_junction.h"
#include "liner_radiation.h"
#include "wsgg_bordbar2021.h"
#include "materials.h"
#include "orifice.h"
#include "registry.h"
#include "slot_film_htc.h"
#include "solver_interface.h"
#include "tee_junction.h"
#include "vatistas.h"
#include "stagnation.h"
#include "state.h"
#include "thermo.h"
#include "transport.h"
#include "units.h"
#include "utils.h"

namespace py = pybind11;
using namespace combaero;

static std::vector<double>
to_vec(py::array_t<double, py::array::c_style | py::array::forcecast> arr) {
  auto buf = arr.unchecked<1>();
  std::vector<double> v(buf.shape(0));
  for (py::ssize_t i = 0; i < buf.shape(0); ++i)
    v[static_cast<std::size_t>(i)] = buf(i);
  return v;
}

// A Python callable installed as the global warning handler.
//
// Deliberately a named type rather than a lambda: it lets
// get_warning_handler() recognise, via std::function::target, that the
// current handler is one this module wrapped, and hand back the ORIGINAL
// Python object instead of wrapping it a second time.
//
// Without that, save-and-restore -- which is exactly what
// suppress_warnings() does -- grew the chain by one layer every time:
// py::object inside a std::function inside a py::cpp_function inside the
// next std::function. The layers cost nothing until a warning is emitted,
// at which point warn() recurses through all of them. Measured: 6000
// round-trips survived, 12000 overflowed the C stack and killed the
// process with no catchable Python error.
struct PyWarningHandler {
  py::object fn;

  void operator()(const std::string &msg) const {
    py::gil_scoped_acquire gil;
    fn(msg);
  }
};

PYBIND11_MODULE(_core, m) {

  // ---------------------------------------------------------------------
  // Parametrised rib correlations
  // ---------------------------------------------------------------------
  py::enum_<combaero::cooling::RibProvenance>(m, "RibProvenance")
      .value("Extracted", combaero::cooling::RibProvenance::Extracted,
             "An equation printed in a source, transcribed and checked "
             "against it.")
      .value("Fitted", combaero::cooling::RibProvenance::Fitted,
             "A form we chose, fitted to data we digitised from a figure. "
             "The exponents are ours.")
      .value("User", combaero::cooling::RibProvenance::User,
             "Supplied by the caller. Carries no claim.");

  py::enum_<combaero::cooling::RibShape>(
      m, "RibShape",
      "Rib shape beyond its angle. A set binds one shape; a series of "
      "another shape is refused rather than scored.")
      .value("Unspecified", combaero::cooling::RibShape::Unspecified,
             "Binds nothing: the default for a user-built set.")
      .value("Transverse", combaero::cooling::RibShape::Transverse)
      .value("Parallel", combaero::cooling::RibShape::Parallel)
      .value("Crossed", combaero::cooling::RibShape::Crossed)
      .value("V", combaero::cooling::RibShape::V,
             "Forward V, pointing downstream.")
      .value("Lambda", combaero::cooling::RibShape::Lambda,
             "Pointing upstream: a V with the flow reversed.");

  py::enum_<combaero::cooling::RibCorrelationSet::RAlphaShape>(
      m, "RAlphaShape")
      .value("PowerLaw",
             combaero::cooling::RibCorrelationSet::RAlphaShape::PowerLaw,
             "R is a power-law product, R_alpha's exponent applied to "
             "alpha/90 like any other geometry term.")
      .value("QuadraticAlpha",
             combaero::cooling::RibCorrelationSet::RAlphaShape::QuadraticAlpha,
             "R/[(p/e/10)^0.35 (W/H)^m] is a quadratic in alpha/90 (the "
             "R_quad_* fields), with m switching on alpha == 90 deg. "
             "R_eD, R_WH and R_alpha are ignored in this shape.")
      .value("QuadraticAlphaTwoBand",
             combaero::cooling::RibCorrelationSet::RAlphaShape::
                 QuadraticAlphaTwoBand,
             "R is a quadratic in alpha/90 with no p/e or W/H power-law "
             "term; W/H alone selects which of two quadratics applies "
             "(R_quad_* above R_WH_band_boundary, R_quad_narrow_* below). "
             "Fitted by this project to a drawn figure line, not an "
             "extraction -- see han_ribbed.md item 40.");

  py::enum_<combaero::cooling::RibCorrelationSet::GShapeModel>(
      m, "GShapeModel")
      .value("Fixed", combaero::cooling::RibCorrelationSet::GShapeModel::Fixed,
             "G_alpha and G_pe are fixed power-law terms, as for a set with "
             "no channel-shape switch.")
      .value("SquareVsRectangular",
             combaero::cooling::RibCorrelationSet::GShapeModel::
                 SquareVsRectangular,
             "G's alpha and p/e exponents switch on whether W/H == 1 "
             "(square, the G_shape_*_square fields) or not (rectangular, "
             "the G_shape_*_rect fields). G_alpha and G_pe are ignored.")
      .value("NarrowChannelAlphaSwitch",
             combaero::cooling::RibCorrelationSet::GShapeModel::
                 NarrowChannelAlphaSwitch,
             "G = C(e+)^n, with C switching on alpha == 90 deg "
             "(G_narrow_C_alpha90/off_axis) and, below "
             "G_narrow_WH_band_boundary, a further (W/H) correction to "
             "both C and n. C_G, G_eD, G_pe, G_WH and G_alpha are ignored.");

  py::class_<combaero::cooling::RibTerm>(m, "RibTerm")
      .def(py::init<>())
      .def(py::init([](double exponent, double reference) {
             return combaero::cooling::RibTerm{exponent, reference};
           }),
           py::arg("exponent") = 0.0, py::arg("reference") = 1.0)
      .def_readwrite("exponent", &combaero::cooling::RibTerm::exponent)
      .def_readwrite("reference", &combaero::cooling::RibTerm::reference,
                     "What the variable is divided by before the exponent is "
                     "applied. Data, not a convention: (p/e/10)^b and (p/e)^b "
                     "are the same function with different constants.");

  py::class_<combaero::cooling::RibRange>(m, "RibRange")
      .def(py::init<>())
      .def(py::init([](double lo, double hi) {
             return combaero::cooling::RibRange{lo, hi};
           }),
           py::arg("lo") = 0.0, py::arg("hi") = 0.0)
      .def_readwrite("lo", &combaero::cooling::RibRange::lo)
      .def_readwrite("hi", &combaero::cooling::RibRange::hi);

  py::class_<combaero::cooling::RibGeometry>(m, "RibGeometry")
      .def(py::init<>())
      .def(py::init([](double e_D, double p_e, double W_H, double alpha_deg) {
             return combaero::cooling::RibGeometry{e_D, p_e, W_H, alpha_deg};
           }),
           py::arg("e_D") = 0.0, py::arg("p_e") = 0.0, py::arg("W_H") = 1.0,
           py::arg("alpha_deg") = 90.0)
      .def_readwrite("e_D", &combaero::cooling::RibGeometry::e_D)
      .def_readwrite("p_e", &combaero::cooling::RibGeometry::p_e)
      .def_readwrite("W_H", &combaero::cooling::RibGeometry::W_H)
      .def_readwrite("alpha_deg", &combaero::cooling::RibGeometry::alpha_deg);

  py::class_<combaero::cooling::RibResult>(m, "RibResult")
      .def_readonly("R", &combaero::cooling::RibResult::R)
      .def_readonly("f", &combaero::cooling::RibResult::f)
      .def_readonly("e_plus", &combaero::cooling::RibResult::e_plus)
      .def_readonly("G", &combaero::cooling::RibResult::G)
      .def_readonly("St_r", &combaero::cooling::RibResult::St_r,
                    "Ribbed-side Stanton number. Combining it with the smooth "
                    "walls is the caller's job.")
      .def_readonly("dSt_dRe", &combaero::cooling::RibResult::dSt_dRe)
      .def_readonly("G_bar", &combaero::cooling::RibResult::G_bar,
                    "Four-wall average G, when the set prints its own; 0 "
                    "otherwise -- check has_G_bar.")
      .def_readonly("has_G_bar", &combaero::cooling::RibResult::has_G_bar)
      .def_readonly("extrapolated",
                    &combaero::cooling::RibResult::extrapolated,
                    "Outside the set's advisory validity. Advisory because a "
                    "band belongs to the source's rig, not the caller's.");

  py::enum_<combaero::cooling::AccuracyProvenance>(
      m, "AccuracyProvenance",
      "Where a correlation's stated accuracy came from. Only Stated may be "
      "used as a band to judge a model against: a Measured figure IS that "
      "model's own error, and scoring a model inside its own error answers "
      "nothing. See docs/VALIDATION_POLICY.md and issue #389.")
      .value("Unstated", combaero::cooling::AccuracyProvenance::Unstated,
             "No figure available; value is NaN and must not be read.")
      .value("Stated", combaero::cooling::AccuracyProvenance::Stated,
             "The source states it. The only judgeable kind.")
      .value("Measured", combaero::cooling::AccuracyProvenance::Measured,
             "This project measured it against digitised data. Reportable, "
             "never a judging band.");

  py::class_<combaero::cooling::StatedAccuracy>(
      m, "StatedAccuracy",
      "A correlation accuracy as a fraction, carried with its provenance.")
      .def(py::init<>())
      .def_readwrite("value", &combaero::cooling::StatedAccuracy::value,
                     "Fraction (0.06 for 'within 6%'); NaN when Unstated.")
      .def_readwrite("provenance",
                     &combaero::cooling::StatedAccuracy::provenance)
      .def("usable_as_band",
           &combaero::cooling::StatedAccuracy::usable_as_band,
           "True only for Stated. A Measured figure must never be the band "
           "a model is judged against.")
      .def_static("stated", &combaero::cooling::StatedAccuracy::stated,
                  py::arg("value"))
      .def_static("measured", &combaero::cooling::StatedAccuracy::measured,
                  py::arg("value"))
      .def_static("unstated", &combaero::cooling::StatedAccuracy::unstated);

  py::class_<combaero::cooling::RibCorrelationSet>(m, "RibCorrelationSet")
      .def(py::init<>())
      .def_readwrite("name", &combaero::cooling::RibCorrelationSet::name)
      .def_readwrite("source", &combaero::cooling::RibCorrelationSet::source)
      .def_readwrite("validity_source",
                     &combaero::cooling::RibCorrelationSet::validity_source)
      .def_readwrite("provenance",
                     &combaero::cooling::RibCorrelationSet::provenance)
      .def_readwrite("symmetric",
                     &combaero::cooling::RibCorrelationSet::symmetric)
      .def_readwrite("shape", &combaero::cooling::RibCorrelationSet::shape)
      .def_readwrite("R_alpha_shape",
                     &combaero::cooling::RibCorrelationSet::R_alpha_shape)
      .def_readwrite("C_R", &combaero::cooling::RibCorrelationSet::C_R)
      .def_readwrite("R_eD", &combaero::cooling::RibCorrelationSet::R_eD)
      .def_readwrite("R_pe", &combaero::cooling::RibCorrelationSet::R_pe)
      .def_readwrite("R_WH", &combaero::cooling::RibCorrelationSet::R_WH)
      .def_readwrite("R_alpha", &combaero::cooling::RibCorrelationSet::R_alpha)
      .def_readwrite("R_quad_c0",
                     &combaero::cooling::RibCorrelationSet::R_quad_c0)
      .def_readwrite("R_quad_c1",
                     &combaero::cooling::RibCorrelationSet::R_quad_c1)
      .def_readwrite("R_quad_c2",
                     &combaero::cooling::RibCorrelationSet::R_quad_c2)
      .def_readwrite(
          "R_quad_WH_exponent_at_90",
          &combaero::cooling::RibCorrelationSet::R_quad_WH_exponent_at_90)
      .def_readwrite(
          "R_quad_WH_exponent_off_90",
          &combaero::cooling::RibCorrelationSet::R_quad_WH_exponent_off_90)
      .def_readwrite("R_quad_WH_cap",
                     &combaero::cooling::RibCorrelationSet::R_quad_WH_cap)
      .def_readwrite("R_quad_narrow_c0",
                     &combaero::cooling::RibCorrelationSet::R_quad_narrow_c0)
      .def_readwrite("R_quad_narrow_c1",
                     &combaero::cooling::RibCorrelationSet::R_quad_narrow_c1)
      .def_readwrite("R_quad_narrow_c2",
                     &combaero::cooling::RibCorrelationSet::R_quad_narrow_c2)
      .def_readwrite("R_WH_band_boundary",
                     &combaero::cooling::RibCorrelationSet::R_WH_band_boundary)
      .def_readwrite("G_shape_model",
                     &combaero::cooling::RibCorrelationSet::G_shape_model)
      .def_readwrite("C_G", &combaero::cooling::RibCorrelationSet::C_G)
      .def_readwrite("G_eD", &combaero::cooling::RibCorrelationSet::G_eD)
      .def_readwrite("G_pe", &combaero::cooling::RibCorrelationSet::G_pe)
      .def_readwrite("G_WH", &combaero::cooling::RibCorrelationSet::G_WH)
      .def_readwrite("G_alpha", &combaero::cooling::RibCorrelationSet::G_alpha)
      .def_readwrite("G_eplus_exponent",
                     &combaero::cooling::RibCorrelationSet::G_eplus_exponent)
      .def_readwrite(
          "G_shape_alpha_exponent_square",
          &combaero::cooling::RibCorrelationSet::G_shape_alpha_exponent_square)
      .def_readwrite(
          "G_shape_alpha_exponent_rect",
          &combaero::cooling::RibCorrelationSet::G_shape_alpha_exponent_rect)
      .def_readwrite(
          "G_shape_pe_exponent_square",
          &combaero::cooling::RibCorrelationSet::G_shape_pe_exponent_square)
      .def_readwrite(
          "G_shape_pe_exponent_rect",
          &combaero::cooling::RibCorrelationSet::G_shape_pe_exponent_rect)
      .def_readwrite(
          "G_narrow_C_alpha90",
          &combaero::cooling::RibCorrelationSet::G_narrow_C_alpha90)
      .def_readwrite(
          "G_narrow_C_off_axis",
          &combaero::cooling::RibCorrelationSet::G_narrow_C_off_axis)
      .def_readwrite(
          "G_narrow_WH_band_boundary",
          &combaero::cooling::RibCorrelationSet::G_narrow_WH_band_boundary)
      .def_readwrite(
          "G_narrow_WH_C_exponent",
          &combaero::cooling::RibCorrelationSet::G_narrow_WH_C_exponent)
      .def_readwrite(
          "G_narrow_WH_n_exponent",
          &combaero::cooling::RibCorrelationSet::G_narrow_WH_n_exponent)
      .def_readwrite("C_Gbar", &combaero::cooling::RibCorrelationSet::C_Gbar,
                     "Printed four-wall G_bar coefficient; 0 = not printed.")
      .def_readwrite(
          "Gbar_eplus_exponent",
          &combaero::cooling::RibCorrelationSet::Gbar_eplus_exponent)
      .def_readwrite("valid_Re",
                     &combaero::cooling::RibCorrelationSet::valid_Re)
      .def_readwrite("valid_eD",
                     &combaero::cooling::RibCorrelationSet::valid_eD)
      .def_readwrite("valid_pe",
                     &combaero::cooling::RibCorrelationSet::valid_pe)
      .def_readwrite("valid_WH",
                     &combaero::cooling::RibCorrelationSet::valid_WH)
      .def_readwrite("valid_alpha",
                     &combaero::cooling::RibCorrelationSet::valid_alpha)
      .def_readwrite("valid_eplus",
                     &combaero::cooling::RibCorrelationSet::valid_eplus)
      .def_readwrite("valid_Pr",
                     &combaero::cooling::RibCorrelationSet::valid_Pr)
      .def_readwrite("accuracy_R",
                     &combaero::cooling::RibCorrelationSet::accuracy_R)
      .def_readwrite("accuracy_G",
                     &combaero::cooling::RibCorrelationSet::accuracy_G);

  m.def("han_1988_orthogonal", &combaero::cooling::han_1988_orthogonal,
        "Han (1988) 90 deg orthogonal ribs. Extracted and confirmed; see "
        "validation/cooling/extractions/han_ribbed.md.");
  m.def("rallabandi_2009_high_re", &combaero::cooling::rallabandi_2009_high_re,
        "Rallabandi, Yang and Han (2009) 45 deg sharp-edged ribs at high "
        "Reynolds number. Extracted and confirmed; see "
        "validation/cooling/extractions/han_ribbed_high_re.md.");
  m.def("han_park_1988_angled", &combaero::cooling::han_park_1988_angled,
        "Han and Park (1988) angled ribs in broad-aspect-ratio rectangular "
        "ducts. Extracted and confirmed; see "
        "validation/cooling/extractions/han_ribbed.md.");
  m.def("han_1989_narrow_channel", &combaero::cooling::han_1989_narrow_channel,
        "Han et al. (1989), Eq. 4.19, narrow-aspect-ratio channels "
        "(1/4 < W/H < 1). G is extracted from the text; R has no printed "
        "equation and is this project's own fit to Fig. 4.48a's drawn "
        "line -- see validation/cooling/extractions/han_ribbed.md, item 40.");
  m.def("han_zhang_lee_1991", &combaero::cooling::han_zhang_lee_1991,
        py::arg("shape"), py::arg("alpha_deg"),
        "Han, Zhang and Lee (1991) JHT 113, 590, Table 2: R, G and G_bar "
        "for Transverse at 90 deg and Parallel/Crossed/V/Lambda at 45 or 60 "
        "deg, one rig (square, e/D 0.0625, P/e 10). Anything else raises. "
        "See validation/cooling/extractions/han_zhang_lee_1991_jht.md.");
  m.def("validate_rib_set", &combaero::cooling::validate_rib_set,
        py::arg("correlation_set"),
        "Reject a set that cannot be evaluated. Hard errors, unlike the "
        "advisory validity ranges: these are mistakes, not operating points.");
  m.def("evaluate_rib", &combaero::cooling::evaluate_rib,
        py::arg("correlation_set"), py::arg("geometry"), py::arg("Re"),
        "Evaluate the chain. Re may be negative or zero; the guards are "
        "smooth through both.");
  m.def("rib_ribbed_wall_fraction", &combaero::cooling::rib_ribbed_wall_fraction,
        py::arg("W_H"), py::arg("n_ribbed_walls"),
        "Ribbed fraction of a W x H duct's wetted perimeter, ribs on the W "
        "walls: W/(W+H) for 2 (two opposite walls, Han's configuration), 1 "
        "for 4, W/(2(W+H)) for 1. Any other count raises.");
  m.def("rib_channel_friction", &combaero::cooling::rib_channel_friction,
        py::arg("correlation_set"), py::arg("geometry"), py::arg("Re"),
        py::arg("n_ribbed_walls"),
        "(f_D, df_D/dRe): the channel's DARCY friction factor, "
        "4 [phi f + (1 - phi) 0.046 Re^-0.2], with f the set's four-sided "
        "Fanning factor and phi the ribbed wall fraction (Han's area "
        "weighting, Han and Park 1988 Eq. 6). Exactly Han's measured "
        "two-wall fbar (x4) for n_ribbed_walls=2. Re may be negative or zero.");

  // -- Ratio-form rib sets (#444): Nu/Nu0 and f/f0 on a declared baseline --
  py::class_<combaero::cooling::RatioBaseline>(
      m, "RatioBaseline",
      "Power-law smooth-duct baseline: coeff Re^re_exponent Pr^pr_exponent. "
      "Dittus-Boelter heating is (0.023, 0.8, 0.4).")
      .def(py::init<>())
      .def(py::init([](double coeff, double re_exponent, double pr_exponent) {
             return combaero::cooling::RatioBaseline{coeff, re_exponent, pr_exponent};
           }),
           py::arg("coeff"), py::arg("re_exponent"), py::arg("pr_exponent") = 0.0)
      .def_readwrite("coeff", &combaero::cooling::RatioBaseline::coeff)
      .def_readwrite("re_exponent", &combaero::cooling::RatioBaseline::re_exponent)
      .def_readwrite("pr_exponent", &combaero::cooling::RatioBaseline::pr_exponent);

  py::enum_<combaero::cooling::RatioBelowFloor>(m, "RatioBelowFloor",
                                  "What Nu hands over to below the fitted Re range.")
      .value("Gnielinski", combaero::cooling::RatioBelowFloor::Gnielinski,
             "Smooth-pipe Gnielinski: laminar-limited as Re -> 0.")
      .value("SourceBaseline", combaero::cooling::RatioBelowFloor::SourceBaseline,
             "Hold the ratio on the source baseline (D-B tends to 0).")
      .value("User", combaero::cooling::RatioBelowFloor::User,
             "The caller's own power-law Nu0 (RibRatioOptions.user_Nu0).");

  py::class_<combaero::cooling::RibRatioOptions>(m, "RibRatioOptions")
      .def(py::init<>())
      .def_readwrite("below_floor", &combaero::cooling::RibRatioOptions::below_floor)
      .def_readwrite("user_Nu0", &combaero::cooling::RibRatioOptions::user_Nu0);

  py::class_<combaero::cooling::RibRatioSet>(
      m, "RibRatioSet",
      "Nu/Nu0 and f/f0 as multipliers on the baseline the source fitted to. "
      "In range Nu = r Nu0_source; below Re_floor a C1 handover to an "
      "extrapolation baseline. See include/rib_ratio_correlation.h.")
      .def(py::init<>())
      .def_readwrite("name", &combaero::cooling::RibRatioSet::name)
      .def_readwrite("source", &combaero::cooling::RibRatioSet::source)
      .def_readwrite("validity_source", &combaero::cooling::RibRatioSet::validity_source)
      .def_readwrite("provenance", &combaero::cooling::RibRatioSet::provenance)
      .def_readwrite("shape", &combaero::cooling::RibRatioSet::shape)
      .def_readwrite("symmetric", &combaero::cooling::RibRatioSet::symmetric)
      .def_readwrite("C_Nu", &combaero::cooling::RibRatioSet::C_Nu)
      .def_readwrite("Nu_Re", &combaero::cooling::RibRatioSet::Nu_Re)
      .def_readwrite("Nu_eD", &combaero::cooling::RibRatioSet::Nu_eD)
      .def_readwrite("Nu_pe", &combaero::cooling::RibRatioSet::Nu_pe)
      .def_readwrite("Nu_WH", &combaero::cooling::RibRatioSet::Nu_WH)
      .def_readwrite("Nu_alpha", &combaero::cooling::RibRatioSet::Nu_alpha)
      .def_readwrite("C_f", &combaero::cooling::RibRatioSet::C_f)
      .def_readwrite("f_Re", &combaero::cooling::RibRatioSet::f_Re)
      .def_readwrite("f_eD", &combaero::cooling::RibRatioSet::f_eD)
      .def_readwrite("f_pe", &combaero::cooling::RibRatioSet::f_pe)
      .def_readwrite("f_WH", &combaero::cooling::RibRatioSet::f_WH)
      .def_readwrite("f_alpha", &combaero::cooling::RibRatioSet::f_alpha)
      .def_readwrite("Nu0_source", &combaero::cooling::RibRatioSet::Nu0_source)
      .def_readwrite("f0_source", &combaero::cooling::RibRatioSet::f0_source)
      .def_readwrite("Re_floor", &combaero::cooling::RibRatioSet::Re_floor)
      .def_readwrite("Re_floor_f", &combaero::cooling::RibRatioSet::Re_floor_f)
      .def_readwrite("valid_Re", &combaero::cooling::RibRatioSet::valid_Re)
      .def_readwrite("valid_eD", &combaero::cooling::RibRatioSet::valid_eD)
      .def_readwrite("valid_pe", &combaero::cooling::RibRatioSet::valid_pe)
      .def_readwrite("valid_WH", &combaero::cooling::RibRatioSet::valid_WH)
      .def_readwrite("valid_alpha", &combaero::cooling::RibRatioSet::valid_alpha)
      .def_readwrite("valid_Pr", &combaero::cooling::RibRatioSet::valid_Pr)
      .def_readwrite("accuracy_Nu", &combaero::cooling::RibRatioSet::accuracy_Nu)
      .def_readwrite("accuracy_f", &combaero::cooling::RibRatioSet::accuracy_f);

  py::class_<combaero::cooling::RibRatioResult>(m, "RibRatioResult")
      .def_readonly("Nu", &combaero::cooling::RibRatioResult::Nu)
      .def_readonly("dNu_dRe", &combaero::cooling::RibRatioResult::dNu_dRe)
      .def_readonly("f", &combaero::cooling::RibRatioResult::f)
      .def_readonly("df_dRe", &combaero::cooling::RibRatioResult::df_dRe)
      .def_readonly("ratio_Nu", &combaero::cooling::RibRatioResult::ratio_Nu)
      .def_readonly("ratio_f", &combaero::cooling::RibRatioResult::ratio_f)
      .def_readonly("below_floor", &combaero::cooling::RibRatioResult::below_floor)
      .def_readonly("extrapolated", &combaero::cooling::RibRatioResult::extrapolated);

  m.def("evaluate_rib_ratio", &combaero::cooling::evaluate_rib_ratio, py::arg("correlation_set"),
        py::arg("geometry"), py::arg("Re"), py::arg("Pr"),
        py::arg("options") = combaero::cooling::RibRatioOptions{},
        "Nu, f and their analytic Re-derivatives. Never raises: Re may be "
        "negative or zero, and every output is finite and non-negative.");
  m.def("validate_rib_ratio_set", &combaero::cooling::validate_rib_ratio_set,
        py::arg("correlation_set"),
        "Reject a ratio set that cannot be evaluated -- mistakes, not "
        "operating points.");
  m.def("validate_rib_ratio_options", &combaero::cooling::validate_rib_ratio_options,
        py::arg("options"));
  m.def("taslim_spring_1987", &combaero::cooling::taslim_spring_1987,
        py::arg("aspect_ratio_taslim"), py::arg("e_D"),
        "Taslim & Spring (1987) two-side transverse turbulators, P/e 10, one "
        "ratio set per tested configuration. aspect_ratio_taslim is HEIGHT/"
        "WIDTH (Han's W/H = 1/AR). Nu/Nu_DB = C (Re/1e4)^-0.2 as the paper "
        "states; C and f digitised from its Figs. 9 and 11 (provenance "
        "Fitted). Untested configurations raise.");

  // ---------------------------------------------------------------------
  // Jet impingement correlations
  // ---------------------------------------------------------------------
  py::enum_<combaero::cooling::ImpingementThermalBC>(m, "ImpingementThermalBC")
      .value("ConstantHeatFlux",
             combaero::cooling::ImpingementThermalBC::ConstantHeatFlux)
      .value("ConstantWallTemperature",
             combaero::cooling::ImpingementThermalBC::ConstantWallTemperature);

  py::class_<combaero::cooling::SingleJetImpingementSet>(
      m, "SingleJetImpingementSet")
      .def(py::init<>())
      .def_readwrite("name",
                     &combaero::cooling::SingleJetImpingementSet::name)
      .def_readwrite("source",
                     &combaero::cooling::SingleJetImpingementSet::source)
      .def_readwrite("A", &combaero::cooling::SingleJetImpingementSet::A)
      .def_readwrite("B", &combaero::cooling::SingleJetImpingementSet::B)
      .def_readwrite("C", &combaero::cooling::SingleJetImpingementSet::C)
      .def_readwrite("Re_exponent",
                     &combaero::cooling::SingleJetImpingementSet::Re_exponent)
      .def_readwrite(
          "n_const_heat_flux",
          &combaero::cooling::SingleJetImpingementSet::n_const_heat_flux)
      .def_readwrite(
          "n_const_wall_temp",
          &combaero::cooling::SingleJetImpingementSet::n_const_wall_temp);

  m.def("goldstein_1986_single_jet",
        &combaero::cooling::goldstein_1986_single_jet,
        "Goldstein, Behbahani and Heppelmann (1986) single round jet on a "
        "flat plate. Extracted and confirmed; see "
        "validation/cooling/extractions/han_impingement.md.");
  m.def("validate_single_jet_set", &combaero::cooling::validate_single_jet_set,
        py::arg("correlation_set"),
        "Reject a set that cannot be evaluated. Hard error, not an "
        "operating-point guard.");
  m.def("single_jet_impingement_nu",
        &combaero::cooling::single_jet_impingement_nu,
        py::arg("correlation_set"), py::arg("bc"), py::arg("Re"),
        py::arg("L_D"), py::arg("R_D"),
        "Nu_bar for a single round jet impinging on a flat plate.");

  py::class_<combaero::cooling::SingleJetImpingementResult>(
      m, "SingleJetImpingementResult")
      .def_readonly("Nu", &combaero::cooling::SingleJetImpingementResult::Nu)
      .def_readonly("dNu_dRe",
                    &combaero::cooling::SingleJetImpingementResult::dNu_dRe);

  m.def("single_jet_impingement", &combaero::cooling::single_jet_impingement,
        py::arg("correlation_set"), py::arg("bc"), py::arg("Re"),
        py::arg("L_D"), py::arg("R_D"),
        "Same as single_jet_impingement_nu, plus its analytic Re-derivative "
        "for a network element's wall-coupling Jacobian.");

  py::enum_<combaero::cooling::JetHolePattern>(m, "JetHolePattern")
      .value("Inline", combaero::cooling::JetHolePattern::Inline)
      .value("Staggered", combaero::cooling::JetHolePattern::Staggered);

  py::class_<combaero::cooling::ImpingementRange>(m, "ImpingementRange")
      .def(py::init<>())
      .def(py::init([](double lo, double hi) {
             return combaero::cooling::ImpingementRange{lo, hi};
           }),
           py::arg("lo") = 0.0, py::arg("hi") = 0.0)
      .def_readwrite("lo", &combaero::cooling::ImpingementRange::lo)
      .def_readwrite("hi", &combaero::cooling::ImpingementRange::hi);

  py::class_<combaero::cooling::JetArrayGeometricFit>(
      m, "JetArrayGeometricFit")
      .def(py::init<>())
      .def(py::init([](double C, double nx, double ny, double nz) {
             return combaero::cooling::JetArrayGeometricFit{C, nx, ny, nz};
           }),
           py::arg("C") = 0.0, py::arg("nx") = 0.0, py::arg("ny") = 0.0,
           py::arg("nz") = 0.0)
      .def_readwrite("C", &combaero::cooling::JetArrayGeometricFit::C)
      .def_readwrite("nx", &combaero::cooling::JetArrayGeometricFit::nx)
      .def_readwrite("ny", &combaero::cooling::JetArrayGeometricFit::ny)
      .def_readwrite("nz", &combaero::cooling::JetArrayGeometricFit::nz);

  py::class_<combaero::cooling::JetArrayCorrelationSet>(
      m, "JetArrayCorrelationSet")
      .def(py::init<>())
      .def_readwrite("name", &combaero::cooling::JetArrayCorrelationSet::name)
      .def_readwrite("source",
                     &combaero::cooling::JetArrayCorrelationSet::source)
      .def_readwrite(
          "validity_source",
          &combaero::cooling::JetArrayCorrelationSet::validity_source)
      .def_readwrite("pattern",
                     &combaero::cooling::JetArrayCorrelationSet::pattern)
      .def_readwrite("A_fit", &combaero::cooling::JetArrayCorrelationSet::A_fit)
      .def_readwrite("m_fit", &combaero::cooling::JetArrayCorrelationSet::m_fit)
      .def_readwrite("B_fit", &combaero::cooling::JetArrayCorrelationSet::B_fit)
      .def_readwrite("n_fit", &combaero::cooling::JetArrayCorrelationSet::n_fit)
      .def_readwrite("valid_Re_j",
                     &combaero::cooling::JetArrayCorrelationSet::valid_Re_j)
      .def_readwrite("valid_Gc_Gj",
                     &combaero::cooling::JetArrayCorrelationSet::valid_Gc_Gj)
      .def_readwrite("valid_xn_d",
                     &combaero::cooling::JetArrayCorrelationSet::valid_xn_d)
      .def_readwrite("valid_yn_d",
                     &combaero::cooling::JetArrayCorrelationSet::valid_yn_d)
      .def_readwrite("valid_z_d",
                     &combaero::cooling::JetArrayCorrelationSet::valid_z_d)
      .def_readwrite(
          "valid_aspect_ratio",
          &combaero::cooling::JetArrayCorrelationSet::valid_aspect_ratio)
      .def_readwrite(
          "standard_error",
          &combaero::cooling::JetArrayCorrelationSet::standard_error);

  py::class_<combaero::cooling::JetArrayImpingementResult>(
      m, "JetArrayImpingementResult")
      .def_readonly("Nu", &combaero::cooling::JetArrayImpingementResult::Nu)
      .def_readonly(
          "dNu_dRe_j",
          &combaero::cooling::JetArrayImpingementResult::dNu_dRe_j)
      .def_readonly(
          "dNu_dGc_Gj",
          &combaero::cooling::JetArrayImpingementResult::dNu_dGc_Gj)
      .def_readonly(
          "extrapolated",
          &combaero::cooling::JetArrayImpingementResult::extrapolated);

  m.attr("FLORSCHUETZ_1981_DEFAULT_CD") =
      combaero::cooling::FLORSCHUETZ_1981_DEFAULT_CD;

  m.def("florschuetz_1981_inline", &combaero::cooling::florschuetz_1981_inline,
        "Florschuetz, Truman and Metzger (1981) jet array with crossflow, "
        "inline hole pattern. Extracted and confirmed against the primary "
        "paper; see validation/cooling/extractions/han_impingement.md.");
  m.def("florschuetz_1981_staggered",
        &combaero::cooling::florschuetz_1981_staggered,
        "Florschuetz, Truman and Metzger (1981) jet array with crossflow, "
        "staggered hole pattern. Extracted and confirmed against the "
        "primary paper; see "
        "validation/cooling/extractions/han_impingement.md.");
  m.def("validate_jet_array_set", &combaero::cooling::validate_jet_array_set,
        py::arg("correlation_set"),
        "Reject a set that cannot be evaluated. Hard error, not an "
        "operating-point guard.");
  m.def("jet_array_impingement_nu", &combaero::cooling::jet_array_impingement_nu,
        py::arg("correlation_set"), py::arg("Re_j"), py::arg("Gc_Gj"),
        py::arg("Pr"), py::arg("xn_d"), py::arg("yn_d"), py::arg("z_d"),
        "Nu at one spanwise row. Re_j and Gc_Gj are that row's own local "
        "values, not an array mean.");
  m.def("crossflow_to_jet_ratio_at_x",
        &combaero::cooling::crossflow_to_jet_ratio_at_x, py::arg("yn_d"),
        py::arg("z_d"), py::arg("C_D"), py::arg("x_over_xn"),
        "Gc/Gj at continuous streamwise position x/xn (Florschuetz Eq. 8). "
        "Depends on (yn/d)(z/d) only, not xn/d.");
  m.def("crossflow_to_jet_ratio_at_row",
        &combaero::cooling::crossflow_to_jet_ratio_at_row, py::arg("yn_d"),
        py::arg("z_d"), py::arg("C_D"), py::arg("row"),
        "Gc/Gj at discrete spanwise row (1-indexed, counting from "
        "upstream). Row 1 is exactly 0.");

  // Jet-plate row from network flows (#465): the generic core the
  // impingement elements call. See impingement_correlation.h.
  py::class_<combaero::cooling::JetRowGeometry>(m, "JetRowGeometry")
      .def(py::init([](double d, double xn_d, double yn_d, double z_d,
                       double n_holes) {
             combaero::cooling::JetRowGeometry g;
             g.d = d;
             g.xn_d = xn_d;
             g.yn_d = yn_d;
             g.z_d = z_d;
             g.n_holes = n_holes;
             return g;
           }),
           py::arg("d"), py::arg("xn_d"), py::arg("yn_d"), py::arg("z_d"),
           py::arg("n_holes"))
      .def_readwrite("d", &combaero::cooling::JetRowGeometry::d)
      .def_readwrite("xn_d", &combaero::cooling::JetRowGeometry::xn_d)
      .def_readwrite("yn_d", &combaero::cooling::JetRowGeometry::yn_d)
      .def_readwrite("z_d", &combaero::cooling::JetRowGeometry::z_d)
      .def_readwrite("n_holes", &combaero::cooling::JetRowGeometry::n_holes);
  py::class_<combaero::cooling::JetRowHeatTransfer>(m, "JetRowHeatTransfer")
      .def_readonly("h", &combaero::cooling::JetRowHeatTransfer::h)
      .def_readonly("Nu", &combaero::cooling::JetRowHeatTransfer::Nu)
      .def_readonly("Re_j", &combaero::cooling::JetRowHeatTransfer::Re_j)
      .def_readonly("Gc_Gj", &combaero::cooling::JetRowHeatTransfer::Gc_Gj)
      .def_readonly("dh_dm_jet",
                    &combaero::cooling::JetRowHeatTransfer::dh_dm_jet)
      .def_readonly("dh_dm_crossflow",
                    &combaero::cooling::JetRowHeatTransfer::dh_dm_crossflow)
      .def_readonly("extrapolated",
                    &combaero::cooling::JetRowHeatTransfer::extrapolated)
      .def_readonly("eta_r", &combaero::cooling::JetRowHeatTransfer::eta_r)
      .def_readonly("deta_r_dm_jet",
                    &combaero::cooling::JetRowHeatTransfer::deta_r_dm_jet)
      .def_readonly("deta_r_dm_crossflow",
                    &combaero::cooling::JetRowHeatTransfer::deta_r_dm_crossflow)
      .def_readonly("eta_r_extrapolated",
                    &combaero::cooling::JetRowHeatTransfer::eta_r_extrapolated);
  py::class_<combaero::cooling::ImpingementEtaR>(m, "ImpingementEtaR")
      .def_readonly("eta", &combaero::cooling::ImpingementEtaR::eta)
      .def_readonly("deta_dGc_Gj",
                    &combaero::cooling::ImpingementEtaR::deta_dGc_Gj)
      .def_readonly("extrapolated",
                    &combaero::cooling::ImpingementEtaR::extrapolated);
  m.attr("FLORSCHUETZ_SU_1985_ETA_A") =
      combaero::cooling::FLORSCHUETZ_SU_1985_ETA_A;
  m.attr("FLORSCHUETZ_SU_1985_ETA_R") =
      combaero::cooling::FLORSCHUETZ_SU_1985_ETA_R;
  m.attr("FLORSCHUETZ_SU_1985_ETA_N") =
      combaero::cooling::FLORSCHUETZ_SU_1985_ETA_N;
  m.def("impingement_eta_r", &combaero::cooling::impingement_eta_r,
        py::arg("Gc_Gj"), py::arg("xn_d"), py::arg("yn_d"), py::arg("z_d"),
        "Crossflow temperature influence factor of one jet row, fitted to\n"
        "Florschuetz and Su (1985) NASA CR-3936 Appendix A (inline arrays):\n"
        "eta_r = X^n/(1+X^n), X = a (Gc/Gj)(z/d)^r, with d(eta)/d(Gc/Gj).\n"
        "A FITTED closure (the source publishes none); extrapolated flags\n"
        "geometry or Gc/Gj outside the inline data.");
  m.def("jet_row_gc_gj_per_flow_ratio",
        &combaero::cooling::jet_row_gc_gj_per_flow_ratio, py::arg("yn_d"),
        py::arg("z_d"),
        "(pi/4)/((yn/d)(z/d)): Gc/Gj per unit m_crossflow/m_jet.");
  m.def("jet_row_heat_transfer", &combaero::cooling::jet_row_heat_transfer,
        py::arg("correlation_set"), py::arg("geometry"), py::arg("m_jet"),
        py::arg("m_crossflow"), py::arg("mu"), py::arg("k"), py::arg("Pr"),
        "Target-side h of one jet-plate row from its own jet flow and the\n"
        "crossflow mass approaching it, with analytic dh/dm_jet and\n"
        "dh/dm_crossflow. Which flows count as crossflow is the caller's\n"
        "configuration.");

  m.attr("ROGERS_HERSH_1976_CAP_BAND") = combaero::ROGERS_HERSH_1976_CAP_BAND;
  m.attr("FLORSCHUETZ_METZGER_SU_1982_CAP_BAND") =
      combaero::FLORSCHUETZ_METZGER_SU_1982_CAP_BAND;
  m.attr("ROGERS_HERSH_1976_INJECTION_A") = combaero::ROGERS_HERSH_1976_INJECTION_A;
  m.attr("ROGERS_HERSH_1976_INJECTION_B") = combaero::ROGERS_HERSH_1976_INJECTION_B;
  m.attr("ROGERS_HERSH_1976_SUCTION_A") = combaero::ROGERS_HERSH_1976_SUCTION_A;
  m.attr("ROGERS_HERSH_1976_SUCTION_B") = combaero::ROGERS_HERSH_1976_SUCTION_B;
  py::class_<combaero::CrossflowCd>(m, "CrossflowCd")
      .def_readonly("Cd", &combaero::CrossflowCd::Cd)
      .def_readonly("dCd_dG_ideal", &combaero::CrossflowCd::dCd_dG_ideal)
      .def_readonly("dCd_dG_c", &combaero::CrossflowCd::dCd_dG_c)
      .def_readonly("capped", &combaero::CrossflowCd::capped);
  m.def("rogers_hersh_1976_cd", &combaero::rogers_hersh_1976_cd, py::arg("Cd0"),
        py::arg("G_ideal"), py::arg("G_c"), py::arg("suction"),
        "Discharge coefficient of a wall hole under a grazing crossflow, Rogers\n"
        "& Hersh (1976) NASA CR-2681 Fig. 11, resolved in q = G_ideal / G_c\n"
        "(ideal hole over crossflow mass velocity) and capped at the hole's own\n"
        "Cd0. suction: the hole draws from the grazing stream (a reversed jet).\n"
        "Analytic dCd/dG_ideal, dCd/dG_c.");
  m.attr("ROGERS_HERSH_1976_VELOCITY_RATIO_MIN") = combaero::ROGERS_HERSH_1976_VELOCITY_RATIO_MIN;
  m.attr("ROGERS_HERSH_1976_VELOCITY_RATIO_MAX") = combaero::ROGERS_HERSH_1976_VELOCITY_RATIO_MAX;
  m.attr("ROGERS_HERSH_1976_L_OVER_D_MIN") = combaero::ROGERS_HERSH_1976_L_OVER_D_MIN;
  m.attr("ROGERS_HERSH_1976_L_OVER_D_MAX") = combaero::ROGERS_HERSH_1976_L_OVER_D_MAX;
  m.attr("ROGERS_HERSH_1976_INJECTION_L_OVER_D_RATIO_MAX") =
      combaero::ROGERS_HERSH_1976_INJECTION_L_OVER_D_RATIO_MAX;
  m.attr("ROGERS_HERSH_1976_SUCTION_L_OVER_D_MIN") = combaero::ROGERS_HERSH_1976_SUCTION_L_OVER_D_MIN;
  m.attr("ROGERS_HERSH_1976_SUCTION_L_OVER_D_MAX") = combaero::ROGERS_HERSH_1976_SUCTION_L_OVER_D_MAX;
  m.attr("ROGERS_HERSH_1976_V_INF_MIN") = combaero::ROGERS_HERSH_1976_V_INF_MIN;
  m.attr("ROGERS_HERSH_1976_V_INF_MAX") = combaero::ROGERS_HERSH_1976_V_INF_MAX;
  py::class_<combaero::RogersHersh1976Envelope>(m, "RogersHersh1976Envelope")
      .def_readonly("velocity_ratio", &combaero::RogersHersh1976Envelope::velocity_ratio)
      .def_readonly("l_over_d", &combaero::RogersHersh1976Envelope::l_over_d)
      .def_readonly("grazing_velocity", &combaero::RogersHersh1976Envelope::grazing_velocity)
      .def("in_range", &combaero::RogersHersh1976Envelope::in_range);
  m.def("rogers_hersh_1976_envelope", &combaero::rogers_hersh_1976_envelope,
        py::arg("Vi_over_V_inf"), py::arg("L_over_D"), py::arg("V_inf"), py::arg("suction"),
        "Which inputs lie outside Rogers & Hersh (1976)'s envelope, by variable\n"
        "and branch: velocity_ratio (Vi/V_inf outside 0.004-0.5, Fig. 11),\n"
        "l_over_d (injection: L/D outside 0.12-2.0, or outside 0.6-1.0 above\n"
        "Vi/V_inf 0.35; suction: outside 0.6-1.0, around Fig. 11's L/D 1),\n"
        "grazing_velocity (V_inf outside 24-85 m/s). Information to report,\n"
        "never a selector: rogers_hersh_1976_cd returns the fit anyway.");
  m.def("florschuetz_metzger_su_1982_has_plate",
        &combaero::florschuetz_metzger_su_1982_has_plate, py::arg("xn_d"), py::arg("yn_d"),
        py::arg("z_d"), "True when CR-3630 Table 4.1 covers the jet plate.");
  m.def("florschuetz_metzger_su_1982_cd", &combaero::florschuetz_metzger_su_1982_cd,
        py::arg("xn_d"), py::arg("yn_d"), py::arg("z_d"), py::arg("G_ideal"), py::arg("G_c"),
        "CR-3630 Table 4.1 / Eq. 4.2 Cd of a B-size impingement plate under its\n"
        "crossflow, resolved in q = G_ideal / G_c; the table's own Cd (no Cd0).\n"
        "Raises for a plate the table does not cover. Analytic derivatives.");
  py::enum_<combaero::Rowbury2001Set>(
      m, "Rowbury2001Set",
      "Rowbury, Oldfield & Lock (2001) Table 2 coefficient set: vane region\n"
      "and hole shape. A flat wall takes the suction-surface set, as the\n"
      "authors do for Gritsch et al. (1998) in their Fig. 9.")
      .value("CylindricalPressureSurface", combaero::Rowbury2001Set::CylindricalPressureSurface)
      .value("CylindricalLeadingEdge", combaero::Rowbury2001Set::CylindricalLeadingEdge)
      .value("CylindricalSuctionSurface", combaero::Rowbury2001Set::CylindricalSuctionSurface)
      .value("FanShapedPressureSurface", combaero::Rowbury2001Set::FanShapedPressureSurface)
      .value("FanShapedSuctionSurface", combaero::Rowbury2001Set::FanShapedSuctionSurface);
  m.attr("ROWBURY_2001_K_PERP_MIN") = combaero::ROWBURY_2001_K_PERP_MIN;
  m.attr("ROWBURY_2001_UNIQUE_MARGIN") = combaero::ROWBURY_2001_UNIQUE_MARGIN;
  m.attr("ROWBURY_2001_FADE_WIDTH") = combaero::ROWBURY_2001_FADE_WIDTH;
  py::class_<combaero::Rowbury2001Coefficients>(m, "Rowbury2001Coefficients")
      .def_readonly("A", &combaero::Rowbury2001Coefficients::A)
      .def_readonly("B", &combaero::Rowbury2001Coefficients::B)
      .def_readonly("C", &combaero::Rowbury2001Coefficients::C);
  py::class_<combaero::Rowbury2001Cd>(m, "Rowbury2001Cd")
      .def_readonly("Cd", &combaero::Rowbury2001Cd::Cd)
      .def_readonly("dCd_dCd0", &combaero::Rowbury2001Cd::dCd_dCd0)
      .def_readonly("dCd_dk", &combaero::Rowbury2001Cd::dCd_dk)
      .def_readonly("delta", &combaero::Rowbury2001Cd::delta)
      .def_readonly("k_used", &combaero::Rowbury2001Cd::k_used)
      .def_readonly("held", &combaero::Rowbury2001Cd::held)
      .def_readonly("faded", &combaero::Rowbury2001Cd::faded);
  py::class_<combaero::Rowbury2001ExitDensity>(m, "Rowbury2001ExitDensity")
      .def_readonly("rho", &combaero::Rowbury2001ExitDensity::rho)
      .def_readonly("T", &combaero::Rowbury2001ExitDensity::T)
      .def_readonly("drho_dPt", &combaero::Rowbury2001ExitDensity::drho_dPt)
      .def_readonly("drho_dTt", &combaero::Rowbury2001ExitDensity::drho_dTt)
      .def_readonly("drho_dP", &combaero::Rowbury2001ExitDensity::drho_dP)
      .def_readonly("drho_dY", &combaero::Rowbury2001ExitDensity::drho_dY);
  py::class_<combaero::Rowbury2001KPerp>(m, "Rowbury2001KPerp")
      .def_readonly("k", &combaero::Rowbury2001KPerp::k)
      .def_readonly("dk_dG", &combaero::Rowbury2001KPerp::dk_dG)
      .def_readonly("dk_drho_h", &combaero::Rowbury2001KPerp::dk_drho_h)
      .def_readonly("dk_dq_m", &combaero::Rowbury2001KPerp::dk_dq_m);
  m.def("rowbury_2001_coefficients", &combaero::rowbury_2001_coefficients, py::arg("set"),
        py::arg("theta_deg"), py::arg("alpha_deg"),
        "Rowbury, Oldfield & Lock (2001) Table 2: (A, B, C) of\n"
        "delta_out = A k_perp^-B - C at hole inclination theta (to the surface)\n"
        "and orientation alpha [deg]. Raises where A <= 0 or B <= 0.");
  m.def("rowbury_2001_k_floor", &combaero::rowbury_2001_k_floor, py::arg("set"),
        py::arg("theta_deg"), py::arg("alpha_deg"),
        "k_perp below which the exit-crossflow loss is held:\n"
        "max(ROWBURY_2001_K_PERP_MIN, k_unique), k_unique (B > 1 only) keeping\n"
        "the hole's flow a unique root.");
  m.def("rowbury_2001_k_fit_max", &combaero::rowbury_2001_k_fit_max, py::arg("set"),
        "k_fit = max(k_data, k_applied): beyond it the offset C fades. k_data,\n"
        "the k_perp the set's Table 1 rows provably reach; k_applied, the\n"
        "largest at which the authors applied the set (Fig. 9, suction).");
  m.def("rowbury_2001_k_data_max", &combaero::rowbury_2001_k_data_max, py::arg("set"),
        "Tier 1 of k_fit: the k_perp the set's Table 1 rows provably reach.");
  m.def("rowbury_2001_in_range", &combaero::rowbury_2001_in_range, py::arg("set"),
        py::arg("theta_deg"), py::arg("alpha_deg"),
        "True inside the set's Table 1 (theta, alpha) span.");
  m.def("rowbury_2001_exit_crossflow_cd", &combaero::rowbury_2001_exit_crossflow_cd,
        py::arg("Cd0"), py::arg("k_perp"), py::arg("set"), py::arg("theta_deg"),
        py::arg("alpha_deg"),
        "Cd of a film hole under an EXIT crossflow, Rowbury, Oldfield & Lock\n"
        "(2001) Eqs. 7-8: 1/Cd^2 = 1/Cd0^2 + A k^-B - C, Cd0 the no-crossflow\n"
        "Cd, k = k_perp held at rowbury_2001_k_floor; beyond\n"
        "rowbury_2001_k_fit_max(set) the offset C fades (Cd -> Cd0 as k -> inf).\n"
        "Analytic dCd/dCd0, dCd/dk.");
  m.def("rowbury_2001_k_perp", &combaero::rowbury_2001_k_perp, py::arg("G_h"),
        py::arg("rho_h"), py::arg("rho_m_u_m2"), py::arg("theta_deg"),
        py::arg("alpha_deg") = 0.0,
        "Perpendicular outlet momentum flux ratio k_perp = G_h^2 f /\n"
        "(rho_h rho_m u_m^2), f = 1/sqrt(cot^2 theta + 1 + tan^2 alpha) (Fig. 3),\n"
        "and its partials.");
  m.def("rowbury_2001_exit_density", &combaero::rowbury_2001_exit_density, py::arg("Pt"),
        py::arg("Tt"), py::arg("P"), py::arg("X"),
        "Coolant exit-plane static density for k_perp: isentropic expansion\n"
        "from (Pt, Tt) to the discharge static P, variable cp; analytic\n"
        "partials, drho_dY over mass fractions.");
  m.attr("STATION_KAPPA_MERGE_NORMAL") = combaero::STATION_KAPPA_MERGE_NORMAL;
  m.attr("STATION_KAPPA_BLEED_BASSETT") = combaero::STATION_KAPPA_BLEED_BASSETT;
  m.attr("WALL_HTC_KNEE") = combaero::WALL_HTC_KNEE;
  m.attr("MIXER_HEAT_MDOT_FLOOR") = combaero::solver::kMixerHeatMdotFloor;
  m.def("mixer_heat_per_unit_flow", &combaero::solver::mixer_heat_per_unit_flow,
        py::arg("mdot_tot"),
        "(phi, dphi/dmdot): the enthalpy rise per unit heat a mixing node\n"
        "applies at total flow mdot_tot. 1/|m| at or above\n"
        "MIXER_HEAT_MDOT_FLOOR; below it the heat taken up, Q |m| phi, fades\n"
        "C1 to none at zero flow.");
  m.attr("SENSIBLE_ENTHALPY_REF_T") = combaero::solver::kSensibleEnthalpyRefT;
  py::class_<combaero::StationHalfDrop>(m, "StationHalfDrop")
      .def_readonly("dP", &combaero::StationHalfDrop::dP)
      .def_readonly("d_dm_a", &combaero::StationHalfDrop::d_dm_a)
      .def_readonly("d_dm_b", &combaero::StationHalfDrop::d_dm_b)
      .def_readonly("d_dP", &combaero::StationHalfDrop::d_dP)
      .def_readonly("d_dT", &combaero::StationHalfDrop::d_dT);
  m.def("station_half_drop", &combaero::station_half_drop, py::arg("m_a"), py::arg("m_b"),
        py::arg("P"), py::arg("T"), py::arg("X"), py::arg("area"), py::arg("kappa"),
        "Half the static-pressure drop over one side-stream station (m_a\n"
        "arriving, m_b leaving along a channel of area A), at the station\n"
        "node's density: (m_b^2 - m_a^2 + kappa m_a (m_a - m_b)) / (2 rho A^2),\n"
        "any flow direction (#815). kappa = 0: normal injection (merge);\n"
        "kappa = 0.75: bleed, Bassett et al. (2001) K2/K5 exactly.");
  m.def("segment_area_change_drop", &combaero::segment_area_change_drop, py::arg("m"),
        py::arg("P"), py::arg("T"), py::arg("X"), py::arg("area_from"), py::arg("area_to"),
        "Static drop over a duct segment from its area change alone (#812):\n"
        "m^2 / (2 rho) (1/A_to^2 - 1/A_from^2), Bernoulli, even in m. Zero for\n"
        "equal areas. Derivative in d_dm_a.");
  m.def("channel_entry_drop", &combaero::channel_entry_drop, py::arg("m"), py::arg("P"),
        py::arg("T"), py::arg("X"), py::arg("area"), py::arg("K_in"),
        "Static drop entering a duct from a reservoir: (1 + K_in) m|m| /\n"
        "(2 rho A^2), at the duct density. Derivative in d_dm_a.");

  m.attr("CROSSFLOW_VI_FLOOR") = combaero::CROSSFLOW_VI_FLOOR;
  py::class_<combaero::CrossflowRatio>(m, "CrossflowRatio")
      .def_readonly("U1_over_Vi", &combaero::CrossflowRatio::U1_over_Vi)
      .def_readonly("U1", &combaero::CrossflowRatio::U1)
      .def_readonly("Vi", &combaero::CrossflowRatio::Vi)
      .def_readonly("d_dm_a", &combaero::CrossflowRatio::d_dm_a)
      .def_readonly("d_dm_b", &combaero::CrossflowRatio::d_dm_b)
      .def_readonly("d_dP", &combaero::CrossflowRatio::d_dP)
      .def_readonly("d_dT", &combaero::CrossflowRatio::d_dT)
      .def_readonly("d_dP_down", &combaero::CrossflowRatio::d_dP_down);
  m.def("crossflow_velocity_ratio", &combaero::crossflow_velocity_ratio, py::arg("m_a"),
        py::arg("m_b"), py::arg("P"), py::arg("T"), py::arg("X"), py::arg("P_down"),
        py::arg("area"),
        "McGreehan-Schotsch's supply-side U1/Vi for a hole in a duct wall: U1\n"
        "the duct's mean velocity at the hole's station, (m_a + m_b)/(2 rho A);\n"
        "Vi the isentropic jet velocity from the duct's STATIC (P, T) to the\n"
        "discharge static P_down. Analytic derivatives.");

  py::class_<combaero::FannoMachState>(m, "FannoMachState")
      .def_readonly("M", &combaero::FannoMachState::M)
      .def_readonly("T", &combaero::FannoMachState::T)
      .def_readonly("P", &combaero::FannoMachState::P)
      .def_readonly("rho", &combaero::FannoMachState::rho)
      .def_readonly("u", &combaero::FannoMachState::u)
      .def_readonly("a", &combaero::FannoMachState::a);
  py::class_<combaero::FannoDuctResult>(m, "FannoDuctResult")
      .def_readonly("choked", &combaero::FannoDuctResult::choked)
      .def_readonly("G", &combaero::FannoDuctResult::G)
      .def_readonly("L_star", &combaero::FannoDuctResult::L_star)
      .def_readonly("inlet", &combaero::FannoDuctResult::inlet)
      .def_readonly("exit", &combaero::FannoDuctResult::exit)
      .def_readonly("Pt_exit", &combaero::FannoDuctResult::Pt_exit);
  m.attr("FANNO_MACH_REL_TOL") = combaero::kFannoMachRelTol;
  m.def("fanno_state_at_mach", &combaero::fanno_state_at_mach, py::arg("G"), py::arg("Tt"),
        py::arg("M"), py::arg("X"),
        "Fanno station state at Mach M for mass flux G and stagnation T (algebraic).");
  m.def("fanno_dx_dmach", &combaero::fanno_dx_dmach, py::arg("G"), py::arg("Tt"), py::arg("M"),
        py::arg("X"), py::arg("D"), py::arg("roughness"), py::arg("friction_model"),
        py::arg("f_multiplier") = 1.0,
        "Fanno dx/dM [m]; friction_model 'fixed' means a constant Darcy f = f_multiplier.");
  m.def("fanno_length_between", &combaero::fanno_length_between, py::arg("G"), py::arg("Tt"),
        py::arg("M1"), py::arg("M2"), py::arg("X"), py::arg("D"), py::arg("roughness"),
        py::arg("friction_model"), py::arg("f_multiplier") = 1.0,
        "Duct length over which Fanno flow goes from M1 to M2 (Mach quadrature).");
  m.def("fanno_inlet_mach", &combaero::fanno_inlet_mach, py::arg("Pt"), py::arg("Tt"),
        py::arg("G"), py::arg("X"),
        "Subsonic Mach of an isentropic inlet from (Pt, Tt) at mass flux G (1 if choked).");
  m.def("fanno_sonic_mass_flux", &combaero::fanno_sonic_mass_flux, py::arg("Pt"), py::arg("Tt"),
        py::arg("X"), "Sonic isentropic mass flux of a stagnation state [kg/(m^2 s)].");
  m.def("fanno_duct", &combaero::fanno_duct, py::arg("Pt"), py::arg("Tt"), py::arg("G"),
        py::arg("X"), py::arg("L"), py::arg("D"), py::arg("roughness"),
        py::arg("friction_model"), py::arg("f_multiplier") = 1.0,
        "Fanno duct fed isentropically from (Pt, Tt) at mass flux G: inlet and exit\n"
        "states, L*, choked flag. Integrated in Mach number, exact to L* = L.");
  m.def("fanno_choked_mass_flux", &combaero::fanno_choked_mass_flux, py::arg("Pt"), py::arg("Tt"),
        py::arg("X"), py::arg("L"), py::arg("D"), py::arg("roughness"),
        py::arg("friction_model"), py::arg("f_multiplier") = 1.0,
        "Mass flux [kg/(m^2 s)] at which a duct of length L chokes exactly at its exit.");

  py::class_<combaero::FannoChannelFlow>(m, "FannoChannelFlow")
      .def_readonly("G", &combaero::FannoChannelFlow::G)
      .def_readonly("choked", &combaero::FannoChannelFlow::choked)
      .def_readonly("M_in", &combaero::FannoChannelFlow::M_in)
      .def_readonly("M_exit", &combaero::FannoChannelFlow::M_exit)
      .def_readonly("M_exit_root", &combaero::FannoChannelFlow::M_exit_root)
      .def_readonly("dG_dPt0", &combaero::FannoChannelFlow::dG_dPt0)
      .def_readonly("dG_dTt0", &combaero::FannoChannelFlow::dG_dTt0)
      .def_readonly("dG_dP_target", &combaero::FannoChannelFlow::dG_dP_target)
      .def_readonly("dG_dTt_exit", &combaero::FannoChannelFlow::dG_dTt_exit)
      .def_readonly("dG_ddir", &combaero::FannoChannelFlow::dG_ddir);
  m.attr("FANNO_FLOW_DRIVE_FLOOR") = combaero::kFannoFlowDriveFloor;
  m.def("fanno_channel_flow", &combaero::fanno_channel_flow, py::arg("Pt0"), py::arg("Tt0"),
        py::arg("X"), py::arg("P_target"), py::arg("exit_total"), py::arg("L"), py::arg("D"),
        py::arg("roughness"), py::arg("friction_model"), py::arg("f_multiplier") = 1.0,
        py::arg("M_exit_guess") = -1.0,
        py::arg("dY_directions") = std::vector<std::vector<double>>{},
        py::arg("with_derivatives") = true,
        "Mass flux a duct passes from (Pt0, Tt0) against P_target (matched to the exit\n"
        "stagnation pressure if exit_total, else the exit static), with derivatives.\n"
        "Saturates at the choked flux; 0 at or below zero drive.");

  py::class_<combaero::FannoRayleighMarch>(m, "FannoRayleighMarch")
      .def_readonly("M_in", &combaero::FannoRayleighMarch::M_in)
      .def_readonly("M_exit", &combaero::FannoRayleighMarch::M_exit)
      .def_readonly("zeta_min", &combaero::FannoRayleighMarch::zeta_min)
      .def_readonly("zeta_peak", &combaero::FannoRayleighMarch::zeta_peak)
      .def_readonly("P_exit", &combaero::FannoRayleighMarch::P_exit)
      .def_readonly("Pt_exit", &combaero::FannoRayleighMarch::Pt_exit)
      .def_readonly("T_exit", &combaero::FannoRayleighMarch::T_exit);
  m.def("fanno_rayleigh_march", &combaero::fanno_rayleigh_march, py::arg("G"), py::arg("Pt0"),
        py::arg("Tt0"), py::arg("Tt_exit"), py::arg("X"), py::arg("L"), py::arg("D"),
        py::arg("roughness"), py::arg("friction_model"), py::arg("f_multiplier") = 1.0,
        "March a duct with friction and a linear Tt rise from Tt0 to Tt_exit at mass flux G\n"
        "(Fanno-Rayleigh, #551): inlet/exit Mach, exit P/Pt/T, zeta_min <= 0 when choked.");
  m.def("fanno_rayleigh_choked_mass_flux", &combaero::fanno_rayleigh_choked_mass_flux,
        py::arg("Pt0"), py::arg("Tt0"), py::arg("Tt_exit"), py::arg("X"), py::arg("L"),
        py::arg("D"), py::arg("roughness"), py::arg("friction_model"),
        py::arg("f_multiplier") = 1.0,
        "Mass flux [kg/(m^2 s)] at which a duct with friction and heat chokes.");
  m.def("fanno_rayleigh_channel_flow", &combaero::fanno_rayleigh_channel_flow, py::arg("Pt0"),
        py::arg("Tt0"), py::arg("Tt_exit"), py::arg("X"), py::arg("P_target"),
        py::arg("exit_total"), py::arg("L"), py::arg("D"), py::arg("roughness"),
        py::arg("friction_model"), py::arg("f_multiplier") = 1.0, py::arg("G_guess") = -1.0,
        py::arg("dY_directions") = std::vector<std::vector<double>>{},
        py::arg("with_derivatives") = true,
        "fanno_channel_flow with the stagnation temperature rising linearly from Tt0 to\n"
        "Tt_exit along the duct (Fanno-Rayleigh, #551); adds dG_dTt_exit.");

  py::class_<combaero::MergeFaceState>(m, "MergeFaceState")
      .def_readonly("P_face", &combaero::MergeFaceState::P_face)
      .def_readonly("Pt_face", &combaero::MergeFaceState::Pt_face)
      .def_readonly("M_face", &combaero::MergeFaceState::M_face)
      .def_readonly("choked", &combaero::MergeFaceState::choked)
      .def_readonly("dPf_dm_main", &combaero::MergeFaceState::dPf_dm_main)
      .def_readonly("dPf_dm_out", &combaero::MergeFaceState::dPf_dm_out)
      .def_readonly("dPf_dJ", &combaero::MergeFaceState::dPf_dJ)
      .def_readonly("dPf_dP", &combaero::MergeFaceState::dPf_dP)
      .def_readonly("dPf_dT", &combaero::MergeFaceState::dPf_dT)
      .def_readonly("dPf_dT_main", &combaero::MergeFaceState::dPf_dT_main)
      .def_readonly("dPtf_dm_main", &combaero::MergeFaceState::dPtf_dm_main)
      .def_readonly("dPtf_dm_out", &combaero::MergeFaceState::dPtf_dm_out)
      .def_readonly("dPtf_dJ", &combaero::MergeFaceState::dPtf_dJ)
      .def_readonly("dPtf_dP", &combaero::MergeFaceState::dPtf_dP)
      .def_readonly("dPtf_dT", &combaero::MergeFaceState::dPtf_dT)
      .def_readonly("dPtf_dT_main", &combaero::MergeFaceState::dPtf_dT_main)
      .def_readonly("dPf_dY", &combaero::MergeFaceState::dPf_dY)
      .def_readonly("dPf_dY_main", &combaero::MergeFaceState::dPf_dY_main)
      .def_readonly("dPtf_dY", &combaero::MergeFaceState::dPtf_dY)
      .def_readonly("dPtf_dY_main", &combaero::MergeFaceState::dPtf_dY_main);
  m.def("chamber_merge_face_state", &combaero::chamber_merge_face_state,
        py::arg("m_main"), py::arg("Tt_main"), py::arg("X_main"), py::arg("m_out"),
        py::arg("P"), py::arg("Tt"), py::arg("X"), py::arg("side_momentum"),
        py::arg("area"), py::arg("incompressible") = false,
        "Main-inlet face of a merge chamber from the exact (compressible)\n"
        "impulse balance over a constant-area control volume: face static P\n"
        "on the subsonic root of the ideal-gas quadratic, face Pt from the\n"
        "same entropy-based stagnation closure as the chamber. P: the\n"
        "chamber (outlet) static pressure; Tt, Tt_main: the STAGNATION\n"
        "temperatures of the chamber and of the main stream (the face and\n"
        "outlet static T are solved inside; the dPf_dT/dPf_dT_main and\n"
        "dPtf_dT/dPtf_dT_main fields are with respect to Tt and Tt_main);\n"
        "side_momentum: sum of the side streams' streamwise impulse.\n"
        "incompressible: the low-Mach limit, densities at the chamber P and\n"
        "the stagnation temperatures; linear in P_face, no choking (#588).");

  py::class_<combaero::DuctStaticState>(m, "DuctStaticState")
      .def_readonly("T", &combaero::DuctStaticState::T)
      .def_readonly("P", &combaero::DuctStaticState::P)
      .def_readonly("rho", &combaero::DuctStaticState::rho)
      .def_readonly("M", &combaero::DuctStaticState::M)
      .def_readonly("dT_dTt", &combaero::DuctStaticState::dT_dTt)
      .def_readonly("dT_dPt", &combaero::DuctStaticState::dT_dPt)
      .def_readonly("dT_dG", &combaero::DuctStaticState::dT_dG)
      .def_readonly("dP_dTt", &combaero::DuctStaticState::dP_dTt)
      .def_readonly("dP_dPt", &combaero::DuctStaticState::dP_dPt)
      .def_readonly("dP_dG", &combaero::DuctStaticState::dP_dG)
      .def_readonly("drho_dTt", &combaero::DuctStaticState::drho_dTt)
      .def_readonly("drho_dPt", &combaero::DuctStaticState::drho_dPt)
      .def_readonly("drho_dG", &combaero::DuctStaticState::drho_dG)
      .def_readonly("choked", &combaero::DuctStaticState::choked)
      .def_readonly("G_eff", &combaero::DuctStaticState::G_eff)
      .def_readonly("dGeff_dTt", &combaero::DuctStaticState::dGeff_dTt)
      .def_readonly("dGeff_dPt", &combaero::DuctStaticState::dGeff_dPt)
      .def_readonly("dGeff_dG", &combaero::DuctStaticState::dGeff_dG);
  m.def("duct_static_state", &combaero::duct_static_state, py::arg("Tt"), py::arg("Pt"),
        py::arg("G"), py::arg("X"),
        "Static state at a duct inlet fed from stagnation (Tt, Pt) at mass flux\n"
        "G = |m_dot| / A: the isentropic subsonic root, or the sonic state when G\n"
        "exceeds the choked flux (choked = True). Partials in Tt, Pt, G (#636).");
  m.attr("DUCT_HTC_SONIC_MARGIN") = combaero::kDuctHtcSonicMargin;
  m.def("duct_htc_static_state", &combaero::duct_htc_static_state, py::arg("Tt"), py::arg("Pt"),
        py::arg("G"), py::arg("X"),
        "The duct state a heat-transfer correlation reads: duct_static_state at a\n"
        "flux bent C1 below the sonic one (G itself below (1 - 2w) G*, capped at\n"
        "(1 - w) G*, w = DUCT_HTC_SONIC_MARGIN), so it is smooth through choke (#551).");

  py::class_<combaero::DuctStaticStateDY>(m, "DuctStaticStateDY")
      .def_readonly("dT_dY", &combaero::DuctStaticStateDY::dT_dY)
      .def_readonly("dP_dY", &combaero::DuctStaticStateDY::dP_dY)
      .def_readonly("drho_dY", &combaero::DuctStaticStateDY::drho_dY)
      .def_readonly("dGeff_dY", &combaero::DuctStaticStateDY::dGeff_dY);
  m.def("duct_htc_static_state_dY", &combaero::duct_htc_static_state_dY, py::arg("Tt"),
        py::arg("Pt"), py::arg("G"), py::arg("X"), py::arg("state"),
        "duct_htc_static_state's gradients in the MASS fractions at fixed Tt, Pt\n"
        "and G (#717): dT_dY, dP_dY, drho_dY and the bent flux's dGeff_dY, one\n"
        "entry per species. Implicit-function derivatives of the isentropic\n"
        "energy and entropy balances. `state` is duct_htc_static_state(Tt, Pt,\n"
        "G, X), already solved: below the bend nothing is solved again.");

  py::class_<combaero::JetImpulse>(m, "JetImpulse")
      .def_readonly("J", &combaero::JetImpulse::J)
      .def_readonly("dJ_dm", &combaero::JetImpulse::dJ_dm)
      .def_readonly("dJ_dPt", &combaero::JetImpulse::dJ_dPt)
      .def_readonly("dJ_dTt", &combaero::JetImpulse::dJ_dTt)
      .def_readonly("dJ_dP", &combaero::JetImpulse::dJ_dP)
      .def_readonly("choked", &combaero::JetImpulse::choked);
  m.def("jet_impulse", &combaero::jet_impulse, py::arg("m"), py::arg("Pt"),
        py::arg("Tt"), py::arg("P"), py::arg("X"),
        "Streamwise impulse of a jet from stagnation (Pt, Tt) into static P,\n"
        "per unit cos(theta): isentropic velocity (variable cp), or the sonic\n"
        "momentum plus pressure thrust when choked. Cd does not enter.");

  // ---------------------------------------------------------------------
  // Pin-fin arrays (#335): provenanced sets, exact basis converters, fin
  // model. Canonical basis: Re_D on D and Vmax at A_min, Nu_D = h D / k,
  // f = dP / (2 rho Vmax^2 N).
  // ---------------------------------------------------------------------
  {
    namespace cc = combaero::cooling;
    m.attr("PIN_FIN_RE_EPS") = cc::PIN_FIN_RE_EPS;
    m.attr("PIN_FIN_F_BLEND") = cc::PIN_FIN_F_BLEND;

    py::enum_<cc::PinArrangement>(m, "PinArrangement")
        .value("Staggered", cc::PinArrangement::Staggered)
        .value("Inline", cc::PinArrangement::Inline);
    py::enum_<cc::PinNuSurface>(m, "PinNuSurface",
                                "Which surface a set's Nu describes.")
        .value("Total", cc::PinNuSurface::Total,
               "Pin and endwall, area-weighted (what an array element uses)")
        .value("Endwall", cc::PinNuSurface::Endwall)
        .value("Pin", cc::PinNuSurface::Pin);
    py::enum_<cc::PinReBasis>(m, "PinReBasis")
        .value("DiameterVmax", cc::PinReBasis::DiameterVmax,
               "Re_D = mdot D / (mu A_min) (canonical)")
        .value("VanFossenDprime", cc::PinReBasis::VanFossenDprime,
               "Re_D' = mdot D' / (mu A'), D' = 4V/A_t, A' = V/L");
    py::enum_<cc::PinFrictionBasis>(m, "PinFrictionBasis")
        .value("PerRowVmax", cc::PinFrictionBasis::PerRowVmax,
               "f = dP / (2 rho Vmax^2 N) (canonical, Metzger)")
        .value("PerRowGapVmax", cc::PinFrictionBasis::PerRowGapVmax,
               "f = dP / (2 rho Vmax^2 (N - 1)) (Damerow 1972)")
        .value("PerRowHalfHeadVmax", cc::PinFrictionBasis::PerRowHalfHeadVmax,
               "f = 2 dP / (rho Vmax^2 N), 4x canonical (Chyu 1990)");

    py::class_<cc::PinFinGeometry>(m, "PinFinGeometry",
                                   "Pin array geometry, pitches and height "
                                   "over the pin diameter.")
        .def(py::init<>())
        .def(py::init([](double S_D, double X_D, double H_D, int N_rows,
                         cc::PinArrangement arrangement) {
               cc::PinFinGeometry g;
               g.S_D = S_D;
               g.X_D = X_D;
               g.H_D = H_D;
               g.N_rows = N_rows;
               g.arrangement = arrangement;
               return g;
             }),
             py::arg("S_D") = 2.5, py::arg("X_D") = 2.5, py::arg("H_D") = 1.0,
             py::arg("N_rows") = 10,
             py::arg("arrangement") = cc::PinArrangement::Staggered)
        .def_readwrite("S_D", &cc::PinFinGeometry::S_D)
        .def_readwrite("X_D", &cc::PinFinGeometry::X_D)
        .def_readwrite("H_D", &cc::PinFinGeometry::H_D)
        .def_readwrite("N_rows", &cc::PinFinGeometry::N_rows)
        .def_readwrite("arrangement", &cc::PinFinGeometry::arrangement);

    m.def("pin_fin_min_gap_D", &cc::pin_fin_min_gap_D, py::arg("geometry"),
          "Minimum gap per transverse pitch / D: min(S - 1, diagonal) "
          "staggered, S - 1 inline.");
    m.def("pin_fin_amin_over_afrontal", &cc::pin_fin_amin_over_afrontal,
          py::arg("geometry"), "A_min / A_frontal = Vchannel / Vmax.");
    m.def("pin_fin_dprime_over_D", &cc::pin_fin_dprime_over_D,
          py::arg("geometry"), "VanFossen D' / D (Armstrong-Winstanley Eq. 8).");
    m.def("pin_fin_aprime_over_amin", &cc::pin_fin_aprime_over_amin,
          py::arg("geometry"), "VanFossen A' / A_min (Eq. 10).");
    m.def("pin_fin_dh_over_D", &cc::pin_fin_dh_over_D, py::arg("geometry"),
          "Tube-bank D_h / D = 4 A_min L / (A_t D) (Eq. 13).");

    py::class_<cc::PinFinAreaFractions>(m, "PinFinAreaFractions")
        .def_readonly("endwall_exposed",
                      &cc::PinFinAreaFractions::endwall_exposed)
        .def_readonly("pin", &cc::PinFinAreaFractions::pin)
        .def_readonly("pin_over_total",
                      &cc::PinFinAreaFractions::pin_over_total);
    m.def("pin_fin_area_fractions", &cc::pin_fin_area_fractions,
          py::arg("geometry"),
          "Per-wall areas over the base (planform) area, pins to mid-height.");

    py::class_<cc::PinFinEfficiency>(m, "PinFinEfficiency")
        .def_readonly("eta_fin", &cc::PinFinEfficiency::eta_fin)
        .def_readonly("eta_t", &cc::PinFinEfficiency::eta_t)
        .def_readonly("deta_fin_dh", &cc::PinFinEfficiency::deta_fin_dh)
        .def_readonly("deta_t_dh", &cc::PinFinEfficiency::deta_t_dh);
    m.def("pin_fin_array_efficiency", &cc::pin_fin_array_efficiency,
          py::arg("h"), py::arg("k_pin"), py::arg("D"), py::arg("H"),
          py::arg("A_f_over_A_t"),
          "Fin and array efficiency, fin length H/2 (Eqs 14-16), with "
          "analytic d/dh.");

    py::class_<cc::PinFinNuSet>(m, "PinFinNuSet",
                                "Pin-fin Nu set in its source's own basis.")
        .def(py::init<>())
        .def_readwrite("name", &cc::PinFinNuSet::name)
        .def_readwrite("source", &cc::PinFinNuSet::source)
        .def_readwrite("validity_source", &cc::PinFinNuSet::validity_source)
        .def_readwrite("provenance", &cc::PinFinNuSet::provenance)
        .def_readwrite("arrangement", &cc::PinFinNuSet::arrangement)
        .def_readwrite("surface", &cc::PinFinNuSet::surface)
        .def_readwrite("re_basis", &cc::PinFinNuSet::re_basis)
        .def_readwrite("C", &cc::PinFinNuSet::C)
        .def_readwrite("Re_exp", &cc::PinFinNuSet::Re_exp)
        .def_readwrite("Pr_exp", &cc::PinFinNuSet::Pr_exp)
        .def_readwrite("term_XD", &cc::PinFinNuSet::term_XD)
        .def_readwrite("term_SD", &cc::PinFinNuSet::term_SD)
        .def_readwrite("term_HD", &cc::PinFinNuSet::term_HD)
        .def_readwrite("valid_Re", &cc::PinFinNuSet::valid_Re)
        .def_readwrite("valid_SD", &cc::PinFinNuSet::valid_SD)
        .def_readwrite("valid_XD", &cc::PinFinNuSet::valid_XD)
        .def_readwrite("valid_HD", &cc::PinFinNuSet::valid_HD)
        .def_readwrite("valid_Nrows", &cc::PinFinNuSet::valid_Nrows)
        .def_readwrite("valid_Pr", &cc::PinFinNuSet::valid_Pr)
        .def_readwrite("accuracy_Nu", &cc::PinFinNuSet::accuracy_Nu);

    py::class_<cc::PinFinFrictionSet>(m, "PinFinFrictionSet",
                                      "Pin-fin friction set, one or two Re "
                                      "segments blended C1 in ln Re.")
        .def(py::init<>())
        .def_readwrite("name", &cc::PinFinFrictionSet::name)
        .def_readwrite("source", &cc::PinFinFrictionSet::source)
        .def_readwrite("validity_source",
                       &cc::PinFinFrictionSet::validity_source)
        .def_readwrite("provenance", &cc::PinFinFrictionSet::provenance)
        .def_readwrite("arrangement", &cc::PinFinFrictionSet::arrangement)
        .def_readwrite("basis", &cc::PinFinFrictionSet::basis)
        .def_readwrite("C1", &cc::PinFinFrictionSet::C1)
        .def_readwrite("Re_exp1", &cc::PinFinFrictionSet::Re_exp1)
        .def_readwrite("C2", &cc::PinFinFrictionSet::C2)
        .def_readwrite("Re_exp2", &cc::PinFinFrictionSet::Re_exp2)
        .def_readwrite("Re_split", &cc::PinFinFrictionSet::Re_split)
        .def_readwrite("term_SD", &cc::PinFinFrictionSet::term_SD)
        .def_readwrite("term_XD", &cc::PinFinFrictionSet::term_XD)
        .def_readwrite("term_HD", &cc::PinFinFrictionSet::term_HD)
        .def_readwrite("valid_Re", &cc::PinFinFrictionSet::valid_Re)
        .def_readwrite("valid_SD", &cc::PinFinFrictionSet::valid_SD)
        .def_readwrite("valid_XD", &cc::PinFinFrictionSet::valid_XD)
        .def_readwrite("valid_HD", &cc::PinFinFrictionSet::valid_HD)
        .def_readwrite("valid_Nrows", &cc::PinFinFrictionSet::valid_Nrows)
        .def_readwrite("accuracy_f", &cc::PinFinFrictionSet::accuracy_f);

    py::class_<cc::PinFinRatioModifier>(
        m, "PinFinRatioModifier",
        "A data-backed ratio applied to another set's output.")
        .def(py::init<>())
        .def_readwrite("name", &cc::PinFinRatioModifier::name)
        .def_readwrite("source", &cc::PinFinRatioModifier::source)
        .def_readwrite("provenance", &cc::PinFinRatioModifier::provenance)
        .def_readwrite("from_arrangement",
                       &cc::PinFinRatioModifier::from_arrangement)
        .def_readwrite("to_arrangement",
                       &cc::PinFinRatioModifier::to_arrangement)
        .def_readwrite("surface", &cc::PinFinRatioModifier::surface)
        .def_readwrite("C_Nu", &cc::PinFinRatioModifier::C_Nu)
        .def_readwrite("Nu_Re_exp", &cc::PinFinRatioModifier::Nu_Re_exp)
        .def_readwrite("has_f", &cc::PinFinRatioModifier::has_f)
        .def_readwrite("C_f", &cc::PinFinRatioModifier::C_f)
        .def_readwrite("f_Re_exp", &cc::PinFinRatioModifier::f_Re_exp)
        .def_readwrite("valid_Re", &cc::PinFinRatioModifier::valid_Re)
        .def_readwrite("valid_SD", &cc::PinFinRatioModifier::valid_SD)
        .def_readwrite("valid_XD", &cc::PinFinRatioModifier::valid_XD)
        .def_readwrite("valid_HD", &cc::PinFinRatioModifier::valid_HD);

    py::class_<cc::PinFinNuResult>(m, "PinFinNuResult")
        .def_readonly("Nu", &cc::PinFinNuResult::Nu)
        .def_readonly("dNu_dRe", &cc::PinFinNuResult::dNu_dRe)
        .def_readonly("Re_native", &cc::PinFinNuResult::Re_native)
        .def_readonly("extrapolated", &cc::PinFinNuResult::extrapolated);
    py::class_<cc::PinFinFrictionResult>(m, "PinFinFrictionResult")
        .def_readonly("f", &cc::PinFinFrictionResult::f)
        .def_readonly("df_dRe", &cc::PinFinFrictionResult::df_dRe)
        .def_readonly("extrapolated",
                      &cc::PinFinFrictionResult::extrapolated);
    py::class_<cc::PinFinModifierResult>(m, "PinFinModifierResult")
        .def_readonly("ratio_Nu", &cc::PinFinModifierResult::ratio_Nu)
        .def_readonly("dratio_Nu_dRe",
                      &cc::PinFinModifierResult::dratio_Nu_dRe)
        .def_readonly("ratio_f", &cc::PinFinModifierResult::ratio_f)
        .def_readonly("dratio_f_dRe", &cc::PinFinModifierResult::dratio_f_dRe)
        .def_readonly("has_f", &cc::PinFinModifierResult::has_f)
        .def_readonly("extrapolated",
                      &cc::PinFinModifierResult::extrapolated);

    m.def("evaluate_pin_fin_nu", &cc::evaluate_pin_fin_nu,
          py::arg("correlation_set"), py::arg("geometry"), py::arg("Re_D"),
          py::arg("Pr"),
          "Canonical Nu_D (and analytic dNu/dRe_D) from a set in any basis. "
          "Never throws.");
    m.def("evaluate_pin_fin_friction", &cc::evaluate_pin_fin_friction,
          py::arg("correlation_set"), py::arg("geometry"), py::arg("Re_D"),
          "Canonical f = dP / (2 rho Vmax^2 N) and df/dRe_D. Never throws.");
    m.def("evaluate_pin_fin_modifier", &cc::evaluate_pin_fin_modifier,
          py::arg("modifier"), py::arg("geometry"), py::arg("Re_D"));
    m.def("validate_pin_fin_geometry", &cc::validate_pin_fin_geometry,
          py::arg("geometry"));
    m.def("validate_pin_fin_nu_set", &cc::validate_pin_fin_nu_set,
          py::arg("correlation_set"));
    m.def("validate_pin_fin_friction_set", &cc::validate_pin_fin_friction_set,
          py::arg("correlation_set"));
    m.def("validate_pin_fin_modifier", &cc::validate_pin_fin_modifier,
          py::arg("modifier"));

    m.def("metzger_1986_staggered_nu", &cc::metzger_1986_staggered_nu,
          "Metzger et al. (1986) via Armstrong-Winstanley Eq. 2. Default.");
    m.def("metzger_1982_staggered_friction",
          &cc::metzger_1982_staggered_friction,
          "Metzger, Fan & Shepard (1982) via Armstrong-Winstanley Eqs 20-21. "
          "Default.");
    m.def("vanfossen_1982_staggered_nu", &cc::vanfossen_1982_staggered_nu,
          "VanFossen (1982), NASA TM-81696 Eq. 16, D' basis.");
    m.def("damerow_1972_staggered_friction",
          &cc::damerow_1972_staggered_friction,
          "Damerow et al. (1972), NASA CR-120883 Eq. 18.");
    m.def("chyu_1998_nu", &cc::chyu_1998_nu, py::arg("arrangement"),
          py::arg("surface"),
          "Chyu et al. (1998), ASME 98-GT-175, via Han Table 4.7.");
    m.def("chyu_1990_nu", &cc::chyu_1990_nu, py::arg("arrangement"),
          py::arg("fillet") = false,
          "Chyu (1990) Table 2, pin surface, straight or fillet pins.");
    m.def("chyu_1990_friction", &cc::chyu_1990_friction,
          py::arg("arrangement"), py::arg("fillet") = false,
          "Chyu (1990) Fig. 6, digitised and fitted (provenance Fitted). The "
          "inline straight set is the default inline friction.");
    m.def("chyu_1990_fillet_over_straight",
          &cc::chyu_1990_fillet_over_straight, py::arg("arrangement"),
          "Chyu (1990) fillet/straight Nu ratio; no friction ratio.");
    m.def("chyu_1998_inline_over_staggered",
          &cc::chyu_1998_inline_over_staggered,
          "Chyu (1998) inline/staggered Total Nu ratio; no friction ratio.");
  }

  m.doc() = "Python bindings for combaero core";

  // Expose the CorrelationValidity enum
  py::enum_<solver::CorrelationValidity>(
      m, "CorrelationValidity", "Status flag for correlation evaluations.")
      .value("VALID", solver::CorrelationValidity::VALID)
      .value("EXTRAPOLATED", solver::CorrelationValidity::EXTRAPOLATED)
      .value("INVALID", solver::CorrelationValidity::INVALID)
      .export_values();

  // Expose the CorrelationResult struct for (f, J) tuples
  py::class_<solver::CorrelationResult<std::tuple<double, double>>>(
      m, "CorrelationResult",
      "Wrapper for correlation results, returning the (value, "
      "jacobian) tuple alongside validity status.")
      .def_readonly(
          "result",
          &solver::CorrelationResult<std::tuple<double, double>>::result)
      .def_readonly(
          "status",
          &solver::CorrelationResult<std::tuple<double, double>>::status)
      .def_readonly(
          "message",
          &solver::CorrelationResult<std::tuple<double, double>>::message)
      .def("__repr__",
           [](const solver::CorrelationResult<std::tuple<double, double>> &r) {
             std::string status_str = "VALID";
             if (r.status == solver::CorrelationValidity::EXTRAPOLATED)
               status_str = "EXTRAPOLATED";
             else if (r.status == solver::CorrelationValidity::INVALID)
               status_str = "INVALID";
             return "<CorrelationResult: (" +
                    std::to_string(std::get<0>(r.result)) + ", " +
                    std::to_string(std::get<1>(r.result)) +
                    "), status=" + status_str + ">";
           });

  py::class_<solver::Stream>(m, "MassStream", "Mass stream descriptor")
      .def(py::init<double, double, double, std::vector<double>>())
      .def_readwrite("m_dot", &solver::Stream::m_dot)
      .def_readwrite("Tt", &solver::Stream::Tt)
      .def_readwrite("Pt", &solver::Stream::Pt)
      .def_readwrite("Y", &solver::Stream::Y);

  py::class_<solver::StreamJacobian>(m, "StreamJacobian",
                                     "Derivatives w.r.t a stream")
      .def_readonly("d_mdot", &solver::StreamJacobian::d_mdot)
      .def_readonly("d_Tt", &solver::StreamJacobian::d_Tt)
      .def_readonly("d_Pt", &solver::StreamJacobian::d_Pt)
      .def_readonly("d_Y", &solver::StreamJacobian::d_Y);

  py::class_<MixtureState>(m, "MixtureState",
                           "High-level state for network components")
      .def(py::init([](double P, double Pt, double T, double Tt,
                       double m_dot, const std::vector<double> &Y) {
             MixtureState s;
             s.P = P;
             s.Pt = Pt;
             s.T = T;
             s.Tt = Tt;
             s.m_dot = m_dot;
             if (!Y.empty()) {
               s.Y = Y;
               double sum = 0.0;
               for (double y : Y)
                 sum += y;
               if (sum > 0.0) {
                 s.X = mass_to_mole(Y);
               } else {
                 s.X.assign(Y.size(), 0.0);
               }
             } else {
               std::size_t n = num_species();
               s.X.assign(n, 0.0);
               s.Y.assign(n, 0.0);
               if (n > 0) {
                 s.X[0] =
                     1.0; // Default to first species (usually N2 or similar)
                 s.Y[0] = 1.0;
               }
             }
             return s;
           }),
           py::arg("P") = 101325.0, py::arg("Pt") = 101325.0,
           py::arg("T") = 298.15, py::arg("Tt") = 298.15,
           py::arg("m_dot") = 0.0, py::arg("Y") = std::vector<double>{})
      .def_readwrite("P", &MixtureState::P)
      .def_readwrite("Pt", &MixtureState::Pt)
      .def_readwrite("T", &MixtureState::T)
      .def_readwrite("Tt", &MixtureState::Tt)
      .def_readwrite("m_dot", &MixtureState::m_dot)
      .def_property(
          "X", [](const MixtureState &s) { return s.X; },
          [](MixtureState &s, const std::vector<double> &v) {
            s.X = v;
            if (!v.empty())
              s.Y = mole_to_mass(v);
          })
      .def_property(
          "Y", [](const MixtureState &s) { return s.Y; },
          [](MixtureState &s, const std::vector<double> &v) {
            s.Y = v;
            if (!v.empty())
              s.X = mass_to_mole(v);
          })
      .def("density",
           [](const MixtureState &s) { return s.static_state().rho(); })
      .def("cp", [](const MixtureState &s) { return cp_mass(s.T, s.X); })
      .def("specific_heat",
           [](const MixtureState &s) { return cp_mass(s.T, s.X); })
      .def("h", [](const MixtureState &s) { return h_mass(s.T, s.X); })
      .def("enthalpy", [](const MixtureState &s) { return h_mass(s.T, s.X); })
      .def("mw", [](const MixtureState &s) { return s.static_state().mw(); })
      .def("molecular_weight",
           [](const MixtureState &s) { return s.static_state().mw(); })
      .def("static_state", &MixtureState::static_state)
      .def("total_state", &MixtureState::total_state);

  py::class_<solver::ChamberResult>(m, "ChamberResult",
                                    "Analytical node residuals")
      .def(py::init<>())
      .def_readonly("residuals", &solver::ChamberResult::residuals)
      .def_readonly("local_jacobian", &solver::ChamberResult::local_jacobian)
      .def_readonly("stream_jacobian", &solver::ChamberResult::stream_jacobian);

  py::class_<solver::MixerResult>(
      m, "MixerResult", "Output and internal gradients of stream mixes")
      .def_readonly("Tt_mix", &solver::MixerResult::Tt_mix)
      .def_readonly("Pt_mix", &solver::MixerResult::Pt_mix)
      .def_readonly("Y_mix", &solver::MixerResult::Y_mix)
      .def_readonly("dTt_mix_d_delta_h", &solver::MixerResult::dTt_mix_d_delta_h)
      .def_readonly("dTt_mix_dQ", &solver::MixerResult::dTt_mix_dQ)
      .def_readonly("dTt_mix_d_stream", &solver::MixerResult::dTt_mix_d_stream)
      .def_readonly("dPt_mix_d_stream",
                    &solver::MixerResult::dPt_mix_d_stream)
      .def_readonly("dY_mix_d_stream", &solver::MixerResult::dY_mix_d_stream);

  py::class_<solver::EntropyGenerationResult>(
      m, "EntropyGenerationResult",
      "Internal entropy generation of a control volume (#545)")
      .def_readonly("sigma", &solver::EntropyGenerationResult::sigma,
                    "Entropy generation [W/K]")
      .def_readonly("sigma_norm", &solver::EntropyGenerationResult::sigma_norm,
                    "sigma / (R m_in) [-]; -ln(Pt_out/Pt_in) for an adiabatic two-port");
  m.def("entropy_generation", &solver::entropy_generation, py::arg("inflows"),
        py::arg("outflows"), py::arg("S_heat") = 0.0,
        R"doc(Internal entropy generation of a control volume (#545).

sigma = sum_out m s - sum_in m s - S_heat, with each MassStream at its
stagnation state (Tt, Pt) and S_heat [W/K] the
entropy the heat brought in (heat_entropy_flux). Returns
EntropyGenerationResult(sigma [W/K], sigma_norm [-]); sigma_norm =
sigma / (R m_in), which for an adiabatic two-port of fixed composition is
-ln(Pt_out / Pt_in).)doc");
  m.def("heat_entropy_flux", &solver::heat_entropy_flux, py::arg("m_dot"),
        py::arg("Tt_in"), py::arg("Tt_out"), py::arg("Y"), py::arg("M_in") = 0.0,
        py::arg("M_out") = 0.0,
        R"doc(Entropy [W/K] heat added to a stream brings in (#545).

Integral of dQ / T as the stream's total temperature runs from Tt_in to
Tt_out [K], dQ = m_dot dh(Tt), at the static temperature T(Tt, M); the Mach
number runs linearly in Tt from M_in to M_out. Variable cp throughout.)doc");

  m.def("mixer_from_streams_and_jacobians",
        &solver::mixer_from_streams_and_jacobians, py::arg("streams"),
        py::arg("Q") = 0.0, py::arg("fraction") = 0.0);
  m.def("adiabatic_T_complete_and_jacobian_T_from_streams",
        &solver::adiabatic_T_complete_and_jacobian_T_from_streams,
        py::arg("streams"), py::arg("P"), py::arg("Q") = 0.0,
        py::arg("fraction") = 0.0);
  m.def("adiabatic_T_equilibrium_and_jacobians_from_streams",
        &solver::adiabatic_T_equilibrium_and_jacobians_from_streams,
        py::arg("streams"), py::arg("P"), py::arg("Q") = 0.0,
        py::arg("fraction") = 0.0);
  m.def(
      "combustor_residuals_and_jacobians",
      [](const std::vector<solver::Stream> &streams, double P, double Q,
         double fraction, py::function py_loss, bool use_equilibrium) {
        PressureLossCorrelation cpp_loss =
            [py_loss](const PressureLossContext &ctx) -> std::tuple<double, double> {
          py::object res = py_loss(ctx);
          try {
            return res.cast<std::tuple<double, double>>();
          } catch (...) {
            return {res.cast<double>(), 0.0};
          }
        };
        return solver::combustor_residuals_and_jacobians(
            streams, P, Q, fraction, cpp_loss, use_equilibrium);
      },
      py::arg("streams"), py::arg("P"), py::arg("Q") = 0.0,
      py::arg("fraction") = 0.0, py::arg("pressure_loss"),
      py::arg("use_equilibrium") = false);

  py::class_<solver::OrificeResult>(
      m, "OrificeResult", "Result of orifice residual and Jacobian evaluation")
      .def_readonly("m_dot_calc", &solver::OrificeResult::m_dot_calc)
      .def_readonly("d_mdot_dP_total_up",
                    &solver::OrificeResult::d_mdot_dP_total_up)
      .def_readonly("d_mdot_dP_static_down",
                    &solver::OrificeResult::d_mdot_dP_static_down)
      .def_readonly("d_mdot_dP_static_up",
                    &solver::OrificeResult::d_mdot_dP_static_up)
      .def_readonly("d_mdot_dP_rho_ref", &solver::OrificeResult::d_mdot_dP_rho_ref)
      .def_readonly("d_mdot_dTt_up", &solver::OrificeResult::d_mdot_dTt_up)
      .def_readonly("d_mdot_dY_up", &solver::OrificeResult::d_mdot_dY_up);

  py::class_<solver::ChannelResult>(
      m, "ChannelSolverResult",
      "Result of channel residual and Jacobian evaluation")
      .def(py::init<>())
      .def_readonly("dP_calc", &solver::ChannelResult::dP_calc)
      .def_readonly("d_dP_d_mdot", &solver::ChannelResult::d_dP_d_mdot)
      .def_readonly("d_dP_dP_static_up",
                    &solver::ChannelResult::d_dP_dP_static_up)
      .def_readonly("d_dP_dTt_up", &solver::ChannelResult::d_dP_dTt_up)
      .def_readonly("d_dP_dY_up", &solver::ChannelResult::d_dP_dY_up);

  py::class_<solver::MomentumChamberResult>(
      m, "MomentumChamberResult",
      "Result of momentum chamber residual and Jacobian evaluation")
      .def_readonly("residual", &solver::MomentumChamberResult::residual)
      .def_readonly("d_res_dP", &solver::MomentumChamberResult::d_res_dP)
      .def_readonly("d_res_dP_total",
                    &solver::MomentumChamberResult::d_res_dP_total)
      .def_readonly("d_res_dmdot", &solver::MomentumChamberResult::d_res_dmdot)
      .def_readonly("d_res_dT", &solver::MomentumChamberResult::d_res_dT)
      .def_readonly("d_res_dY", &solver::MomentumChamberResult::d_res_dY)
      .def_readonly("T_static", &solver::MomentumChamberResult::T_static)
      .def_readonly("dT_static_dP", &solver::MomentumChamberResult::dT_static_dP)
      .def_readonly("dT_static_dT", &solver::MomentumChamberResult::dT_static_dT)
      .def_readonly("dT_static_dmdot",
                    &solver::MomentumChamberResult::dT_static_dmdot)
      .def_readonly("M", &solver::MomentumChamberResult::M)
      .def_readonly("q_dyn", &solver::MomentumChamberResult::q_dyn)
      .def_readonly("rho_static", &solver::MomentumChamberResult::rho_static)
      .def_readonly("drho_dP", &solver::MomentumChamberResult::drho_dP)
      .def_readonly("drho_dT", &solver::MomentumChamberResult::drho_dT)
      .def_readonly("drho_dmdot", &solver::MomentumChamberResult::drho_dmdot)
      .def_readonly("drho_dY", &solver::MomentumChamberResult::drho_dY)
      .def_readonly("dT_static_dY", &solver::MomentumChamberResult::dT_static_dY)
      .def_readonly("bounded", &solver::MomentumChamberResult::bounded)
      .def_readonly("P_sonic", &solver::MomentumChamberResult::P_sonic)
      .def_readonly("Pt_sonic", &solver::MomentumChamberResult::Pt_sonic);

  m.def("orifice_residuals_and_jacobian",
        &solver::orifice_residuals_and_jacobian, py::arg("m_dot"),
        py::arg("P_total_up"), py::arg("P_static_up"), py::arg("Tt_up"),
        py::arg("Y_up"), py::arg("P_static_down"), py::arg("Cd"),
        py::arg("area"), py::arg("beta") = 0.0, py::arg("P_rho_ref") = -1.0,
        "Incompressible orifice with all derivatives. With E applied\n"
        "(0 < beta < 1) the drive is P_static_up - P_static_down, else\n"
        "P_total_up - P_static_down. Density at P_static_up, or at P_rho_ref\n"
        "when > 0 (its derivative is then d_mdot_dP_rho_ref).");

  m.def("channel_residuals_and_jacobian",
        &solver::channel_residuals_and_jacobian, py::arg("m_dot"),
        py::arg("P_total_up"), py::arg("P_static_up"), py::arg("Tt_up"),
        py::arg("Y_up"), py::arg("P_static_down"), py::arg("L"), py::arg("D"),
        py::arg("roughness"), py::arg("friction_model"),
        py::arg("f_multiplier") = 1.0, py::arg("exit_head_lost") = false);

  // Compressible flow elements
  m.def("orifice_compressible_mdot_and_jacobian",
        &solver::orifice_compressible_mdot_and_jacobian, py::arg("Tt"),
        py::arg("Pt"), py::arg("P_back"), py::arg("X"), py::arg("Cd"),
        py::arg("area"), py::arg("beta") = 0.0,
        "Compressible orifice flow, analytic Jacobians exact on both sides of\n"
        "choke (the flow is C1 through P*, #772).\n\n"
        "Returns: (mdot, d_mdot_dPt, d_mdot_dP_back, d_mdot_dTt)");

  m.def("orifice_compressible_residuals_and_jacobian",
        &solver::orifice_compressible_residuals_and_jacobian, py::arg("m_dot"),
        py::arg("P_total_up"), py::arg("P_static_up"), py::arg("Tt_up"),
        py::arg("Y_up"), py::arg("P_static_down"), py::arg("Cd"), py::arg("area"),
        py::arg("beta") = 0.0,
        "Compressible orifice for network solver with all derivatives. With E\n"
        "applied (0 < beta < 1) the nozzle expands from P_static_up, else from\n"
        "P_total_up.");

  py::class_<solver::OrificeModelResult, solver::OrificeResult>(
      m, "OrificeModelResult",
      "OrificeResult of an opt-in compressible orifice model (#854), plus\n"
      "d_mdot_dCd and the model's state.")
      .def_readonly("d_mdot_dCd", &solver::OrificeModelResult::d_mdot_dCd)
      .def_readonly("Ci", &solver::OrificeModelResult::Ci)
      .def_readonly("contraction", &solver::OrificeModelResult::contraction)
      .def_readonly("expansibility", &solver::OrificeModelResult::expansibility)
      .def_readonly("choked", &solver::OrificeModelResult::choked)
      .def_readonly("in_range", &solver::OrificeModelResult::in_range);

  m.def("orifice_jobson_residuals_and_jacobian",
        &solver::orifice_jobson_residuals_and_jacobian, py::arg("m_dot"),
        py::arg("P_total_up"), py::arg("P_static_up"), py::arg("Tt_up"),
        py::arg("Y_up"), py::arg("P_static_down"), py::arg("Cd"), py::arg("area"),
        py::arg("beta") = 0.0, py::arg("Cv") = 1.0,
        "Jobson (1955) compressible orifice with the approach-velocity\n"
        "extension (#854), sub- and supercritical, all derivatives analytic\n"
        "(composition by a forward difference). Cd is the incompressible\n"
        "coefficient in ISO's form; Ci from jobson_ci_from_cd(Cd, beta, Cv).");

  m.def("orifice_iso5167_residuals_and_jacobian",
        &solver::orifice_iso5167_residuals_and_jacobian, py::arg("m_dot"),
        py::arg("P_total_up"), py::arg("P_static_up"), py::arg("Tt_up"),
        py::arg("Y_up"), py::arg("P_static_down"), py::arg("Cd"), py::arg("area"),
        py::arg("beta") = 0.0,
        "ISO 5167-2:2003 Eq. (1) with Eq. (5)'s expansibility (#854).\n"
        "in_range is False below p2/p1 = 0.75, where Eq. (5) is extrapolated.");

  py::class_<solver::SonicStatic>(
      m, "SonicStatic", "Sonic static state of a chamber's mass flux (#637)")
      .def_readonly("P", &solver::SonicStatic::P)
      .def_readonly("T", &solver::SonicStatic::T)
      .def_readonly("dP_dTt", &solver::SonicStatic::dP_dTt)
      .def_readonly("dP_dG", &solver::SonicStatic::dP_dG)
      .def_readonly("dP_dY", &solver::SonicStatic::dP_dY);
  m.def("sonic_static_pressure", &solver::sonic_static_pressure, py::arg("Tt"), py::arg("G"),
        py::arg("Y"),
        "Static pressure P* at which mass flux G is sonic for stagnation Tt\n"
        "(variable cp): the floor of a chamber discharging into a reservoir\n"
        "(#637). P* = G sqrt(R T* / gamma*).");
  py::class_<solver::ApproachStatic>(
      m, "ApproachStatic", "Static pressure of an orifice's approach pipe (#651)")
      .def_readonly("P", &solver::ApproachStatic::P)
      .def_readonly("dP_dPt", &solver::ApproachStatic::dP_dPt)
      .def_readonly("dP_dTt", &solver::ApproachStatic::dP_dTt)
      .def_readonly("dP_dG", &solver::ApproachStatic::dP_dG)
      .def_readonly("dP_dY", &solver::ApproachStatic::dP_dY)
      .def_readonly("choked", &solver::ApproachStatic::choked);
  m.def("approach_static_pressure", &solver::approach_static_pressure,
        py::arg("Pt"), py::arg("Tt"), py::arg("G"), py::arg("Y"),
        py::arg("incompressible"),
        "Static pressure of the flow at mass flux G from stagnation (Pt, Tt):\n"
        "duct_static_state compressible, Pt = P + G^2 R Tt / (2 P)\n"
        "incompressible. The approach state of an orifice with E (#651).");
  py::class_<solver::OrificeTapPressure>(
      m, "OrificeTapPressure",
      "Downstream tap pressure of an orifice plate between pipes (#645)")
      .def_readonly("p_tap", &solver::OrificeTapPressure::p_tap)
      .def_readonly("loss_ratio", &solver::OrificeTapPressure::loss_ratio)
      .def_readonly("d_dP_static_up", &solver::OrificeTapPressure::d_dP_static_up)
      .def_readonly("d_dP_total_up", &solver::OrificeTapPressure::d_dP_total_up)
      .def_readonly("d_dP_total_down", &solver::OrificeTapPressure::d_dP_total_down)
      .def_readonly("d_dC", &solver::OrificeTapPressure::d_dC)
      .def_readonly("floored", &solver::OrificeTapPressure::floored)
      .def_readonly("in_range", &solver::OrificeTapPressure::in_range);
  m.def("orifice_tap_pressure_and_jacobian", &solver::orifice_tap_pressure_and_jacobian,
        py::arg("P_static_up"), py::arg("Pt_up"), py::arg("Pt_down"),
        py::arg("beta"), py::arg("C"), py::arg("compressible"),
        "p_tap = P_static_up - (Pt_up - Pt_down) / k(beta, C): the\n"
        "tap pressure that a recovered drop between two pipes implies, with\n"
        "its derivatives (d_dP_total_up/down are with respect to Pt_up/down).\n"
        "k is ISO 5167-2's permanent-loss ratio.\n"
        "Compressible, p_tap is floored below any critical ratio.");


  m.def("momentum_chamber_residual_and_jacobian",
        &solver::momentum_chamber_residual_and_jacobian, py::arg("P"),
        py::arg("Pt"), py::arg("m_dot"), py::arg("Tt"), py::arg("Y"),
        py::arg("area"), py::arg("incompressible") = false,
        "Momentum chamber residual: Pt = Pt(P, T); Tt is the chamber "
        "stagnation temperature, the static T is solved inside (#576). "
        "The d_res_dT, dT_static_dT and drho_dT fields are with respect to "
        "Tt, d_res_dP_total with respect to Pt. incompressible: the low-Mach "
        "limit, T_static = Tt and Pt = P + (m/A)^2 / (2 rho) (#588).");

  // --- Momentum-CV junction: whole-element (f, J) -------------------------
  py::class_<solver::MpceGeometry>(
      m, "MpceGeometry",
      "Fixed geometry and tuning for the momentum-CV junction kernel")
      .def(py::init<>())
      .def_readwrite("area", &solver::MpceGeometry::area)
      .def_readwrite("theta_rad", &solver::MpceGeometry::theta_rad)
      .def_readwrite("port_sign", &solver::MpceGeometry::port_sign)
      .def_readwrite("joining_etransfer_alpha",
                     &solver::MpceGeometry::joining_etransfer_alpha)
      .def_readwrite("eta_scale", &solver::MpceGeometry::eta_scale)
      .def_readwrite("gamma", &solver::MpceGeometry::gamma);

  py::class_<solver::MpceResidualJacobian>(
      m, "MpceResidualJacobian",
      "Whole-element residuals and Jacobian for the momentum-CV junction")
      .def_readonly("residual", &solver::MpceResidualJacobian::residual)
      .def_readonly("jacobian", &solver::MpceResidualJacobian::jacobian)
      .def_readonly("valid", &solver::MpceResidualJacobian::valid)
      .def_readonly("common_port", &solver::MpceResidualJacobian::common_port)
      .def_readonly("k_term_sign", &solver::MpceResidualJacobian::k_term_sign)
      .def_readonly("k_per_port", &solver::MpceResidualJacobian::k_per_port);

  m.def("mpce_residuals_and_jacobian",
        &solver::mpce_residuals_and_jacobian, py::arg("p_static"),
        py::arg("p_total"), py::arg("rho"), py::arg("drho_dp"),
        py::arg("outer_mdot"), py::arg("pt_jct"), py::arg("geom"),
        "Momentum-CV junction whole-element (f, J), three ports.\n\n"
        "Physics only: the caller owns the degenerate-state guards and the\n"
        "wrong-direction soft barrier (see include/mpce_junction.h).\n\n"
        "Jacobian seed order, matching the returned rows:\n"
        "  0..2  p_static[i]     3..5  p_total[i]\n"
        "  6..8  outer_mdot[i]   9     pt_jct\n\n"
        "``outer_mdot`` is the CONNECTING element's mass flow, not the\n"
        "junction-convention one; ``geom.port_sign`` maps between them, so\n"
        "the Jacobian comes back already in the solver's own unknowns.\n\n"
        "Returns: MpceResidualJacobian (valid=False when the flow pattern\n"
        "is not a junction at all).");

  // --- Ideal momentum junction: whole-element (f, J), #517 ----------------
  py::class_<solver::IdealJunctionGeometry>(
      m, "IdealJunctionGeometry",
      "Port geometry for the ideal momentum junction kernel")
      .def(py::init<>())
      .def_readwrite("area", &solver::IdealJunctionGeometry::area)
      .def_readwrite("axis_rad", &solver::IdealJunctionGeometry::axis_rad)
      .def_readwrite("port_sign", &solver::IdealJunctionGeometry::port_sign)
      .def_readwrite("loss_factor", &solver::IdealJunctionGeometry::loss_factor,
                     "Receiver loss factor f_i >= 0 per port (#526); empty = all 1, "
                     "the parameter-free closure.");

  py::class_<solver::IdealJunctionResult>(
      m, "IdealJunctionResult",
      "Whole-element residuals and Jacobian for the ideal momentum junction")
      .def_readonly("residual", &solver::IdealJunctionResult::residual)
      .def_readonly("jacobian", &solver::IdealJunctionResult::jacobian)
      .def_readonly("velocity", &solver::IdealJunctionResult::velocity)
      .def_readonly("V_x", &solver::IdealJunctionResult::V_x)
      .def_readonly("V_y", &solver::IdealJunctionResult::V_y);

  m.def("ideal_junction_residuals_and_jacobian",
        &solver::ideal_junction_residuals_and_jacobian, py::arg("p_static"),
        py::arg("rho"), py::arg("drho_dp"), py::arg("outer_mdot"),
        py::arg("p_jn"), py::arg("geom"),
        "Parameter-free momentum junction (#517), whole-element (f, J), 2-6\n"
        "ports. Suppliers sit at the junction static p_jn; each receiver j\n"
        "satisfies P_j - p_jn + (m_j/A_j)(u_j - V . n_j) = 0, V the supplied\n"
        "mass-weighted velocity vector. Merge = Blaisdell & Manson theory;\n"
        "split = vector Borda-Carnot loss (Idelchik Eq 7-6, A' = 1).\n"
        "See include/ideal_junction.h.\n\n"
        "Jacobian columns: p_static[0..N-1], outer_mdot[0..N-1], p_jn.");

  // Ejector (supersonic ejector on the MultiPortChamberBase topology,
  // spanning the critical / subcritical / unchoked jet-pump regimes). Physics
  // + analytic Jacobians in include/ejector.h. The critical-mode scalar
  // closures expose Jacobians w.r.t. the 4 thermodynamic inputs (p_g, t_g,
  // p_e, t_e); the operating-regime closures add the P_py / omega / m_dot
  // partials, and ejector_element_residuals_and_jacobian returns the whole
  // 4-row element system w.r.t. all nine unknowns. gamma, r_gas, geometry and
  // recovery_efficiency are frozen coefficients here. Area ratios are passed
  // as plain doubles (a lambda builds the internal EjectorGeometry), keeping
  // the Python call flat like the other solver bindings.
  py::class_<solver::EjectorEntrainmentResult>(
      m, "EjectorEntrainmentResult",
      "Value outputs of ejector_entrainment_ratio (Huang 1999 Eqs. 1-8)")
      .def_readonly("omega", &solver::EjectorEntrainmentResult::omega)
      .def_readonly("mach_nozzle_exit",
                    &solver::EjectorEntrainmentResult::mach_nozzle_exit)
      .def_readonly("mach_hypothetical_throat",
                    &solver::EjectorEntrainmentResult::mach_hypothetical_throat)
      .def_readonly("area_ratio_primary_core",
                    &solver::EjectorEntrainmentResult::area_ratio_primary_core)
      .def_readonly("area_ratio_entrained",
                    &solver::EjectorEntrainmentResult::area_ratio_entrained)
      .def_readonly("p_mixing_pa",
                    &solver::EjectorEntrainmentResult::p_mixing_pa);

  py::class_<solver::EjectorEntrainmentJacobian>(
      m, "EjectorEntrainmentJacobian",
      "Entrainment ratio value plus d(omega)/d(p_g, t_g, p_e, t_e)")
      .def_readonly("value", &solver::EjectorEntrainmentJacobian::value)
      .def_readonly("domega_dp_g",
                    &solver::EjectorEntrainmentJacobian::domega_dp_g)
      .def_readonly("domega_dt_g",
                    &solver::EjectorEntrainmentJacobian::domega_dt_g)
      .def_readonly("domega_dp_e",
                    &solver::EjectorEntrainmentJacobian::domega_dp_e)
      .def_readonly("domega_dt_e",
                    &solver::EjectorEntrainmentJacobian::domega_dt_e);

  py::class_<solver::EjectorCriticalPressureResult>(
      m, "EjectorCriticalPressureResult",
      "Value outputs of ejector_critical_back_pressure (Kracik & Dvorak "
      "2016 Eqs. 7-13)")
      .def_readonly("p_c_pa", &solver::EjectorCriticalPressureResult::p_c_pa)
      .def_readonly(
          "p_mixed_stagnation_pa",
          &solver::EjectorCriticalPressureResult::p_mixed_stagnation_pa)
      .def_readonly(
          "temp_mixed_stagnation_k",
          &solver::EjectorCriticalPressureResult::temp_mixed_stagnation_k)
      .def_readonly("lambda_mixed",
                    &solver::EjectorCriticalPressureResult::lambda_mixed)
      .def_readonly("mach_mixed",
                    &solver::EjectorCriticalPressureResult::mach_mixed);

  py::class_<solver::EjectorCriticalPressureJacobian>(
      m, "EjectorCriticalPressureJacobian",
      "Critical back pressure value plus d(P_c*)/d(p_g, t_g, p_e, t_e)")
      .def_readonly("value", &solver::EjectorCriticalPressureJacobian::value)
      .def_readonly("dpc_dp_g",
                    &solver::EjectorCriticalPressureJacobian::dpc_dp_g)
      .def_readonly("dpc_dt_g",
                    &solver::EjectorCriticalPressureJacobian::dpc_dt_g)
      .def_readonly("dpc_dp_e",
                    &solver::EjectorCriticalPressureJacobian::dpc_dp_e)
      .def_readonly("dpc_dt_e",
                    &solver::EjectorCriticalPressureJacobian::dpc_dt_e);

  py::class_<solver::EjectorCDNozzleJacobian>(
      m, "EjectorCDNozzleJacobian",
      "C-D nozzle mass flow plus d(mdot)/d(p0, t0, P_py)")
      .def_readonly("mdot", &solver::EjectorCDNozzleJacobian::mdot)
      .def_readonly("dmdot_dp0", &solver::EjectorCDNozzleJacobian::dmdot_dp0)
      .def_readonly("dmdot_dt0", &solver::EjectorCDNozzleJacobian::dmdot_dt0)
      .def_readonly("dmdot_dp_py", &solver::EjectorCDNozzleJacobian::dmdot_dp_py);

  py::class_<solver::EjectorCDExitStaticJacobian>(
      m, "EjectorCDExitStaticJacobian",
      "C-D nozzle exit static P_py plus d(P_py)/d(m_dot, p0, t0)")
      .def_readonly("p_py", &solver::EjectorCDExitStaticJacobian::p_py)
      .def_readonly("dp_py_dm_dot", &solver::EjectorCDExitStaticJacobian::dp_py_dm_dot)
      .def_readonly("dp_py_dp0", &solver::EjectorCDExitStaticJacobian::dp_py_dp0)
      .def_readonly("dp_py_dt0", &solver::EjectorCDExitStaticJacobian::dp_py_dt0)
      .def_readonly("out_of_range", &solver::EjectorCDExitStaticJacobian::out_of_range,
                    "+1: m_dot above the flow at P* (P_py clamped to P*); -1: below the\n"
                    "flow at Pt (1 - 1e-12) (P_py clamped there); 0: in range. A clamped\n"
                    "result's partials are the clamped function's (0, k, 0).");

  py::class_<solver::EjectorJetPumpDischargeJacobian>(
      m, "EjectorJetPumpDischargeJacobian",
      "Jet-pump discharge stagnation plus d(p03)/d(p_g, t_g, p_e, t_e, P_py, omega)")
      .def_readonly("p03", &solver::EjectorJetPumpDischargeJacobian::p03)
      .def_readonly("dp03_dp_g", &solver::EjectorJetPumpDischargeJacobian::dp03_dp_g)
      .def_readonly("dp03_dt_g", &solver::EjectorJetPumpDischargeJacobian::dp03_dt_g)
      .def_readonly("dp03_dp_e", &solver::EjectorJetPumpDischargeJacobian::dp03_dp_e)
      .def_readonly("dp03_dt_e", &solver::EjectorJetPumpDischargeJacobian::dp03_dt_e)
      .def_readonly("dp03_dp_py", &solver::EjectorJetPumpDischargeJacobian::dp03_dp_py)
      .def_readonly("dp03_domega", &solver::EjectorJetPumpDischargeJacobian::dp03_domega);

  py::class_<solver::EjectorNozzleExitMach>(
      m, "EjectorNozzleExitMach",
      "Supersonic nozzle-exit Mach M_p1 from A/A* plus its out-of-range flag")
      .def_readonly("mach", &solver::EjectorNozzleExitMach::mach)
      .def_readonly("out_of_range", &solver::EjectorNozzleExitMach::out_of_range);

  m.def("ejector_nozzle_exit_mach", &solver::ejector_nozzle_exit_mach, py::arg("area_ratio"),
        py::arg("gamma"),
        "Supersonic root M of A/A*(M) = area_ratio on the bisection bracket\n"
        "[1 + 1e-12, 50]. An area ratio outside A/A*(bracket) returns the\n"
        "bracket end (clamp) with out_of_range = -1 (below) or +1 (above).\n\n"
        "Returns: EjectorNozzleExitMach (mach, out_of_range)");

  m.def("ejector_choked_mass_flow_and_jacobian",
        &solver::ejector_choked_mass_flow_and_jacobian, py::arg("Pt"),
        py::arg("Tt"), py::arg("area_throat"), py::arg("gamma"),
        py::arg("r_gas"), py::arg("eta"),
        "Choked (sonic-throat) mass flow (Huang 1999 Eqs. 1, 7) from the\n"
        "primary stagnation state (Pt, Tt).\n\n"
        "Returns: (mdot, d_mdot_dPt, d_mdot_dTt)");

  m.def("ejector_cd_nozzle_mass_flow_and_jacobian",
        &solver::ejector_cd_nozzle_mass_flow_and_jacobian, py::arg("Pt"),
        py::arg("Tt"), py::arg("p_static_down"), py::arg("area_throat"),
        py::arg("area_exit"), py::arg("gamma"), py::arg("r_gas"), py::arg("eta"),
        py::arg("eps_frac"),
        "C-D nozzle mass flow (choked/unchoked) from the primary stagnation\n"
        "state (Pt, Tt) with analytic Jacobian.\n\n"
        "Returns: EjectorCDNozzleJacobian (mdot + dmdot_d{p0,t0,p_py}, the\n"
        "p0/t0 suffixes being d/dPt and d/dTt)");

  m.def("ejector_cd_nozzle_exit_static_and_jacobian",
        &solver::ejector_cd_nozzle_exit_static_and_jacobian, py::arg("m_dot"),
        py::arg("Pt"), py::arg("Tt"), py::arg("area_throat"), py::arg("area_exit"),
        py::arg("gamma"), py::arg("r_gas"), py::arg("eta"), py::arg("eps_frac"),
        "C-D nozzle exit static P_py that passes m_dot from the primary\n"
        "stagnation state (Pt, Tt), with implicit Jacobian.\n\n"
        "Returns: EjectorCDExitStaticJacobian (p_py + dp_py_d{m_dot,p0,t0}, the\n"
        "p0/t0 suffixes being d/dPt and d/dTt)");

  m.def("ejector_jetpump_discharge_and_jacobian",
        &solver::ejector_jetpump_discharge_and_jacobian, py::arg("p_g"),
        py::arg("t_g"), py::arg("p_e"), py::arg("t_e"), py::arg("p_py"),
        py::arg("omega"), py::arg("gamma"), py::arg("recovery_efficiency"),
        py::arg("eta_primary") = 1.0, py::arg("eta_secondary") = 1.0,
        "Jet-pump mixed-flow discharge stagnation with analytic Jacobian.\n"
        "eta_primary / eta_secondary: the nozzles' kinetic-energy efficiencies;\n"
        "each stream mixes at sqrt(eta) of its isentropic velocity (#588).\n\n"
        "Returns: EjectorJetPumpDischargeJacobian (p03 + dp03_d{p_g,t_g,p_e,t_e,p_py,omega})");

  m.def(
      "ejector_entrainment_ratio_and_jacobian",
      [](double p_g, double t_g, double p_e, double t_e,
         double area_ratio_nozzle, double area_ratio_mix, double gamma) {
        return solver::ejector_entrainment_ratio_and_jacobian(
            p_g, t_g, p_e, t_e,
            solver::EjectorGeometry{area_ratio_nozzle, area_ratio_mix}, gamma);
      },
      py::arg("p_g"), py::arg("t_g"), py::arg("p_e"), py::arg("t_e"),
      py::arg("area_ratio_nozzle"), py::arg("area_ratio_mix"), py::arg("gamma"),
      "Critical-mode entrainment ratio omega with analytic Jacobian.\n\n"
      "Returns: EjectorEntrainmentJacobian (value + domega_d{p_g,t_g,p_e,t_e})");

  m.def(
      "ejector_critical_back_pressure_and_jacobian",
      [](double p_g, double t_g, double p_e, double t_e,
         double area_ratio_nozzle, double area_ratio_mix, double gamma,
         double r_gas, double recovery_efficiency) {
        return solver::ejector_critical_back_pressure_and_jacobian(
            p_g, t_g, p_e, t_e,
            solver::EjectorGeometry{area_ratio_nozzle, area_ratio_mix}, gamma,
            r_gas, recovery_efficiency);
      },
      py::arg("p_g"), py::arg("t_g"), py::arg("p_e"), py::arg("t_e"),
      py::arg("area_ratio_nozzle"), py::arg("area_ratio_mix"), py::arg("gamma"),
      py::arg("r_gas"), py::arg("recovery_efficiency"),
      "Critical back pressure P_c* with analytic Jacobian.\n\n"
      "Returns: EjectorCriticalPressureJacobian (value + dpc_d{p_g,t_g,p_e,t_e})");

  m.def(
      "ejector_dead_head_back_pressure_and_jacobian",
      [](double p_g, double t_g, double p_e, double t_e,
         double area_ratio_nozzle, double area_ratio_mix, double gamma,
         double r_gas, double recovery_efficiency) {
        return solver::ejector_dead_head_back_pressure_and_jacobian(
            p_g, t_g, p_e, t_e,
            solver::EjectorGeometry{area_ratio_nozzle, area_ratio_mix}, gamma,
            r_gas, recovery_efficiency);
      },
      py::arg("p_g"), py::arg("t_g"), py::arg("p_e"), py::arg("t_e"),
      py::arg("area_ratio_nozzle"), py::arg("area_ratio_mix"), py::arg("gamma"),
      py::arg("r_gas"), py::arg("recovery_efficiency"),
      "Dead-head back pressure P_b0 (the mixing closure at zero entrainment)\n"
      "with analytic Jacobian; value.p_c_pa carries P_b0.\n\n"
      "Returns: EjectorCriticalPressureJacobian (value + dpc_d{p_g,t_g,p_e,t_e})");

  m.def("ejector_huang_mixing_coefficient", &solver::ejector_huang_mixing_coefficient,
        py::arg("area_ratio_mix"),
        "phi_m of Huang 1999 Eq. 19: the mixing-loss coefficient, stepped in A_3/A_t.");

  m.def(
      "ejector_huang_critical_back_pressure_and_jacobian",
      [](double p_g, double t_g, double p_e, double t_e,
         double area_ratio_nozzle, double area_ratio_mix, double gamma,
         double r_gas, double recovery_efficiency) {
        return solver::ejector_huang_critical_back_pressure_and_jacobian(
            p_g, t_g, p_e, t_e,
            solver::EjectorGeometry{area_ratio_nozzle, area_ratio_mix}, gamma,
            r_gas, recovery_efficiency);
      },
      py::arg("p_g"), py::arg("t_g"), py::arg("p_e"), py::arg("t_e"),
      py::arg("area_ratio_nozzle"), py::arg("area_ratio_mix"), py::arg("gamma"),
      py::arg("r_gas"), py::arg("recovery_efficiency"),
      "Critical back pressure P_c* by Huang 1999's printed chain (phi_m mixing,\n"
      "normal shock, isentropic diffuser) with analytic Jacobian.\n\n"
      "Returns: EjectorCriticalPressureJacobian (value + dpc_d{p_g,t_g,p_e,t_e})");

  m.def(
      "ejector_sudden_expansion_dead_head_and_jacobian",
      [](double p_g, double t_g, double p_e, double t_e,
         double area_ratio_nozzle, double area_ratio_mix, double gamma,
         double r_gas, double recovery_efficiency) {
        return solver::ejector_sudden_expansion_dead_head_and_jacobian(
            p_g, t_g, p_e, t_e,
            solver::EjectorGeometry{area_ratio_nozzle, area_ratio_mix}, gamma,
            r_gas, recovery_efficiency);
      },
      py::arg("p_g"), py::arg("t_g"), py::arg("p_e"), py::arg("t_e"),
      py::arg("area_ratio_nozzle"), py::arg("area_ratio_mix"), py::arg("gamma"),
      py::arg("r_gas"), py::arg("recovery_efficiency"),
      "Dead-head back pressure P_b0 as a sudden expansion over A_3 (zero\n"
      "secondary flow) with analytic Jacobian; value.p_c_pa carries P_b0.\n\n"
      "Returns: EjectorCriticalPressureJacobian (value + dpc_d{p_g,t_g,p_e,t_e})");

  py::class_<solver::EjectorElementResidualJacobian>(
      m, "EjectorElementResidualJacobian",
      "4-row ejector element residual system + Jacobian. `residuals` is [R0, "
      "R1, R2, R3]; `jacobian` is jacobian[row][seed] over the 9 unknowns "
      "(mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out, p_py).")
      .def_readonly("residuals", &solver::EjectorElementResidualJacobian::residuals)
      .def_readonly("jacobian", &solver::EjectorElementResidualJacobian::jacobian)
      .def_readonly("nozzle_exit_mach_out_of_range",
                    &solver::EjectorElementResidualJacobian::nozzle_exit_mach_out_of_range,
                    "Nonzero when the nozzle-exit Mach M_p1 is a clamped bracket end:\n"
                    "+1 area_ratio_nozzle above A/A*(M = 50), -1 below A/A*(M = 1).");

  m.def(
      "ejector_element_residuals_and_jacobian",
      [](double mp, double ms, double mdot_out, double p_g, double t_g,
         double p_e, double t_e, double p_out, double p_py,
         double area_ratio_nozzle, double area_ratio_mix, double area_throat,
         double area_nozzle_exit, double area_secondary, double gamma,
         double r_gas, double eta_primary, double eta_secondary,
         double recovery_efficiency, double eps_frac, double s_choke_lo,
         double s_choke_hi, double s_sub_lo, double s_sub_hi,
         const std::string& anchor_closure) {
        solver::EjectorAnchorClosure closure;
        if (anchor_closure == "kracik_dvorak") {
          closure = solver::EjectorAnchorClosure::kracik_dvorak;
        } else if (anchor_closure == "huang1999") {
          closure = solver::EjectorAnchorClosure::huang1999;
        } else {
          throw std::invalid_argument("ejector_element_residuals_and_jacobian: anchor_closure "
                                      "must be 'kracik_dvorak' or 'huang1999', got '" +
                                      anchor_closure + "'");
        }
        return solver::ejector_element_residuals_and_jacobian(
            mp, ms, mdot_out, p_g, t_g, p_e, t_e, p_out, p_py,
            solver::EjectorGeometry{area_ratio_nozzle, area_ratio_mix},
            area_throat, area_nozzle_exit, area_secondary, gamma, r_gas,
            eta_primary, eta_secondary, recovery_efficiency, eps_frac,
            s_choke_lo, s_choke_hi, s_sub_lo, s_sub_hi, closure);
      },
      py::arg("mp"), py::arg("ms"), py::arg("mdot_out"), py::arg("p_g"),
      py::arg("t_g"), py::arg("p_e"), py::arg("t_e"), py::arg("p_out"),
      py::arg("p_py"), py::arg("area_ratio_nozzle"), py::arg("area_ratio_mix"),
      py::arg("area_throat"), py::arg("area_nozzle_exit"),
      py::arg("area_secondary"), py::arg("gamma"), py::arg("r_gas"),
      py::arg("eta_primary"), py::arg("eta_secondary"),
      py::arg("recovery_efficiency"), py::arg("eps_frac"), py::arg("s_choke_lo"),
      py::arg("s_choke_hi"), py::arg("s_sub_lo"), py::arg("s_sub_hi"),
      py::arg("anchor_closure") = "huang1999",
      "Whole-element 4-row ejector (f, J) across all three operating regimes.\n"
      "anchor_closure: 'huang1999' (default) or 'kracik_dvorak', the (P_c*, P_b0)\n"
      "pair of the plateau end and the droop's far end.\n\n"
      "Returns: EjectorElementResidualJacobian (residuals[4] + jacobian[4][9])");

  py::class_<solver::EjectorJetPumpIncompressibleResult>(
      m, "EjectorJetPumpIncompressibleResult",
      "Incompressible jet-pump residuals (#588): `residuals` [R0..R3], "
      "`jacobian`[row][seed] over (mp, ms, mo, Pt_p, v_p, Pt_s, v_s, Pt_o, p3), "
      "and the mixing-duct exit p4, q4.")
      .def_readonly("residuals", &solver::EjectorJetPumpIncompressibleResult::residuals)
      .def_readonly("jacobian", &solver::EjectorJetPumpIncompressibleResult::jacobian)
      .def_readonly("p4", &solver::EjectorJetPumpIncompressibleResult::p4)
      .def_readonly("q4", &solver::EjectorJetPumpIncompressibleResult::q4);

  m.def("ejector_jetpump_incompressible_residuals_and_jacobian",
        &solver::ejector_jetpump_incompressible_residuals_and_jacobian, py::arg("mp"),
        py::arg("ms"), py::arg("mo"), py::arg("pt_p"), py::arg("v_p"), py::arg("pt_s"),
        py::arg("v_s"), py::arg("pt_o"), py::arg("p3"), py::arg("area_nozzle_exit"),
        py::arg("area_secondary"), py::arg("area_mixing"), py::arg("k_nozzle") = 0.0,
        py::arg("k_suction") = 0.0, py::arg("k_throat") = 0.0, py::arg("k_diffuser") = 0.0,
        py::arg("recovery") = 1.0,
        "Incompressible jet pump (#588): Sanger's 1-D theory (NASA TN D-4445,\n"
        "App. B) with each stream at its own constant density, the low-Mach\n"
        "limit of the ejector. Same four rows and owned unknown p3 as the\n"
        "compressible element; v = 1/rho.");

  // Area change (#731): Idelchik closures and the compressible element kernel.
  m.def("sudden_expansion_zeta",
        [](double Re, double area_ratio) {
          const auto z = combaero::sudden_expansion_zeta(Re, area_ratio);
          return std::make_tuple(z.zeta, z.dzeta_dRe);
        },
        py::arg("Re"), py::arg("area_ratio"),
        "Sudden expansion loss coefficient, Idelchik (1966) diagram 4-1, on the\n"
        "small section's dynamic head: 26/Re below Re 8, a C1 interpolation of\n"
        "the printed table to Re 3500, Borda-Carnot (1 - n)^2 above.\n"
        "Returns (zeta, dzeta_dRe).");

  m.def("sudden_contraction_zeta",
        [](double Re, double area_ratio) {
          const auto z = combaero::sudden_contraction_zeta(Re, area_ratio);
          return std::make_tuple(z.zeta, z.dzeta_dRe);
        },
        py::arg("Re"), py::arg("area_ratio"),
        "Sudden contraction loss coefficient, Idelchik (1966) diagrams 3-10 and\n"
        "3-9, on the small section's dynamic head: 27/Re below Re 8, a C1\n"
        "interpolation of Karev's table to Re 1e4, 0.5 (1 - n) above.\n"
        "Returns (zeta, dzeta_dRe).");

  m.def("conical_diffuser_zeta", &combaero::conical_diffuser_zeta,
        py::arg("full_angle_deg"), py::arg("area_ratio"),
        "Conical diffuser expansion loss zeta_exp, Idelchik (1966) diagram 5-2\n"
        "(bilinear through the printed table), on the small section's dynamic\n"
        "head. full_angle_deg is the included angle; friction not included.");

  m.def("conical_nozzle_zeta", &combaero::conical_nozzle_zeta,
        py::arg("full_angle_deg"), py::arg("length_over_D"), py::arg("area_ratio"),
        "Conical contraction loss coefficient zeta' (1 - n), Idelchik (1966)\n"
        "diagram 3-7, zeta'(full angle, l/D) from diagram 3-6 (bilinear through\n"
        "the printed table); l is the cone's axial length, D the small diameter.\n"
        "Re > 1e4; friction not included.");

  py::class_<combaero::BordaCarnotStep>(m, "BordaCarnotStep",
                                        "Compressible Borda-Carnot step")
      .def_readonly("M2", &combaero::BordaCarnotStep::M2)
      .def_readonly("pt_ratio", &combaero::BordaCarnotStep::pt_ratio)
      .def_readonly("dpt_dM1", &combaero::BordaCarnotStep::dpt_dM1)
      .def_readonly("dpt_dgamma", &combaero::BordaCarnotStep::dpt_dgamma);

  m.def("compressible_borda_carnot", &combaero::compressible_borda_carnot,
        py::arg("M1"), py::arg("area_ratio"), py::arg("gamma"),
        "Sudden expansion from M1 by area ratio n (small/large), calorically\n"
        "perfect gas: mass, momentum (base pressure on the whole large face) and\n"
        "total enthalpy conserved. Returns BordaCarnotStep (M2, Pt2/Pt1 and its\n"
        "partials in M1 and gamma).");

  py::class_<combaero::AreaChangeLoss>(m, "AreaChangeLoss",
                                       "Area change loss at a given flow")
      .def_readonly("dPt", &combaero::AreaChangeLoss::dPt)
      .def_readonly("zeta", &combaero::AreaChangeLoss::zeta)
      .def_readonly("mach", &combaero::AreaChangeLoss::mach)
      .def_readonly("Re", &combaero::AreaChangeLoss::Re)
      .def_readonly("expansion", &combaero::AreaChangeLoss::expansion)
      .def_readonly("choked", &combaero::AreaChangeLoss::choked);

  m.def("area_change_loss", &combaero::area_change_loss,
        py::arg("m_dot_abs"), py::arg("Pt0"), py::arg("Tt0"), py::arg("Y"),
        py::arg("F_in"), py::arg("F_out"), py::arg("model") = "sharp",
        py::arg("length") = 0.0, py::arg("D_h") = 0.0, py::arg("incompressible") = false,
        "Stagnation-pressure loss of an area change from F_in to F_out at a given\n"
        "flow, fed from (Pt0, Tt0, Y). model 'sharp' or 'conical' (length sets\n"
        "the cone). Returns AreaChangeLoss.");

  py::class_<combaero::AreaChangeFlow>(m, "AreaChangeFlow",
                                       "Area change flow at a given drive")
      .def_readonly("m_dot", &combaero::AreaChangeFlow::m_dot)
      .def_readonly("choked", &combaero::AreaChangeFlow::choked)
      .def_readonly("mach", &combaero::AreaChangeFlow::mach)
      .def_readonly("Re", &combaero::AreaChangeFlow::Re)
      .def_readonly("zeta", &combaero::AreaChangeFlow::zeta)
      .def_readonly("dm_dPt0", &combaero::AreaChangeFlow::dm_dPt0)
      .def_readonly("dm_dTt0", &combaero::AreaChangeFlow::dm_dTt0)
      .def_readonly("dm_dPt_target", &combaero::AreaChangeFlow::dm_dPt_target)
      .def_readonly("dm_dY", &combaero::AreaChangeFlow::dm_dY);

  m.def("area_change_flow", &combaero::area_change_flow,
        py::arg("Pt0"), py::arg("Tt0"), py::arg("Y"), py::arg("Pt_target"),
        py::arg("F_in"), py::arg("F_out"), py::arg("model") = "sharp",
        py::arg("length") = 0.0, py::arg("D_h") = 0.0, py::arg("incompressible") = false,
        py::arg("with_derivatives") = true,
        "Mass flow of an area change from F_in to F_out whose loss is the drive\n"
        "Pt0 - Pt_target (odd in it), saturating at the sonic flow, with its\n"
        "analytic partials. Returns AreaChangeFlow.");

  py::class_<combaero::AreaChangeElementFlow>(m, "AreaChangeElementFlow",
                                              "Area change element flow between two nodes")
      .def_readonly("m_dot", &combaero::AreaChangeElementFlow::m_dot)
      .def_readonly("choked", &combaero::AreaChangeElementFlow::choked)
      .def_readonly("mach", &combaero::AreaChangeElementFlow::mach)
      .def_readonly("Re", &combaero::AreaChangeElementFlow::Re)
      .def_readonly("zeta", &combaero::AreaChangeElementFlow::zeta)
      .def_readonly("dm_dPt_a", &combaero::AreaChangeElementFlow::dm_dPt_a)
      .def_readonly("dm_dTt_a", &combaero::AreaChangeElementFlow::dm_dTt_a)
      .def_readonly("dm_dPt_b", &combaero::AreaChangeElementFlow::dm_dPt_b)
      .def_readonly("dm_dTt_b", &combaero::AreaChangeElementFlow::dm_dTt_b)
      .def_readonly("dm_dY_a", &combaero::AreaChangeElementFlow::dm_dY_a)
      .def_readonly("dm_dY_b", &combaero::AreaChangeElementFlow::dm_dY_b);

  m.def("area_change_element_flow", &combaero::area_change_element_flow,
        py::arg("Pt_a"), py::arg("Tt_a"), py::arg("Y_a"), py::arg("Pt_b"), py::arg("Tt_b"),
        py::arg("Y_b"), py::arg("F_a"), py::arg("F_b"), py::arg("model") = "sharp",
        py::arg("length") = 0.0, py::arg("D_h") = 0.0, py::arg("incompressible") = false,
        py::arg("with_derivatives") = true,
        "The area change element's flow from node a (face F_a) to node b, both\n"
        "directions blended C1 through zero drive. Returns AreaChangeElementFlow.");

  // Compression and expansion along a path of given efficiency (#893).
  py::class_<combaero::TurbomachinePath>(
      m, "TurbomachinePath",
      "Outlet state of a compression or expansion at a total-to-total pressure\n"
      "ratio, with its Jacobian (#893).")
      .def_readonly("T2", &combaero::TurbomachinePath::T2)
      .def_readonly("dh", &combaero::TurbomachinePath::dh)
      .def_readonly("eta_poly", &combaero::TurbomachinePath::eta_poly)
      .def_readonly("eta_is", &combaero::TurbomachinePath::eta_is)
      .def_readonly("T2s", &combaero::TurbomachinePath::T2s)
      .def_readonly("dT2_dT1", &combaero::TurbomachinePath::dT2_dT1)
      .def_readonly("dT2_dPR", &combaero::TurbomachinePath::dT2_dPR)
      .def_readonly("dT2_deta", &combaero::TurbomachinePath::dT2_deta)
      .def_readonly("dT2_dY", &combaero::TurbomachinePath::dT2_dY)
      .def_readonly("ddh_dT1", &combaero::TurbomachinePath::ddh_dT1)
      .def_readonly("ddh_dPR", &combaero::TurbomachinePath::ddh_dPR)
      .def_readonly("ddh_deta", &combaero::TurbomachinePath::ddh_deta)
      .def_readonly("ddh_dY", &combaero::TurbomachinePath::ddh_dY);

  {
    // Exactly one efficiency: the polytropic one is the path's own, the
    // isentropic one the textbook end-state form.
    auto path = [](double T1, const std::vector<double>& X, double PR,
                   std::optional<double> eta_poly, std::optional<double> eta_is,
                   bool expansion) {
      if (eta_poly.has_value() == eta_is.has_value()) {
        throw py::value_error(std::string(expansion ? "expansion" : "compression") +
                              ": give exactly one of eta_poly and eta_is");
      }
      const bool poly = eta_poly.has_value();
      return combaero::turbomachine_path(T1, X, PR, poly ? *eta_poly : *eta_is, poly,
                                         expansion);
    };
    m.def(
        "compression",
        [path](double T1, const std::vector<double>& X, double PR,
               std::optional<double> eta_poly, std::optional<double> eta_is) {
          return path(T1, X, PR, eta_poly, eta_is, false);
        },
        py::arg("T1"), py::arg("X"), py::arg("PR"), py::kw_only(),
        py::arg("eta_poly") = py::none(), py::arg("eta_is") = py::none(),
        "Compression from total temperature T1 [K], mole fractions X, at the\n"
        "total pressure ratio PR = P2/P1 >= 1, with exactly one of eta_poly or\n"
        "eta_is in (0, 1]. Polytropic: S(T2) - S(T1) = (R/eta_poly) ln PR,\n"
        "exact for a variable cp. Isentropic: h2 = h1 + (h2s - h1)/eta_is. The\n"
        "composition is frozen. Returns TurbomachinePath (T2, dh > 0, both\n"
        "efficiencies, T2s, and the Jacobian of T2 and dh in T1, PR, eta, Y).");
    m.def(
        "expansion",
        [path](double T1, const std::vector<double>& X, double PR,
               std::optional<double> eta_poly, std::optional<double> eta_is) {
          return path(T1, X, PR, eta_poly, eta_is, true);
        },
        py::arg("T1"), py::arg("X"), py::arg("PR"), py::kw_only(),
        py::arg("eta_poly") = py::none(), py::arg("eta_is") = py::none(),
        "Expansion from total temperature T1 [K], mole fractions X, at the\n"
        "total pressure ratio PR = P1/P2 >= 1, with exactly one of eta_poly or\n"
        "eta_is in (0, 1]. Polytropic: S(T2) - S(T1) = -eta_poly R ln PR, exact\n"
        "for a variable cp. Isentropic: h2 = h1 - eta_is (h1 - h2s). The\n"
        "composition is frozen. Returns TurbomachinePath (T2, dh < 0, both\n"
        "efficiencies, T2s, and the Jacobian of T2 and dh in T1, PR, eta, Y).");
  }

  py::class_<combaero::AreaChangeResult>(
      m, "AreaChangeResult",
      "Result of sharp_area_change: dP and Jacobian partials")
      .def_readonly("dP", &combaero::AreaChangeResult::dP)
      .def_readonly("dS_dm", &combaero::AreaChangeResult::dS_dm)
      .def_readonly("dS_drho", &combaero::AreaChangeResult::dS_drho)
      .def_readonly("dS_dmu", &combaero::AreaChangeResult::dS_dmu);

  m.def("sharp_area_change", &combaero::sharp_area_change,
        py::arg("m_dot"), py::arg("rho"), py::arg("mu"),
        py::arg("F0"), py::arg("F1"), py::arg("m_scale") = 1e-4,
        py::arg("D_h") = 0.0,
        "Sharp-edge area change at a given density (incompressible): pressure\n"
        "drop and Jacobian partials, the Idelchik closures above.\n\n"
        "D_h overrides the circular hydraulic diameter when > 0.\n"
        "Returns AreaChangeResult with dP, dS_dm, dS_drho, dS_dmu.");

  m.def("conical_area_change", &combaero::conical_area_change,
        py::arg("m_dot"), py::arg("rho"), py::arg("mu"),
        py::arg("F0"), py::arg("F1"), py::arg("length"), py::arg("m_scale") = 1e-4,
        "Conical area change at a given density (incompressible): pressure drop\n"
        "and Jacobian partials.\n\n"
        "Returns AreaChangeResult with dP, dS_dm, dS_drho, dS_dmu.");

  // Tee junction (Bassett 2001)
  py::class_<solver::TeeJunctionResult>(
      m, "TeeJunctionResult",
      "Result of merging/branching_tee_residuals_and_jacobian for network solver")
      .def(py::init<>())
      .def_readonly("R_straight", &solver::TeeJunctionResult::R_straight)
      .def_readonly("R_branch", &solver::TeeJunctionResult::R_branch)
      .def_readonly("dR_straight_d_mdot_com",
                    &solver::TeeJunctionResult::dR_straight_d_mdot_com)
      .def_readonly("dR_straight_d_mdot_branch",
                    &solver::TeeJunctionResult::dR_straight_d_mdot_branch)
      .def_readonly("dR_straight_dP_static_com",
                    &solver::TeeJunctionResult::dR_straight_dP_static_com)
      .def_readonly("dR_straight_dT_com",
                    &solver::TeeJunctionResult::dR_straight_dT_com)
      .def_readonly("dR_straight_dY_com",
                    &solver::TeeJunctionResult::dR_straight_dY_com)
      .def_readonly("dR_branch_d_mdot_com",
                    &solver::TeeJunctionResult::dR_branch_d_mdot_com)
      .def_readonly("dR_branch_d_mdot_branch",
                    &solver::TeeJunctionResult::dR_branch_d_mdot_branch)
      .def_readonly("dR_branch_dP_static_com",
                    &solver::TeeJunctionResult::dR_branch_dP_static_com)
      .def_readonly("dR_branch_dT_com",
                    &solver::TeeJunctionResult::dR_branch_dT_com)
      .def_readonly("dR_branch_dY_com",
                    &solver::TeeJunctionResult::dR_branch_dY_com)
      .def_readonly("K_straight", &solver::TeeJunctionResult::K_straight)
      .def_readonly("K_branch", &solver::TeeJunctionResult::K_branch)
      .def_readonly("q", &solver::TeeJunctionResult::q)
      .def_readonly("blend_w", &solver::TeeJunctionResult::blend_w)
      .def_readonly("topology_valid", &solver::TeeJunctionResult::topology_valid)
      .def_readonly("status", &solver::TeeJunctionResult::status);

  m.def("merging_tee_residuals_and_jacobian",
        &solver::merging_tee_residuals_and_jacobian,
        py::arg("m_dot_com"), py::arg("m_dot_branch"),
        py::arg("dP0_straight"), py::arg("dP0_branch"),
        py::arg("P_static_com"), py::arg("T_com"), py::arg("Y_com"),
        py::arg("theta"), py::arg("psi"), py::arg("F_C"),
        py::arg("blend_k") = 30.0,
        "MergingTee (joining flow, Bassett 2001 type 6) residuals and Jacobians.\n\n"
        "Parameters:\n"
        "  m_dot_com    : total outlet flow C [kg/s]\n"
        "  m_dot_branch : branch inlet flow B [kg/s]\n"
        "  dP0_straight : p0_A - p0_C [Pa]\n"
        "  dP0_branch   : p0_B - p0_C [Pa]\n"
        "  P_static_com : static pressure at port C [Pa] (for density)\n"
        "  T_com        : temperature at port C [K]\n"
        "  Y_com        : mass fractions at port C [-]\n"
        "  theta        : branch angle [rad] (valid: (0, pi/2])\n"
        "  psi          : area ratio F_C/F_B [-] (valid: >= 0.05)\n"
        "  F_C          : main duct area [m^2]\n"
        "  blend_k      : tanh blend sharpness (default 30.0)\n\n"
        "Returns: TeeJunctionResult");

  m.def("branching_tee_residuals_and_jacobian",
        &solver::branching_tee_residuals_and_jacobian,
        py::arg("m_dot_com"), py::arg("m_dot_branch"),
        py::arg("dP0_straight"), py::arg("dP0_branch"),
        py::arg("P_static_com"), py::arg("T_com"), py::arg("Y_com"),
        py::arg("theta"), py::arg("psi"), py::arg("F_C"),
        py::arg("blend_k") = 30.0,
        "BranchingTee (separating flow, Bassett 2001 type 3) residuals and Jacobians.\n\n"
        "Parameters:\n"
        "  m_dot_com    : total inlet flow C [kg/s]\n"
        "  m_dot_branch : branch outlet flow B [kg/s]\n"
        "  dP0_straight : p0_C - p0_A [Pa]\n"
        "  dP0_branch   : p0_C - p0_B [Pa]\n"
        "  P_static_com : static pressure at port C [Pa] (for density)\n"
        "  T_com        : temperature at port C [K]\n"
        "  Y_com        : mass fractions at port C [-]\n"
        "  theta        : branch angle [rad] (valid: (0, pi/2])\n"
        "  psi          : area ratio F_C/F_B [-] (valid: >= 0.05)\n"
        "  F_C          : main duct area [m^2]\n"
        "  blend_k      : tanh blend sharpness (default 30.0)\n\n"
        "Returns: TeeJunctionResult");

  // Tee junction pure K functions and helpers (for testing and diagnostics)
  m.def("tee_K1", &combaero::K1, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K1: branch loss, separating flow type 1 (reverse of K6).");
  m.def("tee_K2", &combaero::K5, py::arg("q"),
        "Bassett K2: straight-through loss, separating flow types 1/2. Same expression as K5.");
  m.def("tee_K3", &combaero::K3, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K3: branch-to-downstream, separating flow type 2 (reverse of K4).");
  m.def("tee_K4", &combaero::K4, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K4: branch-to-downstream, separating flow type 2.");
  m.def("tee_K5", &combaero::K5, py::arg("q"),
        "Bassett K5: straight-through loss, separating flow. K5(q) = q^2-1.5q+0.5");
  m.def("tee_K6", &combaero::K6, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K6: branch loss, separating flow type 3.");
  m.def("tee_K7", &combaero::K7, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K7: straight path, joining flow type 4.");
  m.def("tee_K8", &combaero::K8, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K8: branch path, joining flow type 4.");
  m.def("tee_K9", &combaero::K9, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K9: branch A path, joining flow type 5.");
  m.def("tee_K10", &combaero::K10, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K10: branch C path, joining flow type 5.");
  m.def("tee_K11", &combaero::K11, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K11: straight-to-outlet, joining flow type 6.");
  m.def("tee_K12", &combaero::K12, py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Bassett K12: branch-to-outlet, joining flow type 6.");
  m.def("tee_blend_weight", &combaero::blend_weight,
        py::arg("r"), py::arg("k") = 30.0,
        "Smooth topology blend weight: 0.5*(1+tanh(k*r)).");
  m.def("tee_check_inputs",
        [](double q, double psi, double theta) {
            auto s = combaero::tee_check_inputs(q, psi, theta);
            return py::dict(py::arg("q_in_range") = s.q_in_range,
                            py::arg("psi_valid") = s.psi_valid,
                            py::arg("theta_valid") = s.theta_valid,
                            py::arg("valid") = s.valid());
        },
        py::arg("q"), py::arg("psi"), py::arg("theta"),
        "Check if tee junction inputs are in the Bassett 2001 validated range.\n"
        "Returns dict with q_in_range, psi_valid, theta_valid, valid keys.");
  m.def("merging_tee_K_straight", &combaero::merging_tee_K_straight,
        py::arg("q"), py::arg("psi"), py::arg("theta"), py::arg("blend_k") = 30.0,
        "MergingTee effective straight-path K (blended, with soft bounds).");
  m.def("merging_tee_K_branch", &combaero::merging_tee_K_branch,
        py::arg("q"), py::arg("psi"), py::arg("theta"), py::arg("blend_k") = 30.0,
        "MergingTee effective branch-path K (blended, with soft bounds).");
  m.def("branching_tee_K_straight", &combaero::branching_tee_K_straight,
        py::arg("q"), py::arg("psi"), py::arg("theta"), py::arg("blend_k") = 30.0,
        "BranchingTee effective straight-path K (blended, with soft bounds).");
  m.def("branching_tee_K_branch", &combaero::branching_tee_K_branch,
        py::arg("q"), py::arg("psi"), py::arg("theta"), py::arg("blend_k") = 30.0,
        "BranchingTee effective branch-path K (blended, with soft bounds).");

  // Compressible Unified0D tee junction (Mynard & Valen-Sendstad 2015, O(M^2) extension)
  py::class_<solver::BranchInput>(m, "BranchInput",
      "Per-branch face state for compressible_branching/merging_tee_rj.")
      .def(py::init<>())
      .def_readwrite("P_static",  &solver::BranchInput::P_static)
      .def_readwrite("Pt",        &solver::BranchInput::Pt)
      .def_readwrite("T",         &solver::BranchInput::T)
      .def_readwrite("m_dot",     &solver::BranchInput::m_dot)
      .def_readwrite("A",         &solver::BranchInput::A)
      .def_readwrite("theta",     &solver::BranchInput::theta)
      .def_readwrite("gamma_eff", &solver::BranchInput::gamma_eff)
      .def_readwrite("R_gas",     &solver::BranchInput::R_gas);

  py::class_<solver::CompressibleTeeResult>(m, "CompressibleTeeResult",
      "Result of compressible_branching/merging_tee_rj.")
      .def(py::init<>())
      .def_readonly("R_0",              &solver::CompressibleTeeResult::R_0)
      .def_readonly("R_1",              &solver::CompressibleTeeResult::R_1)
      .def_readonly("dR0_dPt_com",      &solver::CompressibleTeeResult::dR0_dPt_com)
      .def_readonly("dR0_dP_com",       &solver::CompressibleTeeResult::dR0_dP_com)
      .def_readonly("dR0_dT_com",       &solver::CompressibleTeeResult::dR0_dT_com)
      .def_readonly("dR0_dmdot_com",    &solver::CompressibleTeeResult::dR0_dmdot_com)
      .def_readonly("dR0_dPt_str",      &solver::CompressibleTeeResult::dR0_dPt_str)
      .def_readonly("dR0_dP_str",       &solver::CompressibleTeeResult::dR0_dP_str)
      .def_readonly("dR0_dT_str",       &solver::CompressibleTeeResult::dR0_dT_str)
      .def_readonly("dR0_dPt_bra",      &solver::CompressibleTeeResult::dR0_dPt_bra)
      .def_readonly("dR0_dP_bra",       &solver::CompressibleTeeResult::dR0_dP_bra)
      .def_readonly("dR0_dT_bra",       &solver::CompressibleTeeResult::dR0_dT_bra)
      .def_readonly("dR0_dmdot_branch", &solver::CompressibleTeeResult::dR0_dmdot_branch)
      .def_readonly("dR1_dPt_com",      &solver::CompressibleTeeResult::dR1_dPt_com)
      .def_readonly("dR1_dP_com",       &solver::CompressibleTeeResult::dR1_dP_com)
      .def_readonly("dR1_dT_com",       &solver::CompressibleTeeResult::dR1_dT_com)
      .def_readonly("dR1_dmdot_com",    &solver::CompressibleTeeResult::dR1_dmdot_com)
      .def_readonly("dR1_dPt_str",      &solver::CompressibleTeeResult::dR1_dPt_str)
      .def_readonly("dR1_dP_str",       &solver::CompressibleTeeResult::dR1_dP_str)
      .def_readonly("dR1_dT_str",       &solver::CompressibleTeeResult::dR1_dT_str)
      .def_readonly("dR1_dPt_bra",      &solver::CompressibleTeeResult::dR1_dPt_bra)
      .def_readonly("dR1_dP_bra",       &solver::CompressibleTeeResult::dR1_dP_bra)
      .def_readonly("dR1_dT_bra",       &solver::CompressibleTeeResult::dR1_dT_bra)
      .def_readonly("dR1_dmdot_branch", &solver::CompressibleTeeResult::dR1_dmdot_branch);

  m.def("compressible_branching_tee_rj",
        &solver::compressible_branching_tee_rj,
        py::arg("com"), py::arg("str"), py::arg("bra"),
        "Compressible branching tee (Unified0D, O(M^2)).\n"
        "common=supplier, straight+branch=collectors.\n"
        "  R_0 = Pt_com - Pt_str - K_str * q_ref\n"
        "  R_1 = Pt_com - Pt_bra - K_bra * q_ref\n"
        "Returns: CompressibleTeeResult");

  m.def("compressible_merging_tee_rj",
        &solver::compressible_merging_tee_rj,
        py::arg("com"), py::arg("str"), py::arg("bra"),
        "Compressible merging tee (Unified0D, O(M^2)).\n"
        "straight+branch=suppliers, common=collector.\n"
        "  R_0 = P_static_str - P_static_bra  (supplier pressure equality)\n"
        "  R_1 = p0_dat - Pt_com - K_com * q_ref\n"
        "Returns: CompressibleTeeResult");

  // Species metadata helpers
  m.def("num_species", &num_species,
        "Number of thermo species in the internal tables.");

  m.def("species_name", &species_name, py::arg("index"),
        "Return the canonical species name for a given index.");

  m.def("species_index_from_name", &species_index_from_name, py::arg("name"),
        "Return the species index for a given canonical name.");

  m.def("species_molar_mass", &species_molar_mass, py::arg("index"),
        "Return the molar mass [g/mol] for a given species index.");

  m.def("species_molar_mass_from_name", &species_molar_mass_from_name,
        py::arg("name"),
        "Return the molar mass [g/mol] for a given species name.");

  // Common names for species symbols
  m.def(
      "formula_to_name", []() { return combaero::formula_to_name; },
      "Mapping from canonical species symbols (e.g. 'CH4') to human-readable "
      "common names (e.g. 'Methane').");

  m.def(
      "name_to_formula", []() { return combaero::name_to_formula; },
      "Mapping from human-readable common names (e.g. 'Methane') to canonical "
      "species symbols (e.g. 'CH4').");

  m.def("common_name", &combaero::common_name, py::arg("formula"),
        "Get common name from formula (e.g., 'CH4' -> 'Methane').");

  m.def("formula", &combaero::formula, py::arg("name"),
        "Get formula from common name (e.g., 'Methane' -> 'CH4').");

  m.def(
      "mass_to_mole",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_arr) {
        return mass_to_mole(to_vec(Y_arr));
      },
      py::arg("Y"), "Convert mass fractions Y to mole fractions X");

  m.def(
      "mixture_h",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return h(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mixture enthalpy at temperature T and mole fractions X.");

  // Thermodynamic mixture properties
  m.def(
      "cp",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return cp(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mixture isobaric heat capacity Cp(T, X) [J/(mol*K)].");

  m.def(
      "h",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return h(T, X);
      },
      py::arg("T"), py::arg("X"), "Mixture enthalpy h(T, X) [J/mol].");

  m.def(
      "s",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P, double P_ref) {
        auto X = to_vec(X_arr);
        return s(T, X, P, P_ref);
      },
      py::arg("T"), py::arg("X"), py::arg("P"), py::arg("P_ref") = 101325.0,
      "Mixture entropy s(T, X, P, P_ref) [J/(mol*K)].");

  m.def(
      "cv",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return cv(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mixture isochoric heat capacity Cv(T, X) [J/(mol*K)].");

  m.def(
      "u",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return u(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mixture internal energy u(T, X) [J/mol].\n\n"
      "For ideal gas: u = h - RT\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  X : mole fractions [mol/mol]\n\n"
      "Returns: internal energy [J/mol]");

  m.def(
      "density",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return density(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Mixture density rho(T, P, X) [kg/m^3].");

  m.def(
      "specific_gas_constant",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return specific_gas_constant(X);
      },
      py::arg("X"), "Mixture specific gas constant R_s(X) [J/(kg*K)].");

  m.def(
      "isentropic_expansion_coefficient",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return isentropic_expansion_coefficient(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Isentropic expansion coefficient gamma(T, X) = Cp/Cv.");

  m.def(
      "speed_of_sound",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return speed_of_sound(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Speed of sound [m/s].\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  X : mole fractions [mol/mol]\n\n"
      "Returns: speed of sound [m/s]");

  m.def(
      "cp_secant_mass_and_jacobian",
      [](double T1, double T2,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        return cp_secant_mass_and_jacobian(T1, T2, to_vec(X_arr));
      },
      py::arg("T1"), py::arg("T2"), py::arg("X"),
      "Secant cp over [T1, T2], (h(T1) - h(T2)) / (T1 - T2) [J/(kg*K)], and\n"
      "its partials in T1 and T2 (cp(T1) in the limit T2 -> T1). The capacity\n"
      "rate of an epsilon-NTU wall exchange (#627).");

  m.def(
      "cp_secant_mass_dY",
      [](double T1, double T2,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        return cp_secant_mass_dY(T1, T2, to_vec(X_arr));
      },
      py::arg("T1"), py::arg("T2"), py::arg("X"),
      "The secant cp's gradient in the MASS fractions [J/(kg*K)], one entry\n"
      "per species: a wall side's capacity rate |m| cp moves with its\n"
      "composition (#717).");

  m.def(
      "cp_mass",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return cp_mass(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass-specific heat capacity at constant pressure.\n\n"
      "Cp_mass = Cp_molar / mwmix * 1000\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  X : mole fractions [mol/mol]\n\n"
      "Returns: heat capacity [J/(kg*K)]");

  m.def(
      "cv_mass",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return cv_mass(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass-specific heat capacity at constant volume.\n\n"
      "Cv_mass = Cv_molar / mwmix * 1000\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  X : mole fractions [mol/mol]\n\n"
      "Returns: heat capacity [J/(kg*K)]");

  m.def(
      "h_mass",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return h_mass(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass-specific enthalpy.\n\n"
      "h_mass = h_molar / mwmix * 1000\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  X : mole fractions [mol/mol]\n\n"
      "Returns: enthalpy [J/kg]");

  m.def(
      "s_mass",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P, double P_ref = 101325.0) {
        auto X = to_vec(X_arr);
        return s_mass(T, X, P, P_ref);
      },
      py::arg("T"), py::arg("X"), py::arg("P"), py::arg("P_ref") = 101325.0,
      "Mass-specific entropy.\n\n"
      "s_mass = s_molar / mwmix * 1000\n\n"
      "Parameters:\n"
      "  T     : temperature [K]\n"
      "  X     : mole fractions [mol/mol]\n"
      "  P     : pressure [Pa]\n"
      "  P_ref : reference pressure [Pa] (default: 101325)\n\n"
      "Returns: entropy [J/(kg*K)]");

  m.def(
      "u_mass",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return u_mass(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass-specific internal energy.\n\n"
      "u_mass = u_molar / mwmix * 1000\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  X : mole fractions [mol/mol]\n\n"
      "Returns: internal energy [J/kg]");

  m.def("molar_volume", &molar_volume, py::arg("T"), py::arg("P"),
        "Ideal gas molar volume V_m = R*T/P [m^3/mol].");

  // Transport properties
  py::class_<TransportDerivatives>(m, "TransportDerivatives")
      .def_readonly("mu", &TransportDerivatives::mu)
      .def_readonly("dmu_dT", &TransportDerivatives::dmu_dT)
      .def_readonly("k", &TransportDerivatives::k)
      .def_readonly("dk_dT", &TransportDerivatives::dk_dT);
  m.def(
      "transport_and_dT",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return transport_and_dT(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Mixture viscosity and thermal conductivity with analytic d/dT (C2 in T).");
  py::class_<Omega22Result>(m, "Omega22Result")
      .def_readonly("omega", &Omega22Result::omega)
      .def_readonly("domega_dTstar", &Omega22Result::domega_dTstar);
  m.def("omega22_and_derivative", &omega22_and_derivative, py::arg("T_star"),
        py::arg("delta_star"),
        "Reduced collision integral Omega*(2,2) and dOmega*/dT*: quintic Hermite\n"
        "in (ln T*, ln Omega*) with tabulated slopes and curvatures (C2);\n"
        "delta* = 0 exact Lennard-Jones.");
  m.def("dcp_R_dT", &dcp_R_dT, py::arg("species_idx"), py::arg("T"),
        "d(cp/R)/dT of one species [1/K].");
  m.def(
      "set_thermo_model",
      [](const std::string& name) {
        if (name == "nasa9") {
          set_thermo_model(ThermoModel::NASA9);
        } else if (name == "nasa9_c2") {
          set_thermo_model(ThermoModel::NASA9_C2);
        } else {
          throw std::invalid_argument(
              "set_thermo_model: expected 'nasa9' or 'nasa9_c2', got '" + name + "'");
        }
      },
      py::arg("name"),
      "Select the species thermo every property evaluates (#549): 'nasa9_c2'\n"
      "(the default: each interval join bridged by a C2 cp; h and s identical\n"
      "to NASA-9 outside the |T - T_join| < 0.05 T_join windows) or 'nasa9'\n"
      "(the published polynomials). Process-wide.");
  m.def(
      "thermo_model",
      []() {
        return thermo_model() == ThermoModel::NASA9_C2 ? std::string("nasa9_c2")
                                                       : std::string("nasa9");
      },
      "The species thermo in use: 'nasa9' or 'nasa9_c2' (#549).");
  m.def(
      "viscosity",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return viscosity(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Dynamic viscosity mu(T, P, X) [Pa*s].");

  m.def(
      "thermal_conductivity",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return thermal_conductivity(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Thermal conductivity lambda(T, P, X) [W/(m*K)].");

  m.def(
      "kinematic_viscosity",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return kinematic_viscosity(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Kinematic viscosity nu(T, P, X) [m^2/s].");

  m.def(
      "prandtl",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return prandtl(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Prandtl number Pr(T, P, X) [-].");

  m.def(
      "thermal_diffusivity",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return thermal_diffusivity(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Thermal diffusivity alpha(T, P, X) [m^2/s].");

  m.def(
      "reynolds",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double V, double L) {
        auto X = to_vec(X_arr);
        return reynolds(T, P, X, V, L);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("V"), py::arg("L"),
      "Reynolds number Re = rho*V*L/mu [-].");

  m.def(
      "peclet",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double V, double L) {
        auto X = to_vec(X_arr);
        return peclet(T, P, X, V, L);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("V"), py::arg("L"),
      "Peclet number Pe = V*L/alpha (thermal) [-].");

  // TransportState struct binding - Bundle of transport properties
  py::class_<TransportState>(
      m, "TransportState",
      "Bundle of transport properties for a gas mixture.\n\n"
      "All properties are read-only attributes computed in a single call.\n"
      "Provides IDE autocomplete and type safety.\n\n"
      "Attributes:\n"
      "  T     : Temperature [K] (input, echoed back)\n"
      "  P     : Pressure [Pa] (input, echoed back)\n"
      "  rho   : Density [kg/m^3]\n"
      "  mu    : Dynamic viscosity [Pa*s]\n"
      "  k     : Thermal conductivity [W/(m*K)]\n"
      "  nu    : Kinematic viscosity [m^2/s]\n"
      "  alpha : Thermal diffusivity [m^2/s]\n"
      "  Pr    : Prandtl number [-]\n"
      "  cp    : Specific heat at constant pressure [J/(kg*K)]\n"
      "  cv    : Specific heat at constant volume [J/(kg*K)]\n"
      "  gamma : Heat capacity ratio cp/cv [-]\n"
      "  a     : Speed of sound [m/s]")
      .def_readonly("T", &TransportState::T, "Temperature [K]")
      .def_readonly("P", &TransportState::P, "Pressure [Pa]")
      .def_readonly("rho", &TransportState::rho, "Density [kg/m^3]")
      .def_readonly("mu", &TransportState::mu, "Dynamic viscosity [Pa*s]")
      .def_readonly("k", &TransportState::k, "Thermal conductivity [W/(m*K)]")
      .def_readonly("nu", &TransportState::nu, "Kinematic viscosity [m^2/s]")
      .def_readonly("alpha", &TransportState::alpha,
                    "Thermal diffusivity [m^2/s]")
      .def_readonly("Pr", &TransportState::Pr, "Prandtl number [-]")
      .def_readonly("cp", &TransportState::cp,
                    "Specific heat at constant pressure [J/(kg*K)]")
      .def_readonly("cv", &TransportState::cv,
                    "Specific heat at constant volume [J/(kg*K)]")
      .def_readonly("gamma", &TransportState::gamma,
                    "Heat capacity ratio cp/cv [-]")
      .def_readonly("a", &TransportState::a, "Speed of sound [m/s]")
      .def("__repr__", [](const TransportState &s) {
        return "<TransportState: T=" + std::to_string(s.T) +
               " K, "
               "P=" +
               std::to_string(s.P) +
               " Pa, "
               "mu=" +
               std::to_string(s.mu) +
               " Pa*s, "
               "k=" +
               std::to_string(s.k) +
               " W/(m*K), "
               "Pr=" +
               std::to_string(s.Pr) +
               ", "
               "gamma=" +
               std::to_string(s.gamma) +
               ", "
               "a=" +
               std::to_string(s.a) + " m/s>";
      });

  // EquilibriumResult struct binding
  py::class_<EquilibriumResult>(
      m, "EquilibriumResult",
      "Result of an equilibrium solve.\n\n"
      "Carries the final state plus solver diagnostics, suitable for use\n"
      "in network solvers where non-convergence must be detectable.\n\n"
      "Attributes:\n"
      "  state     : equilibrium State (T, P, X)\n"
      "  T_in      : input temperature [K]\n"
      "  delta_T   : state.T - T_in [K] (negative = endothermic)\n"
      "  converged : solver convergence flag")
      .def_readonly("state", &EquilibriumResult::state,
                    "Equilibrium State (T, P, X)")
      .def_readonly("T_in", &EquilibriumResult::T_in, "Input temperature [K]")
      .def_readonly("delta_T", &EquilibriumResult::delta_T,
                    "state.T - T_in [K]")
      .def_readonly("converged", &EquilibriumResult::converged,
                    "Solver convergence flag")
      .def("__repr__", [](const EquilibriumResult &r) {
        return "<EquilibriumResult: T=" + std::to_string(r.state.T) +
               " K, delta_T=" + std::to_string(r.delta_T) +
               " K, converged=" + (r.converged ? "True" : "False") + ">";
      });

  m.def(
      "transport_state",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        if (T <= 0.0) throw std::invalid_argument("Temperature must be strictly positive");
        if (P <= 0.0) throw std::invalid_argument("Pressure must be strictly positive");
        auto X = to_vec(X_arr);
        return transport_state(T, P, X);
      },
      py::arg("T"), py::arg("P"), py::arg("X"),
      "Compute all transport properties at once.\n\n"
      "Convenience function that computes all transport properties\n"
      "for a gas mixture in a single call. Returns TransportState struct\n"
      "with read-only attributes for IDE autocomplete support.\n\n"
      "Parameters:\n"
      "  T : temperature [K]\n"
      "  P : pressure [Pa]\n"
      "  X : mole fractions [-]\n\n"
      "Returns: TransportState object with attributes:\n"
      "  T, P, rho, mu, k, nu, alpha, Pr, cp\n\n"
      "Example:\n"
      "  >>> X = cb.species.dry_air()\n"
      "  >>> state = cb.transport_state(T=300, P=101325, X=X)\n"
      "  >>> print(state.mu)  # IDE autocomplete works!\n"
      "  1.85e-05\n"
      "  >>> print(state.Pr)\n"
      "  0.707");

  // CompleteState struct binding - Bundle of thermo + transport properties
  py::class_<CompleteState>(
      m, "CompleteState",
      "Bundle of thermodynamic and transport properties for a gas mixture.\n\n"
      "All properties are read-only attributes computed in a single call.\n"
      "Provides IDE autocomplete and type safety.\n\n"
      "Attributes:\n"
      "  thermo    : ThermoState with all 16 thermodynamic properties\n"
      "  transport : TransportState with all 9 transport properties")
      .def_readonly("thermo", &CompleteState::thermo,
                    "ThermoState with all thermodynamic properties")
      .def_readonly("transport", &CompleteState::transport,
                    "TransportState with all transport properties")
      .def("__repr__", [](const CompleteState &s) {
        return "<CompleteState: T=" + std::to_string(s.thermo.T) +
               " K, "
               "P=" +
               std::to_string(s.thermo.P) +
               " Pa, "
               "rho=" +
               std::to_string(s.thermo.rho) +
               " kg/m^3, "
               "mu=" +
               std::to_string(s.transport.mu) + " Pa*s>";
      });

  m.def(
      "complete_state",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P_ref) {
        auto X = to_vec(X_arr);
        return complete_state(T, P, X, P_ref);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("P_ref") = 101325.0,
      "Compute all thermodynamic and transport properties at once.\n\n"
      "Convenience function that computes all thermodynamic and transport\n"
      "properties for a gas mixture in a single call. Returns CompleteState\n"
      "struct with nested thermo and transport states.\n\n"
      "Parameters:\n"
      "  T     : temperature [K]\n"
      "  P     : pressure [Pa]\n"
      "  X     : mole fractions [-]\n"
      "  P_ref : reference pressure for entropy [Pa] (default: 101325.0)\n\n"
      "Returns: CompleteState object with nested attributes:\n"
      "  thermo.T, thermo.P, thermo.rho, thermo.cp, thermo.h, etc. (16 "
      "properties)\n"
      "  transport.mu, transport.k, transport.Pr, etc. (9 properties)\n\n"
      "Example:\n"
      "  >>> X = cb.species.dry_air()\n"
      "  >>> state = cb.complete_state(T=300, P=101325, X=X)\n"
      "  >>> print(state.thermo.h)      # Thermodynamic property\n"
      "  -103.6\n"
      "  >>> print(state.transport.mu)  # Transport property\n"
      "  1.68e-05");

  // Humid air utilities
  m.def("dry_air", &dry_air,
        "Get standard dry air composition (mole fractions) as a vector of size "
        "num_species().\n\n"
        "Orders determined by NASA-9 thermo data index.");

  m.def("humid_air_composition", &humid_air_composition, py::arg("T"),
        py::arg("P"), py::arg("RH"),
        "Humid air mole-fraction composition for temperature T [K], pressure P "
        "[Pa], and relative humidity RH (0-1).");

  m.def("dewpoint", &dewpoint, py::arg("T"), py::arg("P"), py::arg("RH"),
        "Dewpoint temperature [K] for ambient T [K], pressure P [Pa], and RH "
        "(0-1).\n\n"
        "RH in (0, 1]; RH = 1 returns T (saturated). RH > 1 (supersaturated), RH <= 0 "
        "or non-finite RH raises ValueError.");

  m.def("relative_humidity_from_dewpoint", &relative_humidity_from_dewpoint,
        py::arg("T"), py::arg("Tdp"), py::arg("P"),
        "Relative humidity (0-1) from dry-bulb T [K], dewpoint Tdp [K], and "
        "pressure P [Pa].\n\n"
        "P is accepted for API symmetry but not used (dewpoint is "
        "pressure-independent\n"
        "at atmospheric conditions).");

  m.def("saturation_vapor_pressure", &saturation_vapor_pressure, py::arg("T"),
        "Water vapor saturation pressure [Pa] using Huang (2018) formula.\n\n"
        "Validated range: -100\u00b0C to +100\u00b0C (173.15\u2013373.15 K).\n"
        "Extrapolation up to 5% of range span (~10 K) is allowed silently.\n"
        "Raises IndexError beyond that. Over ice at T <= 273.15 K, so an RH\n"
        "below 0 C is relative to ice.\n\n"
        "T : temperature [K]\n\n"
        "Returns: saturation vapor pressure [Pa]");

  m.def("vapor_pressure", &vapor_pressure, py::arg("T"), py::arg("RH"),
        "Actual vapor pressure [Pa] from temperature and relative humidity.\n\n"
        "T  : temperature [K]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: vapor pressure [Pa]");

  m.def("humidity_ratio", &humidity_ratio, py::arg("T"), py::arg("P"),
        py::arg("RH"),
        "Humidity ratio (kg water vapor / kg dry air).\n\n"
        "T  : temperature [K]\n"
        "P  : total pressure [Pa]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: humidity ratio W [-]");

  m.def("water_vapor_mole_fraction", &water_vapor_mole_fraction, py::arg("T"),
        py::arg("P"), py::arg("RH"),
        "Mole fraction of water vapor in humid air.\n\n"
        "T  : temperature [K]\n"
        "P  : total pressure [Pa]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: mole fraction of H2O [-]");

  m.def("wet_bulb_temperature", &wet_bulb_temperature, py::arg("T"),
        py::arg("P"), py::arg("RH"),
        "Wet-bulb temperature [K] via Newton-Raphson on the psychrometric "
        "equation.\n\n"
        "Solves: cp_air*(T - Twb) = hfg(Twb)*(Ws(Twb) - W)\n"
        "where hfg(Twb) = 2501000 - 2381*(Twb - 273.15) J/kg (ASHRAE linear "
        "fit).\n\n"
        "T  : dry-bulb temperature [K]\n"
        "P  : total pressure [Pa]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: wet-bulb temperature [K]");

  m.def("humid_air_enthalpy", &humid_air_enthalpy, py::arg("T"), py::arg("P"),
        py::arg("RH"),
        "Enthalpy of humid air [J/kg dry air] using ASHRAE psychrometric "
        "constants.\n\n"
        "h = cp_air*(T - 273.15) + W*(2501000 + 1860*(T - 273.15))\n"
        "where cp_air=1005, h_fg=2501000, cp_wv=1860 J/(kg\u00b7K).\n\n"
        "T  : temperature [K]\n"
        "P  : total pressure [Pa]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: enthalpy [J/kg dry air]");

  m.def("humid_air_enthalpy_nasa9", &humid_air_enthalpy_nasa9, py::arg("T"),
        py::arg("P"), py::arg("RH"),
        "Enthalpy of humid air [J/kg dry air] using NASA-9 thermodynamic "
        "data.\n\n"
        "Computes h(T, X) from the NASA-9 polynomial database and converts\n"
        "to J/kg dry air basis.\n\n"
        "NOT on the psychrometric datum of humid_air_enthalpy(): the NASA-9\n"
        "datum is the elements at 298.15 K, so water vapour carries its\n"
        "formation enthalpy and the two differ by hundreds of kJ/kg at the\n"
        "same state. Only differences at fixed humidity ratio compare.\n\n"
        "T  : temperature [K] (bounded by saturation_vapor_pressure's range,\n"
        "     163.15-383.15 K, because RH is evaluated at T)\n"
        "P  : total pressure [Pa]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: enthalpy [J/kg dry air]");

  m.def("humid_air_density", &humid_air_density, py::arg("T"), py::arg("P"),
        py::arg("RH"),
        "Density of humid air [kg/m\u00b3] using ideal gas law.\n\n"
        "R_ha = R_da * (1 + W/0.621945) / (1 + W), rho = P / (R_ha * T)\n\n"
        "T  : temperature [K]\n"
        "P  : total pressure [Pa]\n"
        "RH : relative humidity (0-1)\n\n"
        "Returns: density [kg/m\u00b3]");

  // Composition conversion helpers
  m.def(
      "mole_to_mass",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return mole_to_mass(X);
      },
      py::arg("X"), "Convert mole fractions X to mass fractions Y.");

  m.def(
      "mwmix", [](py::array_t<double> X) { return mwmix(to_vec(X)); },
      py::arg("X"),
      "Calculate mixture molecular weight [g/mol] from mole fractions.");

  // Fraction utilities
  m.def(
      "normalize_fractions",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> arr) {
        auto v = to_vec(arr);
        return normalize_fractions(v);
      },
      py::arg("fractions"),
      "Normalize a vector of fractions so it sums to 1.0 (or return all zeros "
      "if input is all zeros).");

  m.def(
      "convert_to_dry_fractions",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> arr) {
        auto v = to_vec(arr);
        return convert_to_dry_fractions(v);
      },
      py::arg("mole_fractions"),
      "Convert mole fractions to dry fractions by removing water vapor and "
      "renormalizing.");

  // Equivalence ratio (mole basis)
  m.def(
      "equivalence_ratio_mole",
      [](py::array_t<double, py::array::c_style | py::array::forcecast>
             X_mix_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             X_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             X_ox_arr) {
        auto X_mix = to_vec(X_mix_arr);
        auto X_fuel = to_vec(X_fuel_arr);
        auto X_ox = to_vec(X_ox_arr);
        return equivalence_ratio_mole(X_mix, X_fuel, X_ox);
      },
      py::arg("X_mix"), py::arg("X_fuel"), py::arg("X_ox"),
      "Equivalence ratio phi (mole basis) for mixture X_mix formed from fuel "
      "X_fuel and oxidizer X_ox.");

  m.def(
      "set_equivalence_ratio_mole",
      [](double phi,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             X_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             X_ox_arr) {
        auto X_fuel = to_vec(X_fuel_arr);
        auto X_ox = to_vec(X_ox_arr);
        return set_equivalence_ratio_mole(phi, X_fuel, X_ox);
      },
      py::arg("phi"), py::arg("X_fuel"), py::arg("X_ox"),
      "Construct mole-fraction mixture X_mix with target phi from fuel X_fuel "
      "and oxidizer X_ox.");

  // Equivalence ratio (mass basis)
  m.def(
      "equivalence_ratio_mass",
      [](py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_mix_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_ox_arr) {
        auto Y_mix = to_vec(Y_mix_arr);
        auto Y_fuel = to_vec(Y_fuel_arr);
        auto Y_ox = to_vec(Y_ox_arr);
        return equivalence_ratio_mass(Y_mix, Y_fuel, Y_ox);
      },
      py::arg("Y_mix"), py::arg("Y_fuel"), py::arg("Y_ox"),
      "Equivalence ratio phi (mass basis) for mixture Y_mix formed from fuel "
      "Y_fuel and oxidizer Y_ox.");

  m.def(
      "set_equivalence_ratio_mass",
      [](double phi,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_ox_arr) {
        auto Y_fuel = to_vec(Y_fuel_arr);
        auto Y_ox = to_vec(Y_ox_arr);
        return set_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
      },
      py::arg("phi"), py::arg("Y_fuel"), py::arg("Y_ox"),
      "Construct mass-fraction mixture Y_mix with target phi from fuel Y_fuel "
      "and oxidizer Y_ox.");

  // Bilger-based helpers (mass basis)
  m.def(
      "bilger_stoich_mixture_fraction_mass",
      [](py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_ox_arr) {
        auto Y_fuel = to_vec(Y_fuel_arr);
        auto Y_ox = to_vec(Y_ox_arr);
        return bilger_stoich_mixture_fraction_mass(Y_fuel, Y_ox);
      },
      py::arg("Y_fuel"), py::arg("Y_ox"),
      "Stoichiometric Bilger mixture fraction Z_st (mass basis) for fuel "
      "Y_fuel and oxidizer Y_ox.");

  m.def(
      "bilger_Z_from_equivalence_ratio_mass",
      [](double phi,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_ox_arr) {
        auto Y_fuel = to_vec(Y_fuel_arr);
        auto Y_ox = to_vec(Y_ox_arr);
        return bilger_Z_from_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
      },
      py::arg("phi"), py::arg("Y_fuel"), py::arg("Y_ox"),
      "Convert mass-basis equivalence ratio phi to Bilger mixture fraction Z "
      "for fuel/oxidizer streams.");

  m.def(
      "equivalence_ratio_from_bilger_Z_mass",
      [](double Z,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_fuel_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             Y_ox_arr) {
        auto Y_fuel = to_vec(Y_fuel_arr);
        auto Y_ox = to_vec(Y_ox_arr);
        return equivalence_ratio_from_bilger_Z_mass(Z, Y_fuel, Y_ox);
      },
      py::arg("Z"), py::arg("Y_fuel"), py::arg("Y_ox"),
      "Convert Bilger mixture fraction Z to mass-basis equivalence ratio phi "
      "for fuel/oxidizer streams.");

  // Inverse temperature calculations
  m.def(
      "calc_T_from_h",
      [](double h_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_h(h_target, X, T_guess, tol, max_iter);
      },
      py::arg("h_target"), py::arg("X"), py::arg("T_guess") = 300.0,
      py::arg("tol") = 1.0e-6, py::arg("max_iter") = 50,
      "Solve for temperature T given target molar enthalpy h_target [J/mol] "
      "and composition X.");

  m.def(
      "calc_T_from_h_mass",
      [](double h_mass_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_h_mass(h_mass_target, X, T_guess, tol, max_iter);
      },
      py::arg("h_mass_target"), py::arg("X"), py::arg("T_guess") = 300.0,
      py::arg("tol") = 1.0e-6, py::arg("max_iter") = 50,
      "Solve for temperature T given target mass enthalpy h_mass_target [J/kg] "
      "and composition X.");

  m.def(
      "calc_T_from_s",
      [](double s_target, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_s(s_target, P, X, T_guess, tol, max_iter);
      },
      py::arg("s_target"), py::arg("P"), py::arg("X"),
      py::arg("T_guess") = 300.0, py::arg("tol") = 1.0e-6,
      py::arg("max_iter") = 50,
      "Solve for temperature T given target entropy s_target, pressure P, and "
      "composition X.");

  m.def(
      "calc_T_from_cp",
      [](double cp_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_cp(cp_target, X, T_guess, tol, max_iter);
      },
      py::arg("cp_target"), py::arg("X"), py::arg("T_guess") = 300.0,
      py::arg("tol") = 1.0e-6, py::arg("max_iter") = 50,
      "Solve for temperature T given target Cp cp_target and composition X: the root\n"
      "nearest T_guess in [1, 6000] K (cp is not monotone). Below the 200 K NASA-9\n"
      "data the ideal-gas tail applies; cp - cp(0 K) ~ T^3 there, so a root near\n"
      "1 K is ill-conditioned (~5e-7 relative for CO). The scan steps 10 K, and\n"
      "10 % of T below 100 K, so two roots within one step can be missed;\n"
      "calc_T_from_cp_roots returns every root. RuntimeError when no root lies in\n"
      "the range, ValueError for a non-finite or non-positive input.");

  py::class_<CpRoot>(m, "CpRoot")
      .def_readonly("T", &CpRoot::T)
      .def_readonly("dT_dcp", &CpRoot::dT_dcp)
      .def_readonly("double_root", &CpRoot::double_root)
      .def("__repr__", [](const CpRoot &r) {
        return "CpRoot(T=" + py::repr(py::float_(r.T)).cast<std::string>() +
               ", dT_dcp=" + py::repr(py::float_(r.dT_dcp)).cast<std::string>() +
               ", double_root=" + (r.double_root ? "True" : "False") + ")";
      });

  m.def(
      "calc_T_from_cp_roots",
      [](double cp_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_min, double T_max, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_cp_roots(cp_target, X, T_min, T_max, tol, max_iter);
      },
      py::arg("cp_target"), py::arg("X"), py::arg("T_min") = thermo::CP_INVERSE_T_MIN,
      py::arg("T_max") = thermo::CP_INVERSE_T_MAX, py::arg("tol") = 1.0e-6,
      py::arg("max_iter") = 50,
      "Every root T of cp(T) = cp_target [J/(mol K)] in [T_min, T_max] K, sorted,\n"
      "as a list of CpRoot (T [K], dT_dcp = 1/(dcp/dT) [K/(J/(mol K))],\n"
      "double_root); empty when there is none. The opt-in alternative to\n"
      "calc_T_from_cp, which returns the root nearest T_guess. cp is split at its\n"
      "extrema (sign changes of dcp/dT on calc_T_from_cp's scan nodes, bisected,\n"
      "and a pair inside one step, searched where |dcp/dT| dips between nodes)\n"
      "into monotone pieces, so roots closer together than a scan step are found.\n"
      "A target equal to cp at an extremum (to 1e-13 relative) is one root with\n"
      "double_root=True and dT_dcp NaN. The range must lie in [1, 6000] K.\n"
      "ValueError for a non-finite or non-positive input or an invalid range, and\n"
      "when cp equals the target on an interval (argon: not isolated roots).");

  m.def(
      "calc_T_from_u",
      [](double u_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_u(u_target, X, T_guess, tol, max_iter);
      },
      py::arg("u_target"), py::arg("X"), py::arg("T_guess") = 300.0,
      py::arg("tol") = 1.0e-6, py::arg("max_iter") = 50,
      "Solve for temperature T given target molar internal energy u_target [J/mol] "
      "and composition X.");

  m.def(
      "calc_T_from_s_mass",
      [](double s_mass_target, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_s_mass(s_mass_target, P, X, T_guess, tol, max_iter);
      },
      py::arg("s_mass_target"), py::arg("P"), py::arg("X"),
      py::arg("T_guess") = 300.0, py::arg("tol") = 1.0e-6,
      py::arg("max_iter") = 50,
      "Solve for temperature T given target mass entropy s_mass_target [J/(kg*K)], "
      "pressure P [Pa], and composition X.");

  m.def(
      "calc_T_from_u_mass",
      [](double u_mass_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_u_mass(u_mass_target, X, T_guess, tol, max_iter);
      },
      py::arg("u_mass_target"), py::arg("X"), py::arg("T_guess") = 300.0,
      py::arg("tol") = 1.0e-6, py::arg("max_iter") = 50,
      "Solve for temperature T given target mass internal energy u_mass_target [J/kg] "
      "and composition X.");

  m.def(
      "calc_T_from_sv_mass",
      [](double s_mass_target, double v_mass_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_sv_mass(s_mass_target, v_mass_target, X, T_guess, tol, max_iter);
      },
      py::arg("s_mass_target"), py::arg("v_mass_target"), py::arg("X"),
      py::arg("T_guess") = 300.0, py::arg("tol") = 1.0e-6,
      py::arg("max_iter") = 50,
      "Solve for temperature T given target mass entropy s_mass_target [J/(kg*K)], "
      "mass specific volume v_mass_target [m^3/kg], and composition X.");

  m.def(
      "calc_T_from_sh_mass",
      [](double s_mass_target, double h_mass_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double T_guess, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return calc_T_from_sh_mass(s_mass_target, h_mass_target, X, T_guess, tol, max_iter);
      },
      py::arg("s_mass_target"), py::arg("h_mass_target"), py::arg("X"),
      py::arg("T_guess") = 300.0, py::arg("tol") = 1.0e-6,
      py::arg("max_iter") = 50,
      "Solve for temperature T given target mass entropy s_mass_target [J/(kg*K)], "
      "mass enthalpy h_mass_target [J/kg], and composition X.");

  // Thermodynamic derivatives
  m.def(
      "dh_dT",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return dh_dT(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass enthalpy derivative dh/dT [J/(kg*K)] at temperature T and composition X.");

  m.def(
      "ds_dT",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return ds_dT(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass entropy derivative ds/dT [J/(kg*K^2)] at temperature T and composition X.");

  m.def(
      "dcp_dT",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return dcp_dT(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Mass specific heat derivative dcp/dT [J/(kg*K^2)] at temperature T and composition X.");

  m.def(
      "dg_over_RT_dT",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return dg_over_RT_dT(T, X);
      },
      py::arg("T"), py::arg("X"),
      "Derivative of dimensionless Gibbs energy g/RT with respect to T [1/K] "
      "at temperature T and composition X.");

  // Oxygen requirement helpers (scalar fuel index and mixtures)
  m.def("oxygen_required_per_mol_fuel", &oxygen_required_per_mol_fuel,
        py::arg("fuel_index"),
        "Moles O2 required per mole of pure fuel species (by index).");

  m.def("oxygen_required_per_kg_fuel", &oxygen_required_per_kg_fuel,
        py::arg("fuel_index"),
        "Mass of O2 [kg] required per kg of pure fuel species (by index).");

  m.def(
      "oxygen_required_per_mol_mixture",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return oxygen_required_per_mol_mixture(X);
      },
      py::arg("X"),
      "Moles O2 required per mole of fuel mixture with mole fractions X.");

  m.def(
      "oxygen_required_per_kg_mixture",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        auto X = to_vec(X_arr);
        return oxygen_required_per_kg_mixture(X);
      },
      py::arg("X"),
      "Mass of O2 [kg] required per kg of fuel mixture with mole fractions X.");

  m.def(
      "fuel_lhv_molar",
      [](py::array_t<double, py::array::c_style | py::array::forcecast>
             X_fuel_arr,
         double reference_temperature) {
        auto X_fuel = to_vec(X_fuel_arr);
        return fuel_lhv_molar(X_fuel, reference_temperature);
      },
      py::arg("X_fuel"), py::arg("reference_temperature") = 298.15,
      "Fuel lower heating value (LHV) on molar basis [J/mol fuel mixture].\n\n"
      "Computed from complete combustion to CO2 + H2O(g) at reference "
      "temperature.");

  m.def(
      "fuel_lhv_mass",
      [](py::array_t<double, py::array::c_style | py::array::forcecast>
             X_fuel_arr,
         double reference_temperature) {
        auto X_fuel = to_vec(X_fuel_arr);
        return fuel_lhv_mass(X_fuel, reference_temperature);
      },
      py::arg("X_fuel"), py::arg("reference_temperature") = 298.15,
      "Fuel lower heating value (LHV) on mass basis [J/kg fuel mixture].\n\n"
      "Computed from complete combustion to CO2 + H2O(g) at reference "
      "temperature.");

  // Complete combustion helper
  m.def(
      "complete_combustion_to_CO2_H2O",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         bool smooth) {
        auto X = to_vec(X_arr);
        return complete_combustion_to_CO2_H2O(X, smooth);
      },
      py::arg("X"), py::arg("smooth") = false,
      "Complete combustion of mixture X to CO2 and H2O (returns product mole "
      "fractions).");

  m.def(
      "complete_combustion_to_CO2_H2O_with_fraction",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         bool smooth) {
        auto X = to_vec(X_arr);
        double f = 0.0;
        auto X_out = complete_combustion_to_CO2_H2O(X, f, smooth);
        return py::make_tuple(X_out, f);
      },
      py::arg("X"), py::arg("smooth") = false,
      "Complete combustion of mixture X to CO2 and H2O, returning (product "
      "mole fractions, fuel burn fraction f).");

  m.def(
      "adiabatic_T_wgs",
      [](double T_in,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             X_in_arr,
         double P) {
        auto X_in = to_vec(X_in_arr);
        if (X_in.empty())
          throw std::invalid_argument("X_in must be non-empty");

        State in;
        in.T = T_in;
        in.P = P;
        in.set_X(X_in);

        EquilibriumResult out = wgs_equilibrium_adiabatic(in);
        return out.state.T;
      },
      py::arg("T_in"), py::arg("X_in"), py::arg("P") = 101325.0,
      "Adiabatic flame temperature with WGS limited equilibrium.\n\n"
      "T_in : inlet temperature [K]\n"
      "X_in : inlet mole fractions (1D NumPy array)\n"
      "P    : pressure [Pa] (default 101325)");

  // State struct binding - Pythonic property-based API
  py::class_<State>(m, "State")
      .def(py::init<>())
      // Mutable state variables (read/write properties)
      .def_readwrite("T", &State::T, "Temperature [K]")
      .def_readwrite("P", &State::P, "Pressure [Pa]")
      .def_property(
          "X", [](const State &s) { return s.X; },
          [](State &s,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 X_arr) { s.set_X(to_vec(X_arr)); },
          "Mole fractions [-]")
      .def_property(
          "Y", [](const State &s) { return s.Y; },
          [](State &s,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 Y_arr) { s.set_Y(to_vec(Y_arr)); },
          "Mass fractions [-]")
      // Computed thermodynamic properties (read-only)
      .def_property_readonly("mw", &State::mw, "Molecular weight [g/mol]")
      .def_property_readonly("cp", &State::cp,
                             "Specific heat at constant pressure [J/(kg*K)]")
      .def_property_readonly("h", &State::h, "Specific enthalpy [J/kg]")
      .def_property_readonly("s", &State::s, "Specific entropy [J/(kg*K)]")
      .def_property_readonly("cv", &State::cv,
                             "Specific heat at constant volume [J/(kg*K)]")
      .def_property_readonly("u", &State::u, "Specific internal energy [J/kg]")
      .def_property_readonly("rho", &State::rho, "Density [kg/m^3]")
      .def_property_readonly("R", &State::R, "Specific gas constant [J/(kg*K)]")
      .def_property_readonly("gamma", &State::gamma,
                             "Isentropic expansion coefficient [-]")
      .def_property_readonly("a", &State::a, "Speed of sound [m/s]")
      // Transport properties (read-only)
      .def_property_readonly("mu", &State::mu, "Dynamic viscosity [Pa*s]")
      .def_property_readonly("k", &State::k, "Thermal conductivity [W/(m*K)]")
      .def_property_readonly("nu", &State::nu, "Kinematic viscosity [m^2/s]")
      .def_property_readonly("Pr", &State::Pr, "Prandtl number [-]")
      .def_property_readonly("alpha", &State::alpha,
                             "Thermal diffusivity [m^2/s]")
      // Fluent setters (return self for chaining)
      .def("set_T", &State::set_T, py::arg("T"),
           "Set temperature [K], returns self")
      .def("set_P", &State::set_P, py::arg("P"),
           "Set pressure [Pa], returns self")
      .def(
          "set_X",
          [](State &s,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 X_arr) -> State & { return s.set_X(to_vec(X_arr)); },
          py::arg("X"), "Set mole fractions [-], returns self")
      .def(
          "set_Y",
          [](State &s,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 Y_arr) -> State & { return s.set_Y(to_vec(Y_arr)); },
          py::arg("Y"), "Set mass fractions [-], returns self")
      .def(
          "set_TPX",
          [](State &s, double T, double P,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 X_arr) -> State & { return s.set_TPX(T, P, to_vec(X_arr)); },
          py::arg("T"), py::arg("P"), py::arg("X"),
          "Set temperature, pressure, and mole fractions, returns self")
      .def(
          "set_TPY",
          [](State &s, double T, double P,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 Y_arr) -> State & { return s.set_TPY(T, P, to_vec(Y_arr)); },
          py::arg("T"), py::arg("P"), py::arg("Y"),
          "Set temperature, pressure, and mass fractions, returns self")
      .def("set_TP", &State::set_TP, py::arg("T"), py::arg("P"),
           "Set temperature and pressure, returns self")
      .def("set_DP_mass", &State::set_DP_mass, py::arg("rho"), py::arg("P"),
           "Set density and pressure, returns self")
      .def("set_HP_mass", &State::set_HP_mass, py::arg("h"), py::arg("P"),
           "Set mass enthalpy and pressure, returns self")
      .def("set_SP_mass", &State::set_SP_mass, py::arg("s"), py::arg("P"),
           "Set mass entropy and pressure, returns self")
      .def("set_UV_mass", &State::set_UV_mass, py::arg("u"), py::arg("v"),
           "Set mass internal energy and mass specific volume, returns self")
      .def("set_SV_mass", &State::set_SV_mass, py::arg("s"), py::arg("v"),
           "Set mass entropy and mass specific volume, returns self")
      .def("set_PV_mass", &State::set_PV_mass, py::arg("P"), py::arg("v"),
           "Set pressure and mass specific volume, returns self")
      .def("set_UP_mass", &State::set_UP_mass, py::arg("u"), py::arg("P"),
           "Set mass internal energy and pressure, returns self")
      .def("set_VH_mass", &State::set_VH_mass, py::arg("v"), py::arg("h"),
           "Set mass specific volume and mass enthalpy, returns self")
      .def("set_SH_mass", &State::set_SH_mass, py::arg("s"), py::arg("h"),
           "Set mass entropy and mass enthalpy, returns self")
      // Pythonic tuple properties
      .def_property(
          "TPX", [](const State &s) { return py::make_tuple(s.T, s.P, s.X); },
          [](State &s, py::tuple t) {
            if (t.size() != 3)
              throw std::invalid_argument("TPX must be a 3-tuple (T, P, X)");
            auto X_arr =
                t[2].cast<py::array_t<double, py::array::c_style |
                                                  py::array::forcecast>>();
            s.set_TPX(t[0].cast<double>(), t[1].cast<double>(), to_vec(X_arr));
          },
          "Get/Set Temperature [K], Pressure [Pa], and Mole Fractions as a "
          "tuple (T, P, X)")
      .def_property(
          "TPY",
          [](const State &s) {
            return py::make_tuple(s.T, s.P, mole_to_mass(s.X));
          },
          [](State &s, py::tuple t) {
            if (t.size() != 3)
              throw std::invalid_argument("TPY must be a 3-tuple (T, P, Y)");
            auto Y_arr =
                t[2].cast<py::array_t<double, py::array::c_style |
                                                  py::array::forcecast>>();
            s.set_TPY(t[0].cast<double>(), t[1].cast<double>(), to_vec(Y_arr));
          },
          "Get/Set Temperature [K], Pressure [Pa], and Mass Fractions as a "
          "tuple (T, P, Y)")
      .def_property(
          "TP", [](const State &s) { return py::make_tuple(s.T, s.P); },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("TP must be a 2-tuple (T, P)");
            s.set_TP(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Temperature [K] and Pressure [Pa] as a tuple (T, P)")
      .def_property(
          "DP", [](const State &s) { return py::make_tuple(s.rho(), s.P); },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("DP must be a 2-tuple (rho, P)");
            s.set_DP_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Density [kg/m^3] and Pressure [Pa] as a tuple (rho, P)")
      .def_property(
          "HP", [](const State &s) { return py::make_tuple(s.h(), s.P); },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("HP must be a 2-tuple (h, P)");
            s.set_HP_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Enthalpy [J/kg] and Pressure [Pa] as a tuple (h, P)")
      .def_property(
          "SP", [](const State &s) { return py::make_tuple(s.s(), s.P); },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("SP must be a 2-tuple (s, P)");
            s.set_SP_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Entropy [J/(kg K)] and Pressure [Pa] as a tuple (s, P)")
      .def_property(
          "UV",
          [](const State &s) {
            double v_mass = 1.0 / s.rho();
            return py::make_tuple(s.u(), v_mass);
          },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("UV must be a 2-tuple (u, v)");
            s.set_UV_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Internal Energy [J/kg] and Mass Specific Volume "
          "[m^3/kg] as a tuple (u, v)")
      .def_property(
          "SV",
          [](const State &s) {
            double v_mass = 1.0 / s.rho();
            return py::make_tuple(s.s(), v_mass);
          },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("SV must be a 2-tuple (s, v)");
            s.set_SV_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Entropy [J/(kg K)] and Mass Specific Volume [m^3/kg] as "
          "a tuple (s, v)")
      .def_property(
          "PV",
          [](const State &s) {
            double v_mass = 1.0 / s.rho();
            return py::make_tuple(s.P, v_mass);
          },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("PV must be a 2-tuple (P, v)");
            s.set_PV_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Pressure [Pa] and Mass Specific Volume [m^3/kg] as a tuple "
          "(P, v)")
      .def_property(
          "UP", [](const State &s) { return py::make_tuple(s.u(), s.P); },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("UP must be a 2-tuple (u, P)");
            s.set_UP_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Internal Energy [J/kg] and Pressure [Pa] as a tuple "
          "(u, P)")
      .def_property(
          "VH",
          [](const State &s) {
            double v_mass = 1.0 / s.rho();
            return py::make_tuple(v_mass, s.h());
          },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("VH must be a 2-tuple (v, h)");
            s.set_VH_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Specific Volume [m^3/kg] and Mass Enthalpy [J/kg] as a "
          "tuple (v, h)")
      .def_property(
          "SH", [](const State &s) { return py::make_tuple(s.s(), s.h()); },
          [](State &s, py::tuple t) {
            if (t.size() != 2)
              throw std::invalid_argument("SH must be a 2-tuple (s, h)");
            s.set_SH_mass(t[0].cast<double>(), t[1].cast<double>());
          },
          "Get/Set Mass Entropy [J/(kg K)] and Mass Enthalpy [J/kg] as a tuple "
          "(s, h)");

  // AirProperties is a type alias for TransportState (same C++ type).
  // Expose as a module attribute for backward compatibility.
  m.attr("AirProperties") = m.attr("TransportState");

  m.def("air_properties", &air_properties, py::arg("T"), py::arg("P"),
        py::arg("humidity") = 0.0,
        "Compute all air properties at once.\n\n"
        "Convenience function that computes all thermodynamic and transport\n"
        "properties of air in a single call. Returns a TransportState\n"
        "(also accessible as AirProperties) with read-only attributes.\n\n"
        "Parameters:\n"
        "  T        : temperature [K]\n"
        "  P        : pressure [Pa]\n"
        "  humidity : relative humidity [0-1] (default: 0.0 = dry air)\n\n"
        "Returns: TransportState with attributes:\n"
        "  T, P, rho, mu, k, nu, alpha, Pr, cp, cv, gamma, a\n\n"
        "Example:\n"
        "  >>> props = cb.air_properties(T=300, P=101325, humidity=0.5)\n"
        "  >>> print(props.rho)  # IDE autocomplete works!\n"
        "  1.177\n"
        "  >>> print(props.gamma)\n"
        "  1.4");

  // ThermoState struct binding - Bundle of thermodynamic properties
  py::class_<ThermoState>(
      m, "ThermoState",
      "Bundle of thermodynamic properties for a gas mixture.\n\n"
      "All properties are read-only attributes computed in a single call.\n"
      "Provides IDE autocomplete and type safety.\n\n"
      "Attributes:\n"
      "  T        : Temperature [K] (input, echoed back)\n"
      "  P        : Pressure [Pa] (input, echoed back)\n"
      "  rho      : Density [kg/m^3]\n"
      "  cp       : Specific heat at constant pressure [J/(kg*K)]\n"
      "  cv       : Specific heat at constant volume [J/(kg*K)]\n"
      "  h        : Specific enthalpy [J/kg]\n"
      "  s        : Specific entropy [J/(kg*K)]\n"
      "  u        : Specific internal energy [J/kg]\n"
      "  gamma    : Isentropic expansion coefficient [-]\n"
      "  a        : Speed of sound [m/s]\n"
      "  cp_mole  : Molar cp [J/(mol*K)]\n"
      "  cv_mole  : Molar cv [J/(mol*K)]\n"
      "  h_mole   : Molar enthalpy [J/mol]\n"
      "  s_mole   : Molar entropy [J/(mol*K)]\n"
      "  u_mole   : Molar internal energy [J/mol]\n"
      "  mw       : Molecular weight [g/mol]")
      .def_readonly("T", &ThermoState::T, "Temperature [K]")
      .def_readonly("P", &ThermoState::P, "Pressure [Pa]")
      .def_readonly("rho", &ThermoState::rho, "Density [kg/m^3]")
      .def_readonly("cp", &ThermoState::cp,
                    "Specific heat at constant pressure [J/(kg*K)]")
      .def_readonly("cv", &ThermoState::cv,
                    "Specific heat at constant volume [J/(kg*K)]")
      .def_readonly("h", &ThermoState::h, "Specific enthalpy [J/kg]")
      .def_readonly("s", &ThermoState::s, "Specific entropy [J/(kg*K)]")
      .def_readonly("u", &ThermoState::u, "Specific internal energy [J/kg]")
      .def_readonly("gamma", &ThermoState::gamma,
                    "Isentropic expansion coefficient [-]")
      .def_readonly("a", &ThermoState::a, "Speed of sound [m/s]")
      .def_readonly("cp_mole", &ThermoState::cp_mole, "Molar cp [J/(mol*K)]")
      .def_readonly("cv_mole", &ThermoState::cv_mole, "Molar cv [J/(mol*K)]")
      .def_readonly("h_mole", &ThermoState::h_mole, "Molar enthalpy [J/mol]")
      .def_readonly("s_mole", &ThermoState::s_mole, "Molar entropy [J/(mol*K)]")
      .def_readonly("u_mole", &ThermoState::u_mole,
                    "Molar internal energy [J/mol]")
      .def_readonly("mw", &ThermoState::mw, "Molecular weight [g/mol]")
      .def("__repr__", [](const ThermoState &s) {
        return "<ThermoState: T=" + std::to_string(s.T) +
               " K, "
               "P=" +
               std::to_string(s.P) +
               " Pa, "
               "rho=" +
               std::to_string(s.rho) + " kg/m^3>";
      });

  m.def("thermo_state", &thermo_state, py::arg("T"), py::arg("P"), py::arg("X"),
        py::arg("P_ref") = 101325.0,
        "Compute all thermodynamic properties at once.\n\n"
        "Convenience function that computes all thermodynamic properties\n"
        "for a gas mixture in a single call. Returns ThermoState struct\n"
        "with read-only attributes for IDE autocomplete support.\n\n"
        "Parameters:\n"
        "  T     : temperature [K]\n"
        "  P     : pressure [Pa]\n"
        "  X     : mole fractions [-]\n"
        "  P_ref : reference pressure for entropy [Pa] (default: 101325.0)\n\n"
        "Returns: ThermoState object with attributes:\n"
        "  T, P, rho, cp, cv, h, s, u, gamma, a,\n"
        "  cp_mole, cv_mole, h_mole, s_mole, u_mole, mw\n\n"
        "Example:\n"
        "  >>> X = cb.species.dry_air()\n"
        "  >>> state = cb.thermo_state(T=300, P=101325, X=X)\n"
        "  >>> print(state.h)  # IDE autocomplete works!\n"
        "  8682.5\n"
        "  >>> print(state.gamma)\n"
        "  1.400");

  // Stream struct binding - Pythonic property-based API
  py::class_<Stream>(m, "Stream")
      .def(py::init<>())
      .def_readwrite("state", &Stream::state, "Thermodynamic state")
      // Mutable state (read/write properties)
      .def_property(
          "T", [](const Stream &s) { return s.T(); },
          [](Stream &s, double T) { s.set_T(T); }, "Temperature [K]")
      .def_property(
          "P", [](const Stream &s) { return s.P(); },
          [](Stream &s, double P) { s.set_P(P); }, "Pressure [Pa]")
      .def_property(
          "X", [](const Stream &s) { return s.X(); },
          [](Stream &s,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 X_arr) { s.set_X(to_vec(X_arr)); },
          "Mole fractions [-]")
      .def_readwrite("mdot", &Stream::mdot, "Mass flow rate [kg/s]")
      // Computed properties (read-only)
      .def_property_readonly("mw", &Stream::mw, "Molecular weight [g/mol]")
      .def_property_readonly("cp", &Stream::cp,
                             "Specific heat at constant pressure [J/(kg*K)]")
      .def_property_readonly("h", &Stream::h, "Specific enthalpy [J/kg]")
      .def_property_readonly("s", &Stream::s, "Specific entropy [J/(kg*K)]")
      .def_property_readonly("rho", &Stream::rho, "Density [kg/m^3]")
      // Fluent setters (return self for chaining)
      .def("set_T", &Stream::set_T, py::arg("T"),
           "Set temperature [K], returns self")
      .def("set_P", &Stream::set_P, py::arg("P"),
           "Set pressure [Pa], returns self")
      .def(
          "set_X",
          [](Stream &s,
             py::array_t<double, py::array::c_style | py::array::forcecast>
                 X_arr) -> Stream & { return s.set_X(to_vec(X_arr)); },
          py::arg("X"), "Set mole fractions [-], returns self")
      .def("set_mdot", &Stream::set_mdot, py::arg("mdot"),
           "Set mass flow rate [kg/s], returns self");

  // -------------------------------------------------------------
  // HumidAir State Object
  // -------------------------------------------------------------
  py::class_<HumidAir>(m, "HumidAir")
      .def(py::init<>())
      .def_readwrite("state", &HumidAir::state,
                     "Core CombAero thermodynamic state")
      // Properties
      .def_property(
          "rh", &HumidAir::RH,
          [](HumidAir &h, double rh) { h.set_TP_RH(h.state.T, h.state.P, rh); },
          "Relative humidity [0-1]")
      .def_property_readonly("dewpoint", &HumidAir::dewpoint,
                             "Dewpoint temperature [K]")
      .def_property_readonly("wet_bulb", &HumidAir::wet_bulb,
                             "Wet-bulb temperature [K]")
      .def_property_readonly("humidity_ratio", &HumidAir::humidity_ratio,
                             "Humidity ratio [-]")
      .def_property_readonly(
          "h_mass", &HumidAir::h_mass,
          "ASHRAE psychrometric specific enthalpy [J/kg dry air]")
      // Tuple-like properties
      .def_property(
          "TP_RH", [](const HumidAir &h) { return h.TP_RH(); },
          [](HumidAir &h, std::tuple<double, double, double> tp_rh) {
            h.set_TP_RH(std::get<0>(tp_rh), std::get<1>(tp_rh),
                        std::get<2>(tp_rh));
          },
          "Get/Set T, P, RH tuple")
      .def_property(
          "TP_dewpoint", [](const HumidAir &h) { return h.TP_dewpoint(); },
          [](HumidAir &h, std::tuple<double, double, double> tp_tdp) {
            h.set_TP_dewpoint(std::get<0>(tp_tdp), std::get<1>(tp_tdp),
                              std::get<2>(tp_tdp));
          },
          "Get/Set T, P, dewpoint tuple")
      // Fluent setters (return self for chaining)
      .def("set_TP_RH", &HumidAir::set_TP_RH, py::arg("T"), py::arg("P"),
           py::arg("RH"), "Set T, P, RH, returns self")
      .def("set_TP_dewpoint", &HumidAir::set_TP_dewpoint, py::arg("T"),
           py::arg("P"), py::arg("Tdp"), "Set T, P, dewpoint, returns self")
      .def("set_TP_omega", &HumidAir::set_TP_omega, py::arg("T"), py::arg("P"),
           py::arg("w"), "Set T, P, humidity ratio, returns self");

  // Stream mixing function
  m.def(
      "mix",
      [](const std::vector<Stream> &streams, double P_out) {
        return mix(streams, P_out);
      },
      py::arg("streams"), py::arg("P_out") = -1.0,
      "Mix multiple streams, solving mass and enthalpy balance.\n\n"
      "streams : list of Stream objects\n"
      "P_out   : output pressure [Pa] (default: minimum inlet pressure)\n\n"
      "Returns mixed Stream.");

  // -------------------------------------------------------------
  // Inverse solvers for fuel/oxidizer streams (complete combustion only)
  // -------------------------------------------------------------

  // --- Find fuel stream (oxidizer mdot fixed) ---

  m.def("set_fuel_stream_for_phi", &set_fuel_stream_for_phi, py::arg("phi"),
        py::arg("fuel"), py::arg("oxidizer"),
        "Find fuel stream mass flow rate for target equivalence ratio phi.\n\n"
        "This is a direct algebraic helper (no iteration):\n"
        "  fuel.mdot = phi * fuel_stoich_mdot(oxidizer stream)\n\n"
        "Parameters:\n"
        "  phi      : target equivalence ratio [-]\n"
        "  fuel     : fuel Stream (T, P, X set; mdot ignored)\n"
        "  oxidizer : oxidizer Stream (with mdot set)\n\n"
        "Returns: fuel Stream with mdot set for requested phi.");

  m.def(
      "set_fuel_stream_for_Tad",
      [](double T_ad_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter, bool lean, double phi_max) {
        return set_fuel_stream_for_Tad(T_ad_target, fuel, oxidizer, tol,
                                       max_iter, lean, phi_max);
      },
      py::arg("T_ad_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1.0, py::arg("max_iter") = 100, py::arg("lean") = true,
      py::arg("phi_max") = 10.0,
      "Find fuel mass flow rate to achieve target adiabatic flame "
      "temperature.\n\n"
      "Uses complete combustion to CO2/H2O (no WGS).\n\n"
      "T_ad_target : target adiabatic flame temperature [K] (must be > "
      "oxidizer T)\n"
      "fuel        : fuel Stream (T, P, X set; mdot ignored)\n"
      "oxidizer    : oxidizer Stream (with mdot set)\n"
      "tol         : temperature tolerance [K] (default: 1.0)\n"
      "max_iter    : maximum iterations (default: 100)\n"
      "lean        : if True (default), search on lean side; if False, search "
      "on rich side\n"
      "phi_max     : maximum equivalence ratio for rich side search (default: "
      "10.0)\n\n"
      "Returns fuel Stream with mdot set to achieve target T_ad.");

  m.def(
      "set_fuel_stream_for_O2",
      [](double X_O2_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_fuel_stream_for_O2(X_O2_target, fuel, oxidizer, tol,
                                      max_iter);
      },
      py::arg("X_O2_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find fuel mass flow rate to achieve target O2 mole fraction in burned "
      "products (wet basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.");

  m.def(
      "set_fuel_stream_for_O2_dry",
      [](double X_O2_dry_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_fuel_stream_for_O2_dry(X_O2_dry_target, fuel, oxidizer, tol,
                                          max_iter);
      },
      py::arg("X_O2_dry_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find fuel mass flow rate to achieve target O2 mole fraction in burned "
      "products (dry basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
      "Dry basis: water vapor removed from products before computing mole "
      "fraction.");

  m.def(
      "set_fuel_stream_for_CO2",
      [](double X_CO2_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_fuel_stream_for_CO2(X_CO2_target, fuel, oxidizer, tol,
                                       max_iter);
      },
      py::arg("X_CO2_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find fuel mass flow rate to achieve target CO2 mole fraction in burned "
      "products (wet basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.");

  m.def(
      "set_fuel_stream_for_CO2_dry",
      [](double X_CO2_dry_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_fuel_stream_for_CO2_dry(X_CO2_dry_target, fuel, oxidizer,
                                           tol, max_iter);
      },
      py::arg("X_CO2_dry_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find fuel mass flow rate to achieve target CO2 mole fraction in burned "
      "products (dry basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
      "Dry basis: water vapor removed from products before computing mole "
      "fraction.");

  // --- Find oxidizer stream (fuel mdot fixed) ---

  m.def(
      "set_oxidizer_stream_for_Tad",
      [](double T_ad_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter, bool lean, double phi_max) {
        return set_oxidizer_stream_for_Tad(T_ad_target, fuel, oxidizer, tol,
                                           max_iter, lean, phi_max);
      },
      py::arg("T_ad_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1.0, py::arg("max_iter") = 100, py::arg("lean") = true,
      py::arg("phi_max") = 10.0,
      "Find oxidizer mass flow rate to achieve target adiabatic flame "
      "temperature.\n\n"
      "Uses complete combustion to CO2/H2O (no WGS).\n\n"
      "T_ad_target : target adiabatic flame temperature [K] (must be > fuel "
      "T)\n"
      "fuel        : fuel Stream (with mdot set)\n"
      "oxidizer    : oxidizer Stream (T, P, X set; mdot ignored)\n"
      "tol         : temperature tolerance [K] (default: 1.0)\n"
      "max_iter    : maximum iterations (default: 100)\n"
      "lean        : if True (default), search on lean side; if False, search "
      "on rich side\n"
      "phi_max     : maximum equivalence ratio for search range (default: "
      "10.0)\n\n"
      "Returns oxidizer Stream with mdot set to achieve target T_ad.");

  m.def(
      "set_oxidizer_stream_for_O2",
      [](double X_O2_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_oxidizer_stream_for_O2(X_O2_target, fuel, oxidizer, tol,
                                          max_iter);
      },
      py::arg("X_O2_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find oxidizer mass flow rate to achieve target O2 mole fraction in "
      "burned products (wet basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.");

  m.def(
      "set_oxidizer_stream_for_O2_dry",
      [](double X_O2_dry_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_oxidizer_stream_for_O2_dry(X_O2_dry_target, fuel, oxidizer,
                                              tol, max_iter);
      },
      py::arg("X_O2_dry_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find oxidizer mass flow rate to achieve target O2 mole fraction in "
      "burned products (dry basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
      "Dry basis: water vapor removed from products before computing mole "
      "fraction.");

  m.def(
      "set_oxidizer_stream_for_CO2",
      [](double X_CO2_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_oxidizer_stream_for_CO2(X_CO2_target, fuel, oxidizer, tol,
                                           max_iter);
      },
      py::arg("X_CO2_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find oxidizer mass flow rate to achieve target CO2 mole fraction in "
      "burned products (wet basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.");

  m.def(
      "set_oxidizer_stream_for_CO2_dry",
      [](double X_CO2_dry_target, const Stream &fuel, const Stream &oxidizer,
         double tol, std::size_t max_iter) {
        return set_oxidizer_stream_for_CO2_dry(X_CO2_dry_target, fuel, oxidizer,
                                               tol, max_iter);
      },
      py::arg("X_CO2_dry_target"), py::arg("fuel"), py::arg("oxidizer"),
      py::arg("tol") = 1e-6, py::arg("max_iter") = 100,
      "Find oxidizer mass flow rate to achieve target CO2 mole fraction in "
      "burned products (dry basis).\n\n"
      "Uses complete combustion to CO2/H2O (no WGS). Lean combustion only.\n"
      "Dry basis: water vapor removed from products before computing mole "
      "fraction.");

  // State-based combustion functions
  m.def(
      "complete_combustion",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P, bool smooth_phi0, bool smooth_phi1, double k0, double k1) {
        State in;
        in.T = T;
        in.P = P;
        in.set_X(to_vec(X_arr));
        return complete_combustion(in, smooth_phi0, smooth_phi1, k0, k1);
      },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      py::arg("smooth_phi0") = true, py::arg("smooth_phi1") = true,
      py::arg("k0") = combaero::SMOOTHING_K_PHI0,
      py::arg("k1") = combaero::SMOOTHING_K_PHI1,
      "Adiabatic complete combustion to CO2 and H2O.\n\n"
      "Returns State with adiabatic flame temperature and burned composition.");

  m.def(
      "complete_combustion_isothermal",
      [](double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P, bool smooth_phi0, bool smooth_phi1, double k0, double k1) {
        State in;
        in.T = T;
        in.P = P;
        in.set_X(to_vec(X_arr));
        return complete_combustion_isothermal(in, smooth_phi0, smooth_phi1, k0,
                                              k1);
      },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      py::arg("smooth_phi0") = true, py::arg("smooth_phi1") = true,
      py::arg("k0") = combaero::SMOOTHING_K_PHI0,
      py::arg("k1") = combaero::SMOOTHING_K_PHI1,
      "Isothermal complete combustion to CO2 and H2O.\n\n"
      "Returns State with reactant temperature and burned composition.");

  // Helper: construct State from (T, P, X_arr) -- avoids 4-line boilerplate in
  // every binding
  auto make_state =
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        State s;
        s.T = T;
        s.P = P;
        s.set_X(to_vec(X_arr));
        return s;
      };

  // State-based WGS equilibrium functions
  m.def(
      "wgs_equilibrium",
      [make_state](
          double T,
          py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
          double P) { return wgs_equilibrium(make_state(T, P, X_arr)); },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      "Isothermal WGS equilibrium (CO + H2O <-> CO2 + H2).\n\n"
      "Returns EquilibriumResult with equilibrium composition at input "
      "temperature.");

  m.def(
      "wgs_equilibrium_adiabatic",
      [make_state](
          double T,
          py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
          double P) {
        return wgs_equilibrium_adiabatic(make_state(T, P, X_arr));
      },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      "Adiabatic WGS equilibrium (CO + H2O <-> CO2 + H2).\n\n"
      "Returns EquilibriumResult with equilibrium temperature and "
      "composition.");

  // General reforming + WGS equilibrium (all hydrocarbons)
  m.def(
      "reforming_equilibrium",
      [make_state](
          double T,
          py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
          double P) { return reforming_equilibrium(make_state(T, P, X_arr)); },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      "Isothermal steam reforming + WGS equilibrium for all hydrocarbons.\n\n"
      "General reforming: CnHm + n*H2O <-> n*CO + (n + m/2)*H2\n"
      "Water-Gas Shift: CO + H2O <-> CO2 + H2\n\n"
      "Handles all hydrocarbons: CH4, C2H6, C3H8, iC4H10, nC5H12, etc.\n"
      "Returns EquilibriumResult with equilibrium composition at input "
      "temperature.");

  m.def(
      "reforming_equilibrium_adiabatic",
      [make_state](
          double T,
          py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
          double P) {
        return reforming_equilibrium_adiabatic(make_state(T, P, X_arr));
      },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      "Adiabatic steam reforming + WGS equilibrium for all hydrocarbons.\n\n"
      "General reforming: CnHm + n*H2O <-> n*CO + (n + m/2)*H2\n"
      "Water-Gas Shift: CO + H2O <-> CO2 + H2\n\n"
      "Handles all hydrocarbons: CH4, C2H6, C3H8, iC4H10, nC5H12, etc.\n"
      "Returns EquilibriumResult with equilibrium temperature and "
      "composition.");

  // Analytic Jacobians of combustion (#716).
  py::class_<combaero::CompleteCombustionJacobian>(m, "CompleteCombustionJacobian",
                                                   "Complete combustion products and dX_b/dX")
      .def_readonly("X_b", &combaero::CompleteCombustionJacobian::X_b)
      .def_readonly("dXb_dX", &combaero::CompleteCombustionJacobian::dXb_dX);
  m.def(
      "complete_combustion_products_and_jacobian",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr, bool smooth_phi0,
         bool smooth_phi1, double k0, double k1) {
        return combaero::complete_combustion_products_and_jacobian(to_vec(X_arr), smooth_phi0,
                                                                   smooth_phi1, k0, k1);
      },
      py::arg("X"), py::arg("smooth_phi0") = true, py::arg("smooth_phi1") = true,
      py::arg("k0") = combaero::SMOOTHING_K_PHI0, py::arg("k1") = combaero::SMOOTHING_K_PHI1,
      "Complete combustion to CO2 and H2O with the analytic Jacobian of the\n"
      "product mole fractions in the reactant ones, dXb_dX[j][k].");
  py::class_<combaero::EquilibriumJacobian>(m, "EquilibriumJacobian",
                                            "Combustion equilibrium with analytic derivatives")
      .def_readonly("state", &combaero::EquilibriumJacobian::state)
      .def_readonly("dT_dT_in", &combaero::EquilibriumJacobian::dT_dT_in)
      .def_readonly("dT_dP", &combaero::EquilibriumJacobian::dT_dP)
      .def_readonly("dT_dX", &combaero::EquilibriumJacobian::dT_dX)
      .def_readonly("dXout_dX", &combaero::EquilibriumJacobian::dXout_dX)
      .def_readonly("dXout_dT_in", &combaero::EquilibriumJacobian::dXout_dT_in)
      .def_readonly("dXout_dP", &combaero::EquilibriumJacobian::dXout_dP);
  m.def(
      "combustion_equilibrium_and_jacobian",
      [make_state](double T, py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
                   double P, bool smooth, double k0) {
        return combaero::combustion_equilibrium_and_jacobian(make_state(T, P, X_arr), smooth, k0);
      },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0, py::arg("smooth") = true,
      py::arg("k0") = combaero::SMOOTHING_K_PHI0,
      "combustion_equilibrium with analytic derivatives of T and the product\n"
      "mole fractions in the reactant T, P and X (implicit in the extents).");

  // Combustion + Equilibrium (convenience function)
  m.def(
      "combustion_equilibrium",
      [make_state](
          double T,
          py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
          double P, bool smooth, double k0) {
        return combustion_equilibrium(make_state(T, P, X_arr), smooth, k0);
      },
      py::arg("T"), py::arg("X"), py::arg("P") = 101325.0,
      py::arg("smooth") = true, py::arg("k0") = combaero::SMOOTHING_K_PHI0,
      "Combustion + equilibrium in one step.\n\n"
      "Combines complete combustion with reforming + WGS equilibrium.\n"
      "Use this when starting from an unburned fuel+air mixture.\n\n"
      "Workflow:\n"
      "  1. Complete combustion (fuel + O2 -> CO2 + H2O)\n"
      "  2. Reforming + WGS equilibrium on products\n\n"
      "Returns State with equilibrium temperature and composition.");

  // -------------------------------------------------------------
  // Friction factor correlations
  // -------------------------------------------------------------

  m.def("friction_haaland", &friction_haaland, py::arg("Re"), py::arg("e_D"),
        "Haaland friction factor correlation (1983).\n\n"
        "Explicit approximation to Colebrook-White equation.\n"
        "Accuracy: ~2-3% vs Colebrook-White.\n\n"
        "f = 1 / [-1.8 * log10((e_D/3.7)^1.11 + 6.9/Re)]^2\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-] (must be > 0)\n"
        "  e_D : relative roughness eps/D [-] (must be >= 0)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Valid for turbulent flow (Re > ~2300).\n"
        "For laminar flow, use f = 64/Re.");

  m.def("friction_serghides", &friction_serghides, py::arg("Re"),
        py::arg("e_D"),
        "Serghides friction factor correlation (1984).\n\n"
        "Explicit approximation using Steffensen acceleration.\n"
        "Accuracy: <0.3% vs Colebrook-White (very accurate).\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-] (must be > 0)\n"
        "  e_D : relative roughness eps/D [-] (must be >= 0)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Valid for turbulent flow (Re > ~2300).");

  m.def("friction_colebrook", &friction_colebrook, py::arg("Re"),
        py::arg("e_D"), py::arg("tol") = 1e-10, py::arg("max_iter") = 20,
        "Colebrook-White friction factor equation (1939).\n\n"
        "Implicit equation solved iteratively using Newton-Raphson.\n"
        "The reference standard for turbulent friction factor.\n\n"
        "1/sqrt(f) = -2 * log10(eps/D/3.7 + 2.51/(Re*sqrt(f)))\n\n"
        "Parameters:\n"
        "  Re       : Reynolds number [-] (must be > 0)\n"
        "  e_D      : relative roughness eps/D [-] (must be >= 0)\n"
        "  tol      : convergence tolerance [-] (default: 1e-10)\n"
        "  max_iter : maximum iterations (default: 20)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Uses Haaland correlation as initial guess.");

  m.def("friction_petukhov", &friction_petukhov, py::arg("Re"),
        "Petukhov friction factor correlation (1970).\n\n"
        "For smooth channels only (zero roughness).\n"
        "f = [0.790 * ln(Re) - 1.64]^(-2)\n\n"
        "Parameters:\n"
        "  Re : Reynolds number [-] (must be > 0)\n\n"
        "Returns: Darcy friction factor f [-] (dimensionless)\n\n"
        "Valid for: 3000 < Re < 5x10^6\n"
        "Often used with Gnielinski/Petukhov heat transfer correlations.");

  // -------------------------------------------------------------
  // Vatistas n-vortex model
  // Reference: Vatistas, Kozel, Mih (1991), Exp. Fluids 11, 73-76.
  // -------------------------------------------------------------

  m.def("vatistas_v0_bar", &solver::vatistas_v0_bar, py::arg("r_bar"),
        py::arg("n") = vatistas::n_default,
        "Normalised tangential velocity profile V0_bar(r_bar; n).\n\n"
        "V0_bar = r_bar / (1 + r_bar^{2n})^{1/n}  [Eq. 2]\n\n"
        "Peak at r_bar = 1 with value 2^{-1/n} for all n >= 1.\n\n"
        "Parameters:\n"
        "  r_bar : normalised radius r/r_c [-]\n"
        "  n     : shape parameter [-] (must be >= 1, default 2)\n\n"
        "Returns: V0_bar [-] (dimensionless)");

  m.def("vatistas_dv0_bar_drbar", &solver::vatistas_dv0_bar_drbar,
        py::arg("r_bar"), py::arg("n") = vatistas::n_default,
        "Analytical d(V0_bar)/d(r_bar).\n\n"
        "= (1 - r_bar^{2n}) / (1 + r_bar^{2n})^{1 + 1/n}\n\n"
        "Positive for r_bar < 1, zero at r_bar = 1, negative for r_bar > 1.\n\n"
        "Parameters:\n"
        "  r_bar : normalised radius r/r_c [-]\n"
        "  n     : shape parameter [-] (default 2)\n\n"
        "Returns: d(V0_bar)/d(r_bar) [-]");

  m.def("vatistas_vr_bar", &solver::vatistas_vr_bar, py::arg("r_bar"),
        py::arg("n") = vatistas::n_default,
        "Normalised radial velocity Vr_bar = V_r * r_c / nu  [Eq. 3].\n\n"
        "= -2*(1+n)*r_bar^{2n-1} / (1 + r_bar^{2n})\n\n"
        "Negative by definition (inward flow toward vortex axis).\n\n"
        "Parameters:\n"
        "  r_bar : normalised radius r/r_c [-]\n"
        "  n     : shape parameter [-] (default 2)\n\n"
        "Returns: Vr_bar [-] (dimensionless, negative)");

  m.def("vatistas_dvr_bar_drbar", &solver::vatistas_dvr_bar_drbar,
        py::arg("r_bar"), py::arg("n") = vatistas::n_default,
        "Analytical d(Vr_bar)/d(r_bar).\n\n"
        "Parameters:\n"
        "  r_bar : normalised radius r/r_c [-]\n"
        "  n     : shape parameter [-] (default 2)\n\n"
        "Returns: d(Vr_bar)/d(r_bar) [-]");

  m.def("vatistas_pressure_integral", &solver::vatistas_pressure_integral,
        py::arg("r_bar"), py::arg("n") = vatistas::n_default,
        "Raw pressure integral I(r_bar; n) = integral_0^{r_bar} V0_bar^2/r' dr'.\n\n"
        "Closed form for n=1 and n=2; composite Simpson quadrature otherwise.\n"
        "  n = 1: r_bar^2 / (2*(1 + r_bar^2))\n"
        "  n = 2: arctan(r_bar^2) / 2\n\n"
        "Dimensional delta_P = (Gamma/(2*pi*r_c))^2 * rho * I(r/r_c; n)\n\n"
        "For the normalised form (0 to 1 as r_bar -> inf):\n"
        "  n=2: delta_P_bar = I / (pi/4)  =  (2/pi)*arctan(r_bar^2)  [Eq. 11]\n\n"
        "Parameters:\n"
        "  r_bar : normalised radius r/r_c [-]\n"
        "  n     : shape parameter [-] (default 2)\n\n"
        "Returns: I(r_bar; n) [-]");

  m.def("vatistas_d_pressure_integral_drbar",
        &solver::vatistas_d_pressure_integral_drbar, py::arg("r_bar"),
        py::arg("n") = vatistas::n_default,
        "Analytical d(I)/d(r_bar) = V0_bar(r_bar)^2 / r_bar  [FTC, all n].\n\n"
        "Parameters:\n"
        "  r_bar : normalised radius r/r_c [-]\n"
        "  n     : shape parameter [-] (default 2)\n\n"
        "Returns: d(I)/d(r_bar) [-]");

  m.def("vatistas_v_theta",
        py::overload_cast<double, double, double, double>(
            &solver::vatistas_v_theta),
        py::arg("r"), py::arg("Gamma"), py::arg("r_c"),
        py::arg("n") = vatistas::n_default,
        "Tangential velocity V_theta(r; Gamma, r_c, n) [m/s].\n\n"
        "V_theta = [Gamma/(2*pi*r_c)] * V0_bar(r/r_c; n)\n\n"
        "Parameters:\n"
        "  r     : radius [m] (>= 0)\n"
        "  Gamma : vortex circulation [m^2/s] (must be > 0)\n"
        "  r_c   : core radius [m] (must be > 0)\n"
        "  n     : shape parameter [-] (must be >= 1, default 2)\n\n"
        "Returns: V_theta [m/s]");

  m.def("vatistas_v_theta_and_jacobians",
        py::overload_cast<double, double, double, double>(
            &solver::vatistas_v_theta_and_jacobians),
        py::arg("r"), py::arg("Gamma"), py::arg("r_c"),
        py::arg("n") = vatistas::n_default,
        "V_theta and its partial derivatives (all analytical).\n\n"
        "Parameters:\n"
        "  r, Gamma, r_c, n : same as vatistas_v_theta\n\n"
        "Returns: tuple (V_theta [m/s], dV/dr [1/s],\n"
        "                dV/dGamma [(m/s)/(m^2/s)], dV/dr_c [(m/s)/m])");

  m.def("vatistas_delta_p",
        py::overload_cast<double, double, double, double, double>(
            &solver::vatistas_delta_p),
        py::arg("r"), py::arg("rho"), py::arg("Gamma"), py::arg("r_c"),
        py::arg("n") = vatistas::n_default,
        "Static pressure rise P(r) - P(0) [Pa].\n\n"
        "delta_P = (Gamma/(2*pi*r_c))^2 * rho * I(r/r_c; n)\n\n"
        "Parameters:\n"
        "  r     : radius [m]\n"
        "  rho   : fluid density [kg/m^3]\n"
        "  Gamma : vortex circulation [m^2/s] (must be > 0)\n"
        "  r_c   : core radius [m] (must be > 0)\n"
        "  n     : shape parameter [-] (must be >= 1, default 2)\n\n"
        "Returns: delta_P [Pa]");

  m.def("vatistas_delta_p_and_jacobians",
        py::overload_cast<double, double, double, double, double>(
            &solver::vatistas_delta_p_and_jacobians),
        py::arg("r"), py::arg("rho"), py::arg("Gamma"), py::arg("r_c"),
        py::arg("n") = vatistas::n_default,
        "delta_P and its partial derivatives (all analytical).\n\n"
        "Parameters:\n"
        "  r, rho, Gamma, r_c, n : same as vatistas_delta_p\n\n"
        "Returns: tuple (delta_P [Pa], d_delta_P/dr [Pa/m],\n"
        "                d_delta_P/dGamma [Pa/(m^2/s)],\n"
        "                d_delta_P/dr_c [Pa/m])");

  // -------------------------------------------------------------
  // Heat transfer correlations
  // -------------------------------------------------------------

  m.def(
      "nusselt_dittus_boelter",
      [](double Re, double Pr, bool heating) {
        return nusselt_dittus_boelter(Re, Pr, heating, nullptr);
      },
      py::arg("Re"), py::arg("Pr"), py::arg("heating") = true,
      "Dittus-Boelter Nusselt number correlation (1930).\n\n"
      "Nu = 0.023 * Re^0.8 * Pr^n\n"
      "where n = 0.4 for heating, n = 0.3 for cooling.\n\n"
      "Parameters:\n"
      "  Re      : Reynolds number [-] (validated: Re > 10,000; warns and "
      "extrapolates outside)\n"
      "  Pr      : Prandtl number [-] (validated: 0.6 < Pr < 160; warns and "
      "extrapolates outside)\n"
      "  heating : True if fluid is being heated, False if cooled (default: "
      "True)\n\n"
      "Returns: Nusselt number Nu [-] (dimensionless)\n\n"
      "Valid for fully developed turbulent flow with L/D > 10.\n"
      "Heat transfer coefficient: h = Nu * k / L [W/(m^2*K)]");

  m.def(
      "nusselt_gnielinski",
      [](double Re, double Pr, double f) {
        if (f < 0.0) {
          return nusselt_gnielinski(Re, Pr, nullptr);
        }
        return nusselt_gnielinski(Re, Pr, f, nullptr);
      },
      py::arg("Re"), py::arg("Pr"), py::arg("f") = -1.0,
      "Gnielinski Nusselt number correlation (1976).\n\n"
      "Nu = (f/8) * (Re - 1000) * Pr / (1 + 12.7 * sqrt(f/8) * (Pr^(2/3) - "
      "1))\n\n"
      "More accurate than Dittus-Boelter, especially in transition region.\n"
      "Below Re=2300 uses a C1-smooth Hermite blend to laminar Nu (3.66).\n\n"
      "Parameters:\n"
      "  Re : Reynolds number [-] (validated: 2300-5e6; warns and extrapolates "
      "outside)\n"
      "  Pr : Prandtl number [-] (validated: 0.5-2000; warns and extrapolates "
      "outside)\n"
      "  f  : Darcy friction factor [-] (if < 0, uses Petukhov "
      "auto-friction)\n\n"
      "Returns: Nusselt number Nu [-] (dimensionless)");

  m.def(
      "nusselt_sieder_tate",
      [](double Re, double Pr, double mu_ratio) {
        return nusselt_sieder_tate(Re, Pr, mu_ratio, nullptr);
      },
      py::arg("Re"), py::arg("Pr"), py::arg("mu_ratio") = 1.0,
      "Sieder-Tate Nusselt number correlation (1936).\n\n"
      "Nu = 0.027 * Re^0.8 * Pr^(1/3) * (mu_bulk / mu_wall)^0.14\n\n"
      "Accounts for viscosity variation between bulk and wall temperatures.\n\n"
      "Parameters:\n"
      "  Re       : Reynolds number at bulk temperature [-] (validated: Re > "
      "10,000; warns outside)\n"
      "  Pr       : Prandtl number at bulk temperature [-] (validated: "
      "0.7-16700; warns outside)\n"
      "  mu_ratio : mu_bulk / mu_wall [-] (default: 1.0, typically 0.5-2.0)\n\n"
      "Returns: Nusselt number Nu [-] (dimensionless)\n\n"
      "Valid for L/D > 10.");

  m.def("htc_from_nusselt", &htc_from_nusselt, py::arg("Nu"), py::arg("k"),
        py::arg("L"),
        "Heat transfer coefficient from Nusselt number.\n\n"
        "h = Nu * k / L\n\n"
        "Parameters:\n"
        "  Nu : Nusselt number [-] (dimensionless)\n"
        "  k  : thermal conductivity [W/(m*K)]\n"
        "  L  : characteristic length [m] (diameter for circular channel "
        "flow)\n\n"
        "Returns: heat transfer coefficient h [W/(m^2*K)]");

  m.def(
      "nusselt_circular_channel",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double velocity, double diameter, const std::string &correlation,
         bool heating, double mu_ratio, double roughness) {
        auto X = to_vec(X_arr);
        auto [h, Nu, Re] =
            htc_circular_channel(T, P, X, velocity, diameter, correlation,
                                 heating, mu_ratio, roughness);
        return Nu;
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("velocity"),
      py::arg("diameter"), py::arg("correlation") = "gnielinski",
      py::arg("heating") = true, py::arg("mu_ratio") = 1.0,
      py::arg("roughness") = 0.0,
      "Composite Nusselt number for circular channel flow [-].");

  m.def(
      "htc_circular_channel",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double velocity, double diameter, const std::string &correlation,
         bool heating, double mu_ratio, double roughness) {
        auto X = to_vec(X_arr);
        return htc_circular_channel(T, P, X, velocity, diameter, correlation,
                                    heating, mu_ratio, roughness);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("velocity"),
      py::arg("diameter"), py::arg("correlation") = "gnielinski",
      py::arg("heating") = true, py::arg("mu_ratio") = 1.0,
      py::arg("roughness") = 0.0,
      "Composite heat transfer coefficient for circular channel flow.\n\n"
      "Computes HTC, Nusselt number, and Reynolds number from thermodynamic "
      "state.\n"
      "Automatically calculates density, viscosity, thermal conductivity, and "
      "Prandtl number.\n\n"
      "Parameters:\n"
      "  T           : static temperature [K]\n"
      "  P           : static pressure [Pa]\n"
      "  X           : mole fractions [mol/mol]\n"
      "  velocity    : flow velocity [m/s]\n"
      "  diameter    : channel diameter [m]\n"
      "  correlation : 'gnielinski' (default), 'dittus_boelter', "
      "'sieder_tate', 'petukhov'\n"
      "  heating     : True for heating, False for cooling (affects "
      "Dittus-Boelter)\n"
      "  mu_ratio    : mu_bulk / mu_wall for Sieder-Tate viscosity correction "
      "(default: 1.0)\n"
      "  roughness   : absolute roughness [m] (default: 0.0 = smooth "
      "channel)\n\n"
      "Returns: tuple (h, Nu, Re)\n"
      "  h  : heat transfer coefficient [W/(m^2*K)]\n"
      "  Nu : Nusselt number [-]\n"
      "  Re : Reynolds number [-]\n\n"
      "Correlations:\n"
      "  - gnielinski    : 2300 < Re < 5e6, 0.5 < Pr < 2000 (best "
      "general-purpose)\n"
      "  - dittus_boelter: Re > 10000, 0.6 < Pr < 160 (simple, fast)\n"
      "  - sieder_tate   : Re > 10000, 0.7 < Pr < 16700 (viscosity "
      "correction)\n"
      "  - petukhov      : Re > 10000, 0.5 < Pr < 2000 (high accuracy)\n\n"
      "Automatically handles laminar flow (Re < 2300) with constant Nu.");

  m.def("lmtd", &lmtd, py::arg("dT1"), py::arg("dT2"),
        "Log mean temperature difference (LMTD) for heat exchangers.\n\n"
        "LMTD = (dT1 - dT2) / ln(dT1 / dT2)\n\n"
        "For counter-flow:\n"
        "  dT1 = T_hot_in - T_cold_out\n"
        "  dT2 = T_hot_out - T_cold_in\n\n"
        "For parallel-flow:\n"
        "  dT1 = T_hot_in - T_cold_in\n"
        "  dT2 = T_hot_out - T_cold_out\n\n"
        "Parameters:\n"
        "  dT1 : temperature difference at one end [K]\n"
        "  dT2 : temperature difference at other end [K]\n\n"
        "Returns: log mean temperature difference [K]\n\n"
        "Note: If dT1 ~ dT2, returns arithmetic mean to avoid 0/0.");

  // -------------------------------------------------------------
  // Geometry utilities
  // -------------------------------------------------------------

  m.def("hydraulic_diameter", &hydraulic_diameter, py::arg("A"),
        py::arg("P_wetted"),
        "Hydraulic diameter for arbitrary cross-section.\n\n"
        "Dh = 4 * A / P_wetted\n\n"
        "Parameters:\n"
        "  A        : cross-sectional area [m^2]\n"
        "  P_wetted : wetted perimeter [m]\n\n"
        "Returns: hydraulic diameter Dh [m]");

  m.def("hydraulic_diameter_rect", &hydraulic_diameter_rect, py::arg("a"),
        py::arg("b"),
        "Hydraulic diameter for rectangular duct.\n\n"
        "Dh = 2*a*b / (a + b)\n\n"
        "Parameters:\n"
        "  a : side length [m]\n"
        "  b : side length [m]\n\n"
        "Returns: hydraulic diameter Dh [m]");

  m.def("hydraulic_diameter_annulus", &hydraulic_diameter_annulus,
        py::arg("D_outer"), py::arg("D_inner"),
        "Hydraulic diameter for annular duct (concentric cylinders).\n\n"
        "Dh = D_outer - D_inner\n\n"
        "Parameters:\n"
        "  D_outer : outer diameter [m]\n"
        "  D_inner : inner diameter [m]\n\n"
        "Returns: hydraulic diameter Dh [m]");

  m.def("residence_time", py::overload_cast<double, double>(&residence_time),
        py::arg("V"), py::arg("Q"),
        "Residence time from volumetric flow rate.\n\n"
        "tau = V / Q_dot\n\n"
        "Time for fluid to pass through a volume.\n"
        "Used in reactor design (Damkoehler number), combustor sizing.\n\n"
        "Parameters:\n"
        "  V : volume [m^3]\n"
        "  Q : volumetric flow rate [m^3/s]\n\n"
        "Returns: residence time tau [s]");

  m.def("residence_time_mdot",
        py::overload_cast<double, double, double>(&residence_time_mdot),
        py::arg("V"), py::arg("mdot"), py::arg("rho"),
        "Residence time from mass flow rate.\n\n"
        "tau = V*rho / mdot\n\n"
        "Parameters:\n"
        "  V    : volume [m^3]\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  rho  : density [kg/m^3]\n\n"
        "Returns: residence time tau [s]");

  m.def("channel_mdot", &channel_mdot, py::arg("v"), py::arg("D"),
        py::arg("rho"),
        "Channel mass flow rate from velocity.\n\n"
        "mdot = rho * v * pi * D^2 / 4\n\n"
        "Parameters:\n"
        "  v   : velocity [m/s]\n"
        "  D   : diameter [m]\n"
        "  rho : density [kg/m^3]\n\n"
        "Returns: mass flow rate [kg/s]");

  // Composite channel flow function
  m.def(
      "pressure_drop_channel",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double v, double D, double L, double roughness,
         const std::string &correlation) {
        auto X = to_vec(X_arr);
        return channel_pressure_drop(T, P, X, v, D, L, roughness, correlation);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("v"), py::arg("D"),
      py::arg("L"), py::arg("roughness") = 0.0,
      py::arg("correlation") = "haaland",
      "Composite pressure drop calculation for channel flow.\n\n"
      "Combines thermodynamic properties, Reynolds number, friction factor,\n"
      "and Darcy-Weisbach equation in one call.\n\n"
      "Steps:\n"
      "  1. Calculate density and viscosity from (T, P, X)\n"
      "  2. Calculate Reynolds number: Re = rho * v * D / mu\n"
      "  3. Calculate friction factor using specified correlation\n"
      "  4. Calculate pressure drop: dP = f * (L/D) * (rho * v^2 / 2)\n\n"
      "Parameters:\n"
      "  T           : temperature [K]\n"
      "  P           : pressure [Pa]\n"
      "  X           : mole fractions [mol/mol]\n"
      "  v           : velocity [m/s]\n"
      "  D           : diameter [m]\n"
      "  L           : length [m]\n"
      "  roughness   : absolute roughness [m] (default: 0.0 = smooth)\n"
      "  correlation : friction correlation (default: 'haaland')\n"
      "                Options: 'haaland', 'serghides', 'colebrook', "
      "'petukhov'\n\n"
      "Returns: tuple of (dP [Pa], Re [-], f [-])\n\n"
      "Example:\n"
      "  >>> dP, Re, f = channel_pressure_drop(300, 101325, X_air, 10.0, 0.1, "
      "100.0)\n"
      "  >>> print(f'Pressure drop: {dP:.1f} Pa, Re: {Re:.0f}, f: {f:.4f}')");

  m.def("space_velocity", &space_velocity, py::arg("Q"), py::arg("V"),
        "Space velocity (inverse of residence time).\n\n"
        "SV = Q_dot / V = 1/tau\n\n"
        "Parameters:\n"
        "  Q : volumetric flow rate [m^3/s]\n"
        "  V : volume [m^3]\n\n"
        "Returns: space velocity SV [1/s]");

  m.def("circular_area", &circular_area, py::arg("D"),
        "Circular cross-sectional area.\n\n"
        "A = pi * (D/2)^2\n\n"
        "Parameters:\n"
        "  D : diameter [m]\n\n"
        "Returns: area [m^2]");

  m.def("annular_area", &annular_area, py::arg("D_outer"), py::arg("D_inner"),
        "Annular cross-sectional area.\n\n"
        "A = pi * ((D_outer/2)^2 - (D_inner/2)^2)\n\n"
        "Parameters:\n"
        "  D_outer : outer diameter [m]\n"
        "  D_inner : inner diameter [m]\n\n"
        "Returns: area [m^2]");

  m.def("cylinder_volume", &cylinder_volume, py::arg("D"), py::arg("L"),
        "Cylindrical volume.\n\n"
        "V = pi * (D/2)^2 * L\n\n"
        "Parameters:\n"
        "  D : diameter [m]\n"
        "  L : length [m]\n\n"
        "Returns: volume [m^3]");

  // -------------------------------------------------------------
  // Compressible flow
  // -------------------------------------------------------------

  // Bind structures
  py::class_<IncompressibleStation>(
      m, "IncompressibleStation",
      "Spatial state data point for incompressible channel flow")
      .def(py::init<>())
      .def_readwrite("x", &IncompressibleStation::x, "Distance from inlet [m]")
      .def_readwrite("P", &IncompressibleStation::P, "Static Pressure [Pa]")
      .def_readwrite("T", &IncompressibleStation::T, "Static Temperature [K]")
      .def_readwrite("rho", &IncompressibleStation::rho, "Density [kg/m^3]")
      .def_readwrite("u", &IncompressibleStation::v, "Velocity [m/s]")
      .def_readwrite("M", &IncompressibleStation::M, "Mach Number [-]")
      .def_readwrite("h", &IncompressibleStation::h, "Enthalpy [J/kg]")
      .def("__repr__", [](const IncompressibleStation &s) {
        std::ostringstream oss;
        oss << "<IncompressibleStation x=" << s.x << "m P=" << s.P
            << "Pa T=" << s.T << "K>";
        return oss.str();
      });

  py::class_<IncompressibleFlowSolution>(m, "IncompressibleFlowSolution",
                                         "Result of incompressible flow solver")
      .def(py::init<>())
      .def_readonly("mdot", &IncompressibleFlowSolution::mdot,
                    "Mass flow rate [kg/s]")
      .def_readonly("v", &IncompressibleFlowSolution::v, "Velocity [m/s]")
      .def_readonly("dP", &IncompressibleFlowSolution::dP, "Pressure drop [Pa]")
      .def_readonly("Re", &IncompressibleFlowSolution::Re,
                    "Reynolds number [-]")
      .def_readwrite("rho", &IncompressibleFlowSolution::rho, "Density [kg/m^3]")
      .def_readwrite("f", &IncompressibleFlowSolution::f,
                     "Friction factor / Discharge coefficient [-]")
      .def_readwrite("profile", &IncompressibleFlowSolution::profile,
                     "Optional spatial profile along the channel (props(L))")
      .def("__repr__", [](const IncompressibleFlowSolution &s) {
        std::ostringstream oss;
        oss << "IncompressibleFlowSolution(mdot=" << s.mdot << ", v=" << s.v
            << ", dP=" << s.dP << ", Re=" << s.Re << ")";
        return oss.str();
      });

  // CompressibleFlowSolution struct
  py::class_<CompressibleFlowSolution>(m, "CompressibleFlowSolution")
      .def(py::init<>())
      .def_readonly("stagnation", &CompressibleFlowSolution::stagnation,
                    "Stagnation state")
      .def_readonly("outlet", &CompressibleFlowSolution::outlet,
                    "Outlet static state")
      .def_readonly("v", &CompressibleFlowSolution::v, "Outlet velocity [m/s]")
      .def_readonly("M", &CompressibleFlowSolution::M, "Outlet Mach number [-]")
      .def_readonly("mdot", &CompressibleFlowSolution::mdot,
                    "Mass flow rate [kg/s]")
      .def_readonly("choked", &CompressibleFlowSolution::choked,
                    "True if flow is choked");

  // NozzleStation struct
  py::class_<NozzleStation>(m, "NozzleStation")
      .def(py::init<>())
      .def_readonly("x", &NozzleStation::x, "Axial position [m]")
      .def_readonly("A", &NozzleStation::A, "Area [m^2]")
      .def_readonly("P", &NozzleStation::P, "Static pressure [Pa]")
      .def_readonly("T", &NozzleStation::T, "Static temperature [K]")
      .def_readonly("rho", &NozzleStation::rho, "Density [kg/m^3]")
      .def_readonly("u", &NozzleStation::u, "Velocity [m/s]")
      .def_readonly("M", &NozzleStation::M, "Mach number [-]")
      .def_readonly("h", &NozzleStation::h, "Specific enthalpy [J/kg]");

  // NozzleSolution struct
  py::class_<NozzleSolution>(m, "NozzleSolution")
      .def(py::init<>())
      .def_readonly("inlet", &NozzleSolution::inlet, "Inlet state")
      .def_readonly("outlet", &NozzleSolution::outlet, "Outlet state")
      .def_readonly("mdot", &NozzleSolution::mdot, "Mass flow rate [kg/s]")
      .def_readonly("h0", &NozzleSolution::h0, "Stagnation enthalpy [J/kg]")
      .def_readonly("Tt", &NozzleSolution::Tt, "Stagnation temperature [K]")
      .def_readonly("Pt", &NozzleSolution::Pt, "Stagnation pressure [Pa]")
      .def_readonly("choked", &NozzleSolution::choked,
                    "True if throat is sonic")
      .def_readonly("x_throat", &NozzleSolution::x_throat,
                    "Throat position [m]")
      .def_readonly("A_throat", &NozzleSolution::A_throat, "Throat area [m^2]")
      .def_readonly("profile", &NozzleSolution::profile, "Axial profile");

  // ThrustResult struct
  py::class_<ThrustResult>(m, "ThrustResult")
      .def(py::init<>())
      .def_readwrite("thrust", &ThrustResult::thrust, "Thrust force [N]")
      .def_readwrite("specific_impulse", &ThrustResult::specific_impulse,
                     "Specific impulse Isp [s]")
      .def_readwrite("thrust_coefficient", &ThrustResult::thrust_coefficient,
                     "Thrust coefficient C_F [-]")
      .def_readwrite("mdot", &ThrustResult::mdot, "Mass flow rate [kg/s]")
      .def_readwrite("u_exit", &ThrustResult::u_exit, "Exit velocity [m/s]")
      .def_readwrite("P_exit", &ThrustResult::P_exit, "Exit pressure [Pa]");

  // FannoStation struct
  py::class_<FannoStation>(m, "FannoStation")
      .def(py::init<>())
      .def_readonly("x", &FannoStation::x, "Position [m]")
      .def_readonly("P", &FannoStation::P, "Static pressure [Pa]")
      .def_readonly("T", &FannoStation::T, "Static temperature [K]")
      .def_readonly("rho", &FannoStation::rho, "Density [kg/m^3]")
      .def_readonly("u", &FannoStation::u, "Velocity [m/s]")
      .def_readonly("M", &FannoStation::M, "Mach number [-]")
      .def_readonly("h", &FannoStation::h, "Specific enthalpy [J/kg]")
      .def_readonly("s", &FannoStation::s, "Specific entropy [J/(kg*K)]")
      .def_readonly("f", &FannoStation::f, "Local Darcy friction factor [-]")
      .def_readonly("Re", &FannoStation::Re, "Local Reynolds number [-]");

  // FannoSolution struct
  py::class_<FannoSolution>(m, "FannoSolution")
      .def(py::init<>())
      .def_readonly("inlet", &FannoSolution::inlet, "Inlet state")
      .def_readonly("outlet", &FannoSolution::outlet, "Outlet state")
      .def_readonly("mdot", &FannoSolution::mdot, "Mass flow rate [kg/s]")
      .def_readonly("h0", &FannoSolution::h0, "Stagnation enthalpy [J/kg]")
      .def_readonly("L", &FannoSolution::L, "Channel length [m]")
      .def_readonly("D", &FannoSolution::D, "Channel diameter [m]")
      .def_readonly("f", &FannoSolution::f,
                    "Darcy friction factor (constant-f) [-]")
      .def_readonly("f_avg", &FannoSolution::f_avg,
                    "Average Darcy friction factor [-]")
      .def_readonly("Re_in", &FannoSolution::Re_in, "Inlet Reynolds number [-]")
      .def_readonly("choked", &FannoSolution::choked,
                    "True if flow reached M=1")
      .def_readonly("L_choke", &FannoSolution::L_choke, "Length to choking [m]")
      .def_readonly("profile", &FannoSolution::profile, "Axial profile");

  // Isentropic nozzle flow
  m.def(
      "nozzle_flow",
      [](double Tt, double Pt, double P_back, double A_eff,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return nozzle_flow(Tt, Pt, P_back, A_eff, X, tol, max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("P_back"), py::arg("A_eff"),
      py::arg("X"), py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Isentropic nozzle flow.\n\n"
      "Tt     : stagnation temperature [K]\n"
      "Pt     : stagnation pressure [Pa]\n"
      "P_back : back pressure [Pa]\n"
      "A_eff  : effective flow area [m^2]\n"
      "X      : mole fractions\n\n"
      "Returns CompressibleFlowSolution.");

  // Inverse solvers
  m.def(
      "solve_A_eff_from_mdot",
      [](double Tt, double Pt, double P_back, double mdot_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return solve_A_eff_from_mdot(Tt, Pt, P_back, mdot_target, X, tol,
                                     max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("P_back"), py::arg("mdot_target"),
      py::arg("X"), py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Find effective area for given mass flow rate.");

  m.def(
      "solve_P_back_from_mdot",
      [](double Tt, double Pt, double A_eff, double mdot_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return solve_P_back_from_mdot(Tt, Pt, A_eff, mdot_target, X, tol,
                                      max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("A_eff"), py::arg("mdot_target"),
      py::arg("X"), py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Find back pressure for given mass flow rate.");

  m.def(
      "solve_Pt_from_mdot",
      [](double Tt, double P_back, double A_eff, double mdot_target,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return solve_Pt_from_mdot(Tt, P_back, A_eff, mdot_target, X, tol,
                                  max_iter);
      },
      py::arg("Tt"), py::arg("P_back"), py::arg("A_eff"),
      py::arg("mdot_target"), py::arg("X"), py::arg("tol") = 1e-8,
      py::arg("max_iter") = 50,
      "Find stagnation pressure for given mass flow rate.");

  // Utility functions
  m.def(
      "critical_pressure_ratio",
      [](double Tt, double Pt,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return critical_pressure_ratio(Tt, Pt, X, tol, max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("X"), py::arg("tol") = 1e-8,
      py::arg("max_iter") = 50, "Critical (sonic) pressure ratio P*/Pt.");

  m.def(
      "mach_from_pressure_ratio",
      [](double Tt, double Pt, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return mach_from_pressure_ratio(Tt, Pt, P, X, tol, max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("P"), py::arg("X"),
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Mach number from pressure ratio P/Pt (isentropic).");

  m.def(
      "mass_flux_isentropic",
      [](double Tt, double Pt, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return mass_flux_isentropic(Tt, Pt, P, X, tol, max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("P"), py::arg("X"),
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Mass flux G = rho*v [kg/(m^2*s)] at given pressure (isentropic).");

  // Quasi-1D nozzle with C-D geometry
  m.def(
      "nozzle_cd",
      [](double Tt, double Pt, double P_exit, double A_inlet, double A_throat,
         double A_exit, double x_throat, double x_exit,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         std::size_t n_stations, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return nozzle_cd(Tt, Pt, P_exit, A_inlet, A_throat, A_exit, x_throat,
                         x_exit, X, n_stations, tol, max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("P_exit"), py::arg("A_inlet"),
      py::arg("A_throat"), py::arg("A_exit"), py::arg("x_throat"),
      py::arg("x_exit"), py::arg("X"), py::arg("n_stations") = 100,
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Converging-diverging nozzle flow.\n\n"
      "Uses smooth cosine area profile between inlet, throat, and exit.\n\n"
      "Returns NozzleSolution with axial profile.");

  // Fanno flow
  m.def(
      "fanno_channel",
      [](double T_in, double P_in, double u_in, double L, double D, double f,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         std::size_t n_steps, bool store_profile) {
        auto X = to_vec(X_arr);
        return fanno_channel(T_in, P_in, u_in, L, D, f, X, n_steps,
                             store_profile);
      },
      py::arg("T_in"), py::arg("P_in"), py::arg("u_in"), py::arg("L"),
      py::arg("D"), py::arg("f"), py::arg("X"), py::arg("n_steps") = 100,
      py::arg("store_profile") = false,
      "Fanno flow (adiabatic channel with friction).\n\n"
      "T_in : inlet static temperature [K]\n"
      "P_in : inlet static pressure [Pa]\n"
      "u_in : inlet velocity [m/s]\n"
      "L    : channel length [m]\n"
      "D    : channel diameter [m]\n"
      "f    : Darcy friction factor [-]\n\n"
      "Returns FannoSolution.");

  // -------------------------------------------------------------
  // Heat transfer correlations
  // -------------------------------------------------------------

  m.def(
      "nusselt_petukhov",
      [](double Re, double Pr, double f) {
        if (f < 0.0) {
          return nusselt_petukhov(Re, Pr, nullptr);
        }
        return nusselt_petukhov(Re, Pr, f, nullptr);
      },
      py::arg("Re"), py::arg("Pr"), py::arg("f") = -1.0,
      "Petukhov correlation for turbulent channel flow.\n\n"
      "Warns and extrapolates outside validated range (1e4 < Re < 5e6, 0.5 < "
      "Pr < 2000).\n\n"
      "f : friction factor (if < 0, uses Petukhov correlation)");

  // Warning handler and validity utility
  m.def(
      "set_warning_handler",
      [](py::object handler) {
        if (handler.is_none()) {
          combaero::set_warning_handler({});
        } else {
          combaero::set_warning_handler(PyWarningHandler{std::move(handler)});
        }
      },
      py::arg("handler"),
      "Set the global warning handler for out-of-range correlation "
      "extrapolations.\n\n"
      "handler : callable(str) or None\n"
      "  If None, resets to the default (prints to stderr).\n"
      "  Pass a no-op lambda to suppress all warnings.\n\n"
      "Example:\n"
      "  import combaero\n"
      "  combaero.set_warning_handler(lambda msg: None)  # suppress\n"
      "  combaero.set_warning_handler(None)              # restore default");

  m.def(
      "get_warning_handler",
      []() -> py::object {
        // Used by suppress_warnings() to save and restore the handler.
        auto h = combaero::get_warning_handler();
        // If this module installed the current handler, hand back the very
        // object that was passed in. Wrapping it again would make every
        // save-and-restore cycle add a layer -- see PyWarningHandler.
        if (const auto *py_handler = h.target<PyWarningHandler>()) {
          return py_handler->fn;
        }
        // A C++-side handler (including the default) still has to be
        // wrapped; there is no Python object to return. That wrapper is not
        // a PyWarningHandler, so re-installing it takes the branch above's
        // negative path exactly once and then stops growing.
        return py::cpp_function([h](const std::string &msg) { h(msg); });
      },
      "Return the current warning handler as a Python callable.\n\n"
      "Used internally by suppress_warnings() to save and restore the "
      "handler.");

  m.def(
      "is_well_behaved",
      [](double v, double lo, double hi) {
        return combaero::is_well_behaved(v, lo, hi);
      },
      py::arg("v"), py::arg("lo") = 1e-12, py::arg("hi") = 1e15,
      "Check whether a correlation result is safe for use in a Jacobian-based "
      "solver.\n\n"
      "Returns True if v is finite, v > lo, and v < hi.\n\n"
      "v  : scalar result to check\n"
      "lo : lower bound (default 1e-12)\n"
      "hi : upper bound (default 1e15)");

  m.def(
      "overall_htc",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> h_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             tk_arr) {
        auto h_values = to_vec(h_arr);
        auto t_over_k = to_vec(tk_arr);
        return overall_htc(h_values, t_over_k);
      },
      py::arg("h_values"), py::arg("t_over_k") = py::array_t<double>(),
      "Overall HTC from thermal resistance network.\n\n"
      "1/U = sum(1/h_i) + sum(t_j/k_j)\n\n"
      "h_values : convective HTCs [W/(m^2*K)]\n"
      "t_over_k : thickness/conductivity ratios [m^2*K/W]");

  m.def("overall_htc_wall",
        py::overload_cast<double, double, double, double, double>(&overall_htc_wall),
        py::arg("h_inner"), py::arg("h_outer"), py::arg("t_wall"),
        py::arg("k_wall"), py::arg("R_fouling") = 0.0,
        "Overall HTC for single-layer wall.\n\n"
        "1/U = 1/h_inner + t_wall/k_wall + 1/h_outer");

  m.def(
      "overall_htc_wall_multilayer",
      [](double h_inner, double h_outer,
         py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr,
         double R_fouling) {
        auto layers = to_vec(tk_arr);
        if (R_fouling >= 0) {
          return overall_htc_wall(h_inner, h_outer, layers, R_fouling);
        }
        return overall_htc_wall(h_inner, h_outer, layers);
      },
      py::arg("h_inner"), py::arg("h_outer"), py::arg("t_over_k_layers"),
      py::arg("R_fouling") = -1.0,
      "Overall HTC for multi-layer wall.\n\n"
      "1/U = 1/h_inner + sum(t_i/k_i) + 1/h_outer [+ R_fouling]\n\n"
      "Example: steel + insulation\n"
      "  overall_htc_wall_multilayer(500, 10, [0.003/50, 0.05/0.04])");

  m.def("thermal_resistance", &thermal_resistance, py::arg("h"), py::arg("A"),
        "Thermal resistance for convection.\n\n"
        "R = 1/(h*A) [K/W]");

  m.def("thermal_resistance_wall", &thermal_resistance_wall,
        py::arg("thickness"), py::arg("k"), py::arg("A"),
        "Thermal resistance for conduction through wall.\n\n"
        "R = t/(k*A) [K/W]");

  m.def("lmtd_counterflow", &lmtd_counterflow, py::arg("T_hot_in"),
        py::arg("T_hot_out"), py::arg("T_cold_in"), py::arg("T_cold_out"),
        "LMTD for counter-flow heat exchanger.\n\n"
        "dT1 = T_hot_in - T_cold_out\n"
        "dT2 = T_hot_out - T_cold_in");

  m.def("lmtd_parallelflow", &lmtd_parallelflow, py::arg("T_hot_in"),
        py::arg("T_hot_out"), py::arg("T_cold_in"), py::arg("T_cold_out"),
        "LMTD for parallel-flow heat exchanger.\n\n"
        "dT1 = T_hot_in - T_cold_in\n"
        "dT2 = T_hot_out - T_cold_out");

  m.def("heat_rate", &heat_rate, py::arg("U"), py::arg("A"), py::arg("dT"),
        "Heat transfer rate Q = U * A * dT [W]");

  m.def("heat_flux", &heat_flux, py::arg("U"), py::arg("dT"),
        "Heat flux q = U * dT [W/m^2]");

  m.def("heat_transfer_area", &heat_transfer_area, py::arg("Q"), py::arg("U"),
        py::arg("dT"), "Required area A = Q / (U * dT) [m^2]");

  m.def("heat_transfer_dT", &heat_transfer_dT, py::arg("Q"), py::arg("U"),
        py::arg("A"), "Temperature difference dT = Q / (U * A) [K]");

  m.def(
      "wall_temperature_profile",
      [](double T_hot, double T_cold, double h_hot, double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr,
         double R_fouling) {
        auto t_over_k = to_vec(tk_arr);
        double q = 0.0;
        auto temps = wall_temperature_profile(T_hot, T_cold, h_hot, h_cold,
                                              t_over_k, R_fouling, q);
        return py::make_tuple(temps, q);
      },
      py::arg("T_hot"), py::arg("T_cold"), py::arg("h_hot"), py::arg("h_cold"),
      py::arg("t_over_k"), py::arg("R_fouling") = 0.0,
      "Temperature profile through multi-layer wall.\n\n"
      "Returns (temps, q) where:\n"
      "  temps : interface temperatures [K]\n"
      "  q     : heat flux [W/m^2]");

  m.def("ntu", &ntu, py::arg("U"), py::arg("A"), py::arg("C_min"),
        "Number of Transfer Units: NTU = U * A / C_min");

  m.def("capacity_ratio", &capacity_ratio, py::arg("C_min"), py::arg("C_max"),
        "Heat capacity ratio: C_r = C_min / C_max");

  m.def("effectiveness_counterflow", &effectiveness_counterflow, py::arg("NTU"),
        py::arg("C_r"),
        "Effectiveness for counter-flow heat exchanger.\n\n"
        "eps = (1 - exp(-NTU*(1-C_r))) / (1 - C_r*exp(-NTU*(1-C_r)))");

  m.def("effectiveness_parallelflow", &effectiveness_parallelflow,
        py::arg("NTU"), py::arg("C_r"),
        "Effectiveness for parallel-flow heat exchanger.\n\n"
        "eps = (1 - exp(-NTU*(1+C_r))) / (1 + C_r)");

  m.def("heat_rate_from_effectiveness", &heat_rate_from_effectiveness,
        py::arg("epsilon"), py::arg("C_min"), py::arg("T_hot_in"),
        py::arg("T_cold_in"),
        "Heat rate from effectiveness: Q = eps * C_min * (T_hot_in - T_cold_in) "
        "[W]");

  m.def(
      "heat_flux_from_T_at_edge",
      [](double T_measured, std::size_t edge_idx, double T_hot, double T_cold,
         double h_hot, double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             tk_arr) {
        auto t_over_k = to_vec(tk_arr);
        return heat_flux_from_T_at_edge(T_measured, edge_idx, T_hot, T_cold,
                                        h_hot, h_cold, t_over_k);
      },
      py::arg("T_measured"), py::arg("edge_idx"), py::arg("T_hot"),
      py::arg("T_cold"), py::arg("h_hot"), py::arg("h_cold"),
      py::arg("t_over_k"),
      "Heat flux from temperature measured at wall edge.\n\n"
      "Edge indexing: 0 = hot surface, N = cold surface.\n"
      "Returns heat flux q [W/m^2].");

  m.def(
      "heat_flux_from_T_at_depth",
      [](double T_measured, double depth_from_hot, double T_hot, double T_cold,
         double h_hot, double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast> t_arr,
         py::array_t<double, py::array::c_style | py::array::forcecast> k_arr) {
        auto thicknesses = to_vec(t_arr);
        auto conductivities = to_vec(k_arr);
        return heat_flux_from_T_at_depth(T_measured, depth_from_hot, T_hot,
                                         T_cold, h_hot, h_cold, thicknesses,
                                         conductivities);
      },
      py::arg("T_measured"), py::arg("depth_from_hot"), py::arg("T_hot"),
      py::arg("T_cold"), py::arg("h_hot"), py::arg("h_cold"),
      py::arg("thicknesses"), py::arg("conductivities"),
      "Heat flux from temperature measured at depth from hot surface.\n\n"
      "Useful when thermocouple is embedded slightly into wall.\n"
      "Returns heat flux q [W/m^2].");

  m.def(
      "bulk_T_from_edge_T_and_q",
      [](double T_measured, std::size_t edge_idx, double q, double h_hot,
         double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast> tk_arr,
         const std::string &solve_for) {
        auto t_over_k = to_vec(tk_arr);
        return bulk_T_from_edge_T_and_q(T_measured, edge_idx, q, h_hot, h_cold,
                                        t_over_k, solve_for);
      },
      py::arg("T_measured"), py::arg("edge_idx"), py::arg("q"),
      py::arg("h_hot"), py::arg("h_cold"), py::arg("t_over_k"),
      py::arg("solve_for"),
      "Infer bulk fluid temperature from edge temperature and known heat "
      "flux.\n\n"
      "solve_for: 'hot' or 'cold'\n"
      "Returns inferred bulk temperature [K].");

  m.def(
      "dT_edge_dT_hot",
      [](std::size_t edge_idx, double h_hot, double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             tk_arr) {
        auto t_over_k = to_vec(tk_arr);
        return dT_edge_dT_hot(edge_idx, h_hot, h_cold, t_over_k);
      },
      py::arg("edge_idx"), py::arg("h_hot"), py::arg("h_cold"),
      py::arg("t_over_k"),
      "Sensitivity of edge temperature to hot-side bulk temperature.\n\n"
      "Returns dT_edge/dT_hot [-] (0 to 1).");

  m.def(
      "dT_edge_dT_cold",
      [](std::size_t edge_idx, double h_hot, double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             tk_arr) {
        auto t_over_k = to_vec(tk_arr);
        return dT_edge_dT_cold(edge_idx, h_hot, h_cold, t_over_k);
      },
      py::arg("edge_idx"), py::arg("h_hot"), py::arg("h_cold"),
      py::arg("t_over_k"),
      "Sensitivity of edge temperature to cold-side bulk temperature.\n\n"
      "Returns dT_edge/dT_cold [-] (0 to 1).");

  m.def(
      "dT_edge_dT_bulk",
      [](std::size_t edge_idx, double h_hot, double h_cold,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             tk_arr) {
        auto t_over_k = to_vec(tk_arr);
        auto [dT_hot, dT_cold] =
            dT_edge_dT_bulk(edge_idx, h_hot, h_cold, t_over_k);
        return py::make_tuple(dT_hot, dT_cold);
      },
      py::arg("edge_idx"), py::arg("h_hot"), py::arg("h_cold"),
      py::arg("t_over_k"),
      "Both sensitivities at once.\n\n"
      "Returns (dT_edge/dT_hot, dT_edge/dT_cold).");

  m.def(
      "dT_edge_dq",
      [](std::size_t edge_idx, double h_hot,
         py::array_t<double, py::array::c_style | py::array::forcecast>
             tk_arr) {
        auto t_over_k = to_vec(tk_arr);
        return dT_edge_dq(edge_idx, h_hot, t_over_k);
      },
      py::arg("edge_idx"), py::arg("h_hot"), py::arg("t_over_k"),
      "Sensitivity of edge temperature to heat flux.\n\n"
      "Returns dT_edge/dq [K*m^2/W] (negative).");

  // =========================================================================
  // Combustion Method enum + Combustion State
  // =========================================================================

  py::enum_<CombustionMethod>(
      m, "CombustionMethod",
      "Product model for combustion_state / combustion_state_from_streams.\n\n"
      "  Complete    : complete combustion to CO2 + H2O (fast, default)\n"
      "  Equilibrium : PARTIAL equilibrium -- complete combustion, then\n"
      "                reforming + water-gas shift among CO, H2, CO2, H2O\n"
      "                and the fuel. No dissociation (no O, OH, H, NO):\n"
      "                identical to Complete when lean, cooler when rich.\n"
      "                By design: a combustor rarely reaches full\n"
      "                equilibrium, least of all for NO.")
      .value("Complete", CombustionMethod::Complete)
      .value("Equilibrium", CombustionMethod::Equilibrium)
      .export_values();

  py::class_<CombustionState>(
      m, "CombustionState",
      "Bundle of combustion properties with nested CompleteState objects.\n\n"
      "Contains reactant and product states with all thermodynamic and "
      "transport properties.\n\n"
      "Attributes:\n"
      "  phi                : equivalence ratio [-]\n"
      "  fuel_name          : optional fuel label (string)\n"
      "  method             : CombustionMethod used for products\n"
      "  reactants          : CompleteState at reactant conditions\n"
      "  products           : CompleteState at product conditions (T_ad)\n"
      "  mixture_fraction   : Bilger mixture fraction [-]\n"
      "  fuel_burn_fraction : fraction of fuel burned [0-1]")
      .def_readonly("phi", &CombustionState::phi, "Equivalence ratio [-]")
      .def_readonly("fuel_name", &CombustionState::fuel_name,
                    "Fuel label (e.g., 'CH4', 'Natural Gas')")
      .def_readonly("method", &CombustionState::method,
                    "CombustionMethod used for products")
      .def_readonly("reactants", &CombustionState::reactants,
                    "Reactant CompleteState (all thermo + transport)")
      .def_readonly("products", &CombustionState::products,
                    "Product CompleteState at T_ad (all thermo + transport)")
      .def_readonly("mixture_fraction", &CombustionState::mixture_fraction,
                    "Bilger mixture fraction [-]")
      .def_readonly("fuel_burn_fraction", &CombustionState::fuel_burn_fraction,
                    "Fraction of fuel burned [0-1]")
      .def("__repr__", [](const CombustionState &s) {
        std::string method_str = (s.method == CombustionMethod::Equilibrium)
                                     ? "Equilibrium"
                                     : "Complete";
        return "<CombustionState: phi=" + std::to_string(s.phi) +
               ", method=" + method_str +
               ", T_reactants=" + std::to_string(s.reactants.thermo.T) + " K" +
               ", T_products=" + std::to_string(s.products.thermo.T) + " K" +
               (s.fuel_name.empty() ? "" : ", fuel='" + s.fuel_name + "'") +
               ">";
      });

  using CombustionStateBaseFn = CombustionState (*)(
      const std::vector<double> &, const std::vector<double> &, double, double,
      double, const std::string &, CombustionMethod, bool, bool, double,
      double);
  py::class_<combaero::CombustorLossResult>(m, "CombustorLossResult")
      .def_readonly("xi", &combaero::CombustorLossResult::xi)
      .def_readonly("dxi_dtheta", &combaero::CombustorLossResult::dxi_dtheta)
      .def_readonly("dxi_dmdot", &combaero::CombustorLossResult::dxi_dmdot)
      .def_readonly("dxi_dTt_in", &combaero::CombustorLossResult::dxi_dTt_in)
      .def_readonly("dxi_dPt_in", &combaero::CombustorLossResult::dxi_dPt_in)
      .def_readonly("dxi_dY", &combaero::CombustorLossResult::dxi_dY);
  m.def("combustor_fraction_loss", &combaero::combustor_fraction_loss,
        py::arg("theta"), py::arg("mdot"), py::arg("k"), py::arg("xi0"),
        py::arg("m_ref"),
        "Fractional combustor total-pressure loss xi = (k theta + xi0)\n"
        "g(|mdot|/m_ref), g = 1 - (1 - r^2)^3 below the design flow (m_ref\n"
        "<= 0: no fade), with its partials in theta and mdot.");
  m.def("combustor_head_loss", &combaero::combustor_head_loss, py::arg("theta"),
        py::arg("mdot"), py::arg("Tt_in"), py::arg("Pt_in"), py::arg("Y_in"),
        py::arg("k"), py::arg("zeta0"), py::arg("area"),
        "Combustor head loss as a fraction of the inlet total pressure Pt_in:\n"
        "xi = (k theta + zeta0) 0.5 rho v^2 / Pt_in, v = mdot/(rho area),\n"
        "ideal-gas rho at (Tt_in, Pt_in, Y_in) -- Lefebvre's total-total\n"
        "low-Mach pressure-loss factor -- with its partials in theta, mdot,\n"
        "Tt_in, Pt_in and each inlet mass fraction Y_k (the others held fixed).");
  // The composition point of PressureLossElement's loss models (#725).
  m.attr("PRESSURE_LOSS_INPUT_THETA") = combaero::PRESSURE_LOSS_INPUT_THETA;
  py::enum_<combaero::PressureLossKind>(
      m, "PressureLossKind", "Loss model kind of pressure_loss_compose.")
      .value("Fraction", combaero::PressureLossKind::Fraction)
      .value("Head", combaero::PressureLossKind::Head);
  py::class_<combaero::PressureLossModel>(
      m, "PressureLossModel", "Resolved loss model of pressure_loss_compose.")
      .def(py::init([](combaero::PressureLossKind kind, double k, double base,
                       double m_ref, double area) {
             combaero::PressureLossModel mdl;
             mdl.kind = kind;
             mdl.k = k;
             mdl.base = base;
             mdl.m_ref = m_ref;
             mdl.area = area;
             return mdl;
           }),
           py::arg("kind") = combaero::PressureLossKind::Fraction,
           py::arg("k") = 0.0, py::arg("base") = 0.0, py::arg("m_ref") = 0.0,
           py::arg("area") = 0.0)
      .def_readwrite("kind", &combaero::PressureLossModel::kind)
      .def_readwrite("k", &combaero::PressureLossModel::k)
      .def_readwrite("base", &combaero::PressureLossModel::base)
      .def_readwrite("m_ref", &combaero::PressureLossModel::m_ref)
      .def_readwrite("area", &combaero::PressureLossModel::area);
  py::class_<combaero::PressureLossInputs>(
      m, "PressureLossInputs", "State a loss model reads (pressure_loss_compose).")
      .def(py::init([](double mdot, double Tt_in, double Pt_in, std::vector<double> Y_in,
                       std::vector<double> X_in, std::vector<double> x) {
             combaero::PressureLossInputs in;
             in.mdot = mdot;
             in.Tt_in = Tt_in;
             in.Pt_in = Pt_in;
             in.Y_in = std::move(Y_in);
             in.X_in = std::move(X_in);
             in.x = std::move(x);
             return in;
           }),
           py::arg("mdot") = 0.0, py::arg("Tt_in") = 0.0, py::arg("Pt_in") = 0.0,
           py::arg("Y_in") = std::vector<double>{},
           py::arg("X_in") = std::vector<double>{},
           py::arg("x") = std::vector<double>{})
      .def_readwrite("mdot", &combaero::PressureLossInputs::mdot)
      .def_readwrite("Tt_in", &combaero::PressureLossInputs::Tt_in)
      .def_readwrite("Pt_in", &combaero::PressureLossInputs::Pt_in)
      .def_readwrite("Y_in", &combaero::PressureLossInputs::Y_in)
      .def_readwrite("X_in", &combaero::PressureLossInputs::X_in)
      .def_readwrite("x", &combaero::PressureLossInputs::x);
  py::class_<combaero::PressureLossResult>(
      m, "PressureLossResult", "Composed loss fraction xi and its partials.")
      .def_readonly("xi", &combaero::PressureLossResult::xi)
      .def_readonly("dxi_dmdot", &combaero::PressureLossResult::dxi_dmdot)
      .def_readonly("dxi_dTt_in", &combaero::PressureLossResult::dxi_dTt_in)
      .def_readonly("dxi_dPt_in", &combaero::PressureLossResult::dxi_dPt_in)
      .def_readonly("dxi_dY", &combaero::PressureLossResult::dxi_dY)
      .def_readonly("dxi_dx", &combaero::PressureLossResult::dxi_dx)
      .def_property_readonly(
          "dxi_dtheta",
          [](const combaero::PressureLossResult &r) {
            return r.dxi_dx.size() > combaero::PRESSURE_LOSS_INPUT_THETA
                       ? r.dxi_dx[combaero::PRESSURE_LOSS_INPUT_THETA]
                       : 0.0;
          },
          "dxi_dx[PRESSURE_LOSS_INPUT_THETA], 0 for an empty x.");
  m.def("pressure_loss_compose", &combaero::pressure_loss_compose, py::arg("model"),
        py::arg("inputs"),
        "PressureLossElement's loss model in one call: xi (Pt_out = Pt_in (1 -\n"
        "xi)) and its analytic partials in mdot, Tt_in, Pt_in, each inlet mass\n"
        "fraction and each generic model input x (x[0] = theta). Fraction:\n"
        "combustor_fraction_loss; Head: combustor_head_loss, its value at the\n"
        "state's mole fractions X_in when given. Internal: the composition point\n"
        "of PressureLossElement's built-in correlations.");
  m.def("combustion_state",
        static_cast<CombustionStateBaseFn>(&combustion_state),
        py::arg("X_fuel"), py::arg("X_ox"), py::arg("phi"),
        py::arg("T_reactants"), py::arg("P"), py::arg("fuel_name") = "",
        py::arg("method") = CombustionMethod::Complete,
        py::arg("smooth_phi0") = false, py::arg("smooth_phi1") = false,
        py::arg("k0") = combaero::SMOOTHING_K_PHI0,
        py::arg("k1") = combaero::SMOOTHING_K_PHI1,
        "Compute combustion state from equivalence ratio.\n\n"
        "Parameters:\n"
        "  X_fuel       : fuel composition (mole fractions) [-]\n"
        "  X_ox         : oxidizer composition (mole fractions) [-]\n"
        "  phi          : equivalence ratio [-]\n"
        "  T_reactants  : reactant temperature [K]\n"
        "  P            : pressure [Pa]\n"
        "  fuel_name    : optional fuel label (default: '')\n"
        "  method       : CombustionMethod.Complete (default) or "
        "CombustionMethod.Equilibrium\n\n"
        "Returns: CombustionState with reactants, products, phi, "
        "mixture_fraction\n\n"
        "Example:\n"
        "  >>> state = cb.combustion_state(X_CH4, X_air, phi=1.0, "
        "T_reactants=300, P=101325)\n"
        "  >>> state_eq = cb.combustion_state(X_CH4, X_air, phi=1.2, "
        "T_reactants=300, P=101325,\n"
        "  ...                               "
        "method=cb.CombustionMethod.Equilibrium)");

  using CombustionStreamsBaseFn =
      CombustionState (*)(const Stream &, const Stream &, const std::string &,
                          CombustionMethod, bool, bool, double, double);
  m.def("combustion_state_from_streams",
        static_cast<CombustionStreamsBaseFn>(&combustion_state_from_streams),
        py::arg("fuel_stream"), py::arg("ox_stream"), py::arg("fuel_name") = "",
        py::arg("method") = CombustionMethod::Complete,
        py::arg("smooth_phi0") = false, py::arg("smooth_phi1") = false,
        py::arg("k0") = combaero::SMOOTHING_K_PHI0,
        py::arg("k1") = combaero::SMOOTHING_K_PHI1,
        "Compute combustion state from measured streams.\n\n"
        "Phi is COMPUTED from mass flow rates (output, not input).\n\n"
        "Parameters:\n"
        "  fuel_stream : fuel stream with mdot, T, X\n"
        "  ox_stream   : oxidizer stream with mdot, T, X\n"
        "  fuel_name   : optional fuel label (default: '')\n"
        "  method      : CombustionMethod.Complete (default) or "
        "CombustionMethod.Equilibrium\n\n"
        "Returns: CombustionState with phi computed from flow rates");

  // PressureLossContext -- read-only view passed to user pressure-loss callables
  py::class_<PressureLossContext>(
      m, "PressureLossContext",
      "Context passed to user pressure-loss correlations.\n"
      "All fields are loop-free (no dependence on unknown outlet pressure).\n"
      "Return value of the callable: fractional dP/P_in [-].")
      .def_readonly("phi", &PressureLossContext::phi, "Equivalence ratio [-]")
      .def_readonly("T_ad", &PressureLossContext::T_ad,
                    "Adiabatic flame temperature [K]")
      .def_readonly("theta", &PressureLossContext::theta,
                    "Dimensionless temperature rise T_ad/T_in - 1 [-]")
      .def_readonly("mdot_fuel", &PressureLossContext::mdot_fuel,
                    "Fuel mass flow [kg/s] (0 for phi-based calls)")
      .def_readonly("mdot_air", &PressureLossContext::mdot_air,
                    "Air mass flow [kg/s] (0 for phi-based calls)")
      .def_readonly("mdot_total", &PressureLossContext::mdot_total,
                    "Total inlet mass flow [kg/s]")
      .def_property_readonly(
          "Y_products",
          [](const PressureLossContext &c) { return c.Y_products; },
          "Burned gas mass fractions [-]")
      .def_property_readonly(
          "T_in", [](const PressureLossContext &c) { return c.state_in.T; },
          "Inlet temperature [K]")
      .def_property_readonly(
          "P_in", [](const PressureLossContext &c) { return c.state_in.P; },
          "Inlet pressure [Pa]")
      .def_property_readonly(
          "X_in", [](const PressureLossContext &c) { return c.state_in.X; },
          "Inlet mole fractions [-]")
      .def_property_readonly(
          "X_products",
          [](const PressureLossContext &c) { return c.X_products; },
          "Burned gas mole fractions [-]");

  // combustion_state with pressure-loss hook
  // Use explicit function pointer cast to resolve overload ambiguity.
  using CombustionStateHookFn = CombustionState (*)(
      const std::vector<double> &, const std::vector<double> &, double, double,
      double, const std::string &, CombustionMethod,
      const PressureLossCorrelation &, bool, bool, double, double);
  m.def(
      "combustion_state",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_ox,
         double phi, double T_reactants, double P, const std::string &fuel_name,
         CombustionMethod method, py::function pressure_loss, bool smooth_phi0,
         bool smooth_phi1, double k0, double k1) {
        PressureLossCorrelation fn =
            [pressure_loss](
                const PressureLossContext &ctx) -> std::tuple<double, double> {
          py::object res = pressure_loss(ctx);
          try {
            return res.cast<std::tuple<double, double>>();
          } catch (...) {
            return {res.cast<double>(), 0.0};
          }
        };
        CombustionStateHookFn hook_fn = &combustion_state;
        return hook_fn(to_vec(X_fuel), to_vec(X_ox), phi, T_reactants, P,
                       fuel_name, method, fn, smooth_phi0, smooth_phi1, k0, k1);
      },
      py::arg("X_fuel"), py::arg("X_ox"), py::arg("phi"),
      py::arg("T_reactants"), py::arg("P"), py::arg("fuel_name") = "",
      py::arg("method") = CombustionMethod::Complete, py::arg("pressure_loss"),
      py::arg("smooth_phi0") = false, py::arg("smooth_phi1") = false,
      py::arg("k0") = combaero::SMOOTHING_K_PHI0,
      py::arg("k1") = combaero::SMOOTHING_K_PHI1,
      "combustion_state with user-supplied pressure-loss callable.\n\n"
      "pressure_loss(ctx: PressureLossContext) -> float\n"
      "  Returns fractional dP/P_in [-]. All ctx fields are loop-free.\n\n"
      "Example:\n"
      "  result = cb.combustion_state(X_fuel, X_ox, phi=0.8, T_reactants=700,\n"
      "      P=500000, pressure_loss=lambda c: 0.03 + 0.005*c.theta)");

  // combustion_state_from_streams with pressure-loss hook
  using CombustionStreamsHookFn = CombustionState (*)(
      const Stream &, const Stream &, const std::string &, CombustionMethod,
      const PressureLossCorrelation &, bool, bool, double, double);
  m.def(
      "combustion_state_from_streams",
      [](const Stream &fuel_stream, const Stream &ox_stream,
         const std::string &fuel_name, CombustionMethod method,
         py::function pressure_loss, bool smooth_phi0, bool smooth_phi1,
         double k0, double k1) {
        PressureLossCorrelation fn =
            [pressure_loss](
                const PressureLossContext &ctx) -> std::tuple<double, double> {
          py::object res = pressure_loss(ctx);
          try {
            return res.cast<std::tuple<double, double>>();
          } catch (...) {
            return {res.cast<double>(), 0.0};
          }
        };
        CombustionStreamsHookFn hook_fn = &combustion_state_from_streams;
        return hook_fn(fuel_stream, ox_stream, fuel_name, method, fn,
                       smooth_phi0, smooth_phi1, k0, k1);
      },
      py::arg("fuel_stream"), py::arg("ox_stream"), py::arg("fuel_name") = "",
      py::arg("method") = CombustionMethod::Complete, py::arg("pressure_loss"),
      py::arg("smooth_phi0") = false, py::arg("smooth_phi1") = false,
      py::arg("k0") = combaero::SMOOTHING_K_PHI0,
      py::arg("k1") = combaero::SMOOTHING_K_PHI1,
      "combustion_state_from_streams with user-supplied pressure-loss "
      "callable.\n\n"
      "pressure_loss(ctx: PressureLossContext) -> float\n"
      "  Returns fractional dP/P_in [-].\n"
      "  ctx.mdot_fuel and ctx.mdot_air are populated from stream mass flows.");

  m.def("linear_pressure_loss", [](double k, double c) {
    return py::cpp_function(
        [k, c](const PressureLossContext &ctx) -> std::tuple<double, double> {
          return {k * ctx.theta + c, k};
        });
  });

  // =========================================================================
  // Acoustics & Geometry
  // =========================================================================

  py::module geom = m.def_submodule(
      "geometry", "Geometric primitives representing flow areas and volumes.");

  // Tube struct
  py::class_<Tube>(geom, "Tube",
                   "Cylindrical tube geometry for acoustic analysis")
      .def(py::init<double, double>(), py::arg("L"), py::arg("D"),
           "Create tube with length L [m] and diameter D [m]")
      .def_readwrite("L", &Tube::L, "Length [m]")
      .def_readwrite("D", &Tube::D, "Diameter [m]")
      .def("area", &Tube::area, "Cross-sectional area [m^2]")
      .def("volume", &Tube::volume, "Volume [m^3]")
      .def("perimeter", &Tube::perimeter, "Circumference [m]");

  py::enum_<OrificeTap>(
      m, "OrificeTap",
      "Pressure taps of a metering orifice (ISO 5167-2:2003, 5.3.2.1). A\n"
      "measurement choice: it fixes which pressure difference C refers to,\n"
      "not what the plate does to the flow.")
      .value("Corner", OrificeTap::Corner, "L1 = L2' = 0")
      .value("D_D2", OrificeTap::D_D2, "D and D/2 taps: L1 = 1, L2' = 0.47")
      .value("Flange", OrificeTap::Flange, "L1 = L2' = 25.4 mm / D (the default)");

  // MeteringCdCorrelation enum (early binding for default arguments)
  py::enum_<MeteringCdCorrelation>(
      m, "MeteringCdCorrelation",
      "Cd correlations for a metering orifice in a pipe -- every member is a\n"
      "function of beta = d/D. For a hole in a wall that dumps coolant out of\n"
      "its circuit, use DischargeCdCorrelation instead.")
      .value("ReaderHarrisGallagher", MeteringCdCorrelation::ReaderHarrisGallagher,
             "ISO 5167-2 / ASME MFC-3M")
      .value("Constant", MeteringCdCorrelation::Constant,
             "Fixed Cd, 0.61 (orifice::defaults::metering_cd)")
      .value("UserFunction", MeteringCdCorrelation::UserFunction,
             "Caller-supplied Cd function (C++ only)");

  // DischargeCdCorrelation enum
  //
  // Separate from MeteringCdCorrelation because the inputs differ: no beta,
  // and a crossflow ratio OrificeState cannot express. See orifice.h.
  py::enum_<DischargeCdCorrelation>(
      m, "DischargeCdCorrelation",
      "Cd correlations for a discharge hole -- bleed, film, effusion. Inputs\n"
      "are L/d, r/d and the inlet-side crossflow ratio, never beta = d/D.")
      .value("McGreehanSchotsch1988", DischargeCdCorrelation::McGreehanSchotsch1988,
             "McGreehan & Schotsch (1988), the composite chain of Eqs. (8)-(17)")
      .value("Idelchik1966Sharp", DischargeCdCorrelation::Idelchik1966Sharp,
             "Idelchik (1966) diagram 4-17: sharp-edged hole in a large wall,\n"
             "tabulated from Re = 25 to 1e6")
      .value("Idelchik1966Thick", DischargeCdCorrelation::Idelchik1966Thick,
             "Idelchik (1966) diagram 4-18a: thick-walled (deep) hole,\n"
             "zeta'(l/d) plus a Darcy bore-friction term")
      .value("Idelchik1966Beveled", DischargeCdCorrelation::Idelchik1966Beveled,
             "Idelchik (1966) diagram 4-18b: beveled edges, 40-60 deg")
      .value("Idelchik1966Rounded", DischargeCdCorrelation::Idelchik1966Rounded,
             "Idelchik (1966) diagram 4-18c: rounded edges")
      .value("Lichtarowicz1965", DischargeCdCorrelation::Lichtarowicz1965,
             "Lichtarowicz, Duggins & Markland (1965) long orifice, Eqs. (7)\n"
             "and (12). Valid l/d 2-10 and Re 10 to 2e4 -- the low-Reynolds\n"
             "regime a cooling hole actually runs in, where\n"
             "McGreehanSchotsch1988 is floored at Re = 1e4. Refuses below\n"
             "l/d = 1.5, which the source says to avoid for hysteresis.")
      .value("Constant", DischargeCdCorrelation::Constant,
             "Fixed Cd, 0.60 (orifice::defaults::discharge_cd). Pass an\n"
             "explicit value to discharge_cd instead to pin a measured one.");

  // Annulus struct
  py::class_<Annulus>(geom, "Annulus",
                      "Annular duct geometry for acoustic analysis")
      .def(py::init<>())
      .def(py::init<double, double, double>(), py::arg("L"), py::arg("D_inner"),
           py::arg("D_outer"),
           "Create annulus with length L [m], inner diameter D_inner [m], "
           "outer diameter D_outer [m]")
      .def_readwrite("L", &Annulus::L, "Length [m]")
      .def_readwrite("D_inner", &Annulus::D_inner, "Inner diameter [m]")
      .def_readwrite("D_outer", &Annulus::D_outer, "Outer diameter [m]")
      .def_readwrite("n_azimuthal_max", &Annulus::n_azimuthal_max,
                     "Maximum azimuthal mode number to search")
      .def("D_mean", &Annulus::D_mean, "Mean diameter [m]")
      .def("gap", &Annulus::gap, "Annular gap [m]")
      .def("area", &Annulus::area, "Cross-sectional area [m^2]")
      .def("volume", &Annulus::volume, "Volume [m^3]")
      .def("circumference", &Annulus::circumference, "Mean circumference [m]")
      // Backward-compat properties for AnnularDuctGeometry API (settable radius
      // accessors)
      .def_property(
          "length", [](const Annulus &g) { return g.L; },
          [](Annulus &g, double v) { g.L = v; },
          "Axial length [m] (alias for L)")
      .def_property(
          "radius_inner", [](const Annulus &g) { return g.radius_inner(); },
          [](Annulus &g, double v) { g.D_inner = v * 2.0; }, "Inner radius [m]")
      .def_property(
          "radius_outer", [](const Annulus &g) { return g.radius_outer(); },
          [](Annulus &g, double v) { g.D_outer = v * 2.0; },
          "Outer radius [m]");

  py::class_<CanAnnularFlowGeometry>(
      geom, "CanAnnularFlowGeometry",
      "Parameterized can-annular flow geometry with primary, transition, and "
      "annular sections")
      .def(py::init<double, double, double, double, double, double>(),
           py::arg("L"), py::arg("D_inner"), py::arg("D_outer"),
           py::arg("L_primary"), py::arg("D_primary"), py::arg("L_transition"),
           "Create geometry with annular section and dedicated can transition "
           "parameters.")
      .def_readwrite("L", &CanAnnularFlowGeometry::L,
                     "Annular section length [m]")
      .def_readwrite("D_inner", &CanAnnularFlowGeometry::D_inner,
                     "Annular inner diameter [m]")
      .def_readwrite("D_outer", &CanAnnularFlowGeometry::D_outer,
                     "Annular outer diameter [m]")
      .def_readwrite("L_primary", &CanAnnularFlowGeometry::L_primary,
                     "Primary zone length [m]")
      .def_readwrite("D_primary", &CanAnnularFlowGeometry::D_primary,
                     "Primary zone diameter [m]")
      .def_readwrite("L_transition", &CanAnnularFlowGeometry::L_transition,
                     "Transition length [m]")
      .def("D_mean", &CanAnnularFlowGeometry::D_mean,
           "Annular mean diameter [m]")
      .def("gap", &CanAnnularFlowGeometry::gap, "Annular gap [m]")
      .def("area", &CanAnnularFlowGeometry::area,
           "Annular cross-sectional area [m^2]")
      .def("volume", &CanAnnularFlowGeometry::volume,
           "Annular section volume [m^3]")
      .def("circumference", &CanAnnularFlowGeometry::circumference,
           "Annular mean circumference [m]")
      .def("area_primary", &CanAnnularFlowGeometry::area_primary,
           "Primary circular area [m^2]")
      .def("volume_primary", &CanAnnularFlowGeometry::volume_primary,
           "Primary zone volume [m^3]")
      .def("volume_transition", &CanAnnularFlowGeometry::volume_transition,
           "Transition volume [m^3]")
      .def("volume_total", &CanAnnularFlowGeometry::volume_total,
           "Total geometry volume [m^3]")
      .def("length_total", &CanAnnularFlowGeometry::length_total,
           "Total geometry length [m]");

  // BoundaryCondition enum
  py::enum_<BoundaryCondition>(m, "BoundaryCondition",
                               "Acoustic boundary condition")
      .value("Open", BoundaryCondition::Open, "Pressure release (p'=0)")
      .value("Closed", BoundaryCondition::Closed, "Rigid wall (u'=0)");

  // AcousticMode struct
  py::class_<AcousticMode>(m, "AcousticMode",
                           "Acoustic mode with frequency and indices")
      .def_readonly("n_axial", &AcousticMode::n_axial,
                    "Longitudinal mode number")
      .def_readonly("n_azimuthal", &AcousticMode::n_azimuthal,
                    "Azimuthal mode number")
      .def_readonly("frequency", &AcousticMode::frequency,
                    "Natural frequency [Hz]")
      .def("label", &AcousticMode::label,
           "Mode label (e.g., '1L', '2T', '1L1T')")
      .def("__repr__", [](const AcousticMode &m) {
        return "<AcousticMode " + m.label() + " @ " +
               std::to_string(m.frequency) + " Hz>";
      });

  // Tube modes
  m.def("tube_axial_modes", &tube_axial_modes, py::arg("tube"), py::arg("c"),
        py::arg("upstream"), py::arg("downstream"), py::arg("n_max") = 5,
        "Compute axial acoustic modes of a tube.\n\n"
        "tube      : Tube geometry\n"
        "c         : Speed of sound [m/s]\n"
        "upstream  : Boundary condition at x=0\n"
        "downstream: Boundary condition at x=L\n"
        "n_max     : Maximum mode number (default: 5)\n\n"
        "Returns list of AcousticMode sorted by frequency.");

  // Annulus modes
  m.def("annulus_axial_modes", &annulus_axial_modes, py::arg("annulus"),
        py::arg("c"), py::arg("upstream"), py::arg("downstream"),
        py::arg("n_max") = 5, "Compute axial acoustic modes of an annulus.");

  m.def("annulus_azimuthal_modes", &annulus_azimuthal_modes, py::arg("annulus"),
        py::arg("c"), py::arg("m_max") = 4,
        "Compute azimuthal (tangential) acoustic modes of an annulus.\n\n"
        "f_m = m * c / (pi * D_mean)");

  m.def("annulus_modes", &annulus_modes, py::arg("annulus"), py::arg("c"),
        py::arg("upstream"), py::arg("downstream"), py::arg("n_max") = 5,
        py::arg("m_max") = 4,
        "Compute all acoustic modes (axial, azimuthal, combined) of an "
        "annulus.\n\n"
        "Returns list of AcousticMode sorted by frequency.");

  // Utility functions
  m.def("modes_in_range", &modes_in_range, py::arg("modes"), py::arg("f_min"),
        py::arg("f_max"),
        "Filter modes within a frequency range [f_min, f_max] Hz.");

  m.def(
      "closest_mode",
      [](const std::vector<AcousticMode> &modes,
         double f_target) -> py::object {
        const auto *m = closest_mode(modes, f_target);
        if (m)
          return py::cast(*m);
        return py::none();
      },
      py::arg("modes"), py::arg("f_target"),
      "Find the mode closest to a target frequency.\n\n"
      "Returns AcousticMode or None if modes is empty.");

  m.def("min_mode_separation", &min_mode_separation, py::arg("modes"),
        "Minimum frequency separation between adjacent modes [Hz].");

  // Mean flow correction
  m.def("axial_mode_upstream", &axial_mode_upstream, py::arg("f0"),
        py::arg("M"),
        "Upstream-propagating wave frequency with mean flow.\n\n"
        "f+ = f0 / (1 - M)\n\n"
        "f0 : quiescent frequency [Hz]\n"
        "M  : Mach number [-]");

  m.def("axial_mode_downstream", &axial_mode_downstream, py::arg("f0"),
        py::arg("M"),
        "Downstream-propagating wave frequency with mean flow.\n\n"
        "f- = f0 / (1 + M)");

  m.def("axial_mode_split", &axial_mode_split, py::arg("f0"), py::arg("M"),
        "Frequency split due to mean flow.\n\n"
        "Returns (f_upstream, f_downstream).");

  // Helmholtz resonator
  m.def("helmholtz_frequency", &helmholtz_frequency, py::arg("V"),
        py::arg("A_neck"), py::arg("L_neck"), py::arg("c"),
        py::arg("end_correction") = 0.85,
        "Helmholtz resonator frequency.\n\n"
        "f = (c / (2*pi)) * sqrt(A / (V * L_eff))\n\n"
        "V              : cavity volume [m^3]\n"
        "A_neck         : neck area [m^2]\n"
        "L_neck         : neck length [m]\n"
        "c              : speed of sound [m/s]\n"
        "end_correction : end correction factor (default: 0.85 flanged)");

  // Strouhal
  m.def("strouhal", &strouhal, py::arg("f"), py::arg("L"), py::arg("u"),
        "Strouhal number: St = f * L / u");

  m.def("frequency_from_strouhal", &frequency_from_strouhal, py::arg("St"),
        py::arg("L"), py::arg("u"),
        "Frequency from Strouhal number: f = St * u / L");

  // Convenience
  m.def("quarter_wave_frequency", &quarter_wave_frequency, py::arg("L"),
        py::arg("c"), "Quarter-wave resonator fundamental: f = c / (4L)");

  m.def("half_wave_frequency", &half_wave_frequency, py::arg("L"), py::arg("c"),
        "Half-wave resonator fundamental: f = c / (2L)");

  // Viscothermal boundary layers
  m.def("stokes_layer", &stokes_layer, py::arg("nu"), py::arg("f"),
        "Stokes (viscous) boundary layer thickness: delta_nu = sqrt(2*nu/omega)\n\n"
        "nu : kinematic viscosity [m^2/s]\n"
        "f  : frequency [Hz]");

  m.def("thermal_layer", &thermal_layer, py::arg("alpha"), py::arg("f"),
        "Thermal penetration depth: delta_kappa = sqrt(2*alpha/omega)\n\n"
        "alpha : thermal diffusivity [m^2/s]\n"
        "f     : frequency [Hz]");

  m.def("effective_viscothermal_layer", &effective_viscothermal_layer,
        py::arg("delta_nu"), py::arg("delta_kappa"), py::arg("gamma"),
        "Effective viscothermal layer: delta_eff = delta_nu + (gamma-1)*delta_kappa");

  // Quality factor screening
  m.def("helmholtz_Q", &helmholtz_Q, py::arg("V"), py::arg("A_neck"),
        py::arg("L_neck"), py::arg("nu"), py::arg("alpha"), py::arg("gamma"),
        py::arg("f"),
        "Helmholtz resonator Q factor (viscothermal losses).\n\n"
        "Screening estimate - expect +/-50% accuracy.\n"
        "Real designs require experimental validation.");

  m.def("tube_Q", &tube_Q, py::arg("L"), py::arg("D"), py::arg("nu"),
        py::arg("alpha"), py::arg("gamma"), py::arg("f"),
        "Quarter/half-wave tube Q factor (viscothermal losses).\n\n"
        "Screening estimate - expect +/-50% accuracy.");

  m.def("damping_ratio", &damping_ratio, py::arg("Q"),
        "Damping ratio from quality factor: zeta = 1/(2Q)");

  m.def("bandwidth", &bandwidth, py::arg("f0"), py::arg("Q"),
        "Half-power bandwidth: df = f0/Q [Hz]");

  // AcousticProperties struct binding - Bundle of acoustic properties
  py::class_<AcousticProperties>(
      m, "AcousticProperties",
      "Bundle of acoustic properties computed from (f, rho, c, p_rms).\n\n"
      "All properties are read-only attributes computed in a single call.\n"
      "Provides IDE autocomplete and type safety.\n\n"
      "Attributes:\n"
      "  wavelength        : Acoustic wavelength [m]\n"
      "  frequency         : Frequency [Hz]\n"
      "  impedance         : Characteristic acoustic impedance [Pa*s/m]\n"
      "  particle_velocity : Particle velocity amplitude [m/s]\n"
      "  spl               : Sound pressure level [dB]")
      .def_readonly("wavelength", &AcousticProperties::wavelength,
                    "Acoustic wavelength [m]")
      .def_readonly("frequency", &AcousticProperties::frequency,
                    "Frequency [Hz]")
      .def_readonly("impedance", &AcousticProperties::impedance,
                    "Characteristic acoustic impedance [Pa*s/m]")
      .def_readonly("particle_velocity", &AcousticProperties::particle_velocity,
                    "Particle velocity amplitude [m/s]")
      .def_readonly("spl", &AcousticProperties::spl,
                    "Sound pressure level [dB]")
      .def("__repr__", [](const AcousticProperties &p) {
        return "<AcousticProperties: f=" + std::to_string(p.frequency) +
               " Hz, "
               "lambda=" +
               std::to_string(p.wavelength) +
               " m, "
               "SPL=" +
               std::to_string(p.spl) + " dB>";
      });

  m.def(
      "acoustic_properties", &acoustic_properties, py::arg("f"), py::arg("rho"),
      py::arg("c"), py::arg("p_rms") = 20e-6, py::arg("p_ref") = 20e-6,
      "Compute all acoustic properties at once.\n\n"
      "Convenience function that computes all acoustic properties in a\n"
      "single call. Returns AcousticProperties struct with read-only "
      "attributes\n"
      "for IDE autocomplete support.\n\n"
      "Parameters:\n"
      "  f     : frequency [Hz]\n"
      "  rho   : density [kg/m^3]\n"
      "  c     : speed of sound [m/s]\n"
      "  p_rms : RMS pressure amplitude [Pa] (default: 20e-6 Pa = 0 dB "
      "reference)\n"
      "  p_ref : reference pressure for SPL [Pa] (default: 20e-6 Pa in air)\n\n"
      "Returns: AcousticProperties object with attributes:\n"
      "  wavelength, frequency, impedance, particle_velocity, spl\n\n"
      "Example:\n"
      "  >>> props = cb.acoustic_properties(f=1000, rho=1.2, c=340, "
      "p_rms=1.0)\n"
      "  >>> print(f'Wavelength: {props.wavelength:.3f} m')\n"
      "  >>> print(f'Impedance: {props.impedance:.1f} Pa*s/m')\n"
      "  >>> print(f'SPL: {props.spl:.1f} dB')");

  // =========================================================================
  // Transfer Matrix Method - Advanced Thermoacoustics
  // =========================================================================

  py::class_<TransferMatrix>(m, "TransferMatrix")
      .def(py::init<>())
      .def_readwrite("T11", &TransferMatrix::T11)
      .def_readwrite("T12", &TransferMatrix::T12)
      .def_readwrite("T21", &TransferMatrix::T21)
      .def_readwrite("T22", &TransferMatrix::T22)
      .def("__mul__", &TransferMatrix::operator*);

  py::class_<AcousticMedium>(m, "AcousticMedium")
      .def(py::init<>())
      .def_readwrite("rho", &AcousticMedium::rho)
      .def_readwrite("c", &AcousticMedium::c);

  py::class_<LinerOrificeGeometry>(m, "LinerOrificeGeometry")
      .def(py::init<>())
      .def_readwrite("d_orifice", &LinerOrificeGeometry::d_orifice)
      .def_readwrite("l_orifice", &LinerOrificeGeometry::l_orifice)
      .def_readwrite("porosity", &LinerOrificeGeometry::porosity)
      .def_readwrite("Cd", &LinerOrificeGeometry::Cd);

  py::class_<LinerFlowState>(m, "LinerFlowState")
      .def(py::init<>())
      .def_readwrite("u_bias", &LinerFlowState::u_bias)
      .def_readwrite("u_grazing", &LinerFlowState::u_grazing);

  py::class_<LinerCavity>(m, "LinerCavity")
      .def(py::init<>())
      .def_readwrite("depth", &LinerCavity::depth);

  m.def("absorption_from_impedance_norm", &absorption_from_impedance_norm,
        py::arg("z_norm"));

  m.def("liner_sdof_impedance_norm", &liner_sdof_impedance_norm,
        py::arg("freq"), py::arg("orifice"), py::arg("cavity"), py::arg("flow"),
        py::arg("medium"));

  m.def("liner_sdof_absorption", &liner_sdof_absorption, py::arg("freq"),
        py::arg("orifice"), py::arg("cavity"), py::arg("flow"),
        py::arg("medium"));

  m.def("sweep_liner_sdof_absorption", &sweep_liner_sdof_absorption,
        py::arg("freqs"), py::arg("orifice"), py::arg("cavity"),
        py::arg("flow"), py::arg("medium"));

  m.def("liner_2dof_serial_impedance_norm", &liner_2dof_serial_impedance_norm,
        py::arg("freq"), py::arg("face_orifice"), py::arg("septum_orifice"),
        py::arg("depth_1"), py::arg("depth_2"), py::arg("face_flow"),
        py::arg("septum_flow"), py::arg("medium"));

  m.def("liner_2dof_serial_absorption", &liner_2dof_serial_absorption,
        py::arg("freq"), py::arg("face_orifice"), py::arg("septum_orifice"),
        py::arg("depth_1"), py::arg("depth_2"), py::arg("face_flow"),
        py::arg("septum_flow"), py::arg("medium"));

  m.def("sweep_liner_2dof_serial_absorption",
        &sweep_liner_2dof_serial_absorption, py::arg("freqs"),
        py::arg("face_orifice"), py::arg("septum_orifice"), py::arg("depth_1"),
        py::arg("depth_2"), py::arg("face_flow"), py::arg("septum_flow"),
        py::arg("medium"));

  m.def("orifice_impedance_with_flow", &orifice_impedance_with_flow,
        py::arg("freq"), py::arg("u_bias"), py::arg("u_grazing"),
        py::arg("d_orifice"), py::arg("l_orifice"), py::arg("porosity"),
        py::arg("Cd"), py::arg("rho"), py::arg("c"));

  m.def("quarter_wave_resonator_tmm", &quarter_wave_resonator_tmm,
        py::arg("freq"), py::arg("L_tube"), py::arg("A_duct"),
        py::arg("A_tube"), py::arg("d_orifice"), py::arg("l_orifice"),
        py::arg("porosity"), py::arg("Cd"), py::arg("u_bias"),
        py::arg("u_grazing"), py::arg("rho"), py::arg("c"));

  m.def("is_whistling_risk", &is_whistling_risk, py::arg("freq"),
        py::arg("u_bias"), py::arg("d_orifice"));

  // =========================================================================
  // Can-Annular Combustor Acoustics (Bloch-Floquet Theory)
  // =========================================================================

  py::class_<CanAnnularGeometry>(m, "CanAnnularGeometry")
      .def(py::init<>())
      .def_readwrite("n_cans", &CanAnnularGeometry::n_cans)
      .def_readwrite("length_can", &CanAnnularGeometry::length_can)
      .def_readwrite("area_can", &CanAnnularGeometry::area_can)
      .def_readwrite("radius_plenum", &CanAnnularGeometry::radius_plenum)
      .def_readwrite("area_plenum", &CanAnnularGeometry::area_plenum)
      .def("__repr__", [](const CanAnnularGeometry &g) {
        return "CanAnnularGeometry(n_cans=" + std::to_string(g.n_cans) +
               ", length_can=" + std::to_string(g.length_can) +
               ", area_can=" + std::to_string(g.area_can) +
               ", radius_plenum=" + std::to_string(g.radius_plenum) +
               ", area_plenum=" + std::to_string(g.area_plenum) + ")";
      });

  py::class_<BlochMode>(m, "BlochMode")
      .def(py::init<>())
      .def_readwrite("m_azimuthal", &BlochMode::m_azimuthal)
      .def_readwrite("frequency", &BlochMode::frequency)
      .def_readwrite("n_cans", &BlochMode::n_cans)
      .def("symmetry_type", &BlochMode::symmetry_type)
      .def("__repr__", [](const BlochMode &m) {
        return "BlochMode(m=" + std::to_string(m.m_azimuthal) +
               ", f=" + std::to_string(m.frequency) + " Hz, " +
               m.symmetry_type() + ")";
      });

  m.def(
      "can_annular_eigenmodes", &can_annular_eigenmodes, py::arg("geom"),
      py::arg("c_can"), py::arg("c_plenum"), py::arg("rho_can"),
      py::arg("rho_plenum"), py::arg("f_max") = 2000.0,
      py::arg("bc_can_top") = BoundaryCondition::Closed,
      "Find all can-annular combustor eigenmodes using Bloch-Floquet "
      "theory.\n\n"
      "Solves the dispersion relation D(omega) = Y_can + Y_annulus = 0\n"
      "for a system of N cans coupled by an annular plenum.\n\n"
      "Uses magnitude minimization + Muller refinement on the pole-free\n"
      "dispersion relation (the Argument Principle variant is C++ only).\n\n"
      "Parameters:\n"
      "  geom       : CanAnnularGeometry with N cans, lengths, areas\n"
      "  c_can      : speed of sound in can [m/s]\n"
      "  c_plenum   : speed of sound in plenum [m/s] (often different!)\n"
      "  rho_can    : density in can [kg/m^3]\n"
      "  rho_plenum : density in plenum [kg/m^3]\n"
      "  f_max      : maximum frequency to search [Hz] (default: 2000)\n"
      "  bc_can_top : boundary condition at can top (default: Closed)\n\n"
      "Returns: List of BlochMode sorted by frequency\n\n"
      "Method: Continuous sliding window scan + magnitude minimization\n"
      "Accuracy: Finds modes wherever they occur (no harmonic assumptions)\n\n"
      "References:\n"
      "  - Noiray et al. (2011): Unified framework for combustion instability\n"
      "  - Evesque & Polifke (2005): Low-order acoustic modelling\n\n"
      "Example:\n"
      "  >>> geom = cb.CanAnnularGeometry()\n"
      "  >>> geom.n_cans = 24\n"
      "  >>> geom.length_can = 0.5  # m\n"
      "  >>> geom.area_can = 0.01   # m^2\n"
      "  >>> geom.radius_plenum = 0.5  # m\n"
      "  >>> geom.area_plenum = 0.02   # m^2\n"
      "  >>> modes = cb.can_annular_eigenmodes(geom, c_can=500, c_plenum=600,\n"
      "  ...                                    rho_can=1.0, rho_plenum=1.2)\n"
      "  >>> for mode in modes[:5]:\n"
      "  ...     print(f'm={mode.m_azimuthal}, f={mode.frequency:.1f} Hz')");

  // Adapter function to convert CanAnnularFlowGeometry to CanAnnularGeometry
  m.def(
      "to_acoustic_geometry", &to_acoustic_geometry, py::arg("flow_geom"),
      py::arg("n_cans"), py::arg("L_can"), py::arg("D_can"),
      "Convert CanAnnularFlowGeometry to CanAnnularGeometry for acoustic "
      "solvers.\n\n"
      "This adapter bridges flow geometry (with primary/transition/annular "
      "sections)\n"
      "to acoustic geometry (cans + annular plenum).\n\n"
      "Parameters:\n"
      "  flow_geom : CanAnnularFlowGeometry from geometry.h\n"
      "  n_cans    : Number of cans (required for acoustic model)\n"
      "  L_can     : Can length [m] (acoustic model uses single section)\n"
      "  D_can     : Can diameter [m] (for circular cross-section)\n\n"
      "Returns: CanAnnularGeometry for can_annular_eigenmodes()\n\n"
      "Notes:\n"
      "  - Plenum radius derived from flow_geom.D_mean() / 2\n"
      "  - Plenum area from flow_geom.area()\n"
      "  - Can area computed as pi*(D_can/2)^2\n\n"
      "Example:\n"
      "  >>> flow_geom = cb.CanAnnularFlowGeometry(L=0.3, D_inner=0.4, "
      "D_outer=0.6,\n"
      "  ...                                        L_primary=0.2, "
      "D_primary=0.1,\n"
      "  ...                                        L_transition=0.05)\n"
      "  >>> geom = cb.to_acoustic_geometry(flow_geom, n_cans=24, L_can=0.5, "
      "D_can=0.1)\n"
      "  >>> modes = cb.can_annular_eigenmodes(geom, c_can=500, c_plenum=550,\n"
      "  ...                                    rho_can=1.2, rho_plenum=1.0)");

  // Annular Duct Geometry and Modes
  // AnnularDuctGeometry is a Python-level alias for Annulus (backward
  // compatibility)
  m.attr("AnnularDuctGeometry") = geom.attr("Annulus");

  py::class_<AnnularMode>(m, "AnnularMode")
      .def(py::init<>())
      .def_readwrite("m_azimuthal", &AnnularMode::m_azimuthal)
      .def_readwrite("n_axial", &AnnularMode::n_axial)
      .def_readwrite("frequency", &AnnularMode::frequency)
      .def("mode_type", &AnnularMode::mode_type)
      .def("__repr__", [](const AnnularMode &mode) {
        return "AnnularMode(m=" + std::to_string(mode.m_azimuthal) +
               ", n=" + std::to_string(mode.n_axial) +
               ", f=" + std::to_string(mode.frequency) + " Hz, " +
               mode.mode_type() + ")";
      });

  m.def("annular_duct_eigenmodes", &annular_duct_eigenmodes, py::arg("geom"),
        py::arg("c"), py::arg("rho"), py::arg("f_max"),
        py::arg("bc_ends") = BoundaryCondition::Closed,
        "Find all annular duct eigenmodes using Argument Principle.\n\n"
        "Pure annular waveguide (no cans) with Bloch-Floquet periodic BCs.\n"
        "Uses Nyquist contour integration for robust root finding.\n\n"
        "Parameters:\n"
        "  geom    : Annulus geometry (from geometry.h)\n"
        "  c       : speed of sound [m/s]\n"
        "  rho     : density [kg/m^3]\n"
        "  f_max   : maximum frequency [Hz]\n"
        "  bc_ends : boundary condition at the far end x = L (default:\n"
        "            Closed). The end x = 0 is always rigid: Closed gives\n"
        "            closed-closed modes f_ax = n c / (2 L); Open gives\n"
        "            closed-open modes f_ax = (2n - 1) c / (4 L), not open-open.\n\n"
        "Returns: List of AnnularMode sorted by frequency\n\n"
        "Method: Argument Principle (gold standard for finding ALL roots)\n"
        "Guarantees: No missed modes, handles complex frequencies\n\n"
        "Example:\n"
        "  >>> geom = cb.Annulus(L=1.0, D_inner=0.4, D_outer=1.0)\n"
        "  >>> geom.n_azimuthal_max = 5\n"
        "  >>> modes = cb.annular_duct_eigenmodes(geom, c=500, rho=1.2, "
        "f_max=1000)\n"
        "  >>> for mode in modes[:5]:\n"
        "  ...     print(f'm={mode.m_azimuthal}, n={mode.n_axial}, "
        "f={mode.frequency:.1f} Hz')");

  m.def("annular_duct_modes_analytical", &annular_duct_modes_analytical,
        py::arg("geom"), py::arg("c"), py::arg("f_max"),
        py::arg("bc_ends") = BoundaryCondition::Closed);

  // Utility functions
  m.def("wavelength", &wavelength, py::arg("f"), py::arg("c"),
        "Acoustic wavelength from frequency.\n\n"
        "lambda = c / f\n\n"
        "Parameters:\n"
        "  f : frequency [Hz]\n"
        "  c : speed of sound [m/s]\n\n"
        "Returns: wavelength [m]");

  m.def("frequency_from_wavelength", &frequency_from_wavelength,
        py::arg("lambda"), py::arg("c"),
        "Frequency from wavelength.\n\n"
        "f = c / lambda\n\n"
        "Parameters:\n"
        "  lambda : wavelength [m]\n"
        "  c      : speed of sound [m/s]\n\n"
        "Returns: frequency [Hz]");

  m.def("acoustic_impedance", &acoustic_impedance, py::arg("rho"), py::arg("c"),
        "Characteristic acoustic impedance.\n\n"
        "Z = rho * c\n\n"
        "Parameters:\n"
        "  rho : density [kg/m^3]\n"
        "  c   : speed of sound [m/s]\n\n"
        "Returns: impedance [Pa*s/m]");

  m.def("sound_pressure_level", &sound_pressure_level, py::arg("p_rms"),
        py::arg("p_ref") = 20e-6,
        "Sound pressure level in dB.\n\n"
        "SPL = 20 * log10(p_rms / p_ref)\n\n"
        "Parameters:\n"
        "  p_rms : RMS pressure amplitude [Pa]\n"
        "  p_ref : reference pressure [Pa] (default: 20e-6 Pa in air)\n\n"
        "Returns: SPL [dB]");

  m.def("particle_velocity", &particle_velocity, py::arg("p"), py::arg("rho"),
        py::arg("c"),
        "Particle velocity from pressure amplitude.\n\n"
        "u = p / (rho * c)\n\n"
        "Parameters:\n"
        "  p   : pressure amplitude [Pa]\n"
        "  rho : density [kg/m^3]\n"
        "  c   : speed of sound [m/s]\n\n"
        "Returns: particle velocity [m/s]");

  // Residence time
  m.def("residence_time_tube",
        py::overload_cast<const Tube &, double>(&residence_time),
        py::arg("tube"), py::arg("Q"), "Residence time for a tube [s].");

  m.def("residence_time_annulus",
        py::overload_cast<const Annulus &, double>(&residence_time),
        py::arg("annulus"), py::arg("Q"), "Residence time for an annulus [s].");

  m.def("residence_time_can_annular",
        py::overload_cast<const CanAnnularFlowGeometry &, double>(
            &residence_time),
        py::arg("geom"), py::arg("Q"),
        "Residence time for CanAnnularFlowGeometry based on total volume [s].");

  m.def("residence_time_mdot_can_annular",
        py::overload_cast<const CanAnnularFlowGeometry &, double, double>(
            &residence_time_mdot),
        py::arg("geom"), py::arg("mdot"), py::arg("rho"),
        "Residence time for CanAnnularFlowGeometry from mass flow [s].");

  // Nozzle thrust (from NozzleSolution)
  m.def("nozzle_thrust",
        py::overload_cast<const NozzleSolution &, double>(&nozzle_thrust),
        py::arg("sol"), py::arg("P_amb"),
        "Calculate rocket nozzle thrust from nozzle solution.\n\n"
        "F = mdot * u_e + (P_e - P_amb) * A_e\n\n"
        "sol   : NozzleSolution from nozzle_cd()\n"
        "P_amb : Ambient pressure [Pa]\n\n"
        "Returns ThrustResult with thrust, Isp, and C_F.");

  // Nozzle thrust (convenience, from parameters)
  m.def(
      "nozzle_thrust_cd",
      [](double Tt, double Pt, double P_design, double P_amb, double A_inlet,
         double A_throat, double A_exit, double x_throat, double x_exit,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         std::size_t n_stations, double tol, std::size_t max_iter) {
        auto X = to_vec(X_arr);
        return nozzle_thrust(Tt, Pt, P_design, P_amb, A_inlet, A_throat, A_exit,
                             x_throat, x_exit, X, n_stations, tol, max_iter);
      },
      py::arg("Tt"), py::arg("Pt"), py::arg("P_design"), py::arg("P_amb"),
      py::arg("A_inlet"), py::arg("A_throat"), py::arg("A_exit"),
      py::arg("x_throat"), py::arg("x_exit"), py::arg("X"),
      py::arg("n_stations") = 100, py::arg("tol") = 1e-8,
      py::arg("max_iter") = 50,
      "Calculate rocket nozzle thrust directly from geometry.\n\n"
      "P_design: nozzle exit design pressure for Mach profile (must be > 0).\n"
      "P_amb: ambient pressure for thrust equation (may be 0 for vacuum).\n"
      "Convenience function that calls nozzle_cd() then nozzle_thrust().\n\n"
      "Returns ThrustResult with thrust, Isp, and C_F.");

  // =========================================================================
  // Incompressible flow
  // =========================================================================

  // Bernoulli equation
  m.def("bernoulli_P2", &bernoulli_P2, py::arg("P1"), py::arg("v1"),
        py::arg("v2"), py::arg("rho"), py::arg("dz") = 0.0,
        py::arg("g") = 9.80665,
        "Bernoulli equation: solve for downstream pressure P2.\n\n"
        "P1 + 0.5*rho*v1^2 + rho*g*h1 = P2 + 0.5*rho*v2^2 + rho*g*h2\n\n"
        "P1  : upstream pressure [Pa]\n"
        "v1  : upstream velocity [m/s]\n"
        "v2  : downstream velocity [m/s]\n"
        "rho : fluid density [kg/m^3]\n"
        "dz  : elevation change z2-z1 [m] (positive=up)\n"
        "g   : gravitational acceleration [m/s^2]");

  m.def("bernoulli_v2", &bernoulli_v2, py::arg("P1"), py::arg("P2"),
        py::arg("v1"), py::arg("rho"), py::arg("dz") = 0.0,
        py::arg("g") = 9.80665,
        "Bernoulli equation: solve for downstream velocity v2.\n\n"
        "P1  : upstream pressure [Pa]\n"
        "P2  : downstream pressure [Pa]\n"
        "v1  : upstream velocity [m/s]\n"
        "rho : fluid density [kg/m^3]\n"
        "dz  : elevation change z2-z1 [m]\n"
        "g   : gravitational acceleration [m/s^2]");

  // Orifice flow
  m.def("orifice_mdot",
        static_cast<double (*)(double, double, double, double, double)>(
            &orifice_mdot),
        py::arg("P1"), py::arg("P2"), py::arg("A"), py::arg("Cd"),
        py::arg("rho"),
        "Orifice mass flow rate (incompressible).\n\n"
        "mdot = Cd * A * sqrt(2 * rho * dP)\n\n"
        "P1  : upstream pressure [Pa]\n"
        "P2  : downstream pressure [Pa]\n"
        "A   : orifice area [m^2]\n"
        "Cd  : discharge coefficient [-] (typically 0.6-0.65)\n"
        "rho : fluid density [kg/m^3]\n\n"
        "Returns: mass flow rate [kg/s]");

  m.def("orifice_Q", &orifice_Q, py::arg("P1"), py::arg("P2"), py::arg("A"),
        py::arg("Cd"), py::arg("rho"),
        "Orifice volumetric flow rate.\n\n"
        "Q = Cd * A * sqrt(2 * dP / rho)\n\n"
        "Returns: volumetric flow rate [m^3/s]");

  m.def("orifice_velocity", &orifice_velocity, py::arg("P1"), py::arg("P2"),
        py::arg("rho"),
        "Ideal orifice velocity (no losses).\n\n"
        "v = sqrt(2 * dP / rho)\n\n"
        "Returns: velocity [m/s]");

  m.def("orifice_area", &orifice_area, py::arg("mdot"), py::arg("P1"),
        py::arg("P2"), py::arg("Cd"), py::arg("rho"),
        "Solve for required orifice area given mass flow.\n\n"
        "A = mdot / (Cd * sqrt(2 * rho * dP))\n\n"
        "Returns: orifice area [m^2]");

  m.def("orifice_dP",
        static_cast<double (*)(double, double, double, double)>(&orifice_dP),
        py::arg("mdot"), py::arg("A"), py::arg("Cd"), py::arg("rho"),
        "Solve for pressure drop given mass flow and area.\n\n"
        "dP = (mdot / (Cd * A))^2 / (2 * rho)\n\n"
        "Returns: pressure drop P1-P2 [Pa]");

  // Channel flow
  m.def("channel_dP", &channel_dP, py::arg("v"), py::arg("L"), py::arg("D"),
        py::arg("f"), py::arg("rho"),
        "Channel pressure drop (Darcy-Weisbach).\n\n"
        "dP = f * (L/D) * (rho * v^2 / 2)\n\n"
        "v   : flow velocity [m/s]\n"
        "L   : channel length [m]\n"
        "D   : channel diameter [m]\n"
        "f   : Darcy friction factor [-]\n"
        "rho : fluid density [kg/m^3]\n\n"
        "Returns: pressure drop [Pa]");

  m.def("channel_dP_mdot", &channel_dP_mdot, py::arg("mdot"), py::arg("L"),
        py::arg("D"), py::arg("f"), py::arg("rho"),
        "Channel pressure drop from mass flow rate.\n\n"
        "Returns: pressure drop [Pa]");

  m.def("channel_velocity", &channel_velocity, py::arg("mdot"), py::arg("D"),
        py::arg("rho"),
        "Channel velocity from mass flow rate.\n\n"
        "v = mdot / (rho * pi * D^2 / 4)\n\n"
        "Returns: velocity [m/s]");

  // Hydraulic utilities
  m.def("dynamic_pressure", &dynamic_pressure, py::arg("v"), py::arg("rho"),
        "Dynamic pressure (velocity head).\n\n"
        "q = 0.5 * rho * v^2\n\n"
        "Returns: dynamic pressure [Pa]");

  m.def("velocity_from_q", &velocity_from_q, py::arg("q"), py::arg("rho"),
        "Velocity from dynamic pressure.\n\n"
        "v = sqrt(2 * q / rho)\n\n"
        "Returns: velocity [m/s]");

  // =========================================================================
  // Orifice Cd Correlations
  // =========================================================================

  // OrificeGeometry struct
  py::class_<OrificeGeometry>(m, "OrificeGeometry",
                              "Orifice geometry for Cd correlations.\n\n"
                              "Attributes:\n"
                              "  d : orifice bore diameter [m]\n"
                              "  D : channel diameter [m]\n"
                              "  t : plate thickness [m] (for thick plate)\n"
                              "  r : inlet edge radius [m] (for rounded entry)")
      .def(py::init<>())
      .def_readwrite("d", &OrificeGeometry::d, "Orifice bore diameter [m]")
      .def_readwrite("D", &OrificeGeometry::D, "Channel diameter [m]")
      .def_readwrite("tap", &OrificeGeometry::tap,
                     "Pressure taps the metering Cd refers to (default Flange)")
      .def_readwrite("t", &OrificeGeometry::t, "Plate thickness [m]")
      .def_readwrite("r", &OrificeGeometry::r, "Inlet edge radius [m]")
      .def_readwrite("bevel", &OrificeGeometry::bevel, "Bevel angle [rad]")
      .def_property_readonly("beta", &OrificeGeometry::beta,
                             "Diameter ratio d/D [-]")
      .def_property_readonly("area", &OrificeGeometry::area,
                             "Orifice area [m^2]")
      .def_property_readonly("t_over_d", &OrificeGeometry::t_over_d,
                             "Thickness ratio t/d [-]")
      .def_property_readonly("r_over_d", &OrificeGeometry::r_over_d,
                             "Radius ratio r/d [-]")
      .def("is_valid", &OrificeGeometry::is_valid,
           "Check if geometry is valid");

  // OrificeState struct
  py::class_<OrificeState>(m, "OrificeState",
                           "Flow state at orifice for Cd correlations.\n\n"
                           "Attributes:\n"
                           "  Re_D : channel Reynolds number (based on D) [-]\n"
                           "  dP   : differential pressure [Pa]\n"
                           "  rho  : fluid density [kg/m^3]\n"
                           "  mu   : dynamic viscosity [Pa*s]")
      .def(py::init<>())
      .def_readwrite("Re_D", &OrificeState::Re_D, "Channel Reynolds number [-]")
      .def_readwrite("dP", &OrificeState::dP, "Differential pressure [Pa]")
      .def_readwrite("rho", &OrificeState::rho, "Fluid density [kg/m^3]")
      .def_readwrite("mu", &OrificeState::mu, "Dynamic viscosity [Pa*s]")
      .def("Re_d", &OrificeState::Re_d, py::arg("beta"),
           "Orifice Reynolds number (based on d)");

  // Cd correlation functions
  m.def("Cd_sharp_thin_plate_and_dRe", &Cd_sharp_thin_plate_and_dRe, py::arg("geom"),
        py::arg("state"),
        "(Cd, dCd/dRe_D) for the sharp thin-plate orifice (ISO 5167-2,\n"
        "Reader-Harris/Gallagher): the analytic Re_D derivative.");
  m.def("Cd_sharp_thin_plate", &Cd_sharp_thin_plate, py::arg("geom"),
        py::arg("state"),
        "Discharge coefficient for sharp thin-plate orifice (ISO 5167-2).\n\n"
        "Uses Reader-Harris/Gallagher at the taps geom.tap, with the\n"
        "small-pipe term below D = 71.12 mm. Range: metering_cd_in_range.");
  m.def(
      "metering_cd_and_dRe", &metering_cd_and_dRe,
      py::arg("correlation"), py::arg("geom"), py::arg("Re_D"),
      "(Cd, dCd/dRe_D) of a metering correlation at the taps geom.tap.\n\n"
      "ReaderHarrisGallagher (ISO 5167-2:2003; Stolz was removed in 0.9.0,\n"
      "#876): below the lower Re_D bound Re_b\n"
      "(metering_cd_re_min) Cd is HELD, entered C1 in ln Re_D over 0.1 just\n"
      "below Re_b (#873): the slope ramps linearly to 0, so Cd and\n"
      "dCd/dRe_D stay continuous and values at Re_D >= Re_b are Eq. (4)'s\n"
      "own. metering_cd_in_range is False there. Constant and UserFunction\n"
      "raise ValueError: their Cd is the caller's.");
  m.def("metering_cd_re_min", &metering_cd_re_min, py::arg("correlation"), py::arg("geom"),
        "Lower Re_D bound of Reader-Harris/Gallagher's stated range for this\n"
        "geometry, where metering_cd_and_dRe's hold begins (ISO 5167-2:2003\n"
        "5.3.1): corner/D-D/2 5000 (beta <= 0.56) or 16000 beta^2; flange\n"
        "max(5000, 170 beta^2 D[mm]). ValueError for every other correlation:\n"
        "Constant/UserFunction carry no formula.");
  m.def(
      "metering_cd_in_range", &metering_cd_in_range,
      py::arg("correlation"), py::arg("geom"), py::arg("Re_D"),
      "Is the plate inside the correlation's stated range? A diagnostic,\n"
      "never a selector. ReaderHarrisGallagher: ISO 5167-2:2003 5.3.1 (d,\n"
      "D, beta and Re_D >= metering_cd_re_min).");
  m.def("permanent_loss_ratio_and_dC", &permanent_loss_ratio_and_dC, py::arg("beta"),
        py::arg("C"),
        "(k, dk/dC): ISO 5167-2:2003 5.4.1's permanent loss of an orifice\n"
        "plate between pipes over its tapping differential,\n"
        "k = [s - C b^2] / [s + C b^2], s = sqrt(1 - b^4 (1 - C^2)).\n"
        "Derived for vena-contracta (D-D/2) taps. Raises ValueError outside\n"
        "0 < beta < 1 or C <= 0.");
  m.def("velocity_of_approach_factor", &velocity_of_approach_factor, py::arg("beta"),
        "ISO 5167-1 velocity-of-approach factor E = 1 / sqrt(1 - beta^4).\n"
        "Returns 1 outside 0 < beta < 1 (no approach pipe, or a bore at\n"
        "least as wide as the pipe): the network's no-correction convention.");
  m.def("jobson_force_defect", &jobson_force_defect, py::arg("Ci"), py::arg("beta") = 0.0,
        py::arg("Cv") = 1.0,
        "Jobson (1955) force-defect coefficient with the approach-velocity\n"
        "extension (#854): f = 1/Ci - 1/(2 Ci^2) - Cv beta^2 + Cv^2 beta^4/2.\n"
        "beta = 0 is the paper's Eq. (5).");
  m.def("jobson_ci_from_cd", &jobson_ci_from_cd, py::arg("Cd"), py::arg("beta") = 0.0,
        py::arg("Cv") = 1.0,
        "Incompressible contraction Ci of an orifice whose incompressible\n"
        "flow is Cd E A sqrt(2 rho1 dp), with Jobson's velocity coefficient\n"
        "Cv: Cv Ci / sqrt(1 - (Cv Ci)^2 beta^4) = Cd E. beta = 0: Cd / Cv.");
  m.def("jobson_contraction_coefficient", &jobson_contraction_coefficient, py::arg("r"),
        py::arg("Ci"), py::arg("n"),
        "Jobson (1955) compressible contraction coefficient C of a\n"
        "reservoir-fed sharp-edged orifice: Eq. (11) for r = p/p0 >= r_c,\n"
        "Eq. (15) below. Mass flow Cv C K_N A sqrt(p0 rho0).");
  m.def("jobson_expansibility", &jobson_expansibility, py::arg("r"), py::arg("beta"),
        py::arg("Ci"), py::arg("n"),
        "Expansibility implied by Jobson's theory with the approach-velocity\n"
        "extension (#854): mdot / mdot_inc at the same upstream static state\n"
        "and Ci, r = p2/p1. Comparable with ISO 5167-2 Eq. (5).");
  // The beta extension's checked range: OrificeElement selects jobson1955
  // by default only at or below it (#861).
  m.attr("JOBSON_1955_BETA_MAX") = orifice::jobson1955::beta_max;




  m.def("mcgreehan_schotsch_1988_cd_and_derivatives",
        &orifice::mcgreehan_schotsch::cd_and_derivatives, py::arg("Re"),
        py::arg("r_over_d"), py::arg("L_over_d"), py::arg("U1_over_Vi") = 0.0,
        "Cd and its derivatives: (Cd, dCd/dRe, dCd/d(U1/Vi)).\n\n"
        "The solver-facing form of mcgreehan_schotsch_1988_cd. r/d and L/d are\n"
        "geometry, never solver unknowns, so they carry no partials.\n\n"
        "TWO NUMERICAL TREATMENTS, both solver aids and both stated:\n\n"
        "1. The crossflow input is regularised,\n"
        "   U1/Vi -> sqrt((U1/Vi)^2 + eps^2) with eps = 1e-6. Eq. (17)'s\n"
        "   Rv^0.6 and Rv^0.9 have unbounded slope at Rv = 0, and U1/Vi = 0 is\n"
        "   both the default and the physically right value for a plenum-fed\n"
        "   jet plate, so the singularity sits where the solver lives. This\n"
        "   bounds dCd/d(U1/Vi) at 1.451 instead of 1.4e4 and DOES move Cd, by\n"
        "   5.5e-5 (0.007%). eps sits mid-plateau: the bound is flat at 1.451\n"
        "   for every eps in [1e-7, 1e-4], so there is no peak being fitted.\n\n"
        "2. Below Re = 1e4 the VALUE stays exactly floored -- Eq. (8) is\n"
        "   invalid and divergent there -- but the DERIVATIVE is continued\n"
        "   from the floor rather than reported as the true zero, so a Newton\n"
        "   step below it can climb back and dCd/dRe is continuous across the\n"
        "   floor. This changes no reported Cd.\n\n"
        "Returns:\n"
        "  (Cd, dCd/dRe, dCd/d(U1_over_Vi))");

  m.def("mcgreehan_schotsch_1988_cd", &orifice::Cd_McGreehanSchotsch,
        py::arg("Re"), py::arg("r_over_d"), py::arg("L_over_d"),
        py::arg("U1_over_Vi") = 0.0,
        py::arg("eps") = orifice::mcgreehan_schotsch::rv_smooth_eps,
        "Discharge coefficient for a long orifice with corner radiusing and\n"
        "inlet crossflow (McGreehan and Schotsch 1988, ASME J. Turbomachinery\n"
        "110(2), 213-217, Eqs. 8-17).\n\n"
        "This is the plenum-to-plenum configuration -- a cooling transfer hole\n"
        "or a jet-plate hole -- NOT an ISO 5167 metering orifice in a pipe run;\n"
        "for that use Cd_sharp_thin_plate and its siblings.\n\n"
        "Args:\n"
        "  Re: orifice Reynolds number, V_i * d / nu. Stated valid Re >= 1e4;\n"
        "      below that the correlation is held at that floor, because Eq. 8\n"
        "      diverges (it reaches Cd = 1.0 at Re = 904).\n"
        "  r_over_d: inlet corner radius / bore. 0 for a sharp hole. An ASME\n"
        "      nozzle is reached at 0.82; further radiusing buys nothing.\n"
        "  L_over_d: orifice flat length / bore (plate thickness ratio).\n"
        "  U1_over_Vi: INLET (approach, supply-side) tangential velocity over\n"
        "      ideal through-flow velocity. Defaults to 0, the plenum-fed case.\n\n"
        "      This is NOT a discharge-side crossflow ratio. Do not pass\n"
        "      Florschuetz's Gc/Gj: that is spent-air crossflow in the\n"
        "      impingement channel, on the far face of the plate, and using it\n"
        "      here gives a plausible number that is biased low.\n\n"
        "      V_i must be built from STATIC inlet conditions, not from a total\n"
        "      pressure that already contains the tangential velocity head.\n\n"
        "Not monotonic in U1_over_Vi: Cd rises up to ~5.7% above its\n"
        "zero-crossflow value near U1_over_Vi ~ 0.09 before falling. That is\n"
        "the source's Fig. 4 and its data, not an artifact.\n\n"
        "For a bare plenum-fed jet plate (r/d=0, t/d=1, Re=1e4) this returns\n"
        "0.785, against Florschuetz's recommended 0.79 for the same thing.\n\n"
        "See validation/cooling/extractions/orifice_discharge_coefficient.md.");

  m.def("mcgreehan_schotsch_1988_crossflow_cd",
        &orifice::mcgreehan_schotsch::cd_with_crossflow, py::arg("cd_base"),
        py::arg("U1_over_Vi"),
        py::arg("eps") = orifice::mcgreehan_schotsch::rv_smooth_eps,
        "McGreehan and Schotsch (1988) Eq. (17) alone: the inlet-crossflow\n"
        "correction applied to a baseline discharge coefficient you supply.\n\n"
        "Use when the plate's zero-crossflow Cd is KNOWN -- measured, or a\n"
        "literature value such as Florschuetz's per-configuration Table 1 --\n"
        "and only the crossflow correction is wanted. Otherwise use\n"
        "mcgreehan_schotsch_1988_cd, which computes the baseline too.\n\n"
        "This is how the source uses Eq. (17) in its own validation: Figs. 5\n"
        "and 6 anchor to 'a set baseline point at U1/Vi = 0' from Rohde's\n"
        "measurements (0.64, 0.73, 0.88) rather than to the chain's own\n"
        "prediction, which runs 3-12% higher for the same geometry.\n\n"
        "U1_over_Vi is the INLET (supply-side) tangential velocity ratio; see\n"
        "mcgreehan_schotsch_1988_cd for the trap this must not be fed.\n\n"
        "eps regularises the input, U1/Vi -> sqrt((U1/Vi)^2 + eps^2), so that\n"
        "Eq. (17)'s Rv^0.6 and Rv^0.9 do not give an unbounded derivative at\n"
        "U1/Vi = 0 -- which is the default AND the physically right value for\n"
        "a plenum-fed jet plate. Costs 5.5e-5 on Cd. Pass eps = 0 for the\n"
        "paper exactly. See mcgreehan_schotsch_1988_cd_and_derivatives.");

  m.def("mcgreehan_schotsch_1988_expansion_orifice",
        &orifice::mcgreehan_schotsch::expansion_orifice, py::arg("S"),
        py::arg("gamma"),
        "Eq. (4): orifice adiabatic expansion factor, Y_o = 1 - 0.41(1-S)/gamma.\n"
        "S is the pressure ratio P_s2/P_t1.");

  m.def("mcgreehan_schotsch_1988_expansion_nozzle",
        &orifice::mcgreehan_schotsch::expansion_nozzle, py::arg("S"),
        py::arg("gamma"),
        "Eq. (5): nozzle adiabatic expansion factor -- the exact isentropic\n"
        "one. Reproduces combaero's own nozzle_flow to 0.008% over\n"
        "0.6 <= S <= 0.99, which is what establishes that applying Y on the\n"
        "regime='compressible' path would double-count compressibility.");

  m.def("mcgreehan_schotsch_1988_expansion_factor",
        &orifice::mcgreehan_schotsch::expansion_factor, py::arg("cd"),
        py::arg("S"), py::arg("gamma"),
        py::arg("eps") = orifice::mcgreehan_schotsch::x_smooth_eps,
        "Eqs. (6)/(7): the blended expansion factor Y = (1-X) Y_o + X Y_n,\n"
        "with X = 8.333 (Cd - 0.82) saturated into [0, 1].\n\n"
        "For the INCOMPRESSIBLE orifice equation only:\n"
        "    m_dot = Cd * Y * A * sqrt(2 rho_t1 (P_t1 - P_s2))\n"
        "Do NOT apply it on regime='compressible', which already solves the\n"
        "isentropic nozzle exactly (and Eq. (5) reproduces that solve to\n"
        "0.008%) -- you would correct for compressibility twice.\n\n"
        "eps is the saturation smoothing width. eps=0 gives the paper's exact\n"
        "hard clamp, which has an exactly-zero derivative outside\n"
        "[Cd 0.82, 0.94] and a discontinuous jump of 0.734 in dY/dCd at both\n"
        "knees. The default trades 0.47% in Y for a continuous, non-zero\n"
        "derivative. Y also saturates at the critical pressure ratio, below\n"
        "which the isentropic form predicts decreasing flow.");

  // DischargeHoleGeometry struct
  py::class_<DischargeHoleGeometry>(
      m, "DischargeHoleGeometry",
      "Geometry of one discharge hole. There is deliberately no pipe\n"
      "diameter: a hole in a wall has no beta.")
      .def(py::init<>())
      .def(py::init([](double d, double L, double r, double bevel, double inclination_deg) {
             DischargeHoleGeometry g;
             g.d = d;
             g.L = L;
             g.r = r;
             g.bevel = bevel;
             g.inclination_deg = inclination_deg;
             return g;
           }),
           py::arg("d"), py::arg("L") = 0.0, py::arg("r") = 0.0,
           py::arg("bevel") = 0.0, py::arg("inclination_deg") = 90.0)
      .def_readwrite("d", &DischargeHoleGeometry::d, "Hole diameter [m]")
      .def_readwrite("L", &DischargeHoleGeometry::L,
                     "Hole length along its axis [m] -- the wall thickness for\n"
                     "a normal hole, t/sin(alpha) for an angled one")
      .def_readwrite("r", &DischargeHoleGeometry::r,
                     "Inlet edge radius [m] (0 = sharp)")
      .def_readwrite("bevel", &DischargeHoleGeometry::bevel,
                     "Bevel depth along the hole axis [m] (0 = not beveled).\n"
                     "Idelchik diagram 4-18b's l/Dh at a bevel angle of\n"
                     "40-60 deg -- NOT the wall thickness.")
      .def_readwrite("inclination_deg", &DischargeHoleGeometry::inclination_deg,
                     "Angle between the hole axis and the wall [deg],\n"
                     "Idelchik's delta. 90 (default) = normal hole, no\n"
                     "inclination term. Below 90 opts in to diagram 3-2's\n"
                     "inclined-inlet loss, composed with Idelchik1966Thick\n"
                     "only; every other correlation raises ValueError.")
      .def("is_inclined", &DischargeHoleGeometry::is_inclined,
           "True when inclination_deg < 90: the opt-in term applies")
      .def("L_over_d", &DischargeHoleGeometry::L_over_d, "Length ratio L/d [-]")
      .def("r_over_d", &DischargeHoleGeometry::r_over_d, "Radius ratio r/d [-]")
      .def("bevel_over_d", &DischargeHoleGeometry::bevel_over_d,
           "Bevel depth ratio l/d [-]")
      .def("area", &DischargeHoleGeometry::area, "Hole area [m^2]")
      .def("is_valid", &DischargeHoleGeometry::is_valid,
           "True when d > 0, L, r >= 0 and 0 < inclination_deg <= 90");

  m.def("discharge_inclination_dzeta", &discharge_inclination_dzeta,
        py::arg("inclination_deg"),
        "Idelchik (1966) diagram 3-2's inclined-entrance increment over the\n"
        "normal sharp entrance: dzeta = 0.3 cos(delta) + 0.2 cos^2(delta) on\n"
        "the hole velocity head, delta the hole-axis-to-wall angle [deg];\n"
        "exactly 0 at 90. Geometry only. Raises ValueError outside (0, 90].");

  // DischargeHoleState struct
  py::class_<DischargeHoleState>(m, "DischargeHoleState",
                                 "Flow state at one discharge hole")
      .def(py::init<>())
      .def(py::init([](double Re, double U1_over_Vi) {
             DischargeHoleState f;
             f.Re = Re;
             f.U1_over_Vi = U1_over_Vi;
             return f;
           }),
           py::arg("Re"), py::arg("U1_over_Vi") = 0.0)
      .def_readwrite("Re", &DischargeHoleState::Re,
                     "Hole Reynolds number, based on d [-]")
      .def_readwrite("U1_over_Vi", &DischargeHoleState::U1_over_Vi,
                     "INLET tangential velocity / ideal through-flow velocity\n"
                     "[-]. A plenum-fed hole has 0. Feeding a discharge-side\n"
                     "crossflow ratio here is wrong -- see\n"
                     "mcgreehan_schotsch_1988_cd.");

  m.def(
      "discharge_cd",
      [](DischargeCdCorrelation correlation, const DischargeHoleGeometry& hole,
         const DischargeHoleState& flow, double Cd_constant) {
        auto c = (correlation == DischargeCdCorrelation::Constant && Cd_constant > 0.0)
                     ? make_constant_discharge_correlation(Cd_constant)
                     : make_discharge_correlation(correlation);
        return c->Cd(hole, flow);
      },
      py::arg("correlation"), py::arg("hole"), py::arg("flow"),
      py::arg("Cd_constant") = 0.0,
      "Cd of a discharge hole through the selected correlation.\n\n"
      "Cd_constant is only read for correlation=Constant, where a positive\n"
      "value pins Cd to a measured or literature number instead of the 0.60\n"
      "default. Raises ValueError for a correlation that is declared but not\n"
      "implemented.");

  py::class_<DischargeCdEnvelope>(m, "DischargeCdEnvelope")
      .def_readonly("re", &DischargeCdEnvelope::re)
      .def_readonly("l_over_d", &DischargeCdEnvelope::l_over_d)
      .def_readonly("edge", &DischargeCdEnvelope::edge)
      .def("in_range", &DischargeCdEnvelope::in_range);
  m.def("discharge_cd_envelope", &discharge_cd_envelope, py::arg("correlation"),
        py::arg("hole"), py::arg("flow"),
        "Which inputs lie outside the range the correlation's source fitted or\n"
        "tabulated, by variable: re, l_over_d, edge (bevel or r/d table, or an\n"
        "inclined hole on a correlation that refuses it). Information to\n"
        "report beside the result, never a selector and never a warning:\n"
        "discharge_cd returns a value regardless (held or continued).");
  m.def("discharge_cd_in_range", &discharge_cd_in_range, py::arg("correlation"),
        py::arg("hole"), py::arg("flow"),
        "True when no input lies outside the correlation's source range:\n"
        "discharge_cd_envelope(...).in_range(). A flag to report, never a\n"
        "selector; discharge_cd still returns a value outside it.");
  m.def("discharge_cd_re_floor", &discharge_cd_re_floor, py::arg("correlation"),
        "Hole Reynolds number below which the correlation HOLDS its Cd value\n"
        "(Re clamped up to it), so the value's true dCd/dRe there is 0.\n"
        "discharge_cd_and_derivatives continues the slope below it for\n"
        "McGreehanSchotsch1988 and Idelchik (a Newton aid); a sensitivity\n"
        "built on the held value must use 0 below this Re. 0 for Constant.");
  m.def("discharge_pipe_zeta", &discharge_pipe_zeta, py::arg("correlation"), py::arg("hole"),
        py::arg("beta"),
        "Idelchik (1966) diagrams 4-13..4-16: the hole's permanent loss as an\n"
        "orifice in a straight pipe, zeta1 = DH / (rho w1^2 / 2), referenced to\n"
        "the PIPE velocity. High-Re, incompressible. Raises ValueError outside\n"
        "0 < beta < 1 and for Constant.");
  m.def("discharge_pipe_cd_ratio", &discharge_pipe_cd_ratio, py::arg("correlation"),
        py::arg("hole"), py::arg("beta"),
        "R with C_pipe = Cd_wall R: carries a plenum-measured discharge Cd into a\n"
        "pipe through ISO 5167-2's permanent-loss relation, so the network's\n"
        "recovery row gives discharge_pipe_zeta. R -> 1 as beta -> 0. Geometry\n"
        "only. 1 outside 0 < beta < 1 and for Constant.");

  m.def(
      "discharge_cd_and_derivatives",
      [](DischargeCdCorrelation correlation, const DischargeHoleGeometry& hole,
         const DischargeHoleState& flow, double Cd_constant) {
        auto c = (correlation == DischargeCdCorrelation::Constant && Cd_constant > 0.0)
                     ? make_constant_discharge_correlation(Cd_constant)
                     : make_discharge_correlation(correlation);
        return c->Cd_and_derivatives(hole, flow);
      },
      py::arg("correlation"), py::arg("hole"), py::arg("flow"),
      py::arg("Cd_constant") = 0.0,
      "Solver-facing (f, J) form of discharge_cd:\n"
      "    (Cd, dCd/dRe, dCd/d(U1_over_Vi))\n"
      "Analytic throughout. L/d and r/d are geometry and never enter the\n"
      "Jacobian; a Constant correlation returns exactly (Cd, 0, 0).");

  // The hole-Cd chain: every term of a hole element's Cd in one call
  // (hole_cd.h). Internal composition point of the network's hole elements.
  m.attr("HOLE_CD_FIXED_MIN") = HOLE_CD_FIXED_MIN;
  py::enum_<HoleCdBase>(m, "HoleCdBase", "Base Cd of the hole-Cd chain (hole_cd_compose).")
      .value("Fixed", HoleCdBase::Fixed)
      .value("Metering", HoleCdBase::Metering)
      .value("Discharge", HoleCdBase::Discharge);
  py::enum_<HoleCdJetCrossflow>(m, "HoleCdJetCrossflow",
                                "Impingement jet-crossflow law of the hole-Cd chain.")
      .value("Off", HoleCdJetCrossflow::Off)
      .value("RogersHersh1976", HoleCdJetCrossflow::RogersHersh1976)
      .value("FlorschuetzMetzgerSu1982", HoleCdJetCrossflow::FlorschuetzMetzgerSu1982);
  py::class_<HoleCdModel>(m, "HoleCdModel",
                          "Resolved options of the hole-Cd chain (hole_cd_compose).")
      .def(py::init<>())
      .def_readwrite("base", &HoleCdModel::base)
      .def_readwrite("fixed_cd", &HoleCdModel::fixed_cd)
      .def_readwrite("metering", &HoleCdModel::metering)
      .def_readwrite("discharge", &HoleCdModel::discharge)
      .def_readwrite("continue_re_floor", &HoleCdModel::continue_re_floor)
      .def_readwrite("exit_crossflow", &HoleCdModel::exit_crossflow)
      .def_readwrite("exit_set", &HoleCdModel::exit_set)
      .def_readwrite("exit_theta_deg", &HoleCdModel::exit_theta_deg)
      .def_readwrite("exit_alpha_deg", &HoleCdModel::exit_alpha_deg)
      .def_readwrite("jet", &HoleCdModel::jet)
      .def_readwrite("xn_d", &HoleCdModel::xn_d)
      .def_readwrite("yn_d", &HoleCdModel::yn_d)
      .def_readwrite("z_d", &HoleCdModel::z_d);
  py::class_<HoleCdInputs>(m, "HoleCdInputs", "State the hole-Cd chain reads (hole_cd_compose).")
      .def(py::init<>())
      .def_readwrite("Re", &HoleCdInputs::Re)
      .def_readwrite("U1_over_Vi", &HoleCdInputs::U1_over_Vi)
      .def_readwrite("Re_D", &HoleCdInputs::Re_D)
      .def_readwrite("G_h", &HoleCdInputs::G_h)
      .def_readwrite("rho_h", &HoleCdInputs::rho_h)
      .def_readwrite("rho_m_u_m2", &HoleCdInputs::rho_m_u_m2)
      .def_readwrite("G_ideal", &HoleCdInputs::G_ideal)
      .def_readwrite("G_c", &HoleCdInputs::G_c)
      .def_readwrite("suction", &HoleCdInputs::suction);
  py::class_<HoleCdResult>(m, "HoleCdResult", "Composed Cd, its partials and each term.")
      .def_readonly("Cd", &HoleCdResult::Cd)
      .def_readonly("dCd_dRe", &HoleCdResult::dCd_dRe)
      .def_readonly("dCd_dU1_over_Vi", &HoleCdResult::dCd_dU1_over_Vi)
      .def_readonly("dCd_dRe_D", &HoleCdResult::dCd_dRe_D)
      .def_readonly("dCd_dG_h", &HoleCdResult::dCd_dG_h)
      .def_readonly("dCd_drho_h", &HoleCdResult::dCd_drho_h)
      .def_readonly("dCd_drho_m_u_m2", &HoleCdResult::dCd_drho_m_u_m2)
      .def_readonly("dCd_dG_ideal", &HoleCdResult::dCd_dG_ideal)
      .def_readonly("dCd_dG_c", &HoleCdResult::dCd_dG_c)
      .def_readonly("cd_base", &HoleCdResult::cd_base)
      .def_readonly("cd_hole", &HoleCdResult::cd_hole)
      .def_readonly("dcd_base_dRe", &HoleCdResult::dcd_base_dRe)
      .def_readonly("dcd_base_dU1_over_Vi", &HoleCdResult::dcd_base_dU1_over_Vi)
      .def_readonly("dcd_base_dRe_D", &HoleCdResult::dcd_base_dRe_D)
      .def_readonly("inclination_dzeta", &HoleCdResult::inclination_dzeta)
      .def_readonly("re_held", &HoleCdResult::re_held)
      .def_readonly("exit_active", &HoleCdResult::exit_active)
      .def_readonly("exit_ratio", &HoleCdResult::exit_ratio)
      .def_readonly("k_perp", &HoleCdResult::k_perp)
      .def_readonly("exit", &HoleCdResult::exit)
      .def_readonly("exit_k_floor", &HoleCdResult::exit_k_floor)
      .def_readonly("exit_k_tier", &HoleCdResult::exit_k_tier)
      .def_readonly("exit_floor_above_fit", &HoleCdResult::exit_floor_above_fit)
      .def_readonly("exit_in_range", &HoleCdResult::exit_in_range)
      .def_readonly("jet_acts", &HoleCdResult::jet_acts)
      .def_readonly("jet_partials", &HoleCdResult::jet_partials)
      .def_readonly("jet_capped", &HoleCdResult::jet_capped);
  m.def("hole_cd_compose", &hole_cd_compose, py::arg("model"), py::arg("hole"),
        py::arg("orifice"), py::arg("inputs"),
        "The hole-Cd chain in one call: base Cd (fixed, ISO 5167 metering, or a\n"
        "discharge-hole correlation with its inclination and inlet-crossflow\n"
        "terms), the opt-in exit-crossflow loss (Rowbury et al. 2001), and the\n"
        "impingement jet-crossflow law (Rogers & Hersh 1976 / Florschuetz,\n"
        "Metzger & Su 1982) that replaces it where it acts. Returns the Cd, its\n"
        "analytic partials by input and each term's value. Internal: the\n"
        "composition point of the network's hole elements (_HoleCdChain).");
  py::enum_<HoleSupplyDrive>(m, "HoleSupplyDrive",
                             "The pressure a duct-fed hole's ideal flow is taken from "
                             "(hole_supply_drive).")
      .value("Static", HoleSupplyDrive::Static)
      .value("Total", HoleSupplyDrive::Total);
  py::class_<HoleSupplyDriveResult>(m, "HoleSupplyDriveResult",
                                    "A duct-fed hole's drive pressure and its partials.")
      .def_readonly("P_drive", &HoleSupplyDriveResult::P_drive)
      .def_readonly("q_dyn", &HoleSupplyDriveResult::q_dyn)
      .def_readonly("dP_dP", &HoleSupplyDriveResult::dP_dP)
      .def_readonly("dP_dTt", &HoleSupplyDriveResult::dP_dTt)
      .def_readonly("dP_dm_dot", &HoleSupplyDriveResult::dP_dm_dot)
      .def_readonly("dP_dY", &HoleSupplyDriveResult::dP_dY);
  m.def("hole_supply_drive", &hole_supply_drive, py::arg("drive"), py::arg("P"), py::arg("Tt"),
        py::arg("Y"), py::arg("m_dot"), py::arg("area"), py::arg("incompressible") = false,
        "The supply drive of a duct-fed hole (#1087): the pressure its ideal flow\n"
        "is taken from. Static: the station static P (McGreehan & Schotsch's\n"
        "convention). Total: the duct's stagnation pressure at the station,\n"
        "P + q, q the momentum-chamber closure's dynamic head at the mean duct\n"
        "flow m_dot through `area` (the plenum-fed Cd data's reference). Returns\n"
        "P_drive, q_dyn and the analytic partials of P_drive by P, Tt, m_dot, Y.");
  py::class_<HoleDriveCdRatio>(m, "HoleDriveCdRatio",
                               "Cd_T / Cd_S of a static-referenced hole Cd on the total drive.")
      .def_readonly("ratio", &HoleDriveCdRatio::ratio)
      .def_readonly("d_dP1", &HoleDriveCdRatio::d_dP1)
      .def_readonly("d_dT1", &HoleDriveCdRatio::d_dT1)
      .def_readonly("d_dPt1", &HoleDriveCdRatio::d_dPt1)
      .def_readonly("d_dTt1", &HoleDriveCdRatio::d_dTt1)
      .def_readonly("d_dP2", &HoleDriveCdRatio::d_dP2)
      .def_readonly("d_df", &HoleDriveCdRatio::d_df)
      .def_readonly("d_dY", &HoleDriveCdRatio::d_dY);
  m.def("hole_supply_recovery_cd_ratio", &hole_supply_recovery_cd_ratio, py::arg("P1"),
        py::arg("T1"), py::arg("Pt1"), py::arg("Tt1"), py::arg("P2"), py::arg("Y"), py::arg("f"),
        py::arg("incompressible") = false,
        "A hole Cd referenced to the partly recovered supply drive Pf = P1 + f q1\n"
        "(q1 = Pt1 - P1, f the recovered fraction of the supply head in [0, 1],\n"
        "HoleCd.supply_recovery) onto the total drive: Cd_T / Cd_base =\n"
        "m_id(Pf, Tf -> P2) / m_id(Pt1, Tt1 -> P2), both ideal flows from the\n"
        "network's compressible nozzle kernel; Tf the isentropic T at Pf on the\n"
        "isentrope through (P1, T1) and (Pt1, Tt1). f = 0: the static-referenced\n"
        "conversion; f = 1: ratio 1. incompressible: sqrt((Pf - P2)/(Pt1 - P2)).\n"
        "0 where Pf <= P2. Analytic partials by P1, T1, Pt1, Tt1, P2, f and Y.");
  m.def("hole_static_to_total_cd_ratio", &hole_static_to_total_cd_ratio, py::arg("P1"),
        py::arg("T1"), py::arg("Pt1"), py::arg("Tt1"), py::arg("P2"), py::arg("Y"),
        py::arg("incompressible") = false,
        "The exact conversion of a static-referenced hole Cd onto the total drive:\n"
        "Cd_T / Cd_S = m_id(P1, T1 -> P2) / m_id(Pt1, Tt1 -> P2), both ideal flows\n"
        "from the network's compressible nozzle kernel; (P1, T1) the supply's STATIC\n"
        "state, (Pt1, Tt1) its stagnation state. incompressible: the incompressible\n"
        "kernel's exact ratio sqrt((P1 - P2)/(Pt1 - P2)) (T1, Tt1, Y do not enter).\n"
        "0 where P1 <= P2. Analytic partials\n"
        "by P1, T1, Pt1, Tt1, P2 and the mass fractions Y.");

  m.def("solve_orifice_mdot",
        [](const OrificeGeometry &geom, double dP, double rho, double mu,
           double P_upstream, double kappa, MeteringCdCorrelation correlation,
           double tol, int max_iter) {
          return solve_orifice_mdot(geom, dP, rho, mu, P_upstream, kappa,
                                    correlation, tol, max_iter);
        },
        py::arg("geom"), py::arg("dP"), py::arg("rho"), py::arg("mu"),
        py::arg("P_upstream") = 101325.0, py::arg("kappa") = 0.0,
        py::arg("correlation") = MeteringCdCorrelation::ReaderHarrisGallagher,
        py::arg("tol") = 1e-6, py::arg("max_iter") = 20,
        "Iterative solver for orifice mass flow with Cd-Re coupling.\n\n"
        "Solves the coupled system:\n"
        "  mdot = eps * Cd(Re_D) * A * sqrt(2 * rho * dP)\n"
        "  Re_D = 4 * mdot / (pi * D * mu)\n\n"
        "Iterates until mdot converges (typically 3-5 iterations).\n\n"
        "Parameters:\n"
        "  geom        : OrificeGeometry\n"
        "  dP          : differential pressure [Pa]\n"
        "  rho         : upstream density [kg/m^3]\n"
        "  mu          : dynamic viscosity [Pa*s]\n"
        "  P_upstream  : absolute upstream pressure [Pa] (default: 101325)\n"
        "  kappa       : isentropic exponent cp/cv (0 = incompressible)\n"
        "  correlation : Cd correlation (default: ReaderHarrisGallagher)\n"
        "  tol         : relative convergence tolerance (default: 1e-6)\n"
        "  max_iter    : maximum iterations (default: 20)\n\n"
        "Returns: converged mass flow rate [kg/s]\n\n"
        "Raises: RuntimeError if fails to converge\n\n"
        "Example:\n"
        "  >>> geom = OrificeGeometry(d=0.05, D=0.1)\n"
        "  >>> mdot = solve_orifice_mdot(geom, 10000, 1.2, 1.8e-5)\n"
        "  >>> print(f'Mass flow: {mdot:.4f} kg/s')");

  // Orifice flow calculations with Cd
  m.def("orifice_mdot_Cd",
        static_cast<double (*)(const OrificeGeometry &, double, double, double,
                               double)>(&orifice_mdot),
        py::arg("geom"), py::arg("Cd"), py::arg("dP"), py::arg("rho"),
        py::arg("epsilon") = 1.0,
        "Mass flow through orifice given Cd.\n\n"
        "mdot = Cd * E * epsilon * A * sqrt(2 * rho * dP)\n\n"
        "Returns: mass flow rate [kg/s]");

  m.def("orifice_dP_Cd",
        static_cast<double (*)(const OrificeGeometry &, double, double, double,
                               double)>(&orifice_dP),
        py::arg("geom"), py::arg("Cd"), py::arg("mdot"), py::arg("rho"),
        py::arg("epsilon") = 1.0,
        "Pressure drop for given mass flow and Cd.\n\n"
        "dP = (mdot / (Cd * E * epsilon * A))^2 / (2 * rho)\n\n"
        "Returns: pressure drop [Pa]");

  m.def("orifice_Cd_from_measurement", &orifice_Cd_from_measurement,
        py::arg("geom"), py::arg("mdot"), py::arg("dP"), py::arg("rho"),
        "Solve for discharge coefficient from measurement.\n\n"
        "Cd = mdot / (E * A * sqrt(2 * rho * dP))\n\n"
        "Parameters:\n"
        "  geom : OrificeGeometry\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  dP   : pressure drop [Pa]\n"
        "  rho  : density [kg/m^3]\n\n"
        "Returns: discharge coefficient Cd [-]");

  m.def("expansibility_factor", &expansibility_factor, py::arg("beta"),
        py::arg("dP"), py::arg("P_upstream"), py::arg("kappa"),
        "Expansibility factor for compressible gas flow through orifices.\n\n"
        "Accounts for gas expansion as it accelerates through the orifice.\n"
        "Formula from ISO 5167-2:2003, Section 5.3.2.2.\n\n"
        "eps = 1 - (0.351 + 0.256*beta^4 + 0.93*beta^8) * [1 - (1 - tau)^(1/kappa)]\n\n"
        "where tau = dP/P_upstream, beta = d/D, kappa = cp/cv\n\n"
        "Usage: mdot = eps * Cd * A * sqrt(2 * rho_upstream * dP)\n\n"
        "Parameters:\n"
        "  beta        : diameter ratio d/D [-]\n"
        "  dP          : differential pressure [Pa]\n"
        "  P_upstream  : absolute upstream pressure [Pa]\n"
        "  kappa       : isentropic exponent cp/cv [-]\n\n"
        "Returns: expansibility factor eps [-] (0 < eps <= 1)\n\n"
        "Valid for: 0.1 <= beta <= 0.75, tau <= 0.25, ideal gas");

  // Namespace functions
  m.def("orifice_K_from_Cd", &orifice::K_from_Cd, py::arg("Cd"),
        py::arg("beta"),
        "Loss coefficient K from discharge coefficient Cd.\n\n"
        "K = (1/Cd^2 - 1) * (1 - beta^4)");

  m.def("orifice_Cd_from_K", &orifice::Cd_from_K, py::arg("K"), py::arg("beta"),
        "Discharge coefficient Cd from loss coefficient K.\n\n"
        "Cd = 1 / sqrt(1 + K / (1 - beta^4))");


  // OrificeFlowResult struct binding - Bundle of orifice flow properties
  py::class_<OrificeFlowResult>(
      m, "OrificeFlowResult",
      "Bundle of orifice flow properties computed from (geom, dP, T, P, mu, "
      "Z).\n\n"
      "All properties are read-only attributes computed in a single call.\n"
      "Provides IDE autocomplete and type safety.\n\n"
      "Attributes:\n"
      "  mdot          : Mass flow rate [kg/s]\n"
      "  v             : Velocity through orifice throat [m/s]\n"
      "  Re_D          : Channel Reynolds number (based on D) [-]\n"
      "  Re_d          : Orifice Reynolds number (based on d) [-]\n"
      "  Cd            : Discharge coefficient [-]\n"
      "  epsilon       : Expansibility factor [-] (1.0 for incompressible)\n"
      "  rho_corrected : Density corrected for compressibility [kg/m^3]")
      .def_readonly("mdot", &OrificeFlowResult::mdot, "Mass flow rate [kg/s]")
      .def_readonly("v", &OrificeFlowResult::v,
                    "Velocity through orifice throat [m/s]")
      .def_readonly("Re_D", &OrificeFlowResult::Re_D,
                    "Channel Reynolds number (based on D) [-]")
      .def_readonly("Re_d", &OrificeFlowResult::Re_d,
                    "Orifice Reynolds number (based on d) [-]")
      .def_readonly("Cd", &OrificeFlowResult::Cd, "Discharge coefficient [-]")
      .def_readonly("epsilon", &OrificeFlowResult::epsilon,
                    "Expansibility factor [-]")
      .def_readonly("rho_corrected", &OrificeFlowResult::rho_corrected,
                    "Density corrected for compressibility [kg/m^3]")
      .def("__repr__", [](const OrificeFlowResult &r) {
        return "<OrificeFlowResult: mdot=" + std::to_string(r.mdot) +
               " kg/s, "
               "v=" +
               std::to_string(r.v) +
               " m/s, "
               "Cd=" +
               std::to_string(r.Cd) +
               ", "
               "Re_d=" +
               std::to_string(r.Re_d) + ">";
      });

  m.def(
      "orifice_flow",
      [](const OrificeGeometry &geom, double dP, double T, double P, double mu,
         double Z, const std::vector<double> &X, double kappa,
         MeteringCdCorrelation correlation) {
        return orifice_flow(geom, dP, T, P, mu, Z, X, kappa, correlation);
      },
      py::arg("geom"), py::arg("dP"),
      py::arg("T"), py::arg("P"), py::arg("mu"), py::arg("Z") = 1.0,
      py::arg("X") = std::vector<double>(), py::arg("kappa") = 0.0,
      py::arg("correlation") = MeteringCdCorrelation::ReaderHarrisGallagher,
      "Compute all orifice flow properties at once with real gas "
      "correction.\n\n"
      "Convenience function that computes all orifice flow properties in a\n"
      "single call. Returns OrificeFlowResult struct with read-only "
      "attributes\n"
      "for IDE autocomplete support.\n\n"
      "Parameters:\n"
      "  geom        : OrificeGeometry\n"
      "  dP          : differential pressure [Pa]\n"
      "  T           : upstream static temperature [K]\n"
      "  P           : absolute upstream static pressure [Pa]\n"
      "  mu          : dynamic viscosity [Pa*s]\n"
      "  Z           : compressibility factor [-] (default: 1.0 = ideal gas)\n"
      "  X           : mole fractions [-] (optional, uses air if empty)\n"
      "  kappa       : isentropic exponent cp/cv [-] (default: 0.0 = "
      "incompressible)\n"
      "  correlation : Cd correlation (default: ReaderHarrisGallagher)\n\n"
      "Returns: OrificeFlowResult object with attributes:\n"
      "  mdot, v, Re_D, Re_d, Cd, epsilon, rho_corrected\n\n"
      "Real Gas Support:\n"
      "  The compressibility factor Z corrects density for non-ideal "
      "behavior:\n"
      "    rho_corrected = rho_ideal / Z\n"
      "  For ideal gas: Z = 1.0 (default)\n"
      "  For real gas: Z typically ranges 0.7-1.0 depending on P, T, "
      "composition\n"
      "  Users can obtain Z from REFPROP, Peng-Robinson, or other EOS "
      "libraries\n\n"
      "Example:\n"
      "  >>> geom = cb.OrificeGeometry()\n"
      "  >>> geom.d = 0.05  # 50 mm orifice\n"
      "  >>> geom.D = 0.1   # 100 mm channel\n"
      "  >>> result = cb.orifice_flow(geom, dP=50000, T=300, P=10e6, "
      "mu=1.1e-5, Z=0.85)\n"
      "  >>> print(f'Mass flow: {result.mdot:.3f} kg/s')\n"
      "  >>> print(f'Velocity: {result.v:.2f} m/s')\n"
      "  >>> print(f'Cd: {result.Cd:.3f}')");

  // Convenience wrapper: computes mu from thermodynamic state
  m.def(
      "orifice_flow_state",
      [](const OrificeGeometry &geom, double dP, double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double Z, double kappa, MeteringCdCorrelation correlation) {
        auto X = to_vec(X_arr);
        double mu = viscosity(T, P, X);
        return orifice_flow(geom, dP, T, P, mu, Z, X, kappa, correlation);
      },
      py::arg("geom"), py::arg("dP"), py::arg("T"), py::arg("P"), py::arg("X"),
      py::arg("Z") = 1.0, py::arg("kappa") = 0.0,
      py::arg("correlation") = MeteringCdCorrelation::ReaderHarrisGallagher,
      "Compute all orifice flow properties from thermodynamic state.\n\n"
      "Convenience wrapper around orifice_flow() that computes dynamic "
      "viscosity\n"
      "internally from (T, P, X) using the built-in transport correlations.\n\n"
      "Parameters:\n"
      "  geom        : OrificeGeometry\n"
      "  dP          : differential pressure [Pa]\n"
      "  T           : upstream static temperature [K]\n"
      "  P           : absolute upstream static pressure [Pa]\n"
      "  X           : mole fractions [-]\n"
      "  Z           : compressibility factor [-] (default: 1.0 = ideal gas)\n"
      "  kappa       : isentropic exponent cp/cv [-] (default: 0.0 = "
      "incompressible)\n"
      "  correlation : Cd correlation (default: ReaderHarrisGallagher)\n\n"
      "Returns: OrificeFlowResult with mdot, v, Re_D, Re_d, Cd, epsilon, "
      "rho_corrected\n\n"
      "Example:\n"
      "  >>> geom = cb.OrificeGeometry()\n"
      "  >>> geom.d = 0.05\n"
      "  >>> geom.D = 0.1\n"
      "  >>> X_air = cb.species.dry_air()\n"
      "  >>> result = cb.orifice_flow_state(geom, dP=10000, T=300, P=101325, "
      "X=X_air)\n"
      "  >>> print(f'Mass flow: {result.mdot:.4f} kg/s')");

  // Utility functions
  m.def("orifice_velocity_from_mdot", &orifice_velocity_from_mdot,
        py::arg("mdot"), py::arg("rho"), py::arg("d"), py::arg("Z") = 1.0,
        "Velocity through orifice from mass flow rate with real gas "
        "correction.\n\n"
        "v = mdot / (rho_corrected * A)\n"
        "where rho_corrected = rho / Z, A = pi*d^2/4\n\n"
        "Parameters:\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  rho  : density [kg/m^3]\n"
        "  d    : orifice diameter [m]\n"
        "  Z    : compressibility factor [-] (default: 1.0 = ideal gas)\n\n"
        "Returns: velocity [m/s]");

  // Aliases for backward compatibility in roughness DB
  m.def("orifice_area_from_beta", &orifice_area_from_beta, py::arg("D"),
        py::arg("beta"),
        "Orifice area from beta ratio.\n\n"
        "A = pi*(D*beta/2)^2\n\n"
        "Parameters:\n"
        "  D    : channel diameter [m]\n"
        "  beta : diameter ratio d/D [-]\n\n"
        "Returns: orifice area [m^2]");

  m.def("beta_from_diameters", &beta_from_diameters, py::arg("d"), py::arg("D"),
        "Beta ratio from diameters.\n\n"
        "beta = d / D\n\n"
        "Parameters:\n"
        "  d : orifice diameter [m]\n"
        "  D : channel diameter [m]\n\n"
        "Returns: beta ratio [-]");

  m.def("orifice_Re_d_from_mdot", &orifice_Re_d_from_mdot, py::arg("mdot"),
        py::arg("d"), py::arg("mu"),
        "Orifice Reynolds number from mass flow rate.\n\n"
        "Re_d = 4*mdot / (pi*d*mu)\n\n"
        "Parameters:\n"
        "  mdot : mass flow rate [kg/s]\n"
        "  d    : orifice diameter [m]\n"
        "  mu   : dynamic viscosity [Pa*s]\n\n"
        "Returns: Reynolds number [-]");

  // =========================================================================
  // Channel Roughness Database
  // =========================================================================

  m.def("channel_roughness", &channel_roughness, py::arg("material"),
        "Get absolute roughness for a standard channel material.\n\n"
        "Returns the absolute roughness eps [m] for common channel materials.\n"
        "Material names are case-insensitive.\n\n"
        "Parameters:\n"
        "  material : channel material name (str)\n\n"
        "Returns: absolute roughness eps [m]\n\n"
        "Available materials:\n"
        "  - 'smooth', 'drawn_tubing', 'pvc', 'plastic'\n"
        "  - 'commercial_steel', 'new_steel', 'wrought_iron'\n"
        "  - 'galvanized_iron', 'galvanized_steel', 'rusted_steel'\n"
        "  - 'cast_iron', 'asphalted_cast_iron'\n"
        "  - 'concrete', 'rough_concrete'\n"
        "  - 'riveted_steel', 'wood_stave', 'corrugated_metal'\n\n"
        "Raises: ValueError if material not found\n\n"
        "References:\n"
        "  - Moody (1944), Colebrook (1939), White (2011),\n"
        "    Munson (2013), Crane (2009)\n\n"
        "Example:\n"
        "  >>> eps = channel_roughness('commercial_steel')\n"
        "  >>> print(f'Roughness: {eps*1e6:.1f} um')  # 45.0 um");

  m.def("standard_channel_roughness", &standard_channel_roughness,
        "Get all standard channel roughness values.\n\n"
        "Returns a dictionary mapping material names to absolute\n"
        "roughness values [m].\n\n"
        "Returns: dict[str, float] - material -> roughness [m]\n\n"
        "Example:\n"
        "  >>> roughness_db = standard_channel_roughness()\n"
        "  >>> for material, eps in sorted(roughness_db.items()):\n"
        "  ...     print(f'{material:20s}: {eps*1e6:8.2f} um')");

  // =========================================================================
  // Materials Database
  // =========================================================================

  // Material thermal conductivity functions
  // Range-checked by default (a ValueError outside the source's range);
  // extrapolate=True opts into the database's linear continuation of the end
  // slope (materials::material_conductivity).
  const auto material_k = [](double (*fn)(double), const char *db_name) {
    return [fn, db_name](double T, bool extrapolate) {
      return extrapolate ? combaero::materials::material_conductivity(db_name, T, true) : fn(T);
    };
  };
  m.def("k_inconel718", material_k(&combaero::materials::k_inconel718, "inconel718"), py::arg("T"), py::arg("extrapolate") = false,
        "Thermal conductivity of Inconel 718 [W/(m*K)].\n\n"
        "Ni-based superalloy commonly used in turbine applications.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 300-1200 K\n"
        "Source: Haynes International, Special Metals Corporation\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");

  m.def("k_haynes230", material_k(&combaero::materials::k_haynes230, "haynes230"), py::arg("T"), py::arg("extrapolate") = false,
        "Thermal conductivity of Haynes 230 [W/(m*K)].\n\n"
        "Ni-Cr-W-Mo alloy with high-temperature oxidation resistance.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 298-1273 K (the brochure table ends at 1000 C)\n"
        "Source: HAYNES 230 brochure H-3000O (Haynes International, 2021):\n"
        "k = alpha * cp * rho from its diffusivity and specific-heat columns\n"
        "(#881). Its printed conductivity column is 13% above that at 500 C and\n"
        "outside Gen3 CSP's measured band; the product is inside it.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");

  m.def("k_hastelloy_x", material_k(&combaero::materials::k_hastelloy_x, "hastelloy_x"), py::arg("T"), py::arg("extrapolate") = false,
        "Hastelloy X thermal conductivity [W/(m*K)] from the HASTELLOY X brochure\n"
        "and data sheet (Haynes International), interpolated linearly. Valid 300-1273 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_haynes188", material_k(&combaero::materials::k_haynes188, "haynes188"), py::arg("T"), py::arg("extrapolate") = false,
        "Haynes 188 thermal conductivity [W/(m*K)] from the HAYNES 188 brochure\n"
        "table (Haynes International), interpolated linearly. Valid 300-1273 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_haynes263", material_k(&combaero::materials::k_haynes263, "haynes263"), py::arg("T"), py::arg("extrapolate") = false,
        "Haynes 263 thermal conductivity [W/(m*K)] from the HAYNES 263 brochure\n"
        "table (Haynes International), interpolated linearly. Valid 300-1273 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_haynes25", material_k(&combaero::materials::k_haynes25, "haynes25"), py::arg("T"), py::arg("extrapolate") = false,
        "Haynes 25 (L-605) thermal conductivity [W/(m*K)] from the HAYNES 25\n"
        "brochure table (Haynes International), interpolated linearly. Valid 298-1273 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_inconel625", material_k(&combaero::materials::k_inconel625, "inconel625"), py::arg("T"), py::arg("extrapolate") = false,
        "Inconel 625 thermal conductivity [W/(m*K)] from the Special Metals\n"
        "INCONEL alloy 625 Table 3 (Battelle, annealed), interpolated linearly.\n"
        "Valid 116-1255 K. Gen3 CSP (Georgia Tech) measures 14-25% more at\n"
        "260-700 C, from a higher diffusivity; kept until a third source settles\n"
        "it (#881).\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_nimonic75", material_k(&combaero::materials::k_nimonic75, "nimonic75"), py::arg("T"), py::arg("extrapolate") = false,
        "Nimonic 75 thermal conductivity [W/(m*K)] from the Virgamet Alloy N75 table\n"
        "(20 C value = Special Metals sheet), interpolated linearly. Valid 293-1273 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_copper", material_k(&combaero::materials::k_copper, "copper"), py::arg("T"), py::arg("extrapolate") = false,
        "High-purity copper thermal conductivity [W/(m*K)] from NSRDS-NBS 8\n"
        "Table 2b, interpolated linearly. Valid 100-1200 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_alumina", material_k(&combaero::materials::k_alumina, "alumina"), py::arg("T"), py::arg("extrapolate") = false,
        "Dense Al2O3 (99.5% pure, 98% dense) thermal conductivity [W/(m*K)] from\n"
        "NSRDS-NBS 8 Table 14, interpolated linearly. Valid 100-2000 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_x22crmov12_1", material_k(&combaero::materials::k_x22crmov12_1, "x22crmov12_1"), py::arg("T"), py::arg("extrapolate") = false,
        "X22CrMoV12-1 (1.4923) thermal conductivity [W/(m*K)]: linear between the\n"
        "data sheet's 20 C and 650 C values (Edelstahlwerke Suedwestfalen). Valid 293-923 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_bondcoat_nicocraly", material_k(&combaero::materials::k_bondcoat_nicocraly, "nicocraly"), py::arg("T"), py::arg("extrapolate") = false,
        "NiCoCrAlY bond coat thermal conductivity [W/(m*K)]: Zhu & Miller 1999\n"
        "(NASA/TM-1999-209069) Fig. 2, k = 6.729572 + 0.02344 T_C. Valid 848-1373 K.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");
  m.def("k_stainless_steel_316", material_k(&combaero::materials::k_stainless_steel_316, "ss316"),
        py::arg("T"), py::arg("extrapolate") = false,
        "Thermal conductivity of Stainless Steel 316 [W/(m*K)].\n\n"
        "Austenitic stainless steel with corrosion resistance.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 533-1093 K (260-820 C, the measured range)\n"
        "Source: Georgia Tech Gen3 CSP thermophysical property database, measured\n"
        "k with a 95% band of 0.9-2.7 W/(m*K) (#881), interpolated linearly.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");

  m.def("k_stainless_steel_310", material_k(&combaero::materials::k_stainless_steel_310, "ss310"),
        py::arg("T"), py::arg("extrapolate") = false,
        "Thermal conductivity of Stainless Steel 310 [W/(m*K)].\n\n"
        "Austenitic 25Cr-20Ni heat-resisting stainless steel.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 533-1073 K (260-800 C, the measured range)\n"
        "Source: Georgia Tech Gen3 CSP thermophysical property database, measured\n"
        "k with a 95% band of 0.8-2.5 W/(m*K) (#881), interpolated linearly. Not\n"
        "joined to a room-temperature value: NIST's cryogenic fit (11.6 at 300 K)\n"
        "is a different sample.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");

  m.def("k_stainless_steel_304", material_k(&combaero::materials::k_stainless_steel_304, "ss304"),
        py::arg("T"), py::arg("extrapolate") = false,
        "Thermal conductivity of Stainless Steel 304 (and 304L) [W/(m*K)].\n\n"
        "Austenitic 18-8 stainless steel.\n\n"
        "Parameters:\n"
        "  T : temperature [K]\n\n"
        "Valid range: 200-1672 K (to the solidus, the table's melting range 1672-1727 K)\n"
        "Source: Ho & Chu 1977, CINDAS Report 45 (DTIC ADA129160) Table 2,\n"
        "recommended values (uncertainty about +-5%), interpolated linearly.\n\n"
        "extrapolate : bool, default False. True continues the end slope linearly\n"
        "outside the valid range instead of raising (unsourced; opt-in).");

  m.def(
      "k_tbc_ysz", &combaero::materials::k_tbc_ysz, py::arg("T"),
      py::arg("hours") = 0.0, py::arg("is_ebpvd") = false,
      "Thermal conductivity of a ZrO2-8wt%Y2O3 thermal barrier coating [W/(m*K)].\n\n"
      "Plasma-sprayed (default): Zhu & Miller 2000, NASA/TM-2000-210237. About\n"
      "1.0 W/(m*K) as sprayed, flat in T (p. 6); sintering at T raises it by\n"
      "their Eq. 1, (k - k0)/(k_inf - k0) = A exp(-Q/RT) [1 - exp(-t/tau)],\n"
      "with k0 = 1.0 and k_inf = 1.74. This is the overall conductivity of a\n"
      "254 um coating under a 64 W/cm^2 gradient at surface temperature T: an\n"
      "effective through-thickness value.\n\n"
      "Fidelity on the paper's own Fig. 4 (32 h): matches the 1320 C curve and\n"
      "runs 2% high at 1100 C and 3% high at 990 C. Eq. 1's printed A and Q\n"
      "give a weaker temperature dependence than the figure; the paper's own\n"
      "fitted curves are not a least-squares fit of its data either. Eq. 1 is\n"
      "used as printed.\n\n"
      "EB-PVD (is_ebpvd=True): as deposited only, 1.5 + 0.0002 T_C, uncited and\n"
      "without sintering until a public source is held.\n\n"
      "Parameters:\n"
      "  T        : surface temperature [K]\n"
      "  hours    : time at temperature [h] (default: 0 = as sprayed)\n"
      "  is_ebpvd : True for EB-PVD, False for plasma-sprayed (default)\n\n"
      "Valid range: 300-1700 K\n\n"
      "Example:\n"
      "  >>> k_fresh = cb.k_tbc_ysz(T=1500)  # as sprayed, 1.0\n"
      "  >>> k_aged = cb.k_tbc_ysz(T=1500, hours=32)  # sintered");

  m.def("list_materials", &combaero::materials::list_materials,
        "List all available materials in the database.\n\n"
        "Returns: List of material names (strings)\n\n"
        "Example:\n"
        "  >>> materials = cb.list_materials()\n"
        "  >>> print(materials)\n"
        "  ['alumina', 'copper', 'hastelloy_x', 'haynes188', 'haynes230',\n"
        "   'haynes25', 'haynes263', 'inconel625', 'inconel718', 'nicocraly', 'nimonic75',\n"
        "   'ss304', 'ss310', 'ss316', 'x22crmov12_1', 'ysz']");

  m.def(
      "get_material_conductivity",
      [](const std::string &name, double T, bool extrapolate) {
        // Outside the valid range: the #882 envelope by default, the linear
        // end-slope continuation on request.
        const double k = combaero::materials::material_conductivity(name, T, extrapolate);
        return std::max(1e-3, k);
      },
      py::arg("name"), py::arg("T"), py::arg("extrapolate") = false,
      "Get thermal conductivity [W/(m*K)] for a named material at temperature T "
      "[K].\n"
      "Outside the material's sourced range, extrapolate=False (default) follows\n"
      "the envelope (#882): the end slope ramped linearly to 0 over 100 K, then\n"
      "flat, from 200 K (or T_min, if lower) up to the material's solidus (the\n"
      "generic metal cap 1533.15 K where none is printed; T_max for alumina and\n"
      "ysz), and held beyond. extrapolate=True continues the end slope\n"
      "linearly instead. Neither is sourced: material_in_range(name, T) says\n"
      "which part is data. Raises if k reaches <= 0; floored at 1e-3. Raises\n"
      "ValueError for an unknown material or a non-finite T.");

  m.def(
      "get_material_conductivity_and_slope",
      [](const std::string &name, double T, bool extrapolate) {
        const auto r = combaero::materials::material_conductivity_and_slope(name, T, extrapolate);
        // The same floor as get_material_conductivity; flat where it acts.
        if (!(r.k > 1e-3)) {
          return std::make_tuple(1e-3, 0.0);
        }
        return std::make_tuple(r.k, r.dk_dT);
      },
      py::arg("name"), py::arg("T"), py::arg("extrapolate") = false,
      "(k [W/(m*K)], dk/dT [W/(m*K^2)]) of a named material at T [K] (#932):\n"
      "get_material_conductivity's k, bit for bit, and its analytic slope on\n"
      "every path -- the data's (a table's segment slope, the left one at a\n"
      "printed point), the #882 envelope's fade, 0 where it holds or the\n"
      "1e-3 floor acts, the end slope with extrapolate=True.");

  // The hot end #882 gives a metal whose source prints no solidus: the lowest
  // sourced solidus in the database. An effusion plate with no material
  // named is judged against it (#998).
  m.attr("GENERIC_METAL_SOLIDUS") = combaero::materials::GENERIC_METAL_SOLIDUS;

  m.def("material_in_range", &combaero::materials::material_in_range, py::arg("name"),
        py::arg("T"),
        "Whether T [K] lies in the material's sourced conductivity range\n"
        "[T_min, T_max] (#882). Outside it, get_material_conductivity follows the\n"
        "envelope. Raises ValueError for an unknown material or a non-finite T.");

  m.def(
      "get_material_info",
      [](const std::string &name) {
        const auto &mat = combaero::materials::get_material(name);
        py::dict d;
        d["name"] = mat.name;
        d["T_min"] = mat.T_min;
        d["T_max"] = mat.T_max;
        d["source"] = mat.source;
        d["notes"] = mat.notes;
        d["T_solidus"] = mat.T_solidus;
        d["solidus_source"] = mat.solidus_source;
        d["T_envelope_min"] = combaero::materials::envelope_T_min(mat);
        d["T_envelope_max"] = combaero::materials::envelope_T_max(mat);
        return d;
      },
      py::arg("name"),
      "A database material's record: its sourced conductivity range T_min..T_max\n"
      "[K] and source, notes, and its solidus T_solidus [K] -- the lower end of\n"
      "the melting range its source prints -- with that source (#882). T_solidus\n"
      "is NaN where no held source prints one. T_envelope_min/T_envelope_max [K]\n"
      "are the ends of its k(T) envelope. Raises ValueError for an unknown\n"
      "material.");

  // =========================================================================
  // Advanced Cooling Correlations
  // =========================================================================

  // Rib-enhanced cooling



  m.def("thermal_performance_factor",
        &combaero::cooling::thermal_performance_factor, py::arg("Nu_ratio"),
        py::arg("f_ratio"),
        "Thermal-hydraulic performance factor (Webb & Eckert 1972).\n\n"
        "eta = (Nu/Nu0) / (f/f0)^(1/3)\n"
        "Compares heat transfer gain against pumping power penalty at equal "
        "velocity.\n"
        "eta > 1: net improvement over smooth; eta < 1: pressure-drop "
        "generator.\n"
        "Generic: works for ribs, dimples, pin fins, impingement, or any "
        "surface.\n\n"
        "Parameters:\n"
        "  Nu_ratio: Nu_enhanced / Nu_smooth [-]\n"
        "  f_ratio : f_enhanced  / f_smooth  [-]\n\n"
        "Source: Webb & Eckert (1972) Int. J. Heat Mass Transfer 15(8), "
        "1647-1658\n\n"
        "Returns: thermal-hydraulic performance factor eta [-]");


  // Impingement cooling




  m.def("adiabatic_wall_temperature",
        &combaero::cooling::adiabatic_wall_temperature, py::arg("T_hot"),
        py::arg("T_coolant"), py::arg("eta"),
        "Adiabatic wall temperature from cooling effectiveness.\n\n"
        "Uses: eta = (T_hot - T_aw) / (T_hot - T_coolant).\n\n"
        "Parameters:\n"
        "  T_hot     : mainstream hot-gas temperature [K]\n"
        "  T_coolant : coolant supply temperature [K]\n"
        "  eta       : adiabatic effectiveness [-]\n\n"
        "Returns: adiabatic wall temperature T_aw [K]");

  py::class_<combaero::cooling::EffusionPanelFilm>(m, "EffusionPanelFilm")
      .def_readonly("eta", &combaero::cooling::EffusionPanelFilm::eta)
      .def_readonly("extrapolated", &combaero::cooling::EffusionPanelFilm::extrapolated);
  m.def("effusion_panel_film_effectiveness",
        &combaero::cooling::effusion_panel_film_effectiveness, py::arg("n_rows"),
        py::arg("pitch_x_over_D"), py::arg("s_over_D"), py::arg("M"),
        py::arg("density_ratio"), py::arg("alpha_deg"), py::arg("Tu"),
        "Panel-averaged adiabatic film effectiveness of an effusion plate:\n"
        "Baldauf (2002) per row, Sellers-superposed, averaged over the rows.\n"
        "The alternative closure Andrews 88-GT-290 refuses offered alone --\n"
        "select it explicitly; extrapolated flags Baldauf's envelope.");

  py::class_<combaero::cooling::FilmSpanAverage>(m, "FilmSpanAverage")
      .def_readonly("eta", &combaero::cooling::FilmSpanAverage::eta)
      .def_readonly("deta_dM", &combaero::cooling::FilmSpanAverage::deta_dM)
      .def_readonly("deta_dP", &combaero::cooling::FilmSpanAverage::deta_dP)
      .def_readonly("extrapolated", &combaero::cooling::FilmSpanAverage::extrapolated);
  m.def("film_span_average_baldauf_2002",
        &combaero::cooling::film_span_average_baldauf_2002, py::arg("n_rows"),
        py::arg("row_pitch_over_D"), py::arg("x_start_over_D"),
        py::arg("x_end_over_D"), py::arg("M"), py::arg("density_ratio"),
        py::arg("alpha_deg"), py::arg("s_over_D"), py::arg("Tu"),
        "Film a downstream wall sees from rows of CYLINDRICAL holes: Baldauf\n"
        "(2002) per row, Sellers-superposed, averaged over [x_start, x_end]\n"
        "downstream of the last row (in D). Returns eta, d eta/dM, d eta/dP\n"
        "(density ratio) and Baldauf's envelope flag.");
  py::class_<combaero::cooling::FilmSpanAverageGroups>(m, "FilmSpanAverageGroups")
      .def_readonly("eta", &combaero::cooling::FilmSpanAverageGroups::eta)
      .def_readonly("weight", &combaero::cooling::FilmSpanAverageGroups::weight)
      .def_readonly("deta_dM", &combaero::cooling::FilmSpanAverageGroups::deta_dM)
      .def_readonly("deta_dP", &combaero::cooling::FilmSpanAverageGroups::deta_dP)
      .def_readonly("d_weight_dM", &combaero::cooling::FilmSpanAverageGroups::d_weight_dM)
      .def_readonly("d_weight_dP", &combaero::cooling::FilmSpanAverageGroups::d_weight_dP)
      .def_readonly("extrapolated", &combaero::cooling::FilmSpanAverageGroups::extrapolated);
  m.def("film_span_average_baldauf_2002_groups",
        &combaero::cooling::film_span_average_baldauf_2002_groups, py::arg("n_rows"),
        py::arg("row_pitch_over_D"), py::arg("x_start_over_D"),
        py::arg("x_end_over_D"), py::arg("M"), py::arg("density_ratio"),
        py::arg("alpha_deg"), py::arg("s_over_D"), py::arg("Tu"),
        "Film a downstream wall sees from GROUPS of cylindrical-hole rows (an\n"
        "effusion liner's panels, upstream first), each group at its own M and\n"
        "density ratio: Baldauf (2002) per row, layered row by row so each\n"
        "group keeps its own coolant temperature. T_aw,bar = T_ref (1 - eta)\n"
        "+ sum_g weight[g] T_c,g; derivatives per group.");
  m.def("film_span_average_colban_2011",
        &combaero::cooling::film_span_average_colban_2011, py::arg("n_rows"),
        py::arg("row_pitch_over_D"), py::arg("x_start_over_D"),
        py::arg("x_end_over_D"), py::arg("M"), py::arg("area_ratio"),
        py::arg("P_over_D"), py::arg("t_over_P"), py::arg("alpha_deg") = 30.0,
        py::arg("density_ratio") = 1.7,
        "Film a downstream wall sees from rows of SHAPED holes: Colban (2011)\n"
        "per row, Sellers-superposed, averaged over [x_start, x_end] downstream\n"
        "of the last breakout (in metering D). deta_dP is 0: Eq. (19) has no\n"
        "density ratio.");

  m.def("slot_film_span_average_marek_tacina_1975",
        &combaero::cooling::slot_film_span_average_marek_tacina_1975,
        py::arg("x_start_over_S"), py::arg("x_end_over_S"), py::arg("M"),
        py::arg("Tu"), py::arg("C_m_multiplier") = 1.0,
        "Film a downstream wall sees from a 2-D SLOT: Marek & Tacina (1975)\n"
        "Eqs. (4)-(5), averaged over [x_start, x_end] downstream of the slot\n"
        "lip (in slot heights S), Tu a FRACTION (0.14 = 14%). deta_dP is 0:\n"
        "Eq. (4) has no density ratio.");

  m.def("film_superposition_sellers",
        &combaero::cooling::film_superposition_sellers, py::arg("eta_rows"),
        "Sellers superposition of per-row film effectiveness:\n"
        "    eta = 1 - prod_i (1 - eta_i)\n\n"
        "Assumes the rows are independent. Gao et al. (2025) measure what\n"
        "that costs: the error accumulates with row count and always\n"
        "OVERESTIMATES. Use film_superposition_corrected for many rows.");

  m.def("film_superposition_sellers_and_gradient",
        &combaero::cooling::film_superposition_sellers_and_gradient,
        py::arg("eta_rows"),
        "(eta, [d eta / d eta_i]) for Sellers superposition. The partials\n"
        "are formed from prefix/suffix products rather than by division, so\n"
        "a fully effective row does not produce a NaN.");

  m.def("film_superposition_corrected",
        &combaero::cooling::film_superposition_corrected,
        py::arg("eta_rows"), py::arg("alpha_between_rows"),
        "Gao et al. (2025) Eq. (7): Sellers with a per-row mainstream\n"
        "temperature correction.\n\n"
        "    eta = sum_i [ eta_i prod_{j>=i} alpha_j prod_{k>i} (1 - eta_k) ]\n\n"
        "alpha_between_rows has one FEWER entry than eta_rows: alpha_j is\n"
        "the fraction of the film's temperature deficit that survives\n"
        "mixing between row j and row j+1. All ones reproduces Sellers\n"
        "exactly.\n\n"
        "alpha is the tuning knob for matching a rig. Gao publishes its\n"
        "functional form (see mainstream_temperature_correction) but NOT\n"
        "the fitted coefficients, so there is no defensible default beyond\n"
        "1.0, which is plain Sellers.");

  m.def("film_superposition_corrected_and_gradient",
        &combaero::cooling::film_superposition_corrected_and_gradient,
        py::arg("eta_rows"), py::arg("alpha_between_rows"),
        "(eta, [d eta / d eta_i]) for Gao Eq. (7).");

  m.def("mainstream_temperature_correction",
        &combaero::cooling::mainstream_temperature_correction,
        py::arg("mass_flow_ratio"), py::arg("a"), py::arg("b"),
        "Gao et al. (2025) Eq. (5), the published FORM of the mainstream\n"
        "temperature correction:\n\n"
        "    alpha = a r / (a r + 1) + b,    r = m_coolant / m_mainstream\n\n"
        "Derived from an energy balance on mainstream entrained into the\n"
        "boundary layer at each injection (Eqs. 3-4), so alpha is the\n"
        "fraction of the film's temperature deficit surviving to the next\n"
        "row.\n\n"
        "a and b are REQUIRED, not defaulted. The paper does print them\n"
        "(12 and 0.9465, section 4.3.1), but a r/(a r + 1) >= 0 means that\n"
        "pair gives alpha >= 0.9465 for EVERY r, so it cannot reach the\n"
        "0.69 to 0.85 per-row damping a tighter-pitched plate needs. Eq. (5)\n"
        "carries no streamwise-spacing term, and that is the variable\n"
        "Murray and Ireland (2018) show to dominate. Supply your own\n"
        "fitted values.");

  m.def("mills_entry_length_factor",
        &combaero::cooling::mills_entry_length_factor, py::arg("L_over_D"),
        "Mills' entry-length factor R_Nu for a sharp-edged tube entry.\n\n"
        "Andrews 86-GT-225 Eqs. (13) and (14), two curve fits meeting at\n"
        "L/D = 2. A short hole never develops, so R_Nu > 1; it decays to 1\n"
        "as the hole lengthens.\n\n"
        "  L_over_D : hole length / hole diameter [-]\n\n"
        "Returns: R_Nu [-]");

  m.def("effusion_approach_nusselt",
        &combaero::cooling::effusion_approach_nusselt, py::arg("Re"),
        py::arg("Pr"), py::arg("X_over_L"),
        "Hole-approach contribution to effusion internal heat transfer.\n\n"
        "Sparrow's multi-hole correlation as Andrews 86-GT-225 Eq. (18)\n"
        "rebases it onto the hole's own Nusselt definition:\n\n"
        "  Nu = 0.881 Re^0.476 Pr^(1/3) X/(pi L)\n\n"
        "  Re       : on the HOLE diameter [-]\n"
        "  Pr       : coolant Prandtl number [-]\n"
        "  X_over_L : hole pitch / hole length [-]\n\n"
        "Returns: Nu on the hole diameter and hole internal area [-]");

  m.def("effusion_throat_nusselt",
        &combaero::cooling::effusion_throat_nusselt, py::arg("Re"),
        py::arg("Pr"), py::arg("L_over_D"),
        "Throat contribution to effusion internal heat transfer.\n\n"
        "Mills' short-hole data, Andrews 86-GT-225 Eq. (12):\n\n"
        "  Nu = 0.023 Re^0.8 Pr^(1/3) R_Nu(L/D)\n\n"
        "Returns: Nu on the hole diameter and hole internal area [-]");

  m.def("effusion_internal_nusselt",
        &combaero::cooling::effusion_internal_nusselt, py::arg("Re"),
        py::arg("Pr"), py::arg("X_over_L"), py::arg("L_over_D"),
        "Effusion plate internal Nusselt number, Andrews 86-GT-225\n"
        "Eq. (19): the hole-approach and throat terms summed.\n\n"
        "  Nu = (0.881 (X/pi L) Re^0.476 + 0.023 Re^0.8 R_Nu) Pr^(1/3)\n\n"
        "Referenced to the hole diameter and the HOLE INTERNAL AREA. For a\n"
        "coefficient per unit PLATE area, which is what an effusion element\n"
        "needs, multiply by k/D and then by A_h/A with A_h = pi D L and\n"
        "A = X^2 - pi D^2/4. That factor is 3.46 for Andrews' plate C, so\n"
        "it is not optional.\n\n"
        "Returns: Nu [-]");

  m.def("flat_plate_turbulent_nusselt_and_derivative",
        [](double Re_x, double Pr) {
          return combaero::cooling::flat_plate_turbulent_nusselt_and_derivative(
              Re_x, Pr, nullptr);
        },
        py::arg("Re_x"), py::arg("Pr"),
        "Local turbulent flat-plate Nusselt number on x, and its derivative.\n\n"
        "The unblown reference of the Stanford full-coverage film-cooling\n"
        "data (NASA CR-2648, printed under every table):\n"
        "St_o Pr^0.4 = 0.0295 Re_x^-0.2, so Nu_x = 0.0295 Re_x^0.8 Pr^0.6.\n\n"
        "  Re_x : x-Reynolds number from where the boundary layer starts [-]\n"
        "  Pr   : gas Prandtl number [-]\n\n"
        "Data span Re_x 1.5e5-5e6, air (see flat_plate_turbulent_panel_nusselt\n"
        "for the flag).\n\n"
        "Returns: (Nu_x, dNu_x/dRe_x)");

  py::class_<combaero::cooling::FlatPlatePanelNusselt>(m, "FlatPlatePanelNusselt")
      .def_readonly("Nu", &combaero::cooling::FlatPlatePanelNusselt::Nu)
      .def_readonly("dNu_dRe_L", &combaero::cooling::FlatPlatePanelNusselt::dNu_dRe_L)
      .def_readonly("Re_x_start", &combaero::cooling::FlatPlatePanelNusselt::Re_x_start)
      .def_readonly("Re_x_end", &combaero::cooling::FlatPlatePanelNusselt::Re_x_end)
      .def_readonly("extrapolated", &combaero::cooling::FlatPlatePanelNusselt::extrapolated);

  m.def("flat_plate_turbulent_panel_nusselt",
        &combaero::cooling::flat_plate_turbulent_panel_nusselt, py::arg("Re_L"),
        py::arg("Pr"), py::arg("x_start_over_L"),
        "Panel-averaged turbulent flat-plate Nusselt number over\n"
        "x in [x0, x0 + L], on the panel length L.\n\n"
        "  Nu_L = (0.0295/0.8) Pr^0.6 Re_L^0.8 [(1 + s)^0.8 - s^0.8], s = x0/L\n\n"
        "CR-2648's St_o (the Stanford full-coverage unblown reference)\n"
        "integrated in closed form; no unheated starting length.\n\n"
        "  Re_L           : rho U L / mu on the panel length [-], > 0\n"
        "  Pr             : gas Prandtl number [-], > 0\n"
        "  x_start_over_L : boundary-layer start upstream of the panel's\n"
        "                   leading edge, in panel lengths [-], >= 0\n\n"
        "Returns FlatPlatePanelNusselt: Nu, dNu_dRe_L, Re_x_start, Re_x_end,\n"
        "extrapolated (the panel's Re_x leaves the data's 1.5e5-5e6 span; it\n"
        "is evaluated anyway, never refused). h_bar = Nu k / L.");

  m.def("flat_plate_turbulent_panel_nusselt_and_derivative",
        [](double Re_L, double Pr, double x_start_over_L) {
          return combaero::cooling::flat_plate_turbulent_panel_nusselt_and_derivative(
              Re_L, Pr, x_start_over_L, nullptr);
        },
        py::arg("Re_L"), py::arg("Pr"), py::arg("x_start_over_L"),
        "Solver-facing (f, J) form: (Nu_L, dNu_L/dRe_L) at fixed x0/L.\n"
        "Closed form, 0.8 Nu_L / Re_L.");

  m.def("equivalent_slot_width", &combaero::cooling::equivalent_slot_width,
        py::arg("hole_area"), py::arg("pitch"),
        "Gao Eq. (9): s = A_hole / pitch. The equivalent slot width that\n"
        "lets one row's measured distribution stand in for another spacing\n"
        "through the X/(M s) scaling of Eq. (8).");

  m.def("equivalent_blowing_ratio",
        &combaero::cooling::equivalent_blowing_ratio,
        py::arg("M_baseline"), py::arg("area_baseline"),
        py::arg("area_equivalent"),
        "Gao Eq. (10): M_e = M_0 A_0 / A_e. Normalises a hole count and\n"
        "spacing onto the baseline single-row configuration.");

  m.def("film_effectiveness_baldauf_2002",
        [](double x_over_D, double M, double P, double alpha_deg,
           double s_over_D, double Tu, std::optional<double> b_0) {
          return combaero::cooling::film_effectiveness_baldauf_2002(
              x_over_D, M, P, alpha_deg, s_over_D, Tu, nullptr,
              b_0 ? *b_0 : std::numeric_limits<double>::quiet_NaN());
        },
        py::arg("x_over_D"), py::arg("M"), py::arg("P"), py::arg("alpha_deg"),
        py::arg("s_over_D"), py::arg("Tu"), py::arg("b_0") = py::none(),
        "Laterally averaged adiabatic film-cooling effectiveness downstream\n"
        "of ONE row of cylindrical, streamwise-inclined holes.\n\n"
        "Baldauf, Scheurlen, Schulz & Wittig (2002), ASME J. Turbomachinery\n"
        "124(4), 686-698. Valid from the ejection point to far downstream,\n"
        "and it carries the adjacent jet interaction (lateral hole spacing,\n"
        "jet lift-off) rather than excluding it.\n\n"
        "eta is (T_G - T_AW)/(T_G - T_C), the same convention as\n"
        "adiabatic_wall_temperature(), so the two compose directly.\n\n"
        "Parameters:\n"
        "  x_over_D  : distance downstream of ejection, in hole diameters\n"
        "  M         : blowing rate (rho u)_C/(rho u)_G [-], 0.2 to 2.5\n"
        "  P         : density ratio rho_C/rho_G [-], 1.2 to 1.8\n"
        "  alpha_deg : ejection angle to the SURFACE [deg], 30 to 90\n"
        "  s_over_D  : lateral hole spacing / diameter [-], 2 to 5\n"
        "  Tu        : mainstream turbulence intensity [-], 0.0035 to 0.075\n\n"
        "Alpha is given in DEGREES and converted internally; the paper's\n"
        "trigonometry is in radians.\n\n"
        "NOTE Eq. (31) is implemented as printed and disagrees with the\n"
        "paper's own Table 4 by 36% in b_0. Under 5% effect below M ~ 0.5,\n"
        "up to 50% at M = 2.5. See the extraction record.\n\n"
        "Returns: laterally averaged eta [-]");

  m.def("film_effectiveness_baldauf_2002_and_derivatives",
        [](double x_over_D, double M, double P, double alpha_deg,
           double s_over_D, double Tu, std::optional<double> b_0) {
          return combaero::cooling::film_effectiveness_baldauf_2002_and_derivatives(
              x_over_D, M, P, alpha_deg, s_over_D, Tu, nullptr,
              b_0 ? *b_0 : std::numeric_limits<double>::quiet_NaN());
        },
        py::arg("x_over_D"), py::arg("M"), py::arg("P"), py::arg("alpha_deg"),
        py::arg("s_over_D"), py::arg("Tu"), py::arg("b_0") = py::none(),
        "Solver-facing (f, J) form: (eta, d eta/dM, d eta/dP).\n\n"
        "M and P are what a network solve varies -- M through the coolant\n"
        "mass flow, P through the temperature ratio. Geometry and Tu are\n"
        "fixed per element and carry no partials. Analytic throughout, via\n"
        "forward-mode dual numbers over the same equation chain the value\n"
        "uses, so the two cannot drift apart.");

  m.def("film_effectiveness_colban_2011",
        [](double x_over_D, double M, double area_ratio, double P_over_D,
           double t_over_P, double alpha_deg, double density_ratio) {
          return combaero::cooling::film_effectiveness_colban_2011(
              x_over_D, M, area_ratio, P_over_D, t_over_P, alpha_deg,
              density_ratio, nullptr);
        },
        py::arg("x_over_D"), py::arg("M"), py::arg("area_ratio"),
        py::arg("P_over_D"), py::arg("t_over_P"), py::arg("alpha_deg") = 30.0,
        py::arg("density_ratio") = 1.7,
        "Laterally averaged adiabatic film-cooling effectiveness downstream\n"
        "of ONE row of SHAPED (fan or laidback fan) holes, Eq. (19).\n\n"
        "Colban, Thole & Bogard (2011), ASME J. Turbomachinery 133,\n"
        "011002. eta = 1/(P/t + C1 M^C2 xi^C3), xi = (4/pi)(x/D)(P/D)/(M AR);\n"
        "eta = t/P at the breakout. Same eta convention as\n"
        "adiabatic_wall_temperature().\n\n"
        "Parameters:\n"
        "  x_over_D      : distance downstream of the breakout, in METERING\n"
        "                  diameters, >= 0\n"
        "  M             : blowing ratio at the metering section, 0.2 to 2.5\n"
        "  area_ratio    : A_exit/A_inlet, >= 1 (2.5 to 4.7 fitted)\n"
        "  P_over_D      : lateral pitch / metering diameter (4 to 8 fitted)\n"
        "  t_over_P      : breakout width / pitch, 0.31 to 0.65\n"
        "  alpha_deg     : surface angle [deg]; range check only, fitted at 30\n"
        "  density_ratio : rho_jet/rho_ext; range check only, 1.7 to 2.0\n\n"
        "Outside Table 4's limits (or the fitted angle, density ratio, AR\n"
        "and P/D spans) it warns and extrapolates.\n\n"
        "Returns: laterally averaged eta [-]");

  m.def("film_effectiveness_colban_2011_and_derivatives",
        [](double x_over_D, double M, double area_ratio, double P_over_D,
           double t_over_P, double alpha_deg, double density_ratio) {
          return combaero::cooling::film_effectiveness_colban_2011_and_derivatives(
              x_over_D, M, area_ratio, P_over_D, t_over_P, alpha_deg,
              density_ratio, nullptr);
        },
        py::arg("x_over_D"), py::arg("M"), py::arg("area_ratio"),
        py::arg("P_over_D"), py::arg("t_over_P"), py::arg("alpha_deg") = 30.0,
        py::arg("density_ratio") = 1.7,
        "Solver-facing (f, J) form: (eta, d eta/dM). M is the only input a\n"
        "network solve varies; the derivative is closed form.");

  m.def("slot_film_effectiveness_marek_tacina_1975",
        [](double x_over_S, double M, double Tu, double C_m_multiplier) {
          return combaero::cooling::slot_film_effectiveness_marek_tacina_1975(
              x_over_S, M, Tu, C_m_multiplier, nullptr);
        },
        py::arg("x_over_S"), py::arg("M"), py::arg("Tu"),
        py::arg("C_m_multiplier") = 1.0,
        "Film-cooling effectiveness downstream of a 2-D SLOT (a liner's\n"
        "cooling ring), Juhasz-Marek mixing with C_m from turbulence.\n\n"
        "Marek & Tacina (1975), NASA TN D-7958, Eqs. (4) and (5):\n"
        "eta = 1/(1 + C_m x/(M S)), C_m = 0.25 Tu (Tu a fraction; the report\n"
        "writes 0.25 Tu/100 with Tu in percent). eta = 1 at the lip.\n"
        "Same eta convention as adiabatic_wall_temperature().\n\n"
        "Parameters:\n"
        "  x_over_S       : distance downstream of the slot lip / slot height, >= 0\n"
        "  M              : blowing ratio rho_s U_s/(rho_H U_H), > 0\n"
        "  Tu             : axial free-stream turbulence intensity, a FRACTION\n"
        "                   (0.14 = 14%; 0.07 to 0.35 run; above 1 raises)\n"
        "  C_m_multiplier : the user's rig-matching factor on C_m (1 = the\n"
        "                   report's line; its band is 0.4 to 1.6). Never scored.\n\n"
        "Outside the run's Tu, M or x/S span it warns and extrapolates.\n\n"
        "Returns: eta [-]");

  m.def("slot_film_effectiveness_marek_tacina_1975_and_derivatives",
        [](double x_over_S, double M, double Tu, double C_m_multiplier) {
          return combaero::cooling::slot_film_effectiveness_marek_tacina_1975_and_derivatives(
              x_over_S, M, Tu, C_m_multiplier, nullptr);
        },
        py::arg("x_over_S"), py::arg("M"), py::arg("Tu"),
        py::arg("C_m_multiplier") = 1.0,
        "Solver-facing (f, J) form: (eta, d eta/dM). Closed form.");

  m.def("cooled_wall_heat_flux", &combaero::cooling::cooled_wall_heat_flux,
        py::arg("T_hot"), py::arg("T_coolant"), py::arg("h_hot"),
        py::arg("h_coolant"), py::arg("eta"), py::arg("t_wall"),
        py::arg("k_wall"),
        "Film/effusion-cooled wall heat flux with single wall layer.\n\n"
        "Computes q = U * (T_aw - T_coolant), where T_aw is from cooling\n"
        "effectiveness and U includes both convection sides plus wall "
        "conduction.\n\n"
        "Parameters:\n"
        "  T_hot      : mainstream hot-gas temperature [K]\n"
        "  T_coolant  : coolant supply temperature [K]\n"
        "  h_hot      : hot-side HTC [W/(m^2*K)]\n"
        "  h_coolant  : coolant-side HTC [W/(m^2*K)]\n"
        "  eta        : adiabatic effectiveness [-]\n"
        "  t_wall     : wall thickness [m]\n"
        "  k_wall     : wall thermal conductivity [W/(m*K)]\n\n"
        "Returns: wall heat flux q [W/m^2]");

  // =========================================================================
  // Units System
  // =========================================================================

  py::class_<combaero::units::UnitInfo>(m, "UnitInfo")
      .def_readonly("input", &combaero::units::UnitInfo::input,
                    "Input units description")
      .def_readonly("output", &combaero::units::UnitInfo::output,
                    "Output units description")
      .def("__repr__", [](const combaero::units::UnitInfo &u) {
        return "UnitInfo(input='" + std::string(u.input) + "', output='" +
               std::string(u.output) + "')";
      });

  m.def(
      "get_units",
      [](const std::string &name) -> py::object {
        auto u = combaero::units::get_units(name);
        if (u) {
          return py::cast(combaero::units::UnitInfo{u->input, u->output});
        }
        return py::none();
      },
      py::arg("function_name"),
      "Get input/output units for a function by name.\n\n"
      "Returns UnitInfo with 'input' and 'output' fields, or None if not "
      "found.\n\n"
      "Example:\n"
      "    >>> combaero.get_units('density')\n"
      "    UnitInfo(input='T: K, P: Pa, X: mol/mol', output='kg/m^3')");

  m.def(
      "input_units",
      [](const std::string &name) {
        return std::string(combaero::units::input_units(name));
      },
      py::arg("function_name"),
      "Get input units string for a function (empty string if not found).");

  m.def(
      "output_units",
      [](const std::string &name) {
        return std::string(combaero::units::output_units(name));
      },
      py::arg("function_name"),
      "Get output units string for a function (empty string if not found).");

  m.def(
      "has_units",
      [](const std::string &name) { return combaero::units::has_units(name); },
      py::arg("function_name"), "Check if a function has registered units.");

  m.def(
      "list_functions_with_units",
      []() {
        std::vector<std::string> names;
        for (auto it = combaero::units::begin(); it != combaero::units::end();
             ++it) {
          names.emplace_back(it->name);
        }
        return names;
      },
      "List all function names that have registered units.");

  m.def(
      "all_units",
      []() {
        py::dict result;
        for (auto it = combaero::units::begin(); it != combaero::units::end();
             ++it) {
          py::dict entry;
          entry["input"] = std::string(it->input);
          entry["output"] = std::string(it->output);
          result[py::str(std::string(it->name))] = entry;
        }
        return result;
      },
      "Get all registered units as a dictionary.\n\n"
      "Returns: dict mapping function_name -> {'input': ..., 'output': ...}");

  // =========================================================================
  // Incompressible flow -- scalar primitives
  // =========================================================================

  m.def(
      "channel_area", [](double D) { return circular_area(D); }, py::arg("D"),
      "Channel cross-sectional area [m^2].");

  m.def(
      "channel_volume",
      [](double D, double L) { return cylinder_volume(D, L); }, py::arg("D"),
      py::arg("L"), "Channel volume [m^3].");

  m.def(
      "pressure_loss",
      [](double v, double rho, double zeta) {
        return pressure_loss(v, rho, zeta);
      },
      py::arg("v"), py::arg("rho"), py::arg("zeta"),
      "Pressure drop from loss coefficient: dP = zeta * 0.5 * rho * v^2 [Pa].");

  m.def(
      "velocity_from_pressure_loss",
      [](double dP, double rho, double zeta) {
        return velocity_from_pressure_loss(dP, rho, zeta);
      },
      py::arg("dP"), py::arg("rho"), py::arg("zeta"),
      "Velocity from pressure drop and loss coefficient [m/s].");

  m.def(
      "zeta_from_Cd", [](double Cd) { return zeta_from_Cd(Cd); }, py::arg("Cd"),
      "Loss coefficient from discharge coefficient: zeta = 1/Cd^2 [-].");

  m.def(
      "Cd_from_zeta", [](double zeta) { return Cd_from_zeta(zeta); },
      py::arg("zeta"),
      "Discharge coefficient from loss coefficient: Cd = 1/sqrt(zeta) [-].");

  // =========================================================================
  // Incompressible flow -- thermo-aware high-level API
  // =========================================================================

  m.def(
      "orifice_flow_thermo",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P_back, double A, double Cd) {
        return orifice_flow_thermo(T, P, to_vec(X_arr), P_back, A, Cd);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("P_back"), py::arg("A"),
      py::arg("Cd") = 1.0,
      "Thermo-aware incompressible orifice flow.\n"
      "Evaluates rho from (T, P, X) internally.\n"
      "Returns IncompressibleFlowSolution.");

  m.def(
      "orifice_flow_thermo",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double P_back, double A, py::function cd_fn) {
        auto vec = to_vec(X_arr);
        IncompressibleCdFn fn = [cd_fn, vec](double T_, double P_,
                                             const std::vector<double> &,
                                             double Re) {
          return cd_fn(T_, P_, vec, Re).cast<double>();
        };
        return orifice_flow_thermo(T, P, vec, P_back, A, fn);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("P_back"), py::arg("A"),
      py::arg("cd_fn"),
      "Thermo-aware orifice flow with user-supplied Cd callable.\n\n"
      "cd_fn(T, P, X, Re) -> float  (must return Cd in (0, 1])\n"
      "Re is computed from the ideal throat velocity (loop-free).\n"
      "Returns IncompressibleFlowSolution.");

  m.def(
      "channel_flow",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double u, double L, double D, double f, std::size_t n_steps,
         bool store_profile) {
        return channel_flow(T, P, to_vec(X_arr), u, L, D, f, n_steps,
                            store_profile);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("u"), py::arg("L"),
      py::arg("D"), py::arg("f"), py::arg("n_steps") = 10,
      py::arg("store_profile") = false,
      "Thermo-aware incompressible channel flow with explicit friction "
      "factor.\n"
      "Returns IncompressibleFlowSolution.");

  m.def(
      "channel_flow_rough",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double u, double L, double D, double roughness,
         const std::string &correlation, std::size_t n_steps,
         bool store_profile) {
        return channel_flow_rough(T, P, to_vec(X_arr), u, L, D, roughness,
                                  correlation, n_steps, store_profile);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("u"), py::arg("L"),
      py::arg("D"), py::arg("roughness") = 0.0,
      py::arg("correlation") = "haaland", py::arg("n_steps") = 10,
      py::arg("store_profile") = false,
      "Thermo-aware incompressible channel flow with roughness-based friction "
      "factor.\n"
      "Returns IncompressibleFlowSolution.");

  m.def(
      "channel_flow_rough",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double u, double L, double D, double roughness,
         const std::string &correlation, py::function k_loss_fn) {
        auto vec = to_vec(X_arr);
        IncompressibleKLossFn fn = [k_loss_fn, vec](double T_, double P_,
                                                    const std::vector<double> &,
                                                    double Re) {
          return k_loss_fn(T_, P_, vec, Re).cast<double>();
        };
        return channel_flow_rough(T, P, vec, u, L, D, roughness, correlation,
                                  fn);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("u"), py::arg("L"),
      py::arg("D"), py::arg("roughness") = 0.0,
      py::arg("correlation") = "haaland", py::arg("k_loss_fn"),
      "Channel flow with additional K-loss callable (bends, fittings, "
      "etc.).\n\n"
      "k_loss_fn(T, P, X, Re) -> float  (returns K [-])\n"
      "Extra dP = K * 0.5 * rho * u^2 is added to Darcy-Weisbach friction "
      "loss.\n"
      "Returns IncompressibleFlowSolution.");

  m.def(
      "channel_pressure_drop",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double v, double D, double L, double roughness,
         const std::string &correlation) {
        return channel_pressure_drop(T, P, to_vec(X_arr), v, D, L, roughness,
                                     correlation);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("v"), py::arg("D"),
      py::arg("L"), py::arg("roughness") = 0.0,
      py::arg("correlation") = "haaland",
      "Composite channel pressure drop.\n"
      "Returns tuple (dP [Pa], Re [-], f [-]).");

  // =========================================================================
  // Compressible flow -- new functions
  // =========================================================================

  m.def(
      "fanno_channel_rough",
      [](double T_in, double P_in, double u_in, double L, double D,
         double roughness,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         const std::string &correlation, double f_multiplier,
         std::size_t n_steps, bool store_profile) {
        return fanno_channel_rough(T_in, P_in, u_in, L, D, roughness,
                                   to_vec(X_arr), correlation, f_multiplier,
                                   n_steps, store_profile);
      },
      py::arg("T_in"), py::arg("P_in"), py::arg("u_in"), py::arg("L"),
      py::arg("D"), py::arg("roughness"), py::arg("X"),
      py::arg("correlation") = "haaland", py::arg("f_multiplier") = 1.0,
      py::arg("n_steps") = 100, py::arg("store_profile") = false,
      "Fanno flow with roughness-based variable friction factor f(x).\n"
      "T_in, P_in: the inlet STATIC state [K, Pa]; u_in the inlet velocity.\n"
      "Recomputes f at each RK4 stage from local Re and roughness.\n"
      "Returns FannoSolution with f_avg and Re_in populated.");

  m.def(
      "fanno_max_length",
      [](double T_in, double P_in, double u_in, double D, double f,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return fanno_max_length(T_in, P_in, u_in, D, f, to_vec(X_arr), tol,
                                max_iter);
      },
      py::arg("T_in"), py::arg("P_in"), py::arg("u_in"), py::arg("D"),
      py::arg("f"), py::arg("X"), py::arg("tol") = 1e-6,
      py::arg("max_iter") = 100,
      "Maximum channel length before choking (Fanno L*) [m] from the inlet\n"
      "STATIC state (T_in [K], P_in [Pa]) and velocity u_in [m/s].");

  // =========================================================================
  // Stagnation / static conversion utilities  (stagnation.h)
  // =========================================================================

  m.def("kinetic_energy", &kinetic_energy, py::arg("v"),
        "Specific kinetic energy v^2/2 [J/kg].");

  m.def("h0_from_static", &h0_from_static, py::arg("h_static_J_per_kg"),
        py::arg("v"),
        "Stagnation enthalpy from static enthalpy [J/kg] and velocity [m/s].\n"
        "h0 = h_static + v^2/2  [J/kg].");

  m.def("v_from_h0", &v_from_h0, py::arg("h0"), py::arg("h_static_J_per_kg"),
        "Velocity [m/s] from stagnation and static enthalpy [J/kg].\n"
        "v = sqrt(2*(h0 - h_static)). Returns 0 if h0 <= h_static.");

  m.def(
      "mach_number",
      [](double v, double T,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        return mach_number(v, T, to_vec(X_arr));
      },
      py::arg("v"), py::arg("T"), py::arg("X"),
      "Mach number from velocity [m/s] and static temperature [K].\n"
      "M = v / a(T, X).");

  m.def(
      "Tt_from_static",
      [](double T, double M,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return Tt_from_static(T, M, to_vec(X_arr), tol, max_iter);
      },
      py::arg("T"), py::arg("M"), py::arg("X"), py::arg("tol") = 1e-8,
      py::arg("max_iter") = 50,
      "Stagnation temperature Tt [K] from static T [K] and Mach number.\n"
      "Iterative for variable cp(T) (thermally perfect gas).");

  m.def(
      "Tt_from_static_v",
      [](double T, double v,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return Tt_from_static_v(T, v, to_vec(X_arr), tol, max_iter);
      },
      py::arg("T"), py::arg("v"), py::arg("X"), py::arg("tol") = 1e-8,
      py::arg("max_iter") = 50,
      "Stagnation temperature Tt [K] from static T [K] and velocity v [m/s].\n"
      "Iterative for variable cp(T).");

  m.def(
      "T_from_stagnation",
      [](double Tt, double M,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return T_from_stagnation(Tt, M, to_vec(X_arr), tol, max_iter);
      },
      py::arg("Tt"), py::arg("M"), py::arg("X"), py::arg("tol") = 1e-8,
      py::arg("max_iter") = 50,
      "Static temperature T [K] from stagnation Tt [K] and Mach number.\n"
      "Iterative for variable cp(T).");

  m.def(
      "Pt_from_static",
      [](double P, double T, double M,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return Pt_from_static(P, T, M, to_vec(X_arr), tol, max_iter);
      },
      py::arg("P"), py::arg("T"), py::arg("M"), py::arg("X"),
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Stagnation pressure Pt [Pa] from static P [Pa], T [K] and Mach number.\n"
      "Uses isentropic entropy conservation (variable cp).");

  m.def(
      "P_from_stagnation",
      [](double Pt, double Tt, double M,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return P_from_stagnation(Pt, Tt, M, to_vec(X_arr), tol, max_iter);
      },
      py::arg("Pt"), py::arg("Tt"), py::arg("M"), py::arg("X"),
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Static pressure P [Pa] from stagnation Pt [Pa], Tt [K] and Mach "
      "number.\n"
      "Uses isentropic entropy conservation (variable cp).");

  m.def("recovery_factor", &recovery_factor, py::arg("Pr"),
        py::arg("turbulent") = true,
        "Recovery factor r [-] for adiabatic wall temperature.\n"
        "r = Pr^(1/3) for turbulent (default), Pr^(1/2) for laminar.");

  m.def(
      "T_adiabatic_wall",
      [](double T_static, double v, double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         bool turbulent) {
        return T_adiabatic_wall(T_static, v, T, P, to_vec(X_arr), turbulent);
      },
      py::arg("T_static"), py::arg("v"), py::arg("T"), py::arg("P"),
      py::arg("X"), py::arg("turbulent") = true,
      "Adiabatic wall temperature T_aw [K] from static T [K] and velocity v "
      "[m/s].\n"
      "T_aw = T_static + r * v^2 / (2*cp)  where r = Pr^(1/3) turbulent.\n"
      "Use T_aw (not T_static or Tt) as the hot-side driving temperature\n"
      "in channel evaluations and cooled_wall_heat_flux() for M > 0.3.");

  m.def(
      "T_adiabatic_wall_mach",
      [](double T_static, double M, double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         bool turbulent) {
        return T_adiabatic_wall_mach(T_static, M, T, P, to_vec(X_arr),
                                     turbulent);
      },
      py::arg("T_static"), py::arg("M"), py::arg("T"), py::arg("P"),
      py::arg("X"), py::arg("turbulent") = true,
      "Adiabatic wall temperature T_aw [K] from static T [K] and Mach number.\n"
      "Computes v = M * a(T, X) internally, then applies recovery factor "
      "correction.");

  m.def(
      "stagnation_from_static",
      [](double T, double P, double v,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double tol, std::size_t max_iter) {
        return stagnation_from_static(T, P, v, to_vec(X_arr), tol, max_iter);
      },
      py::arg("T"), py::arg("P"), py::arg("v"), py::arg("X"),
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Combined stagnation state (Tt, Pt) from static (T, P) and velocity v.\n"
      "Returns tuple (Tt [K], Pt [Pa]).\n"
      "Preferred over calling Tt_from_static_v + Pt_from_static separately\n"
      "as it uses a single Newton loop and avoids redundant physics calls.");

  m.def("bulk_velocity", &bulk_velocity, py::arg("m_dot"), py::arg("rho"),
        py::arg("area"),
        "Bulk velocity v = m_dot / (rho * area)  [m/s].");

  m.def("reynolds_from_state", &reynolds_from_state, py::arg("rho"),
        py::arg("v"), py::arg("L"), py::arg("mu"),
        "Reynolds number Re = rho * v * L / mu  [-].\n"
        "Use when rho and mu are already available from CompleteState or\n"
        "TransportState to avoid recomputing mixture properties.");

  // -----------------------------------------------------------------
  // ChannelResult struct
  // -----------------------------------------------------------------
  py::class_<ChannelResult>(
      m, "ChannelResult",
      "Combined convective heat transfer + pressure loss result.\n\n"
      "Returned by channel_smooth, channel_ribbed, channel_dimpled,\n"
      "channel_pin_fin, and channel_impingement.\n\n"
      "All fluid properties are evaluated in a single pass from (T, P, X).\n\n"
      "T_aw is always computed continuously for all M:\n"
      "  T_aw = T_static + r * v^2 / (2*cp)  (r = Pr^1/3 turbulent)\n"
      "  At v=0: T_aw = T_static exactly. No threshold, no discontinuity.\n\n"
      "q = h * (T_aw - T_hot) when T_hot is supplied, else nan.\n\n"
      "f meaning depends on geometry:\n"
      "  smooth/ribbed/dimpled : Darcy friction factor [-]\n"
      "  pin_fin               : pin-array friction coefficient [-]\n"
      "  impingement           : jet-plate loss coefficient 1/Cd^2 [-]")
      .def(py::init<>())
      .def_readonly("h", &ChannelResult::h, "Convective HTC [W/(m^2*K)]")
      .def_readonly("Nu", &ChannelResult::Nu, "Nusselt number [-]")
      .def_readonly("Re", &ChannelResult::Re, "Reynolds number [-]")
      .def_readonly("Pr", &ChannelResult::Pr, "Prandtl number [-]")
      .def_readonly("f", &ChannelResult::f, "Friction / loss coefficient [-]")
      .def_readonly("dP", &ChannelResult::dP, "Pressure drop [Pa]")
      .def_readonly("M", &ChannelResult::M, "Mach number [-]")
      .def_readonly("T_aw", &ChannelResult::T_aw,
                    "Adiabatic wall temperature [K]")
      .def_readonly("q", &ChannelResult::q,
                    "Heat flux [W/m^2] (nan if T_hot not supplied)")
      .def_readonly("dh_dmdot", &ChannelResult::dh_dmdot,
                    "Jacobian: dh/dmdot [W/(m^2*K*kg/s)]")
      .def_readonly("dh_dT", &ChannelResult::dh_dT,
                    "Jacobian: dh/dT [W/(m^2*K^2)]")
      .def_readonly("ddP_dmdot", &ChannelResult::ddP_dmdot,
                    "Jacobian: d(dP)/dmdot [Pa*s/kg]")
      .def_readonly("ddP_dvelocity", &ChannelResult::ddP_dvelocity,
                    "Jacobian: d(dP)/d(velocity) [Pa*s/m]. Chain from this, not "
                    "ddP_dmdot, when converting to a caller's own mass flow: "
                    "ddP_dmdot uses the correlation's internal flow area.")
      .def_readonly("ddP_dT", &ChannelResult::ddP_dT,
                    "Jacobian: d(dP)/dT [Pa/K]")
      .def_readonly("dT_aw_dmdot", &ChannelResult::dT_aw_dmdot,
                    "Jacobian: dT_aw/dmdot [K*s/kg]")
      .def_readonly("dT_aw_dT", &ChannelResult::dT_aw_dT,
                    "Jacobian: dT_aw/dT [-] (approx 1 at low Mach)")
      .def_readonly("dT_aw_dP", &ChannelResult::dT_aw_dP,
                    "Jacobian: dT_aw/dP [K/Pa] at fixed mdot and T")
      .def_readonly("dq_dmdot", &ChannelResult::dq_dmdot,
                    "Jacobian: dq/dmdot [W*s/(m^2*kg)]")
      .def_readonly("dq_dT", &ChannelResult::dq_dT,
                    "Jacobian: dq/dT [W/(m^2*K)]")
      .def_readonly("dq_dT_hot", &ChannelResult::dq_dT_hot,
                    "Jacobian: dq/dT_hot [W/(m^2*K)]")
      .def_readonly("dh_dY", &ChannelResult::dh_dY,
                    "Jacobian: dh/dY [W/(m^2*K)] per species (mass fractions) at fixed\n"
                    "mdot, T, P (#717); empty for a correlation that does not provide it")
      .def_readonly("dT_aw_dY", &ChannelResult::dT_aw_dY,
                    "Jacobian: dT_aw/dY [K] per species (mass fractions) at fixed\n"
                    "mdot, T, P (#717); empty for a correlation that does not provide it")
      .def("__repr__", [](const ChannelResult &r) {
        auto fmt = [](double v) {
          if (std::isnan(v))
            return std::string("nan");
          char buf[32];
          std::snprintf(buf, sizeof(buf), "%.6g", v);
          return std::string(buf);
        };
        return "ChannelResult(h=" + fmt(r.h) + ", Nu=" + fmt(r.Nu) +
               ", Re=" + fmt(r.Re) + ", dP=" + fmt(r.dP) +
               ", T_aw=" + fmt(r.T_aw) + ", q=" + fmt(r.q) + ")";
      });

  // -----------------------------------------------------------------
  // WallCouplingResult struct
  // -----------------------------------------------------------------
  py::class_<WallCouplingResult>(
      m, "WallCouplingResult",
      "Heat transfer through a wall between two elements.\n\n"
      "Returned by wall_coupling_and_jacobian.\n\n"
      "Q is positive when heat flows from side A to side B.\n"
      "T_hot is the hot-side surface temperature.\n"
      "Jacobians enable proper Newton coupling in the network solver.")
      .def_readonly("Q", &WallCouplingResult::Q,
                    "Heat transfer rate [W] (positive A->B)")
      .def_readonly("T_hot", &WallCouplingResult::T_hot,
                    "Hot-side wall surface temperature [K].")
      .def_readonly("dQ_dh_a", &WallCouplingResult::dQ_dh_a,
                    "Jacobian: dQ/dh_a [W*m^2*K/W]")
      .def_readonly("dQ_dh_b", &WallCouplingResult::dQ_dh_b,
                    "Jacobian: dQ/dh_b [W*m^2*K/W]")
      .def_readonly("dQ_dT_aw_a", &WallCouplingResult::dQ_dT_aw_a,
                    "Jacobian: dQ/dT_aw_a [W/K]")
      .def_readonly("dQ_dT_aw_b", &WallCouplingResult::dQ_dT_aw_b,
                    "Jacobian: dQ/dT_aw_b [W/K]")
      .def_readonly("dQ_dC_a", &WallCouplingResult::dQ_dC_a,
                    "Jacobian: dQ/dC_a [K] (wall_exchange_and_jacobian; 0 for the "
                    "linear coupling)")
      .def_readonly("dQ_dC_b", &WallCouplingResult::dQ_dC_b,
                    "Jacobian: dQ/dC_b [K] (wall_exchange_and_jacobian)")
      .def_readonly("UA", &WallCouplingResult::UA,
                    "Overall conductance [W/K] (wall_exchange_and_jacobian)")
      .def_readonly("effectiveness", &WallCouplingResult::effectiveness,
                    "Q / (C_min dT) [-] (wall_exchange_and_jacobian; 0 without "
                    "a finite C_min)")
      .def_readonly("dQ_dR_wall", &WallCouplingResult::dQ_dR_wall,
                    "Jacobian: dQ/d(t/k) of any one layer, = dQ/dR_fouling "
                    "[W^2/(m^2*K)] (#932)")
      .def("__repr__",
           [](const WallCouplingResult &r) {
             auto fmt = [](double v) {
               char buf[32];
               std::snprintf(buf, sizeof(buf), "%.6g", v);
               return std::string(buf);
             };
             return "WallCouplingResult(Q=" + fmt(r.Q) +
                    ", T_hot=" + fmt(r.T_hot) + ")";
           })
      .def("__iter__", [](const WallCouplingResult &r) {
        return py::iter(py::make_tuple(r.Q, r.T_hot));
      });

  // -----------------------------------------------------------------
  // wall_coupling_and_jacobian
  // -----------------------------------------------------------------
  m.def(
      "wall_coupling_and_jacobian",
      py::overload_cast<double, double, double, double, double, double>(
          &wall_coupling_and_jacobian),
      py::arg("h_a"), py::arg("T_aw_a"), py::arg("h_b"), py::arg("T_aw_b"),
      py::arg("t_over_k"), py::arg("A"),
      "Compute heat transfer and Jacobians for single-layer wall coupling.\n\n"
      "Parameters:\n"
      "  h_a, h_b       : convective HTCs on sides A and B [W/(m^2*K)]\n"
      "  T_aw_a, T_aw_b : adiabatic wall temperatures on sides A and B [K]\n"
      "  t_over_k       : wall_thickness / wall_conductivity [m^2*K/W]\n"
      "  A              : effective heat transfer area [m^2]\n\n"
      "Returns: WallCouplingResult with Q and analytical Jacobians");

  m.def("wall_coupling_and_jacobian_multilayer",
        py::overload_cast<double, double, double, double,
                          const std::vector<double> &, double, double>(
            &wall_coupling_and_jacobian),
        py::arg("h_a"), py::arg("T_aw_a"), py::arg("h_b"), py::arg("T_aw_b"),
        py::arg("t_over_k_layers"), py::arg("A"), py::arg("R_fouling") = 0.0,
        "Compute heat transfer and Jacobians for multi-layer wall coupling.\n\n"
        "Parameters:\n"
        "  h_a, h_b       : convective HTCs on sides A and B [W/(m^2*K)]\n"
        "  T_aw_a, T_aw_b : adiabatic wall temperatures on sides A and B [K]\n"
        "  t_over_k_layers: thicknesses / conductivities for each layer "
        "[m^2*K/W]\n"
        "  A              : effective heat transfer area [m^2]\n"
        "  R_fouling      : additional fouling resistance [m^2*K/W]\n\n"
        "Returns: WallCouplingResult with Q and analytical Jacobians");

  m.def("wall_exchange_and_jacobian", &wall_exchange_and_jacobian, py::arg("h_a"),
        py::arg("T_aw_a"), py::arg("C_a"), py::arg("h_b"), py::arg("T_aw_b"),
        py::arg("C_b"), py::arg("t_over_k_layers"), py::arg("A"),
        py::arg("R_fouling") = 0.0, py::arg("counterflow") = false,
        "Heat two streams exchange through a wall, epsilon-NTU (#627).\n\n"
        "Q = eps(NTU, C_r) C_min (T_aw_a - T_aw_b), T_aw at each stream's inlet.\n"
        "  C_a, C_b    : capacity rates |m| cp [W/K]; inf for a node side "
        "(a reservoir)\n"
        "  counterflow : the streams run opposite ways along the wall\n\n"
        "Returns: WallCouplingResult with Q and analytic Jacobians, including\n"
        "dQ_dC_a and dQ_dC_b.");

  m.def("wall_resistance_fixed_point_derivatives",
        &wall_resistance_fixed_point_derivatives, py::arg("T_hot"),
        py::arg("T_cold"), py::arg("h_hot"), py::arg("h_cold"),
        py::arg("t_over_k"), py::arg("R_fouling"), py::arg("dr_dTm"),
        "Total derivatives of a k(T) wall's R_wall = sum t_i/k_i at its fixed\n"
        "point (#932): each layer's t/k at its mean temperature in\n"
        "wall_temperature_profile(T_hot, T_cold, h_hot, h_cold, t_over_k,\n"
        "R_fouling), dr_dTm the layers' d(t/k)/dT [m^2/W] there (0 for a\n"
        "constant k). Returns dR_wall/d(T_hot, T_cold, h_hot, h_cold) by\n"
        "implicit differentiation; analytic.");

  // -----------------------------------------------------------------
  // Liner radiation (#793): Lefebvre & Ballal Ch. 8, Eqs. 8.7-8.17
  // -----------------------------------------------------------------
  m.attr("STEFAN_BOLTZMANN") = STEFAN_BOLTZMANN;
  m.attr("REEVES_EMISSIVITY_COEFF") = REEVES_EMISSIVITY_COEFF;
  m.attr("LUMINOSITY_HYDROGEN_COEFF") = LUMINOSITY_HYDROGEN_COEFF;
  m.attr("BEAM_LENGTH_COEFF_LEFEBVRE") = BEAM_LENGTH_COEFF_LEFEBVRE;
  m.attr("BEAM_LENGTH_COEFF_HOWELL") = BEAM_LENGTH_COEFF_HOWELL;

  py::class_<GasEmissivityResult>(
      m, "GasEmissivityResult",
      "Gas emissivity (Lefebvre & Ballal Eq. 8.11) and its partials.")
      .def_readonly("eps", &GasEmissivityResult::eps, "eps_g [-]")
      .def_readonly("deps_dT", &GasEmissivityResult::deps_dT,
                    "d eps_g / d T_g [1/K]")
      .def_readonly("deps_dP", &GasEmissivityResult::deps_dP,
                    "d eps_g / d P [1/Pa]")
      .def_readonly("deps_dq", &GasEmissivityResult::deps_dq,
                    "d eps_g / d q [-]")
      .def_readonly("deps_dl_b", &GasEmissivityResult::deps_dl_b,
                    "d eps_g / d l_b [1/m]")
      .def_readonly("deps_dL", &GasEmissivityResult::deps_dL,
                    "d eps_g / d L [-]")
      .def("__iter__", [](const GasEmissivityResult &r) {
        return py::iter(py::make_tuple(r.eps, r.deps_dT, r.deps_dP, r.deps_dq,
                                       r.deps_dl_b, r.deps_dL));
      });

  m.def("gas_emissivity_luminous", &gas_emissivity_luminous, py::arg("T_g"),
        py::arg("P"), py::arg("q"), py::arg("l_b"), py::arg("L") = 1.0,
        "Luminous gas emissivity, Lefebvre & Ballal Eq. 8.11 (Reeves Eq. 8.9\n"
        "with L = 1):\n"
        "  eps_g = 1 - exp[-290 (P/1000) L (q l_b)^0.5 T_g^-1.5]\n\n"
        "  T_g : mean gas temperature [K]\n"
        "  P   : gas pressure [Pa] (the book's kPa formed internally)\n"
        "  q   : fuel/air mass ratio [-]\n"
        "  l_b : mean beam length [m]\n"
        "  L   : luminosity factor [-] (1: non-luminous)\n\n"
        "A kerosine-product fit (Reeves 1956): gray, single-zone; not a\n"
        "composition-based emissivity. Returns GasEmissivityResult with the\n"
        "analytic partials in T_g, P, q, l_b and L.");

  m.def("luminosity_factor_hydrogen", &luminosity_factor_hydrogen,
        py::arg("H_mass_percent"),
        "Luminosity factor L = 336 / H^2 (Lefebvre & Ballal Eq. 8.15), H the\n"
        "fuel hydrogen content in mass percent. Returns (L, dL/dH). The\n"
        "alternatives Eqs. 8.12-8.14 (in the fuel C/H ratio) are not shipped;\n"
        "pass L directly to use one.");

  py::class_<LinerGasRadiationResult>(
      m, "LinerGasRadiationResult",
      "Gas-to-liner radiation flux R1 (Lefebvre & Ballal Eq. 8.7) and its\n"
      "partials at fixed eps_g.")
      .def_readonly("R1", &LinerGasRadiationResult::R1,
                    "R1 [W/m^2], positive gas -> wall")
      .def_readonly("dR1_dT_g", &LinerGasRadiationResult::dR1_dT_g,
                    "dR1/dT_g at fixed eps_g [W/(m^2 K)]")
      .def_readonly("dR1_dT_w", &LinerGasRadiationResult::dR1_dT_w,
                    "dR1/dT_w [W/(m^2 K)]")
      .def_readonly("dR1_deps_g", &LinerGasRadiationResult::dR1_deps_g,
                    "dR1/deps_g [W/m^2]")
      .def_readonly("dR1_deps_w", &LinerGasRadiationResult::dR1_deps_w,
                    "dR1/deps_w [W/m^2]")
      .def("__iter__", [](const LinerGasRadiationResult &r) {
        return py::iter(py::make_tuple(r.R1, r.dR1_dT_g, r.dR1_dT_w,
                                       r.dR1_deps_g, r.dR1_deps_w));
      });

  m.def("liner_gas_radiation", &liner_gas_radiation, py::arg("T_g"),
        py::arg("T_w"), py::arg("eps_g"), py::arg("eps_w"),
        "Gas-to-liner radiation flux, Lefebvre & Ballal Eq. 8.7:\n"
        "  R1 = 0.5 sigma (1 + eps_w) eps_g T_g^1.5 (T_g^2.5 - T_w^2.5)\n"
        "[W/m^2]. eps_g from gas_emissivity_luminous. Returns\n"
        "LinerGasRadiationResult (partials at fixed eps_g).");

  m.def("liner_casing_radiation", &liner_casing_radiation, py::arg("T_w"),
        py::arg("T_c"), py::arg("eps_w"), py::arg("eps_c"),
        py::arg("area_ratio"),
        "Liner-to-casing radiation flux per unit liner area, Lefebvre &\n"
        "Ballal Eq. 8.17 (gray two-surface enclosure, view factor 1):\n"
        "  R2 = sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c)\n"
        "       (T_w^4 - T_c^4)\n"
        "area_ratio = A_w/A_c. Returns (R2, dR2/dT_w, dR2/dT_c).");

  m.def("liner_casing_exchange_factor", &liner_casing_exchange_factor,
        py::arg("eps_w"), py::arg("eps_c"), py::arg("area_ratio"),
        "sigma eps_w eps_c / (eps_c + eps_w (1 - eps_c) A_w/A_c)\n"
        "[W/(m^2 K^4)], the exchange factor of Eq. 8.17.");

  m.def("mean_beam_length", &mean_beam_length, py::arg("volume"),
        py::arg("area"), py::arg("coeff") = BEAM_LENGTH_COEFF_LEFEBVRE,
        "Mean beam length l_b = coeff V / A [m]: 3.4 (Lefebvre & Ballal\n"
        "Eq. 8.10, default) or BEAM_LENGTH_COEFF_HOWELL = 3.6 (0.9 x 4 V/A,\n"
        "the textbook value).");

  m.def(
      "elemental_fuel_air_ratio",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> Y) {
        return elemental_fuel_air_ratio(to_vec(Y));
      },
      py::arg("Y"),
      "Fuel/air mass ratio from the elements: mass of C and H atoms over the\n"
      "rest, q = w / (1 - w). The same for reactants and products. Water\n"
      "vapour's H counts as fuel; a fuel's own O or N counts as air.\n"
      "Returns (q, dq/dY) in the mass fractions Y.");

  // Bordbar et al. (2021) pressure-dependent WSGG emissivity (#1071).
  m.attr("WSGG_BORDBAR2021_T_REF") = WSGG_BORDBAR2021_T_REF;
  m.attr("WSGG_BORDBAR2021_T_MIN") = WSGG_BORDBAR2021_T_MIN;
  m.attr("WSGG_BORDBAR2021_T_MAX") = WSGG_BORDBAR2021_T_MAX;
  m.attr("WSGG_BORDBAR2021_MR_MIN") = WSGG_BORDBAR2021_MR_MIN;
  m.attr("WSGG_BORDBAR2021_MR_MAX") = WSGG_BORDBAR2021_MR_MAX;
  m.attr("WSGG_PA_PER_ATM") = WSGG_PA_PER_ATM;

  py::class_<WsggGrayGases>(
      m, "WsggGrayGases",
      "The four gray gases of one Bordbar (2021) pressure set: weights a [-]\n"
      "and pressure absorption coefficients kappa [1/(atm m)].")
      .def_readonly("a", &WsggGrayGases::a, "gray-gas weights a_i [-]")
      .def_readonly("kappa", &WsggGrayGases::kappa,
                    "pressure absorption coefficients kappa_p,i [1/(atm m)]");

  m.def("wsgg_bordbar2021_gray_gases", &wsgg_bordbar2021_gray_gases,
        py::arg("p_set_atm"), py::arg("T"), py::arg("Mr"),
        "Bordbar et al. (2021) Eqs. 8, 10, 11 at one tabulated pressure set\n"
        "(0.1, 0.5, 1, 2, 5, 10, 20, 30, 50 or 80 atm): the gray-gas weights\n"
        "and kappas at temperature T [K] (300-3000) and mole ratio\n"
        "Mr = x_H2O/x_CO2 (0.1-4), for an RTE solve. Unclamped: raises\n"
        "outside the fitted range. Returns WsggGrayGases.");

  py::class_<WsggEmissivityResult>(
      m, "WsggEmissivityResult",
      "Bordbar (2021) WSGG emissivity and its partials.")
      .def_readonly("eps", &WsggEmissivityResult::eps, "eps_g [-]")
      .def_readonly("deps_dT", &WsggEmissivityResult::deps_dT,
                    "d eps_g / d T_g [1/K]")
      .def_readonly("deps_dP", &WsggEmissivityResult::deps_dP,
                    "d eps_g / d P [1/Pa]")
      .def_readonly("deps_dx_h2o", &WsggEmissivityResult::deps_dx_h2o,
                    "d eps_g / d x_H2O [-]")
      .def_readonly("deps_dx_co2", &WsggEmissivityResult::deps_dx_co2,
                    "d eps_g / d x_CO2 [-]")
      .def_readonly("deps_dl_b", &WsggEmissivityResult::deps_dl_b,
                    "d eps_g / d l_b [1/m]");

  m.def("gas_emissivity_bordbar2021", &gas_emissivity_bordbar2021,
        py::arg("T_g"), py::arg("P"), py::arg("x_h2o"), py::arg("x_co2"),
        py::arg("l_b"),
        "Non-luminous H2O-CO2 gas emissivity, Bordbar et al. (2021) pressure-\n"
        "dependent WSGG (Int. J. Heat Mass Transfer 173, 121207, App. A):\n"
        "  eps = sum_i a_i(T, Mr) [1 - exp(-kappa_i(Mr) (x_H2O + x_CO2) P L)]\n"
        "with eps interpolated linearly in ln P between the ten pressure sets\n"
        "(each at the actual pressure path). Clamps T to 300-3000 K, Mr to\n"
        "0.1-4 and eps to [0, 1] (zero partials where clamped).\n\n"
        "  T_g : gas temperature [K]\n"
        "  P   : total pressure [Pa]\n"
        "  x_h2o, x_co2 : mole fractions [-]\n"
        "  l_b : path or mean beam length [m]\n\n"
        "Returns WsggEmissivityResult with the analytic partials.");

  py::class_<WsggMixtureEmissivityResult>(
      m, "WsggMixtureEmissivityResult",
      "Bordbar (2021) WSGG emissivity of a mixture and its partials in the\n"
      "mass fractions.")
      .def_readonly("eps", &WsggMixtureEmissivityResult::eps, "eps_g [-]")
      .def_readonly("deps_dT", &WsggMixtureEmissivityResult::deps_dT,
                    "d eps_g / d T_g [1/K]")
      .def_readonly("deps_dP", &WsggMixtureEmissivityResult::deps_dP,
                    "d eps_g / d P [1/Pa]")
      .def_readonly("deps_dl_b", &WsggMixtureEmissivityResult::deps_dl_b,
                    "d eps_g / d l_b [1/m]")
      .def_readonly("deps_dY", &WsggMixtureEmissivityResult::deps_dY,
                    "d eps_g / d Y_k [-], one per species");

  m.def(
      "gas_emissivity_bordbar2021_mixture",
      [](double T_g, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> Y,
         double l_b) {
        return gas_emissivity_bordbar2021_mixture(T_g, P, to_vec(Y), l_b);
      },
      py::arg("T_g"), py::arg("P"), py::arg("Y"), py::arg("l_b"),
      "gas_emissivity_bordbar2021 from a mixture's mass fractions Y (one per\n"
      "species): x_H2O and x_CO2 formed internally, deps_dY chained through\n"
      "them. Returns WsggMixtureEmissivityResult.");

  py::class_<WallRadiationSides>(
      m, "WallRadiationSides",
      "A wall's local balance with radiation on its faces (#793), folded\n"
      "into effective sides (H, T) per face. d / dq_rad: total derivatives\n"
      "in [T_aw_a, h_a, T_aw_b, h_b, R_wall, T_src_a, k_a, T_src_b, k_b].")
      .def_readonly("H_a", &WallRadiationSides::H_a)
      .def_readonly("T_a", &WallRadiationSides::T_a)
      .def_readonly("H_b", &WallRadiationSides::H_b)
      .def_readonly("T_b", &WallRadiationSides::T_b)
      .def_readonly("T_s_a", &WallRadiationSides::T_s_a)
      .def_readonly("T_s_b", &WallRadiationSides::T_s_b)
      .def_readonly("h_r_a", &WallRadiationSides::h_r_a)
      .def_readonly("h_r_b", &WallRadiationSides::h_r_b)
      .def_readonly("q_rad_a", &WallRadiationSides::q_rad_a)
      .def_readonly("q_rad_b", &WallRadiationSides::q_rad_b)
      .def_readonly("q", &WallRadiationSides::q)
      .def_readonly("iterations", &WallRadiationSides::iterations)
      .def_readonly("d", &WallRadiationSides::d)
      .def_readonly("dq_rad", &WallRadiationSides::dq_rad);

  m.def("wall_radiation_sides", &wall_radiation_sides, py::arg("T_aw_a"),
        py::arg("h_a"), py::arg("T_aw_b"), py::arg("h_b"), py::arg("R_wall"),
        py::arg("kind_a"), py::arg("T_src_a"), py::arg("k_a"),
        py::arg("kind_b"), py::arg("T_src_b"), py::arg("k_b"),
        "Solve a wall's two surface temperatures with radiation on either\n"
        "face (kind 0 none, 1 gas Eq. 8.7, 2 surface Eq. 8.17) and fold it\n"
        "into effective (H, T) per side; analytic total derivatives.");

  py::class_<WallAdiabaticFace>(
      m, "WallAdiabaticFace",
      "The face of a wall whose other side is adiabatic (#970): its surface\n"
      "temperature from h (T_aw - T_s) + h_r (T_src - T_s) = 0. dT_s: in\n"
      "[T_aw, h, T_src, k].")
      .def_readonly("T_s", &WallAdiabaticFace::T_s)
      .def_readonly("h_r", &WallAdiabaticFace::h_r)
      .def_readonly("q_rad", &WallAdiabaticFace::q_rad)
      .def_readonly("iterations", &WallAdiabaticFace::iterations)
      .def_readonly("dT_s", &WallAdiabaticFace::dT_s);

  m.def("wall_adiabatic_face", &wall_adiabatic_face, py::arg("T_aw"),
        py::arg("h"), py::arg("kind"), py::arg("T_src"), py::arg("k"),
        "Surface temperature of a wall face whose other side is adiabatic\n"
        "(#970), with radiation on the face (kind 0 none: T_s = T_aw; 1 gas\n"
        "Eq. 8.7; 2 surface Eq. 8.17); analytic derivatives.");

  // -----------------------------------------------------------------
  // Film-side h downstream of a slot (#962): Seban 1960, Ballal-Lefebvre
  // -----------------------------------------------------------------
  m.attr("SEBAN1960_C_NEAR") = seban1960::C_near;
  m.attr("SEBAN1960_C_NEAR_PRINTED") = seban1960::C_near_printed;
  m.attr("SEBAN1960_C_FAR") = seban1960::C_far;
  m.attr("SEBAN1960_C_FAR_PRINTED") = seban1960::C_far_printed;
  m.attr("SLOT_FILM_HTC_HANDBACK_BAND") = SLOT_FILM_HTC_HANDBACK_BAND;
  m.attr("SLOT_FILM_BRANCH_NEAR") = SLOT_FILM_BRANCH_NEAR;
  m.attr("SLOT_FILM_BRANCH_BLEND") = SLOT_FILM_BRANCH_BLEND;
  m.attr("SLOT_FILM_BRANCH_FAR") = SLOT_FILM_BRANCH_FAR;
  m.attr("SLOT_FILM_BRANCH_WALL_JET") = SLOT_FILM_BRANCH_WALL_JET;
  m.attr("SLOT_FILM_BRANCH_NONE") = SLOT_FILM_BRANCH_NONE;

  py::class_<SlotFilmStanton>(
      m, "SlotFilmStanton",
      "Seban (1960) slot-film Stanton number, span mean over the in-box part.")
      .def_readonly("St", &SlotFilmStanton::St, "mean St_s [-]")
      .def_readonly("dSt_dRe", &SlotFilmStanton::dSt_dRe, "d St / d Re_s")
      .def_readonly("fraction_in_box", &SlotFilmStanton::fraction_in_box,
                    "share of the span at x/s <= x_over_s_max")
      .def_readonly("branch", &SlotFilmStanton::branch,
                    "SLOT_FILM_BRANCH_* code")
      .def_readonly("x_in_box", &SlotFilmStanton::x_in_box)
      .def_readonly("Re_in_box", &SlotFilmStanton::Re_in_box);

  m.def("slot_film_stanton_seban_1960", &slot_film_stanton_seban_1960,
        py::arg("x_start_over_s"), py::arg("x_end_over_s"), py::arg("Re_s"),
        py::arg("C_near") = seban1960::C_near,
        py::arg("C_far") = seban1960::C_far,
        "Seban (1960) Fig. 12 slot-film Stanton number on slot-fluid\n"
        "properties, St_s = h/(rho_s u_s cp_s), span mean over\n"
        "[x_start_over_s, x_end_over_s] (equal: the point value):\n"
        "  curve B: 0.135 Re_s^-0.3 (x/s)^-0.3 (Eq. 13 as drawn; printed 0.10)\n"
        "  curve C: 0.40 Re_s^-0.3 (x/s)^-0.6  (Eq. 14 as drawn; printed 0.41)\n"
        "joined by a C1 blend over +-0.1 decade around their crossing\n"
        "(x/s 37.4). Only the part at x/s <= 300 is averaged; fraction_in_box\n"
        "is its share of the span.");

  py::class_<SlotFilmNusselt>(
      m, "SlotFilmNusselt",
      "Ballal & Lefebvre (1973) wall-jet h s / k_s, span mean over the in-box "
      "part.")
      .def_readonly("Nu_s", &SlotFilmNusselt::Nu_s, "mean h s / k_s [-]")
      .def_readonly("dNu_dRe", &SlotFilmNusselt::dNu_dRe, "d Nu_s / d Re_s")
      .def_readonly("fraction_in_box", &SlotFilmNusselt::fraction_in_box)
      .def_readonly("branch", &SlotFilmNusselt::branch)
      .def_readonly("x_in_box", &SlotFilmNusselt::x_in_box);

  m.def("slot_film_nusselt_ballal_lefebvre_1973",
        &slot_film_nusselt_ballal_lefebvre_1973, py::arg("x_start_over_s"),
        py::arg("x_end_over_s"), py::arg("Re_s"),
        "Ballal & Lefebvre (1973) SME4 Eq. 6 (Lefebvre & Ballal Eq. 8.43),\n"
        "the m > 1.3 wall-jet form Nu_x = 0.10 Re_s^0.8 (x/s)^0.44, as\n"
        "h s / k_s = 0.10 Re_s^0.8 (x/s)^-0.56, span mean over the part at\n"
        "x/s <= 50. Over-predicts the Cranfield rig's own data by 20-67 % at\n"
        "low Re (SME4 Fig. 2): a named option, not the default.");

  py::class_<SlotFilmHtc>(
      m, "SlotFilmHtc",
      "A wall side's h under a slot film (#962): film h over its in-box span,\n"
      "the duct's h past it and below the hand-back mass-velocity ratio.")
      .def_readonly("h", &SlotFilmHtc::h, "the side's h [W/(m^2 K)]")
      .def_readonly("dh_dG_s", &SlotFilmHtc::dh_dG_s,
                    "d h / d G_s at fixed M")
      .def_readonly("dh_dT_s", &SlotFilmHtc::dh_dT_s,
                    "d h / d T_s (slot-fluid properties)")
      .def_readonly("dh_dM", &SlotFilmHtc::dh_dM,
                    "d h / d M (hand-back weight)")
      .def_readonly("dh_dh_duct", &SlotFilmHtc::dh_dh_duct, "1 - w f_x")
      .def_readonly("h_film", &SlotFilmHtc::h_film,
                    "the film's own span-mean h over its in-box part")
      .def_readonly("Re_s", &SlotFilmHtc::Re_s)
      .def_readonly("fraction_in_box", &SlotFilmHtc::fraction_in_box)
      .def_readonly("weight_M", &SlotFilmHtc::weight_M)
      .def_readonly("branch", &SlotFilmHtc::branch)
      .def_readonly("x_in_box", &SlotFilmHtc::x_in_box)
      .def_readonly("Re_in_box", &SlotFilmHtc::Re_in_box)
      .def_readonly("M_in_box", &SlotFilmHtc::M_in_box);

  m.def(
      "slot_film_htc_seban_1960",
      [](double x_start, double x_end, double s, double G_s, double M,
         double h_duct, double T_s, double P_s,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_s) {
        return slot_film_htc_seban_1960(x_start, x_end, s, G_s, M, h_duct,
                                        T_s, P_s, to_vec(X_s));
      },
      py::arg("x_start"), py::arg("x_end"), py::arg("s"), py::arg("G_s"),
      py::arg("M"), py::arg("h_duct"), py::arg("T_s"), py::arg("P_s"),
      py::arg("X_s"),
      "A wall side's h under a slot film, Seban (1960) (the default):\n"
      "  h = h_duct + w(M) f_x (h_film - h_duct)\n"
      "h_film the span mean of St_s G_s cp_s over [x_start, x_end] [m] at\n"
      "x/s <= 300, f_x its share of the span, w(M) a C1 step from 0 at\n"
      "M = 1.0 to 1 at 1.1 (Seban: flat-plate h below m = 1). Slot-fluid\n"
      "properties at T_s [K], P_s [Pa], X_s (mole fractions); G_s the slot\n"
      "mass flux [kg/(m^2 s)], M = G_s / G_gas, h_duct [W/(m^2 K)].");

  m.def(
      "slot_film_htc_ballal_lefebvre_1973",
      [](double x_start, double x_end, double s, double G_s, double M,
         double h_duct, double T_s, double P_s,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_s) {
        return slot_film_htc_ballal_lefebvre_1973(x_start, x_end, s, G_s, M,
                                                  h_duct, T_s, P_s,
                                                  to_vec(X_s));
      },
      py::arg("x_start"), py::arg("x_end"), py::arg("s"), py::arg("G_s"),
      py::arg("M"), py::arg("h_duct"), py::arg("T_s"), py::arg("P_s"),
      py::arg("X_s"),
      "As slot_film_htc_seban_1960 with the Ballal & Lefebvre (1973) wall-jet\n"
      "Eq. 6 (book Eq. 8.44), h_film = 0.10 (k_s/s) Re_s^0.8 (x/s)^-0.56, at\n"
      "x/s <= 50 and handing back below M = 1.3. The named option.");

  // -----------------------------------------------------------------
  // channel_smooth
  // -----------------------------------------------------------------
  m.def(
      "channel_smooth",
      [](double T, double P,
         py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
         double velocity, double diameter, double length, double T_hot,
         const std::string &correlation, bool heating, double mu_ratio,
         double roughness, double Nu_multiplier, double f_multiplier,
         double flow_area) {
        return channel_smooth(T, P, to_vec(X_arr), velocity, diameter, length,
                              T_hot, correlation, heating, mu_ratio, roughness,
                              Nu_multiplier, f_multiplier, flow_area);
      },
      py::arg("T"), py::arg("P"), py::arg("X"), py::arg("velocity"),
      py::arg("diameter"), py::arg("length"),
      py::arg("T_hot") = std::numeric_limits<double>::quiet_NaN(),
      py::arg("correlation") = "gnielinski", py::arg("heating") = true,
      py::arg("mu_ratio") = 1.0, py::arg("roughness") = 0.0,
      py::arg("Nu_multiplier") = 1.0, py::arg("f_multiplier") = 1.0,
      py::arg("flow_area") = std::numeric_limits<double>::quiet_NaN(),
      "Smooth channel/duct: combined HTC + pressure drop.\n\n"
      "Parameters:\n"
      "  T, P, X    : bulk static thermodynamic state\n"
      "  velocity   : bulk velocity [m/s]\n"
      "  diameter   : hydraulic diameter [m]\n"
      "  length     : channel length [m]\n"
      "  T_hot      : wall temperature [K] (nan to skip q)\n"
      "  correlation: 'gnielinski' (default), 'dittus_boelter', 'sieder_tate', "
      "'petukhov'\n"
      "  heating    : True = fluid heated\n"
      "  mu_ratio   : mu_bulk/mu_wall for Sieder-Tate\n"
      "  roughness  : absolute wall roughness [m]\n"
      "  Nu_multiplier : empirical correction factor on Nu (default 1.0)\n"
      "  f_multiplier  : empirical correction factor on f (default 1.0)\n"
      "  flow_area  : flow cross-section [m^2] velocity was taken on (nan = "
      "pi D^2/4); the mass-flow derivatives use it\n\n"
      "Returns: ChannelResult with Jacobian fields populated");





  // -----------------------------------------------------------------
  // Internal Solver Tools (f, J) exact Jacobians
  // -----------------------------------------------------------------
  m.def("orifice_mdot_and_jacobian", &solver::orifice_mdot_and_jacobian,
        py::arg("dP"), py::arg("rho"), py::arg("Cd"), py::arg("area"),
        py::arg("beta") = 0.0,
        "Fast-path orifice flow and analytic derivative (mdot, "
        "d(mdot)/d(dP)).\n\n"
        "Parameters:\n"
        "  dP   : Pressure drop [Pa]\n"
        "  rho  : Density [kg/m^3]\n"
        "  Cd   : Discharge coefficient [-]\n"
        "  area : Orifice area [m^2]\n"
        "  beta : Diameter ratio beta = d/D [-] (default 0.0)\n\n"
        "Returns: tuple(mass_flow [kg/s], derivative [kg/(s*Pa)])");

  m.def("pressure_loss_and_jacobian", &solver::pressure_loss_and_jacobian,
        py::arg("v"), py::arg("rho"), py::arg("K"),
        "Fast-path K-factor pressure loss and analytic derivative (dP, "
        "d(dP)/d(v)).\n\n"
        "Parameters:\n"
        "  v   : Velocity [m/s]\n"
        "  rho : Density [kg/m^3]\n"
        "  K   : Loss coefficient [-]\n\n"
        "Returns: tuple(dP [Pa], derivative [Pa/(m/s)])");

  m.def("lossless_pressure_and_jacobian",
        &solver::lossless_pressure_and_jacobian, py::arg("P_in"),
        py::arg("P_out"),
        "Calculate ideal lossless outlet pressure P_out = P_in and its "
        "Jacobian dP_out/dP_in.");

  // Heat Transfer (f, J) tuples
  m.def("nusselt_and_jacobian_dittus_boelter",
        &solver::nusselt_and_jacobian_dittus_boelter, py::arg("Re"),
        py::arg("Pr"), py::arg("heating") = true,
        "Calculate Nusselt number and its Jacobian w.r.t Re using "
        "Dittus-Boelter.");

  m.def("nusselt_and_jacobian_gnielinski",
        &solver::nusselt_and_jacobian_gnielinski, py::arg("Re"), py::arg("Pr"),
        py::arg("f"),
        "Calculate Nusselt number and its Jacobian w.r.t Re using Gnielinski.");

  m.def(
      "nusselt_and_jacobian_sieder_tate",
      &solver::nusselt_and_jacobian_sieder_tate, py::arg("Re"), py::arg("Pr"),
      py::arg("mu_ratio"),
      "Calculate Nusselt number and its Jacobian w.r.t Re using Sieder-Tate.");

  m.def("nusselt_and_jacobian_petukhov", &solver::nusselt_and_jacobian_petukhov,
        py::arg("Re"), py::arg("Pr"), py::arg("f"),
        "Calculate Nusselt number and its Jacobian w.r.t Re using Petukhov.");

  // Cooling Correlations (f, J) tuples







  // Stagnation (f, J) tuples
  m.def("mach_number_and_jacobian_v", &solver::mach_number_and_jacobian_v,
        py::arg("v"), py::arg("T"), py::arg("X"),
        "Calculate Mach number and its Jacobian w.r.t Velocity.");

  m.def(
      "T_adiabatic_wall_and_jacobian_v",
      &solver::T_adiabatic_wall_and_jacobian_v, py::arg("T_static"),
      py::arg("v"), py::arg("T"), py::arg("P"), py::arg("X"),
      py::arg("turbulent") = true,
      "Calculate Adiabatic Wall Temperature and its Jacobian w.r.t Velocity.");

  m.def("Tt_from_static_and_jacobian_M", &solver::Tt_from_static_and_jacobian_M,
        py::arg("T"), py::arg("M"), py::arg("X"),
        "Calculate Stagnation Temperature from Static and its Jacobian w.r.t "
        "Mach.");

  m.def(
      "Pt_from_static_and_jacobian_M", &solver::Pt_from_static_and_jacobian_M,
      py::arg("P"), py::arg("T"), py::arg("M"), py::arg("X"),
      py::arg("tol") = 1e-8, py::arg("max_iter") = 50,
      "Calculate Stagnation Pressure from Static and its Jacobian w.r.t Mach.");

  m.def("friction_and_jacobian_haaland", &solver::friction_and_jacobian_haaland,
        py::arg("Re"), py::arg("e_D"),
        "Fast-path Haaland channel friction and analytic derivative (f, "
        "d(f)/d(Re)).\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-]\n"
        "  e_D : Relative roughness eps/D [-]\n\n"
        "Returns: tuple(friction_factor [-], derivative [-])");

  m.def("friction_and_jacobian_serghides",
        &solver::friction_and_jacobian_serghides, py::arg("Re"), py::arg("e_D"),
        "Fast-path Serghides channel friction and analytic derivative (f, "
        "d(f)/d(Re)).\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-]\n"
        "  e_D : Relative roughness eps/D [-]\n\n"
        "Returns: tuple(friction_factor [-], derivative [-])");

  m.def("friction_and_jacobian_colebrook",
        &solver::friction_and_jacobian_colebrook, py::arg("Re"), py::arg("e_D"),
        "Fast-path Colebrook-White channel friction and implicit derivative "
        "(f, "
        "d(f)/d(Re)).\n\n"
        "Parameters:\n"
        "  Re  : Reynolds number [-]\n"
        "  e_D : Relative roughness eps/D [-]\n\n"
        "Returns: tuple(friction_factor [-], derivative [-])");

  m.def("friction_and_jacobian_petukhov",
        &solver::friction_and_jacobian_petukhov, py::arg("Re"),
        "Fast-path Petukhov smooth-channel friction and analytic derivative "
        "(f, "
        "d(f)/d(Re)).\n\n"
        "Parameters:\n"
        "  Re : Reynolds number [-]\n\n"
        "Returns: tuple(friction_factor [-], derivative [-])");

  m.def("friction_and_jacobian", &solver::friction_and_jacobian, py::arg("tag"),
        py::arg("Re"), py::arg("e_D"),
        "Dispatcher: selects friction_and_jacobian_* by tag string.\n\n"
        "tag must be one of: 'haaland', 'serghides', 'colebrook', "
        "'petukhov'\n\n"
        "Parameters:\n"
        "  tag : friction model name\n"
        "  Re  : Reynolds number [-]\n"
        "  e_D : Relative roughness eps/D [-]\n\n"
        "Returns: tuple(friction_factor [-], d(f)/d(Re) [-])");

  m.def("density_and_jacobians", &solver::density_and_jacobians, py::arg("T"),
        py::arg("P"), py::arg("X"),
        "Fast-path density and analytic derivatives (rho, d(rho)/d(T), "
        "d(rho)/d(P)).\n\n"
        "Parameters:\n"
        "  T : Temperature [K]\n"
        "  P : Pressure [Pa]\n"
        "  X : Mole fractions [-]\n\n"
        "Returns: tuple(rho [kg/m^3], d_rho_d_T [(kg/m^3)/K], d_rho_d_P "
        "[(kg/m^3)/Pa])");

  m.def("enthalpy_and_jacobian", &solver::enthalpy_and_jacobian, py::arg("T"),
        py::arg("X"),
        "Fast-path enthalpy and analytic derivative (h, d(h)/d(T)).\n\n"
        "Parameters:\n"
        "  T : Temperature [K]\n"
        "  X : Mole fractions [-]\n\n"
        "Returns: tuple(h [J/kg], cp_mass [J/(kg*K)])");

  m.def("viscosity_and_jacobians", &solver::viscosity_and_jacobians,
        py::arg("T"), py::arg("P"), py::arg("X"),
        "Fast-path viscosity and central difference derivatives (mu, "
        "d(mu)/d(T), d(mu)/d(P)).\n\n"
        "Parameters:\n"
        "  T : Temperature [K]\n"
        "  P : Pressure [Pa]\n"
        "  X : Mole fractions [-]\n\n"
        "Returns: tuple(mu [Pa*s], d_mu_d_T [(Pa*s)/K], d_mu_d_P "
        "[(Pa*s)/Pa])");

  m.def("_viscosity_mass_gradient", &solver::viscosity_mass_gradient, py::arg("T"),
        py::arg("P"), py::arg("X"),
        "Internal (#964): d(mu)/dY_k [Pa*s] at fixed (T, P), analytic (Wilke's "
        "rule chained to the mass fractions). X: mole fractions of the state.");

  m.def(
      "equivalence_ratio",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr) {
        return equivalence_ratio(to_vec(X_arr));
      },
      py::arg("X"),
      "Calculate the equivalence ratio (phi) of a mixture based on its "
      "elemental composition.\n\n"
      "Phi = O2_required / O2_available\n\n"
      "Parameters:\n"
      "  X : Mole fractions [-]\n\n"
      "Returns: Phi [-]");

  m.def(
      "equivalence_ratio_mass",
      [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_arr) {
        return equivalence_ratio_mass(to_vec(Y_arr));
      },
      py::arg("Y"),
      "Calculate the equivalence ratio (phi) of a mixture based on its "
      "elemental composition (mass basis).\n\n"
      "Parameters:\n"
      "  Y : Mass fractions [-]\n\n"
      "Returns: Phi [-]");

  m.def("adiabatic_T_complete_and_jacobian_T",
        &solver::adiabatic_T_complete_and_jacobian_T, py::arg("T_in"),
        py::arg("P"), py::arg("X_in"),
        "Fast-path adiabatic flame temperature (complete combustion) and "
        "analytic derivative.\n\n"
        "Parameters:\n"
        "  T_in : Inlet Temperature [K]\n"
        "  P    : Pressure [Pa]\n"
        "  X_in : Inlet Mole fractions [-]\n\n"
        "Returns: tuple(T_ad [K], dT_ad_dT_in [-], X_products [-])");

  m.def("adiabatic_T_equilibrium_and_jacobians",
        &solver::adiabatic_T_equilibrium_and_jacobians, py::arg("T_in"),
        py::arg("P"), py::arg("X_in"),
        "Fast-path adiabatic flame temperature (equilibrium combustion) and "
        "analytic derivatives.\n\n"
        "Parameters:\n"
        "  T_in : Inlet Temperature [K]\n"
        "  P    : Pressure [Pa]\n"
        "  X_in : Inlet Mole fractions [-]\n\n"
        "Returns: tuple(T_ad [K], dT_ad_dT_in [-], dT_ad_dP [K/Pa], X_products "
        "[-])");

  // Registry bindings
  py::class_<Registry>(m, "Registry",
                       "Global function registry for C++ physics models")
      .def_static("get", &Registry::get, py::return_value_policy::reference,
                  "Get the singleton registry instance")
      .def("has_friction_model", &Registry::has_friction_model, py::arg("name"),
           "Check if friction model is registered")
      .def("available_friction_models", &Registry::available_friction_models,
           "List all available friction models");

  // Drop any Python warning handler while the interpreter is still alive.
  //
  // set_warning_handler() stores its argument in a function-local static in
  // correlation_status.cpp. When that argument is a Python callable, the
  // static owns a py::object, and a function-local static is destroyed by
  // __cxa_atexit -- which runs AFTER Py_Finalize(). Releasing a Python
  // reference on a finalised interpreter segfaults, so merely installing a
  // handler and leaving it there crashed the process at exit, including the
  // `set_warning_handler(lambda msg: None)` this module's own docstring
  // recommends.
  //
  // A capsule stored on the module has its destructor run during module
  // teardown, before finalisation completes. Resetting to the default
  // handler there drops the py::object at a point where decref is still
  // legal, and leaves the static holding a plain C++ function.
  m.add_object("_warning_handler_cleanup",
               py::capsule([]() { combaero::set_warning_handler({}); }));
}
