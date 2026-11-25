#ifndef COMBAERO_REGISTRY_H
#define COMBAERO_REGISTRY_H

#include "solver_interface.h"
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace combaero {

class Registry {
public:
  // Unified signature for friction models used by the dispatcher.
  // Note: smooth-pipe models (e.g., Petukhov) do not use e_D, but keep this
  // parameter to avoid special casing in the registry and call sites.
  using FrictionFunc =
      std::function<solver::CorrelationResult<std::tuple<double, double>>(
          double Re, double e_D)>;
  // Regularized F_Re(Re) = f(Re) * Re and dF_Re/dRe without the laminar
  // re_floor (#1216), registered alongside FrictionFunc as one composition
  // point for friction models.
  using FrictionFReFunc =
      std::function<std::tuple<double, double>(double Re, double e_D)>;

  // Get the global singleton instance
  static Registry &get() {
    static Registry instance;
    return instance;
  }

  // Retrieve a friction function by name
  FrictionFunc get_friction_model(const std::string &name) const;

  // Retrieve the F_Re(Re) = f(Re) * Re function by name (#1216)
  FrictionFReFunc get_friction_f_re_model(const std::string &name) const;

  // Check if name exists
  bool has_friction_model(const std::string &name) const;

  // lim_{Re -> 0} f Re of a registered friction model (#1156): the laminar
  // blend's share of Hagen-Poiseuille's f Re = 64 at zero flow
  // (solver::laminar_blend_f_re_at_zero), for every registered model (all four
  // are blended since #1210; Petukhov was the raw fit, f Re -> 0). A duct's drop
  // goes as m f(Re) Re, so this is its slope at exactly zero flow, where the
  // kernel's |m| f is 0 * inf.
  double friction_f_re_at_zero(const std::string &name) const;

  // Get a list of all available friction models
  std::vector<std::string> available_friction_models() const;

private:
  Registry(); // Private constructor builds the mappings

  std::map<std::string, FrictionFunc> friction_models_;
  std::map<std::string, FrictionFReFunc> friction_f_re_models_;
  std::map<std::string, double> friction_f_re_at_zero_;
};

} // namespace combaero

#endif // COMBAERO_REGISTRY_H
