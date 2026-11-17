#include "registry.h"
#include <stdexcept>

namespace combaero {

Registry::Registry() {
  friction_models_["haaland"] = solver::friction_and_jacobian_haaland;
  friction_models_["serghides"] = solver::friction_and_jacobian_serghides;
  friction_models_["colebrook"] = solver::friction_and_jacobian_colebrook;
  friction_models_["petukhov"] = [](double Re, double) {
    return solver::friction_and_jacobian_petukhov(Re);
  };
  const double blended = solver::laminar_blend_f_re_at_zero();
  friction_f_re_at_zero_["haaland"] = blended;
  friction_f_re_at_zero_["serghides"] = blended;
  friction_f_re_at_zero_["colebrook"] = blended;
  friction_f_re_at_zero_["petukhov"] = blended;
}

Registry::FrictionFunc
Registry::get_friction_model(const std::string &name) const {
  auto it = friction_models_.find(name);
  if (it != friction_models_.end()) {
    return it->second;
  }
  throw std::invalid_argument("Unknown friction model registry_id: " + name);
}

double Registry::friction_f_re_at_zero(const std::string &name) const {
  auto it = friction_f_re_at_zero_.find(name);
  if (it != friction_f_re_at_zero_.end()) {
    return it->second;
  }
  throw std::invalid_argument("Unknown friction model registry_id: " + name);
}

bool Registry::has_friction_model(const std::string &name) const {
  return friction_models_.find(name) != friction_models_.end();
}

std::vector<std::string> Registry::available_friction_models() const {
  std::vector<std::string> models;
  models.reserve(friction_models_.size());
  for (const auto &pair : friction_models_) {
    models.push_back(pair.first);
  }
  return models;
}

} // namespace combaero
