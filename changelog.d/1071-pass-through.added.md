- **Radiation now reaches the walls and panels that the assembly helpers build (#1071).**
  - `ImpingementArray.add_wall(..., radiation_hot=..., radiation_target=...)` puts a `GasRadiation` or `CasingRadiation` on each row's `ThermalWall`. It maps the hot face and the target face to sides A and B by `array_is_side_a`.
  - `EffusionLiner(radiation=GasRadiation(...))` hands the spec to every panel's hot face.

  Both default to None, which is unchanged behaviour.
