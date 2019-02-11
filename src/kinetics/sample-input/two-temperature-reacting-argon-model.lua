model = "TwoTemperatureReactingArgon"

TwoTemperatureReactingArgon = {
  ion_tol = 1.0e-15,
  chem_dt = 1.0e-11,
  integration_method = "Backward_Euler",
  -- integration_method = "RK4",
  -- integration_method = "Forward_Euler",
  Newton_Raphson_tol = 1.0e-10
}

