model = "MultiTemperatureGas"
species = {"N2", "O2", "N", "O", "NO"}
energy_modes = {"N2", "O2", "NO"}
mode_components = {
   N2 = {'N2:vib', 'N2:electronic'},
   O2 = {'O2:vib', 'O2:electronic', 'O:electronic'},
   NO = {'NO:vib', 'NO:electronic', 'N:electronic'}
}
