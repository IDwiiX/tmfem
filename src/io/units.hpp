#pragma once

namespace tmfem {

// Case files give moduli in GPa and stresses in MPa; the loaders convert at the input boundary.
inline constexpr double kPaPerMPa = 1e6;
inline constexpr double kPaPerGPa = 1e9;

}  // namespace tmfem
