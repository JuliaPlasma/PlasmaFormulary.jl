"""
This module provides aliases of the most common plasma functionality for user convenience. 

Aliases are denoted with a trailing underscore (e.g., `alias_`).
"""
module Aliases
using ..PlasmaFormularyZ

# Ref: https://docs.plasmapy.org/en/stable/formulary/drifts.html
const vd_ = diamagnetic_drift
const veb_ = ExB_drift
const vfd_ = force_drift
export vd_, veb_, vfd_

# Ref: https://docs.plasmapy.org/en/stable/formulary/lengths.html
const rc_ = gyroradius
const rhoc_ = gyroradius
const cwp_ = inertial_length
export rc_, rhoc_, cwp_

# https://docs.plasmapy.org/en/stable/formulary/speeds.html
const cs_ = ion_sound_speed
const va_ = Alfven_speed
const vth_ = thermal_velocity
export cs_, va_, vth_
end