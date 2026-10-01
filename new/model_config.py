"""Serializable controls for the single-pad model and its static coupling."""
from dataclasses import asdict, dataclass
import math


@dataclass
class ModelConfig:
    mechanics_source: str = 'native'
    ansys_csv: str = 'old/ansys-results/250_500nm_3.5_SiO2.csv'
    ansys_interpolation: str = 'linear'
    temperature_C: float = 250.
    reference_C: float = 22.
    plasticity: bool = True
    duration_s: float = 7000.
    dt_s: float = 1.
    min_dt_s: float = .03125
    max_dt_s: float = 1.
    measure_interval_s: float = 50.
    save_interval_s: float = 100.
    mechanics_nr: int = 80
    mechanics_nz: int = 80
    thermal_steps: int = 40
    phase_nr: int = 100
    phase_nz: int = 220
    diffusivity_factor: float = .1
    interface_width_nm: float = 5.
    wall_energy_J_m2: float = -.651
    mechanical_energy_multiplier: float = 1.
    write_fields: bool = True
    # Freeze the previously studied SiO2 geometry and material set for this bridge.
    # More general geometry/material sweeps remain available in mechanics drivers.

    def validate(self):
        if self.mechanics_source not in ('native', 'ansys', 'none'):
            raise ValueError('mechanics_source must be native, ansys, or none')
        if self.ansys_interpolation not in ('linear', 'nearest'):
            raise ValueError('ansys_interpolation must be linear or nearest')
        for key,value in asdict(self).items():
            if isinstance(value,float) and not math.isfinite(value):
                raise ValueError(f'{key} must be finite')
        if self.duration_s < 0 or not 0 < self.min_dt_s <= self.dt_s <= self.max_dt_s:
            raise ValueError('Require duration >= 0 and 0 < min_dt <= dt <= max_dt')
        if self.temperature_C <= -273.15 or self.reference_C <= -273.15:
            raise ValueError('Temperatures must exceed absolute zero')
        if min(self.measure_interval_s,self.save_interval_s,self.interface_width_nm) <= 0:
            raise ValueError('Intervals and interface width must be positive')
        if self.diffusivity_factor < 0 or self.mechanical_energy_multiplier < 0:
            raise ValueError('Mobility and energy multipliers must be nonnegative')
        if abs(self.wall_energy_J_m2)>1.7:
            raise ValueError('Wall energy requires |Delta gamma| <= gamma_Cu = 1.7 J/m²')
        for key in ('mechanics_nr','mechanics_nz','phase_nr','phase_nz','thermal_steps'):
            value=getattr(self,key)
            if not isinstance(value,int) or isinstance(value,bool) or value<1:
                raise ValueError(f'{key} must be a positive integer')
        if self.mechanics_nr<4 or self.mechanics_nr%2 or self.mechanics_nz<4:
            raise ValueError('Mechanics nr must be even and >=4; nz must be >=4')
        if min(self.phase_nr,self.phase_nz)<4:
            raise ValueError('Phase mesh must have at least four cells per direction')
        if not isinstance(self.plasticity,bool) or not isinstance(self.write_fields,bool):
            raise ValueError('plasticity and write_fields must be booleans')
        return self

    def to_dict(self):
        return asdict(self)


def smoke_config():
    return ModelConfig(duration_s=.1,dt_s=.05,min_dt_s=.001,max_dt_s=.05,
                       measure_interval_s=.05,save_interval_s=.1,mechanics_nr=12,
                       mechanics_nz=12,thermal_steps=4,phase_nr=20,phase_nz=44)
