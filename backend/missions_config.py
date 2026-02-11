# missions_config.py
MISSION_DEFINITIONS = {
    "m1": {
        "title": "Stealth Black",
        "field": "Optical Stealth / NIR Control",
        "parameters": {
            "nir_reflectance": {"label": "NIR反射率", "unit": "%", "target": "< 10"},
            "scattering_albedo": {"label": "散乱アルベド", "unit": "", "target": "0.1 - 0.3"}
        }
    },
    "m2": {
        "title": "Void Black",
        "field": "Deep Absorption",
        "parameters": {
            "absorption_coeff": {"label": "吸収係数", "unit": "cm⁻¹", "target": "> 5000"},
            "aspect_ratio": {"label": "構造アスペクト比", "unit": "", "target": "100 - 500"}
        }
    },
    "m3": {
        "title": "Thermal Black",
        "field": "Thermodynamics / Emissivity",
        "parameters": {
            "emissivity": {"label": "放射率 ε", "unit": "", "target": "< 0.05"},
            "thermal_diffusivity": {"label": "熱拡散率", "unit": "mm²/s", "target": "< 0.1"}
        }
    },
    "m4": {
        "title": "Meta Black",
        "field": "Metamaterials / Wave Control",
        "parameters": {
            "refractive_index": {"label": "負の屈折率 n", "unit": "", "target": "-1.5 - -0.5"},
            "phase_offset": {"label": "位相オフセット", "unit": "rad", "target": "π/2 - π"}
        }
    }
}