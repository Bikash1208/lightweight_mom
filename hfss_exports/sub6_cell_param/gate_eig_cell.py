"""Cell 2 SD eigenmode runner with exact gate-confirm impedance states."""
import cell_param_pyaedt as cp
import cell2_cell as base

EXACT_C_STATES = {"OFF_C013": (29.56, -390.9), "OFF_C018": (15.45, -275.0)}


def configure():
    base.configure()
    cp.SENS_STATES.update(EXACT_C_STATES)


def values(*args, **kwargs):
    return base.values(*args, **kwargs)


if __name__ == "__main__":
    configure()
    base.main()
