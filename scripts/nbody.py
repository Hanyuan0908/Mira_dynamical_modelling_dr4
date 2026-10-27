"""Reader and units for the isolated N-body snapshot (NEMO format, mock/Nbody/t_t0_7).

Snapshot units: length kpc, velocity km/s, G = 1, so the mass unit is 1/G = 232,508 Msun and the time
unit is kpc/(km/s) = 0.9778 Gyr. Stored items (float32): Mass, Position, Velocity, Potential,
Acceleration. Components are identified by particle mass: lightest = disc, middle = bulge,
heaviest = dark matter.
"""
import numpy as np

from common import MOCK_DIR

SNAPSHOT = MOCK_DIR / "Nbody" / "t_t0_7"
G_KPC_KMS_MSUN = 4.300917270e-6
MASS_UNIT_MSUN = 1.0 / G_KPC_KMS_MSUN
TIME_UNIT_GYR = 0.977792221
COMP_CODE = {"disc": 0, "bulge": 1, "dm": 2}
_DTYPE = {b"f": np.float32, b"d": np.float64, b"i": np.int32}


def read_nemo(path=SNAPSHOT, tags=("Mass", "Position", "Velocity", "Potential", "Acceleration")):
    """Read array items of a single-snapshot NEMO file: magic, type, tag, int32 dims ending in 0, data."""
    raw = open(path, "rb").read()
    out = {}
    for tag in tags:
        key = b"\x00" + tag.encode() + b"\x00"
        i = raw.find(key)
        if i < 0:
            raise KeyError(tag)
        typ = raw[i - 1:i]                       # single-letter type code before the tag
        p = i + len(key)
        dims = []
        while True:
            d = int.from_bytes(raw[p:p + 4], "little")
            p += 4
            if d == 0:
                break
            dims.append(d)
        n = int(np.prod(dims))
        out[tag] = np.frombuffer(raw, dtype=_DTYPE[typ], count=n, offset=p).reshape(dims).copy()
    t = raw.find(b"\x00Time\x00")
    out["Time"] = float(np.frombuffer(raw, dtype=_DTYPE[raw[t - 1:t]], count=1, offset=t + 6)[0])
    return out


def components(mass):
    """Component code per particle from the three distinct particle masses."""
    mu = np.unique(mass)
    if len(mu) != 3:
        raise ValueError(f"expected 3 particle masses, found {mu}")
    code = np.empty(len(mass), np.int8)
    for name, m in zip(("disc", "bulge", "dm"), mu):
        code[mass == m] = COMP_CODE[name]
    return code
