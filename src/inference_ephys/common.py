"""Identifiers, deterministic random streams, and file output."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

UNIT = ["Group", "Mouse_ID", "Session_ID", "Region", "Unit_ID", "Stage"]
ANIMAL = ["Group", "Mouse_ID", "Region", "Stage"]
BIN = ["Time_Bin_Start_s", "Time_Bin_End_s"]
TRIAL = ["Trial_ID", "Trial_Order"]

def rng(seed, *labels):
    token = json.dumps([int(seed), *map(str, labels)], ensure_ascii=True).encode()
    return np.random.default_rng(int.from_bytes(hashlib.sha256(token).digest()[:8], "little"))

def checksum(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(2**20), b""):
            h.update(block)
    return h.hexdigest()

def save(frame, root, name):
    dest = Path(root) / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(dest, index=False, encoding="utf-8", float_format="%.12g")
    return frame

def write_json(value, path):
    def convert(x):
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, Path):
            return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=convert), encoding="utf-8")

def correlation(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or len(x) != len(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
        return np.nan
    if np.ptp(x) <= 1e-12 or np.ptp(y) <= 1e-12:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])

def frame(rows, columns=None):
    return pd.DataFrame(rows) if rows else pd.DataFrame(columns=columns or [])
