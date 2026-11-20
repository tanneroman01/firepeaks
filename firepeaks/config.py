"""Settings and local secrets (eg API keys, gee project IDs)

load()      config/default.yaml as a dict, plus every `paths:` entry resolved to an absolute Path
ee_init()   Earth Engine initialised with the project in .env (EE_PROJECT); never hard-code the project id

.env lives at the repo root, is gitignored, and holds KEY=VALUE lines (see .env.example)
"""
from __future__ import annotations
import os
from pathlib import Path
import yaml

REPO = Path(__file__).resolve().parents[1]


def load(path: str | Path | None = None) -> tuple[dict, dict[str, Path]]:
    """(cfg, P): the yaml as a dict, and every `paths:` entry as an absolute Path"""
    path = Path(path) if path else REPO / "config" / "default.yaml"
    cfg = yaml.safe_load(path.read_text())
    P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()}
    return cfg, P


def load_env(path: str | Path | None = None) -> None:
    """read KEY=VALUE lines from .env into os.environ (blank lines and # comments ignored; quotes stripped)"""
    path = Path(path) if path else REPO / ".env"
    if not path.exists():
        return
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip().strip("'\""))


def env(key: str) -> str:
    load_env()
    v = os.environ.get(key)
    if not v:
        raise RuntimeError(f"{key} is not set: add `{key}=...` to {REPO / '.env'} (see .env.example)")
    return v


def ee_init(project: str | None = None) -> None:
    """ee.Initialize with the given project, else EE_PROJECT from .env"""
    import ee
    ee.Initialize(project=project or env("EE_PROJECT"))
