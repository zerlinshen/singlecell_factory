from __future__ import annotations

import hashlib
import logging
import os
import socket
import subprocess
import sys
import threading
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

_ENV_PREFIXES = ("SC_", "SCF_", "CUDA_", "CONDA_")
_ENV_ALLOWLIST = frozenset({"HOSTNAME", "USER", "HOME", "LANG", "LC_ALL", "TZ", "PYTHONPATH"})


def _filter_env() -> dict[str, str]:
    result: dict[str, str] = {}
    for key, val in os.environ.items():
        if any(key.startswith(p) for p in _ENV_PREFIXES) or key in _ENV_ALLOWLIST:
            result[key] = val
    return result


def _git_info() -> tuple[str, str, bool]:
    """Return (sha, branch, dirty). Falls back to empty strings on error."""
    try:
        sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        dirty_out = subprocess.check_output(
            ["git", "status", "--porcelain"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        return sha, branch, bool(dirty_out)
    except Exception:
        return "", "", False


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


class _RSSPoller:
    """Background thread that samples RSS of the current process at 1s intervals."""

    def __init__(self) -> None:
        self._peak: int = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._poll, daemon=True, name="ledger-rss-poller")

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        self._thread.join(timeout=3)

    @property
    def peak_bytes(self) -> int:
        return self._peak

    def _poll(self) -> None:
        try:
            import psutil  # type: ignore
            proc = psutil.Process()
            while not self._stop.wait(1.0):
                try:
                    rss = proc.memory_info().rss
                    if rss > self._peak:
                        self._peak = rss
                except Exception:
                    pass
        except ImportError:
            pass


class RunLedger:
    """Observational-only audit trail for a single pipeline run.

    All public methods are wrapped to be non-raising — a ledger failure must
    never crash the production run.
    """

    def __init__(self, ctx, project_name: str, output_dir: Path) -> None:
        self._ctx = ctx
        self._project = project_name
        self._output_dir = output_dir
        self._start_ts: str = ""
        self._record: dict = {}
        self._module_results: list[dict] = []
        self._rss_poller = _RSSPoller()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def record_start(self) -> None:
        try:
            self._start_ts = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
            sha, branch, dirty = _git_info()
            cfg = getattr(self._ctx, "cfg", None)
            self._record = {
                "schema_version": "run_ledger_v1",
                "timestamp_utc": self._start_ts,
                "project": self._project,
                "cli_args": list(sys.argv),
                "env": _filter_env(),
                "git_sha": sha,
                "git_branch": branch,
                "git_dirty": dirty,
                "conda_env": os.environ.get("CONDA_DEFAULT_ENV", ""),
                "python_version": sys.version,
                "hostname": socket.gethostname(),
                "sample_root": str(cfg.cellranger.sample_root) if cfg and hasattr(cfg, "cellranger") else "",
                "optional_modules": list(cfg.optional_modules) if cfg and hasattr(cfg, "optional_modules") else [],
            }
            self._rss_poller.start()
        except Exception as exc:
            logger.warning("RunLedger.record_start failed: %s", exc)

    def record_module(
        self,
        name: str,
        status: str,
        message: str,
        wall_seconds: float,
        rss_peak_bytes: int,
    ) -> None:
        try:
            self._module_results.append(
                {
                    "name": name,
                    "status": status,
                    "message": message,
                    "wall_seconds": round(float(wall_seconds), 3),
                    "rss_peak_bytes": int(rss_peak_bytes),
                }
            )
        except Exception as exc:
            logger.warning("RunLedger.record_module failed: %s", exc)

    def record_end(self, final_adata_path: Path | None) -> None:
        try:
            self._rss_poller.stop()
            end_ts = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")

            # Compute total wall seconds from timestamps
            try:
                from datetime import datetime as _dt
                t0 = _dt.fromisoformat(self._start_ts.replace("Z", "+00:00"))
                t1 = _dt.fromisoformat(end_ts.replace("Z", "+00:00"))
                total_wall = round((t1 - t0).total_seconds(), 3)
            except Exception:
                total_wall = None

            sha256 = ""
            if final_adata_path is not None and Path(final_adata_path).exists():
                try:
                    sha256 = _sha256_file(Path(final_adata_path))
                except Exception:
                    pass

            self._record.update(
                {
                    "end_timestamp_utc": end_ts,
                    "total_wall_seconds": total_wall,
                    "peak_rss_bytes": self._rss_poller.peak_bytes,
                    "final_adata_sha256": sha256,
                    "module_results": self._module_results,
                }
            )
        except Exception as exc:
            logger.warning("RunLedger.record_end failed: %s", exc)

    def write(self) -> Path | None:
        try:
            import json

            ledger_dir = Path(self._output_dir) / "ops" / "run_ledger"
            ledger_dir.mkdir(parents=True, exist_ok=True)

            ts_tag = self._start_ts.replace(":", "").replace("-", "").replace("Z", "")[:15]
            filename = f"{self._project}_{ts_tag}.json"
            out_path = ledger_dir / filename
            out_path.write_text(json.dumps(self._record, indent=2, ensure_ascii=False), encoding="utf-8")
            logger.info("RunLedger written: %s", out_path)
            return out_path
        except Exception as exc:
            logger.warning("RunLedger.write failed: %s", exc)
            return None
