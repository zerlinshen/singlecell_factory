from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.image as mpimg
import numpy as np
import pandas as pd

from ..context import PipelineContext

logger = logging.getLogger(__name__)

# Parsed by scripts/update_references.py for citation tracking.
__references__ = {
    "reproducible_research_principles": {
        "title": "Ten Simple Rules for Reproducible Computational Research",
        "authors": "Sandve et al.",
        "journal": "PLOS Computational Biology",
        "year": "2013",
        "doi": "10.1371/journal.pcbi.1003285",
        "description": "Guides provenance capture and reproducibility reporting.",
    }
}


class PaperReproModule:
    """Optional module: paper-driven reproduction ledger and figure parity checks.

    The module does not alter adata structure. It focuses on:
    1) tracking provenance for external papers/repos used to evolve modules
    2) validating expected figure reproduction against current pipeline outputs
    """

    name = "paper_repro"
    mutates_structure = False
    requires_keys: dict[str, list[str]] = {}
    provides_keys: dict[str, list[str]] = {"uns": ["paper_repro_report"]}

    def run(self, ctx: PipelineContext) -> None:
        cfg = getattr(ctx.cfg, "paper_repro", None)
        spec_path = getattr(cfg, "spec_json", None)
        strict = bool(getattr(cfg, "strict", False))

        if spec_path is None:
            template_path = ctx.table_dir / "paper_repro_spec.template.json"
            template_path.write_text(
                json.dumps(self._default_template(), indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            msg = (
                "No --paper-spec-json provided. Generated "
                "paper_repro_spec.template.json."
            )
            ctx.metadata["paper_repro_status"] = "skipped_no_spec"
            ctx.metadata["paper_repro_template"] = str(template_path)
            ctx.status(self.name, "skipped", msg)
            return

        spec_path = Path(spec_path)
        if not spec_path.exists():
            raise FileNotFoundError(f"paper_repro spec not found: {spec_path}")

        payload = json.loads(spec_path.read_text(encoding="utf-8"))
        papers = payload.get("papers", [])
        if not isinstance(papers, list) or not papers:
            raise ValueError("paper_repro spec must contain a non-empty 'papers' list.")

        registry_rows: list[dict[str, Any]] = []
        figure_rows: list[dict[str, Any]] = []
        spec_root = spec_path.parent

        for item in papers:
            if not isinstance(item, dict):
                continue
            row, fig_checks = self._process_one_paper(item, spec_root=spec_root, run_dir=ctx.run_dir)
            registry_rows.append(row)
            figure_rows.extend(fig_checks)

        registry_df = pd.DataFrame(registry_rows)
        figures_df = pd.DataFrame(figure_rows)
        registry_csv = ctx.table_dir / "paper_repro_registry.csv"
        figures_csv = ctx.table_dir / "paper_repro_figures.csv"
        registry_df.to_csv(registry_csv, index=False)
        figures_df.to_csv(figures_csv, index=False)

        strict_failures = int((~registry_df["source_complete"]).sum()) if not registry_df.empty else 0
        if not figures_df.empty:
            strict_failures += int((figures_df["status"] != "pass").sum())

        report = {
            "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
            "spec_path": str(spec_path),
            "papers_total": int(len(registry_df)),
            "papers_with_complete_source": int(registry_df["source_complete"].sum()) if not registry_df.empty else 0,
            "figure_checks_total": int(len(figures_df)),
            "figure_checks_passed": int((figures_df["status"] == "pass").sum()) if not figures_df.empty else 0,
            "figure_checks_failed": int((figures_df["status"] != "pass").sum()) if not figures_df.empty else 0,
            "strict_mode": strict,
            "strict_failures": int(strict_failures),
            "next_actions": self._next_actions(registry_df, figures_df),
        }
        report_path = ctx.table_dir / "paper_repro_report.json"
        report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

        ctx.metadata["paper_repro_status"] = "completed"
        ctx.metadata["paper_repro_report"] = {
            "papers_total": report["papers_total"],
            "figure_checks_total": report["figure_checks_total"],
            "figure_checks_passed": report["figure_checks_passed"],
            "strict_failures": report["strict_failures"],
            "report_path": str(report_path),
        }
        if ctx.adata is not None:
            ctx.adata.uns["paper_repro_report"] = report

        if strict and strict_failures > 0:
            raise ValueError(
                f"paper_repro strict mode failed with {strict_failures} unresolved checks."
            )

    @staticmethod
    def _default_template() -> dict[str, Any]:
        return {
            "papers": [
                {
                    "paper_id": "example_2026_method",
                    "title": "Example Paper For Method Reproduction",
                    "source": {
                        "paper_path": "docs/papers/example.pdf",
                        "repo_url": "https://github.com/example/repro-repo",
                        "repo_commit": "abcdef1234567890",
                        "license": "MIT",
                    },
                    "adaptation_targets": ["clustering", "trajectory"],
                    "reproductions": [
                        {
                            "figure_id": "fig2a_umap",
                            "reference_path": "docs/papers/example_fig2a.png",
                            "pipeline_output": "clustering/umap_leiden.png",
                            "metric": "mae",
                            "min_score": 0.75,
                        }
                    ],
                    "notes": "Track what was adapted and what still diverges.",
                }
            ]
        }

    def _process_one_paper(
        self,
        item: dict[str, Any],
        *,
        spec_root: Path,
        run_dir: Path,
    ) -> tuple[dict[str, Any], list[dict[str, Any]]]:
        paper_id = str(item.get("paper_id") or item.get("id") or "unknown_paper")
        title = str(item.get("title") or "")
        source = item.get("source", {}) if isinstance(item.get("source"), dict) else {}
        paper_path = self._resolve_path(source.get("paper_path"), base=spec_root)
        repo_url = str(source.get("repo_url") or "")
        repo_commit = str(source.get("repo_commit") or "")
        license_name = str(source.get("license") or "")

        targets = item.get("adaptation_targets", [])
        if not isinstance(targets, list):
            targets = []
        target_text = ",".join(str(t) for t in targets)

        source_complete = bool(repo_url and repo_commit and license_name)
        paper_exists = bool(paper_path and paper_path.exists())

        row = {
            "paper_id": paper_id,
            "title": title,
            "paper_path": str(paper_path) if paper_path else "",
            "paper_exists": paper_exists,
            "repo_url": repo_url,
            "repo_commit": repo_commit,
            "license": license_name,
            "adaptation_targets": target_text,
            "source_complete": source_complete,
        }

        checks: list[dict[str, Any]] = []
        reproductions = item.get("reproductions", [])
        if not isinstance(reproductions, list):
            reproductions = []

        for rep in reproductions:
            if not isinstance(rep, dict):
                continue
            checks.append(
                self._run_one_figure_check(
                    paper_id=paper_id, rep=rep, spec_root=spec_root, run_dir=run_dir
                )
            )
        return row, checks

    def _run_one_figure_check(
        self,
        *,
        paper_id: str,
        rep: dict[str, Any],
        spec_root: Path,
        run_dir: Path,
    ) -> dict[str, Any]:
        figure_id = str(rep.get("figure_id") or "unknown_figure")
        metric = str(rep.get("metric") or "mae").lower()
        default_min = 0.75 if metric == "mae" else 0.60
        min_score = float(rep.get("min_score", default_min))

        ref_path = self._resolve_path(rep.get("reference_path"), base=spec_root)
        out_path = self._resolve_path(rep.get("pipeline_output"), base=run_dir)
        ref_exists = bool(ref_path and ref_path.exists())
        out_exists = bool(out_path and out_path.exists())

        score: float | None = None
        status = "pending"
        error = ""
        if not ref_exists:
            status = "missing_reference"
        elif not out_exists:
            status = "missing_pipeline_output"
        else:
            try:
                score = self._compare_images(ref_path, out_path, metric=metric)
                status = "pass" if score >= min_score else "below_threshold"
            except Exception as exc:  # defensive: corrupted image or decode issues
                status = "comparison_error"
                error = str(exc)

        return {
            "paper_id": paper_id,
            "figure_id": figure_id,
            "metric": metric,
            "min_score": min_score,
            "score": round(float(score), 4) if score is not None else None,
            "status": status,
            "reference_path": str(ref_path) if ref_path else "",
            "reference_exists": ref_exists,
            "pipeline_output": str(out_path) if out_path else "",
            "pipeline_output_exists": out_exists,
            "error": error,
        }

    @staticmethod
    def _resolve_path(raw: Any, *, base: Path) -> Path | None:
        if raw is None:
            return None
        text = str(raw).strip()
        if not text:
            return None
        path = Path(text)
        return path if path.is_absolute() else (base / path)

    def _compare_images(self, ref_path: Path, out_path: Path, *, metric: str) -> float:
        ref = self._load_image_gray(ref_path)
        out = self._load_image_gray(out_path)
        h = min(ref.shape[0], out.shape[0])
        w = min(ref.shape[1], out.shape[1])
        if h == 0 or w == 0:
            raise ValueError("empty image after shape alignment")
        ref = ref[:h, :w]
        out = out[:h, :w]

        if metric == "pearson":
            return self._pearson_similarity(ref, out)
        return self._mae_similarity(ref, out)

    @staticmethod
    def _load_image_gray(path: Path) -> np.ndarray:
        arr = np.asarray(mpimg.imread(path), dtype=np.float32)
        if arr.ndim == 3:
            if arr.shape[-1] == 4:
                arr = arr[..., :3]
            arr = arr.mean(axis=-1)
        if arr.ndim != 2:
            raise ValueError(f"unsupported image dimensions for {path}: {arr.shape}")
        if arr.size == 0:
            return arr
        vmax = float(np.max(arr))
        if vmax > 1.0:
            arr = arr / 255.0
        return arr

    @staticmethod
    def _pearson_similarity(a: np.ndarray, b: np.ndarray) -> float:
        flat_a = a.ravel()
        flat_b = b.ravel()
        if flat_a.size == 0 or flat_b.size == 0:
            return 0.0
        std_a = float(np.std(flat_a))
        std_b = float(np.std(flat_b))
        if std_a == 0.0 or std_b == 0.0:
            return 1.0 if np.allclose(flat_a, flat_b) else 0.0
        corr = float(np.corrcoef(flat_a, flat_b)[0, 1])
        if np.isnan(corr):
            return 0.0
        return corr

    @staticmethod
    def _mae_similarity(a: np.ndarray, b: np.ndarray) -> float:
        mae = float(np.mean(np.abs(a - b)))
        return float(np.clip(1.0 - mae, 0.0, 1.0))

    @staticmethod
    def _next_actions(registry_df: pd.DataFrame, figures_df: pd.DataFrame) -> list[str]:
        actions: list[str] = []
        if not registry_df.empty:
            missing_source = registry_df.loc[~registry_df["source_complete"], "paper_id"].tolist()
            if missing_source:
                actions.append(
                    "Fill missing source metadata (repo_url/repo_commit/license) for papers: "
                    + ", ".join(str(x) for x in missing_source)
                )
            missing_pdf = registry_df.loc[~registry_df["paper_exists"], "paper_id"].tolist()
            if missing_pdf:
                actions.append(
                    "Add missing local paper files for papers: "
                    + ", ".join(str(x) for x in missing_pdf)
                )
        if not figures_df.empty:
            failed = figures_df.loc[figures_df["status"] != "pass", ["paper_id", "figure_id", "status"]]
            if not failed.empty:
                tags = [f"{r.paper_id}:{r.figure_id}({r.status})" for r in failed.itertuples()]
                actions.append("Investigate failed figure checks: " + ", ".join(tags))
        if not actions:
            actions.append("All configured paper reproduction checks passed.")
        return actions
