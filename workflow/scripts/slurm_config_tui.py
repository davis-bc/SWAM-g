#!/usr/bin/env python3
"""Textual editor for the OHM-g Slurm profile config."""

from __future__ import annotations

import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

if len(sys.argv) > 1 and sys.argv[1] in {"-h", "--help"}:
    print(
        "Usage: python workflow/scripts/slurm_config_tui.py [path/to/config.yaml]\n\n"
        "Defaults to config/slurm/config.yaml."
    )
    sys.exit(0)

try:
    from ruamel.yaml import YAML
    from textual.app import App, ComposeResult
    from textual.containers import Container, Horizontal, ScrollableContainer, Vertical
    from textual.widgets import Button, Header, Input, Label, Select, Static, Switch
except ModuleNotFoundError as exc:  # pragma: no cover - exercised by end users.
    missing = {"ruamel": "ruamel.yaml"}.get(exc.name or "", exc.name or "a required package")
    print(
        f"Missing dependency: {missing}\n\n"
        "Install the OHM-g setup dependencies, then run this command again:\n"
        "  conda install textual ruamel.yaml\n\n"
        "Or install only the TUI dependencies with pip:\n"
        "  python -m pip install -r requirements-tui.txt\n",
        file=sys.stderr,
    )
    sys.exit(1)


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = REPO_ROOT / "config" / "slurm" / "config.yaml"

MODULE_FIELDS = [
    ("run_coverage", "Coverage"),
    ("run_checkm2", "CheckM2 quality"),
    ("run_mlst", "MLST"),
    ("run_txsscan", "TXSScan"),
    ("run_mef", "MobileElementFinder"),
    ("run_mashtree", "Mashtree"),
    ("run_resfinder", "ResFinder"),
    ("run_salmonella_serotyping", "Salmonella serotyping"),
    ("run_ecoli_pathotyping", "E. coli pathotyping"),
    ("pd_lookup", "Pathogen Detection lookup"),
]

PROFILE_FIELDS = [
    ("jobs", "Max submitted jobs"),
    ("cores", "Total available cores"),
    ("restart-times", "Restart attempts"),
    ("latency-wait", "Latency wait seconds"),
]

DEFAULT_RESOURCE_FIELDS = [
    ("slurm_account", "Slurm account"),
    ("slurm_partition", "Slurm partition"),
    ("mem_mb", "Default memory MB"),
    ("runtime", "Default walltime"),
]

BOOLEAN_TOP_LEVEL_FIELDS = [
    ("use-conda", "Use conda"),
    ("keep-incomplete", "Keep incomplete files"),
    ("keep-going", "Keep going after failures"),
    ("rerun-incomplete", "Rerun incomplete jobs"),
    ("show-failed-logs", "Show failed logs"),
]


def load_yaml(path: Path) -> tuple[YAML, Any]:
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.load(handle)
    return yaml, data


def as_text(value: Any) -> str:
    return "" if value is None else str(value)


class FieldRow(Horizontal):
    """A label and one form control."""

    def __init__(self, label: str, control: Any, help_text: str = "") -> None:
        super().__init__(classes="field-row")
        self.label_text = label
        self.help_text = help_text
        self.control = control

    def compose(self) -> ComposeResult:
        yield Label(self.label_text, classes="field-label")
        yield self.control
        if self.help_text:
            yield Static(self.help_text, classes="field-help")


class Section(Vertical):
    """A titled form section."""

    def __init__(self, title: str, *children: Any) -> None:
        super().__init__(classes="section")
        self.title = title
        self.section_children = children

    def compose(self) -> ComposeResult:
        yield Static(self.title, classes="section-title")
        yield from self.section_children


class SlurmConfigTui(App[None]):
    """Lightweight editor for config/slurm/config.yaml."""

    TITLE = "OHM-g Slurm Config"

    CSS = """
    Screen {
        background: $surface;
    }

    #body {
        height: 1fr;
        width: 100%;
        overflow-y: auto;
        padding: 1 2;
    }

    #toolbar {
        dock: bottom;
        height: 5;
        width: 100%;
        padding: 0 2;
        background: $panel;
        border-top: solid $primary;
        align: right middle;
    }

    #status {
        width: 1fr;
        content-align: left middle;
    }

    Button {
        margin: 0 1;
        height: 3;
        min-width: 18;
        color: white;
        text-style: bold;
        content-align: center middle;
    }

    .section {
        height: auto;
        width: 100%;
        margin-bottom: 1;
        padding: 1 2;
        border: solid $primary 30%;
    }

    .section-title {
        text-style: bold;
        color: $accent;
        margin-bottom: 1;
    }

    .field-row {
        height: 3;
        width: 100%;
        margin-bottom: 1;
    }

    .field-label {
        width: 30;
        color: $text;
        content-align: left middle;
    }

    .field-help {
        width: 1fr;
        color: $text-muted;
        content-align: left middle;
        padding-left: 1;
    }

    Input {
        width: 32;
        height: 3;
        color: $text;
    }

    Select {
        width: 32;
        height: 3;
        color: $text;
    }

    Switch {
        margin-top: 0;
    }

    .toggle-grid {
        layout: grid;
        grid-size: 2;
        grid-columns: 1fr 1fr;
        grid-gutter: 1 2;
    }

    .toggle-row {
        height: 3;
    }

    .rule-header {
        height: 1;
        color: $text-muted;
        text-style: bold;
        margin-bottom: 1;
    }

    .rule-row {
        height: 3;
        width: 100%;
        margin-bottom: 0;
    }

    .rule-name {
        width: 24;
        content-align: left middle;
    }

    .rule-input {
        width: 16;
        margin-right: 1;
    }
    """

    BINDINGS = [
        ("ctrl+s", "save", "Save"),
        ("ctrl+r", "reload", "Reload"),
        ("q", "quit", "Quit"),
    ]

    def __init__(self, config_path: Path) -> None:
        super().__init__()
        self.config_path = config_path
        self.yaml, self.data = load_yaml(config_path)
        self.status: Static | None = None

    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        with ScrollableContainer(id="body"):
            yield Section("Slurm Scheduler", *self.profile_fields())
            yield Section("Default Resources", *self.default_resource_fields())
            yield Section("Analysis Modules", self.module_toggles())
            yield Section("Pathogen Detection Lookup", *self.pd_fields())
            yield Section("Per-Rule Resources", self.rule_resource_rows())
        with Horizontal(id="toolbar"):
            self.status = Static(f"Editing {self.config_path}", id="status")
            yield self.status
            yield Button("Reload Config", id="reload", variant="primary")
            yield Button("Save Config", id="save", variant="success")
            yield Button("Quit", id="quit", variant="error")

    def profile_fields(self) -> list[FieldRow]:
        rows: list[FieldRow] = [
            FieldRow("Executor", Input(value=as_text(self.data.get("executor")), id="top__executor")),
            FieldRow(
                "Conda frontend",
                Select(
                    [("conda", "conda"), ("mamba", "mamba")],
                    value=self.data.get("conda-frontend", "conda"),
                    id="top__conda-frontend",
                    allow_blank=False,
                ),
            ),
        ]
        for key, label in PROFILE_FIELDS:
            rows.append(FieldRow(label, Input(value=as_text(self.data.get(key)), id=f"top__{key}")))
        rows.extend(
            FieldRow(label, Switch(value=bool(self.data.get(key)), id=f"top__{key}"))
            for key, label in BOOLEAN_TOP_LEVEL_FIELDS
        )
        return rows

    def default_resource_fields(self) -> list[FieldRow]:
        resources = self.data.setdefault("default-resources", {})
        return [
            FieldRow(label, Input(value=as_text(resources.get(key)), id=f"default__{key}"))
            for key, label in DEFAULT_RESOURCE_FIELDS
        ]

    def module_toggles(self) -> Container:
        config = self.data.setdefault("config", {})
        rows = [
            FieldRow(label, Switch(value=bool(config.get(key)), id=f"config__{key}"))
            for key, label in MODULE_FIELDS
        ]
        return Container(*rows, classes="toggle-grid")

    def pd_fields(self) -> list[FieldRow]:
        config = self.data.setdefault("config", {})
        return [
            FieldRow(
                "PD backend",
                Select(
                    [("ftp", "ftp"), ("local files", "local")],
                    value=config.get("pd_backend", "ftp"),
                    id="config__pd_backend",
                    allow_blank=False,
                ),
                "Use local only when TSV paths below are populated.",
            ),
            FieldRow(
                "Comparator limit",
                Input(value=as_text(config.get("pd_comparator_limit", 10)), id="config__pd_comparator_limit"),
            ),
            FieldRow(
                "Sample metadata TSV",
                Input(value=as_text(config.get("pd_sample_metadata_tsv")), id="config__pd_sample_metadata_tsv"),
            ),
            FieldRow(
                "PD isolates TSV",
                Input(value=as_text(config.get("pd_isolates_tsv")), id="config__pd_isolates_tsv"),
            ),
            FieldRow(
                "PD exceptions TSV",
                Input(value=as_text(config.get("pd_exceptions_tsv")), id="config__pd_exceptions_tsv"),
            ),
        ]

    def rule_resource_rows(self) -> Vertical:
        rule_names = sorted(
            {
                *self.data.get("set-threads", {}).keys(),
                *self.data.get("set-resources", {}).keys(),
                *self.data.get("group-components", {}).keys(),
            }
        )
        rows: list[Any] = [
            Horizontal(
                Static("Rule", classes="rule-name"),
                Static("Threads", classes="rule-input"),
                Static("Memory MB", classes="rule-input"),
                Static("Runtime", classes="rule-input"),
                Static("Group size", classes="rule-input"),
                classes="rule-header",
            )
        ]
        for rule in rule_names:
            rule_resources = self.data.get("set-resources", {}).get(rule, {})
            rows.append(
                Horizontal(
                    Static(rule, classes="rule-name"),
                    Input(value=as_text(self.data.get("set-threads", {}).get(rule, "")), id=f"threads__{rule}", classes="rule-input"),
                    Input(value=as_text(rule_resources.get("mem_mb", "")), id=f"mem__{rule}", classes="rule-input"),
                    Input(value=as_text(rule_resources.get("runtime", "")), id=f"runtime__{rule}", classes="rule-input"),
                    Input(value=as_text(self.data.get("group-components", {}).get(rule, "")), id=f"group_components__{rule}", classes="rule-input"),
                    classes="rule-row",
                )
            )
        return Vertical(*rows)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "save":
            self.action_save()
        elif event.button.id == "reload":
            self.action_reload()
        elif event.button.id == "quit":
            self.exit()

    def action_reload(self) -> None:
        self.yaml, self.data = load_yaml(self.config_path)
        self.refresh(recompose=True)
        self.set_status("Reloaded from disk.")

    def action_save(self) -> None:
        try:
            self.apply_form_to_data()
        except ValueError as exc:
            self.set_status(str(exc), error=True)
            return

        backup_path = self.backup_path()
        shutil.copy2(self.config_path, backup_path)
        with self.config_path.open("w", encoding="utf-8") as handle:
            self.yaml.dump(self.data, handle)
        self.set_status(f"Saved {self.config_path} (backup: {backup_path.name}).")

    def apply_form_to_data(self) -> None:
        for key in ["executor"]:
            self.data[key] = self.input_value(f"top__{key}")
        self.data["conda-frontend"] = self.select_value("top__conda-frontend")
        for key, _label in PROFILE_FIELDS:
            self.data[key] = self.int_value(f"top__{key}", key)
        for key, _label in BOOLEAN_TOP_LEVEL_FIELDS:
            self.data[key] = self.switch_value(f"top__{key}")

        default_resources = self.data.setdefault("default-resources", {})
        for key, _label in DEFAULT_RESOURCE_FIELDS:
            raw_value = self.input_value(f"default__{key}")
            default_resources[key] = self.parse_int(raw_value, "default mem_mb") if key == "mem_mb" else raw_value

        config = self.data.setdefault("config", {})
        for key, _label in MODULE_FIELDS:
            config[key] = self.switch_value(f"config__{key}")
        config["pd_backend"] = self.select_value("config__pd_backend")
        config["pd_comparator_limit"] = self.int_value("config__pd_comparator_limit", "pd_comparator_limit")
        for key in ["pd_sample_metadata_tsv", "pd_isolates_tsv", "pd_exceptions_tsv"]:
            config[key] = self.input_value(f"config__{key}")

        rule_names = sorted(
            {
                *self.data.get("set-threads", {}).keys(),
                *self.data.get("set-resources", {}).keys(),
                *self.data.get("group-components", {}).keys(),
            }
        )
        for rule in rule_names:
            self.data.setdefault("set-threads", {})[rule] = self.int_value(f"threads__{rule}", f"{rule} threads")
            mem_value = self.input_value(f"mem__{rule}")
            runtime_value = self.input_value(f"runtime__{rule}")
            set_resources = self.data.setdefault("set-resources", {})
            if mem_value or runtime_value:
                rule_resources = set_resources.setdefault(rule, {})
                if mem_value:
                    rule_resources["mem_mb"] = self.parse_int(mem_value, f"{rule} mem_mb")
                else:
                    rule_resources.pop("mem_mb", None)
                if runtime_value:
                    rule_resources["runtime"] = runtime_value
                else:
                    rule_resources.pop("runtime", None)
            elif rule in set_resources:
                del set_resources[rule]

            group_value = self.input_value(f"group_components__{rule}")
            if group_value:
                self.data.setdefault("group-components", {})[rule] = self.parse_int(group_value, f"{rule} group size")
            else:
                self.data.setdefault("group-components", {}).pop(rule, None)

    def input_value(self, widget_id: str) -> str:
        return self.query_one(f"#{widget_id}", Input).value.strip()

    def switch_value(self, widget_id: str) -> bool:
        return bool(self.query_one(f"#{widget_id}", Switch).value)

    def select_value(self, widget_id: str) -> str:
        return str(self.query_one(f"#{widget_id}", Select).value)

    def int_value(self, widget_id: str, label: str) -> int:
        return self.parse_int(self.input_value(widget_id), label)

    @staticmethod
    def parse_int(value: str, label: str) -> int:
        try:
            parsed = int(value)
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer.") from exc
        if parsed < 0:
            raise ValueError(f"{label} must be zero or greater.")
        return parsed

    def backup_path(self) -> Path:
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        return self.config_path.with_name(f"{self.config_path.name}.{stamp}.bak")

    def set_status(self, message: str, error: bool = False) -> None:
        if self.status is None:
            return
        self.status.update(message)
        self.status.styles.color = "red" if error else "green"


def main() -> int:
    config_path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_CONFIG_PATH
    if not config_path.exists():
        print(f"Config file not found: {config_path}", file=sys.stderr)
        return 1
    SlurmConfigTui(config_path).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
