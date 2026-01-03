#!/usr/bin/env python3
"""
sync_base_env.py

Compare conda base environment package sets between local machine and a remote machine over SSH.

- Uses: `conda list --json` (for conda + pip packages visible to conda)
- Optionally checks `pip list --format=json` to catch pip-only installs.

Outputs:
- Missing on remote
- Missing on local
- Version mismatches
- Optional suggested install commands for remote

Requirements:
- local: conda available in PATH (or provide --conda)
- remote: conda available in PATH for the SSH session (or provide --remote-conda)
- ssh available locally and key-based auth (or you'll be prompted by ssh)

Example:
  python sync_base_env.py --remote mai@192.168.1.102
  python sync_base_env.py --remote mai@192.168.1.102 --include-pip
  python sync_base_env.py --remote mai@192.168.1.102 --include-pip --suggest-install
  python sync_base_env.py --remote mai@192.168.1.102 --check-project-path --remote-project-path '~/natacha/python'


"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Pkg:
    name: str
    version: str
    channel: str = ""
    build: str = ""
    manager: str = "conda"  # or "pip"


def run(cmd: List[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True)


def run_local(command: str) -> subprocess.CompletedProcess:
    return run(["bash", "-lc", command])


def run_remote(remote: str, command: str) -> subprocess.CompletedProcess:
    return run(["ssh", remote, "bash", "-lc", command])


def parse_json(text: str, context: str) -> object:
    s = (text or "").strip()
    if not s:
        raise RuntimeError(f"[ERROR] Empty output while expecting JSON ({context}).")
    try:
        return json.loads(s)
    except json.JSONDecodeError as e:
        snippet = "\n".join(s.splitlines()[:80])
        raise RuntimeError(
            f"[ERROR] Invalid JSON ({context}): {e}\n"
            f"--- stdout snippet ---\n{snippet}\n--- end ---"
        ) from e


def load_conda_list(output: str) -> Dict[str, Pkg]:
    data = parse_json(output, "conda list --json")
    out: Dict[str, Pkg] = {}
    for item in data:
        name = (item.get("name") or "").strip()
        if not name:
            continue
        out[name.lower()] = Pkg(
            name=name,
            version=(item.get("version") or "").strip(),
            channel=item.get("channel") or "",
            build=item.get("build_string") or item.get("build") or "",
            manager="conda",
        )
    return out


def load_pip_list(output: str) -> Dict[str, Pkg]:
    data = parse_json(output, "pip list --format=json")
    out: Dict[str, Pkg] = {}
    for item in data:
        name = (item.get("name") or "").strip()
        if not name:
            continue
        out[name.lower()] = Pkg(name=name, version=(item.get("version") or "").strip(), manager="pip")
    return out


def merge_pkgs(conda_pkgs: Dict[str, Pkg], pip_pkgs: Dict[str, Pkg]) -> Dict[str, Pkg]:
    merged = dict(pip_pkgs)
    # conda wins on collision
    for k, v in conda_pkgs.items():
        merged[k] = v
    return merged


def inventory_local(include_pip: bool, conda_exe: str = "conda") -> Dict[str, Pkg]:
    cp = run_local(f"{shlex.quote(conda_exe)} list -n base --json")
    if cp.returncode != 0:
        raise RuntimeError(f"[ERROR] local conda failed:\nSTDERR:\n{cp.stderr}\nSTDOUT:\n{cp.stdout}\n")
    conda_pkgs = load_conda_list(cp.stdout)

    if not include_pip:
        return conda_pkgs

    cp2 = run_local(f"{shlex.quote(conda_exe)} run -n base python -m pip list --format=json")
    if cp2.returncode != 0:
        raise RuntimeError(f"[ERROR] local pip failed:\nSTDERR:\n{cp2.stderr}\nSTDOUT:\n{cp2.stdout}\n")
    pip_pkgs = load_pip_list(cp2.stdout)
    return merge_pkgs(conda_pkgs, pip_pkgs)


def inventory_remote(remote: str, include_pip: bool) -> Dict[str, Pkg]:
    # Force your remote initialization path:
    # - source ~/.bashrc to define condainit
    # - run condainit to set PATH / conda shell hook
    # - run conda with --no-plugins for quieter behavior
    prefix = "source ~/.bashrc >/dev/null 2>&1; condainit >/dev/null 2>&1 || true; "

    cmd_list = prefix + "conda --no-plugins list -n base --json"
    cp = run_remote(remote, cmd_list)
    if cp.returncode != 0:
        raise RuntimeError(
            f"[ERROR] remote conda list failed (rc={cp.returncode}).\n"
            f"CMD:\n{cmd_list}\n\nSTDERR:\n{cp.stderr}\n\nSTDOUT:\n{cp.stdout}\n"
        )
    conda_pkgs = load_conda_list(cp.stdout)

    if not include_pip:
        return conda_pkgs

    cmd_pip = prefix + "conda --no-plugins run -n base python -m pip list --format=json"
    cp2 = run_remote(remote, cmd_pip)
    if cp2.returncode != 0:
        raise RuntimeError(
            f"[ERROR] remote pip list failed (rc={cp2.returncode}).\n"
            f"CMD:\n{cmd_pip}\n\nSTDERR:\n{cp2.stderr}\n\nSTDOUT:\n{cp2.stdout}\n"
        )
    pip_pkgs = load_pip_list(cp2.stdout)
    return merge_pkgs(conda_pkgs, pip_pkgs)


def diff_pkgs(local: Dict[str, Pkg], remote: Dict[str, Pkg]) -> Tuple[List[Pkg], List[Pkg], List[Tuple[Pkg, Pkg]]]:
    missing_on_remote: List[Pkg] = []
    missing_on_local: List[Pkg] = []
    mismatches: List[Tuple[Pkg, Pkg]] = []

    for k, lp in sorted(local.items(), key=lambda kv: kv[0]):
        rp = remote.get(k)
        if rp is None:
            missing_on_remote.append(lp)
        elif lp.version != rp.version:
            mismatches.append((lp, rp))

    for k, rp in sorted(remote.items(), key=lambda kv: kv[0]):
        if k not in local:
            missing_on_local.append(rp)

    return missing_on_remote, missing_on_local, mismatches


def suggest_remote_commands(missing_on_remote: List[Pkg]) -> Tuple[str, str]:
    conda_specs = [f"{p.name}=={p.version}" for p in missing_on_remote if p.manager == "conda"]
    pip_specs = [f"{p.name}=={p.version}" for p in missing_on_remote if p.manager == "pip"]

    conda_cmd = "conda install -n base " + " ".join(map(shlex.quote, conda_specs)) if conda_specs else "(none)"
    pip_cmd = "python -m pip install " + " ".join(map(shlex.quote, pip_specs)) if pip_specs else "(none)"
    return conda_cmd, pip_cmd


def fmt(pkgs: List[Pkg]) -> str:
    if not pkgs:
        return "(none)"
    return "\n".join(f"{p.name}=={p.version} [{p.manager}]" for p in pkgs)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--remote", default="mai@192.168.1.102")
    ap.add_argument("--include-pip", action="store_true")
    ap.add_argument("--suggest-install", action="store_true")
    ap.add_argument("--conda", default="conda", help="Local conda executable (default: conda)")
    args = ap.parse_args()

    print(f"[INFO] Remote: {args.remote}")
    print(f"[INFO] Comparing conda base env (include_pip={args.include_pip})")

    local_pkgs = inventory_local(args.include_pip, conda_exe=args.conda)
    remote_pkgs = inventory_remote(args.remote, args.include_pip)

    missing_on_remote, missing_on_local, mismatches = diff_pkgs(local_pkgs, remote_pkgs)

    print("\n=== Missing on REMOTE (present locally) ===")
    print(fmt(missing_on_remote))

    print("\n=== Missing on LOCAL (present remotely) ===")
    print(fmt(missing_on_local))

    print("\n=== Version mismatches (local vs remote) ===")
    if mismatches:
        for lp, rp in mismatches:
            print(f"{lp.name}: local={lp.version} ({lp.manager})  remote={rp.version} ({rp.manager})")
    else:
        print("(none)")

    if args.suggest_install:
        conda_cmd, pip_cmd = suggest_remote_commands(missing_on_remote)
        print("\n=== Suggested commands to run on REMOTE to match LOCAL (best effort) ===")
        print("# Conda:")
        print(conda_cmd)
        print("\n# Pip (run inside base on remote):")
        print(pip_cmd)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())


