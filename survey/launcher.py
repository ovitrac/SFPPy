#!/usr/bin/env python3
"""
SFPPy Survey Applications Launcher
===================================

Launch the Survey web applications (Family Editor and Survey Simulator).

Usage:
    python survey/launcher.py              # Launch both apps
    python survey/launcher.py --app editor # Launch only Family Editor (port 8000)
    python survey/launcher.py --app simulator # Launch only Survey Simulator (port 8001)
    python survey/launcher.py --check      # Check dependencies only

Author: Olivier Vitrac
"""

import os
import sys
import argparse
import subprocess
import time
import signal
import webbrowser
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))


def check_dependencies(verbose: bool = True) -> dict:
    """Check all dependencies and return status."""
    results = {
        "core": {},
        "web": {},
        "optional": {},
        "all_ok": True,
    }

    # Core dependencies
    core_deps = [
        ("numpy", "numpy"),
        ("scipy", "scipy"),
        ("yaml", "pyyaml"),
    ]

    for name, pip_name in core_deps:
        try:
            __import__(name)
            results["core"][name] = True
            if verbose:
                print(f"  [OK] {name}")
        except ImportError:
            results["core"][name] = False
            results["all_ok"] = False
            if verbose:
                print(f"  [MISSING] {name} - install with: pip install {pip_name}")

    # Web dependencies
    web_deps = [
        ("fastapi", "fastapi"),
        ("uvicorn", "uvicorn"),
        ("pydantic", "pydantic"),
        ("jinja2", "jinja2"),
    ]

    for name, pip_name in web_deps:
        try:
            __import__(name)
            results["web"][name] = True
            if verbose:
                print(f"  [OK] {name}")
        except ImportError:
            results["web"][name] = False
            results["all_ok"] = False
            if verbose:
                print(f"  [MISSING] {name} - install with: pip install {pip_name}")

    # Optional dependencies
    optional_deps = [
        ("openpyxl", "openpyxl", "XLSX spreadsheet support"),
        ("odf", "odfpy", "ODS spreadsheet support"),
        ("requests", "requests", "PubChem integration"),
    ]

    for name, pip_name, desc in optional_deps:
        try:
            __import__(name)
            results["optional"][name] = True
            if verbose:
                print(f"  [OK] {name} ({desc})")
        except ImportError:
            results["optional"][name] = False
            if verbose:
                print(f"  [OPTIONAL] {name} ({desc}) - install with: pip install {pip_name}")

    return results


def check_port_available(port: int) -> bool:
    """Check if a port is available."""
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.bind(("127.0.0.1", port))
            return True
        except OSError:
            return False


def wait_for_server(port: int, timeout: float = 10.0) -> bool:
    """Wait for server to become available."""
    import socket
    start = time.time()
    while time.time() - start < timeout:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.connect(("127.0.0.1", port))
                return True
            except (ConnectionRefusedError, OSError):
                time.sleep(0.2)
    return False


def launch_app(app: str, port: int, reload: bool = True) -> subprocess.Popen:
    """Launch a single application."""
    survey_dir = Path(__file__).parent

    if app == "editor":
        module = "survey.app.main:app"
    elif app == "simulator":
        module = "survey.app.simulator:app"
    else:
        raise ValueError(f"Unknown app: {app}")

    cmd = [
        sys.executable, "-m", "uvicorn",
        module,
        "--host", "127.0.0.1",
        "--port", str(port),
    ]

    if reload:
        cmd.append("--reload")

    # Launch from parent directory (so imports work)
    cwd = survey_dir.parent

    proc = subprocess.Popen(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    return proc


def main():
    parser = argparse.ArgumentParser(
        description="SFPPy Survey Applications Launcher",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python survey/launcher.py              # Launch both apps
    python survey/launcher.py --app editor # Family Editor only (port 8000)
    python survey/launcher.py --app simulator # Survey Simulator only (port 8001)
    python survey/launcher.py --check      # Check dependencies
    python survey/launcher.py --no-browser # Don't open browser
        """
    )

    parser.add_argument(
        "--app",
        choices=["editor", "simulator", "both"],
        default="both",
        help="Which app to launch (default: both)"
    )

    parser.add_argument(
        "--check",
        action="store_true",
        help="Check dependencies and exit"
    )

    parser.add_argument(
        "--no-reload",
        action="store_true",
        help="Disable auto-reload on file changes"
    )

    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Don't open browser automatically"
    )

    parser.add_argument(
        "--editor-port",
        type=int,
        default=8000,
        help="Port for Family Editor (default: 8000)"
    )

    parser.add_argument(
        "--simulator-port",
        type=int,
        default=8001,
        help="Port for Survey Simulator (default: 8001)"
    )

    args = parser.parse_args()

    print("=" * 60)
    print("SFPPy Survey Applications Launcher")
    print("=" * 60)

    # Check dependencies
    print("\nChecking dependencies...")
    deps = check_dependencies(verbose=True)

    if args.check:
        print("\n" + "=" * 60)
        if deps["all_ok"]:
            print("All required dependencies are installed.")
        else:
            print("Some required dependencies are missing.")
            print("Install missing packages with pip.")
        return 0 if deps["all_ok"] else 1

    if not deps["all_ok"]:
        print("\nCannot launch: missing required dependencies.")
        print("Run: pip install fastapi uvicorn pydantic jinja2 pyyaml numpy scipy")
        return 1

    # Check ports
    apps_to_launch = []
    if args.app in ["editor", "both"]:
        if not check_port_available(args.editor_port):
            print(f"\nWarning: Port {args.editor_port} is in use (Family Editor may already be running)")
        else:
            apps_to_launch.append(("editor", args.editor_port))

    if args.app in ["simulator", "both"]:
        if not check_port_available(args.simulator_port):
            print(f"\nWarning: Port {args.simulator_port} is in use (Survey Simulator may already be running)")
        else:
            apps_to_launch.append(("simulator", args.simulator_port))

    if not apps_to_launch:
        print("\nNo apps to launch (ports may already be in use)")
        print(f"  Family Editor: http://127.0.0.1:{args.editor_port}")
        print(f"  Survey Simulator: http://127.0.0.1:{args.simulator_port}")
        return 0

    # Launch applications
    print("\nLaunching applications...")
    processes = []
    reload = not args.no_reload

    for app_name, port in apps_to_launch:
        print(f"  Starting {app_name} on port {port}...")
        proc = launch_app(app_name, port, reload=reload)
        processes.append((app_name, port, proc))

    # Wait for servers to start
    print("\nWaiting for servers to start...")
    urls = []
    for app_name, port, proc in processes:
        if wait_for_server(port, timeout=15.0):
            url = f"http://127.0.0.1:{port}"
            urls.append(url)
            print(f"  [OK] {app_name}: {url}")
        else:
            print(f"  [TIMEOUT] {app_name} on port {port}")

    # Open browser
    if urls and not args.no_browser:
        print("\nOpening browser...")
        for url in urls:
            webbrowser.open(url)

    # Print summary
    print("\n" + "=" * 60)
    print("Applications running:")
    print("=" * 60)
    if args.app in ["editor", "both"]:
        print(f"  Family Editor:    http://127.0.0.1:{args.editor_port}")
    if args.app in ["simulator", "both"]:
        print(f"  Survey Simulator:  http://127.0.0.1:{args.simulator_port}")
    print("\nPress Ctrl+C to stop all applications")
    print("=" * 60)

    # Handle shutdown
    def signal_handler(sig, frame):
        print("\n\nShutting down...")
        for app_name, port, proc in processes:
            print(f"  Stopping {app_name}...")
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
        print("Done.")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Keep running and show output
    try:
        while True:
            for app_name, port, proc in processes:
                if proc.poll() is not None:
                    print(f"\n{app_name} exited with code {proc.returncode}")
                    # Read any remaining output
                    output, _ = proc.communicate()
                    if output:
                        print(output)
            time.sleep(1)
    except KeyboardInterrupt:
        signal_handler(None, None)


if __name__ == "__main__":
    sys.exit(main())
