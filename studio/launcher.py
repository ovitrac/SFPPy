#!/usr/bin/env python3
"""
SFPPy Studio Launcher

Launches the comprehensive migration analysis web application.

Usage:
    python launcher.py [--port PORT] [--host HOST] [--reload]

Author: Olivier Vitrac, PhD, HDR
Organization: INRAE + Generative Simulation
"""

import argparse
import subprocess
import sys
import signal
import socket
import time
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from studio.version import __version__, get_version_info


def is_port_in_use(port: int) -> bool:
    """Check if a port is already in use."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('127.0.0.1', port)) == 0


def check_dependencies():
    """Check required dependencies are installed."""
    required = [
        ('numpy', 'numpy'),
        ('scipy', 'scipy'),
        ('yaml', 'pyyaml'),
        ('fastapi', 'fastapi'),
        ('uvicorn', 'uvicorn'),
        ('pydantic', 'pydantic'),
        ('jinja2', 'jinja2'),
    ]

    missing = []
    for module, package in required:
        try:
            __import__(module)
            print(f"  [OK] {module}")
        except ImportError:
            print(f"  [MISSING] {module} - install with: pip install {package}")
            missing.append(package)

    return len(missing) == 0


def main():
    parser = argparse.ArgumentParser(
        description="SFPPy Studio - Food Contact Migration Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python launcher.py                    # Start on default port 8002
    python launcher.py --port 8080        # Start on custom port
    python launcher.py --reload           # Enable auto-reload for development
        """
    )
    parser.add_argument('--port', type=int, default=8002,
                        help='Port to run the server on (default: 8002)')
    parser.add_argument('--host', type=str, default='127.0.0.1',
                        help='Host to bind to (default: 127.0.0.1)')
    parser.add_argument('--reload', action='store_true',
                        help='Enable auto-reload for development')

    args = parser.parse_args()

    # Banner
    print("=" * 60)
    print("SFPPy Studio - Food Contact Migration Analysis")
    print("=" * 60)
    print()

    # Version info
    info = get_version_info()
    print(f"  Studio Version: {info['studio_version']}")
    print(f"  SFPPy Version:  {info['sfppy_version']}")
    print(f"  Python:         {info['python_version']}")
    print(f"  Machine:        {info['machine_name']} ({info['machine_os']})")
    print()

    # Check dependencies
    print("Checking dependencies...")
    if not check_dependencies():
        print("\nError: Missing dependencies. Install them and try again.")
        return 1

    print()

    # Check port availability
    if is_port_in_use(args.port):
        print(f"Warning: Port {args.port} is already in use.")
        print(f"  Either stop the existing service or use --port to specify another port.")
        print(f"  Example: python launcher.py --port 8003")
        return 1

    # Start the server
    print(f"Starting SFPPy Studio on http://{args.host}:{args.port}")
    print()

    studio_dir = Path(__file__).parent
    cmd = [
        sys.executable, "-m", "uvicorn",
        "studio.app.main:app",
        "--host", args.host,
        "--port", str(args.port),
    ]

    if args.reload:
        cmd.append("--reload")
        print("  Auto-reload enabled (development mode)")

    print("=" * 60)
    print(f"  SFPPy Studio: http://{args.host}:{args.port}")
    print()
    print("  Press Ctrl+C to stop")
    print("=" * 60)
    print()

    # Handle shutdown
    process = None

    def signal_handler(sig, frame):
        print("\n\nShutting down SFPPy Studio...")
        if process:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
        print("Done.")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        process = subprocess.Popen(
            cmd,
            cwd=str(studio_dir.parent),  # Run from python/ directory
        )
        process.wait()
    except KeyboardInterrupt:
        signal_handler(None, None)

    return 0


if __name__ == "__main__":
    sys.exit(main())
