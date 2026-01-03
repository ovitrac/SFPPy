#!/bin/bash
# Kill SFPPy Survey Server processes

echo "Killing Survey Server processes..."

# Kill uvicorn/simulator processes
pkill -f "uvicorn.*simulator:app" 2>/dev/null && echo "  Killed uvicorn simulator" || echo "  No uvicorn simulator found"
pkill -f "python.*simulator.py" 2>/dev/null && echo "  Killed simulator.py" || echo "  No simulator.py found"
pkill -f "python.*launcher.py" 2>/dev/null && echo "  Killed launcher.py" || echo "  No launcher.py found"

# Kill any process on port 8001 (simulator default port)
fuser -k 8001/tcp 2>/dev/null && echo "  Freed port 8001" || echo "  Port 8001 was free"

# Kill any process on port 8000 (family editor default port)
fuser -k 8000/tcp 2>/dev/null && echo "  Freed port 8000" || echo "  Port 8000 was free"

echo "Done."
