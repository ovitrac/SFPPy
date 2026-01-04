#!/bin/bash
# =============================================================================
# SFPPy Full - Entrypoint Script
# =============================================================================
# Launches all SFPPy services (Studio + Survey Editor + Survey Simulator)
#
# Author: Olivier Vitrac, PhD, HDR
# Organization: INRAE + Generative Simulation
# =============================================================================

set -e

echo "=============================================="
echo "SFPPy Full - Starting All Services"
echo "=============================================="

# Start Survey Editor (port 8000)
echo "Starting Survey Editor on port 8000..."
python -m uvicorn survey.app.main:app --host 0.0.0.0 --port 8000 &
PID_EDITOR=$!

# Start Survey Simulator (port 8001)
echo "Starting Survey Simulator on port 8001..."
python -m uvicorn survey.app.simulator:app --host 0.0.0.0 --port 8001 &
PID_SIMULATOR=$!

# Start Studio (port 8002)
echo "Starting Studio on port 8002..."
python -m uvicorn studio.app.main:app --host 0.0.0.0 --port 8002 &
PID_STUDIO=$!

echo ""
echo "=============================================="
echo "All services started:"
echo "  Studio:           http://localhost:8002"
echo "  Survey Editor:    http://localhost:8000"
echo "  Survey Simulator: http://localhost:8001"
echo "=============================================="
echo ""

# Handle shutdown
cleanup() {
    echo "Shutting down services..."
    kill $PID_EDITOR $PID_SIMULATOR $PID_STUDIO 2>/dev/null || true
    wait
    echo "All services stopped."
}

trap cleanup SIGTERM SIGINT

# Wait for any process to exit
wait -n

# If any process exits, stop all
cleanup
