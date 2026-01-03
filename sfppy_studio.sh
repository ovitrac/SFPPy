#!/bin/bash
# =============================================================================
# SFPPy Studio Launcher
# Food Contact Migration Analysis Web Application
#
# Author: Olivier Vitrac, PhD, HDR
# Organization: INRAE + Generative Simulation
# =============================================================================

set -e

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STUDIO_DIR="${SCRIPT_DIR}/studio"
DEFAULT_PORT=8002
DEFAULT_HOST="127.0.0.1"

# Print banner
print_banner() {
    echo -e "${GREEN}"
    echo "  ╔═══════════════════════════════════════════════════════════════╗"
    echo "  ║                                                               ║"
    echo "  ║   🍏⏩🍎  SFPPy Studio - Food Contact Migration Analysis      ║"
    echo "  ║                                                               ║"
    echo "  ╚═══════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

# Print usage
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -p, --port PORT      Server port (default: ${DEFAULT_PORT})"
    echo "  -H, --host HOST      Server host (default: ${DEFAULT_HOST})"
    echo "  -d, --debug          Enable debug mode"
    echo "  -r, --reload         Enable auto-reload for development"
    echo "  --check              Check dependencies only"
    echo "  --kill               Kill any running SFPPy Studio server"
    echo "  -h, --help           Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0                   # Start on default port 8002"
    echo "  $0 -p 8080           # Start on port 8080"
    echo "  $0 -d -r             # Start in debug mode with auto-reload"
    echo "  $0 --kill            # Stop running server"
}

# Check Python dependencies
check_dependencies() {
    echo -e "${BLUE}Checking dependencies...${NC}"

    local missing=0

    # Map package names to import names (some differ)
    declare -A pkg_imports=(
        ["numpy"]="numpy"
        ["scipy"]="scipy"
        ["fastapi"]="fastapi"
        ["uvicorn"]="uvicorn"
        ["pydantic"]="pydantic"
        ["jinja2"]="jinja2"
        ["pyyaml"]="yaml"  # pyyaml is imported as 'yaml'
    )

    for pkg in numpy scipy fastapi uvicorn pydantic jinja2 pyyaml; do
        import_name="${pkg_imports[$pkg]}"
        if python3 -c "import ${import_name}" 2>/dev/null; then
            echo -e "  ✓ ${pkg}"
        else
            echo -e "  ${RED}✗ ${pkg} (missing)${NC}"
            missing=1
        fi
    done

    # Check optional dependencies
    for pkg in openpyxl reportlab pillow; do
        if python3 -c "import ${pkg}" 2>/dev/null; then
            echo -e "  ✓ ${pkg} (optional)"
        else
            echo -e "  ${YELLOW}○ ${pkg} (optional, for exports)${NC}"
        fi
    done

    if [ $missing -eq 1 ]; then
        echo -e "\n${RED}Missing required dependencies. Install with:${NC}"
        echo "  pip install numpy scipy fastapi uvicorn pydantic jinja2 pyyaml"
        return 1
    fi

    echo -e "${GREEN}All required dependencies available.${NC}"
    return 0
}

# Kill running server
kill_server() {
    echo -e "${YELLOW}Stopping SFPPy Studio servers...${NC}"

    # Find and kill processes
    local pids=$(pgrep -f "uvicorn.*studio.app.main:app" 2>/dev/null || true)

    if [ -n "$pids" ]; then
        echo "Found processes: $pids"
        kill $pids 2>/dev/null || true
        sleep 1
        # Force kill if still running
        kill -9 $pids 2>/dev/null || true
        echo -e "${GREEN}Server stopped.${NC}"
    else
        echo "No running SFPPy Studio server found."
    fi

    # Also kill by port
    local port_pids=$(lsof -t -i:${DEFAULT_PORT} 2>/dev/null || true)
    if [ -n "$port_pids" ]; then
        echo "Killing process on port ${DEFAULT_PORT}: $port_pids"
        kill $port_pids 2>/dev/null || true
    fi
}

# Check if port is available
check_port() {
    local port=$1
    if lsof -i:${port} > /dev/null 2>&1; then
        echo -e "${RED}Port ${port} is already in use.${NC}"
        echo "Use --kill to stop existing server, or -p to specify another port."
        return 1
    fi
    return 0
}

# Start the server
start_server() {
    local host=$1
    local port=$2
    local debug=$3
    local reload=$4

    print_banner

    echo -e "${BLUE}Starting SFPPy Studio...${NC}"
    echo "  Host: ${host}"
    echo "  Port: ${port}"
    echo "  Debug: ${debug}"
    echo "  Auto-reload: ${reload}"
    echo ""

    # Check port availability
    check_port ${port} || exit 1

    # Change to SFPPy root
    cd "${SCRIPT_DIR}"

    # Build uvicorn command
    local cmd="python3 -m uvicorn studio.app.main:app --host ${host} --port ${port}"

    if [ "${reload}" = "true" ]; then
        cmd="${cmd} --reload --reload-dir studio"
    fi

    if [ "${debug}" = "true" ]; then
        cmd="${cmd} --log-level debug"
    else
        cmd="${cmd} --log-level info"
    fi

    echo -e "${GREEN}Server starting at http://${host}:${port}${NC}"
    echo -e "${YELLOW}Press Ctrl+C to stop${NC}"
    echo ""

    # Run the server
    exec ${cmd}
}

# Parse arguments
HOST="${DEFAULT_HOST}"
PORT="${DEFAULT_PORT}"
DEBUG="false"
RELOAD="false"
CHECK_ONLY="false"
KILL_SERVER="false"

while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--port)
            PORT="$2"
            shift 2
            ;;
        -H|--host)
            HOST="$2"
            shift 2
            ;;
        -d|--debug)
            DEBUG="true"
            shift
            ;;
        -r|--reload)
            RELOAD="true"
            shift
            ;;
        --check)
            CHECK_ONLY="true"
            shift
            ;;
        --kill)
            KILL_SERVER="true"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            usage
            exit 1
            ;;
    esac
done

# Execute requested action
if [ "${KILL_SERVER}" = "true" ]; then
    kill_server
    exit 0
fi

if [ "${CHECK_ONLY}" = "true" ]; then
    check_dependencies
    exit $?
fi

# Normal startup
check_dependencies || exit 1
echo ""
start_server "${HOST}" "${PORT}" "${DEBUG}" "${RELOAD}"
