# SFPPy Docker Configuration 🍏⏩🍎

Docker deployment for SFPPy applications.

## Quick Start

```bash
# Start all services
cd SFPPy
docker-compose -f docker/docker-compose.yml up

# Start specific service
docker-compose -f docker/docker-compose.yml up sfppy-studio
```

## Services

| Service | Port | URL | Description |
|---------|------|-----|-------------|
| **Studio** | 8002 | http://localhost:8002 | Full migration GUI |
| **Survey Editor** | 8000 | http://localhost:8000 | Scenario configuration |
| **Survey Simulator** | 8001 | http://localhost:8001 | Batch processing |

## Docker Images

| Image | Dockerfile | Size | Use Case |
|-------|------------|------|----------|
| `sfppy-core` | `Dockerfile.core` | ~500MB | Scripts, examples |
| `sfppy-studio` | `Dockerfile.studio` | ~600MB | Interactive analysis |
| `sfppy-survey` | `Dockerfile.survey` | ~600MB | Population studies |
| `sfppy-full` | `Dockerfile.full` | ~800MB | All-in-one deployment |

## Build Images

```bash
# Build specific image
docker build -f docker/Dockerfile.studio -t sfppy-studio .

# Build all via compose
docker-compose -f docker/docker-compose.yml build
```

## Run Commands

### Studio (Interactive GUI)

```bash
docker run -p 8002:8002 \
  -v $(pwd)/data:/app/data \
  sfppy-studio
```

### Survey (Both Apps)

```bash
# Editor
docker run -p 8000:8000 sfppy-survey \
  python -m uvicorn survey.app.main:app --host 0.0.0.0 --port 8000

# Simulator
docker run -p 8001:8001 sfppy-survey \
  python -m uvicorn survey.app.simulator:app --host 0.0.0.0 --port 8001
```

### Core (Run Examples)

```bash
docker run -it --rm \
  -v $(pwd)/examples:/app/examples \
  -v $(pwd)/tmp:/app/tmp \
  sfppy-core python examples/example1.py
```

## Volumes

The compose file defines persistent volumes:

| Volume | Purpose |
|--------|---------|
| `sfppy-data` | Shared data files |
| `sfppy-studio-exports` | Studio export files |
| `sfppy-studio-jobs` | Studio job storage |
| `sfppy-survey-cache` | Survey computation cache |
| `sfppy-pubchem-cache` | PubChem data cache |

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `PYTHONUNBUFFERED` | `1` | Immediate log output |

## Health Checks

All services include health checks:

```bash
# Check service health
docker inspect --format='{{.State.Health.Status}}' sfppy-studio
```

## Troubleshooting

### Port Already in Use

```bash
# Find process using port
lsof -i :8002

# Or use different port
docker run -p 9002:8002 sfppy-studio
```

### Permission Issues

```bash
# Run with user mapping
docker run --user $(id -u):$(id -g) sfppy-studio
```

### View Logs

```bash
# All services
docker-compose -f docker/docker-compose.yml logs -f

# Specific service
docker-compose -f docker/docker-compose.yml logs -f sfppy-studio
```

## Production Deployment

For production, consider:

1. Use `restart: always` in compose
2. Set up reverse proxy (nginx/traefik)
3. Configure SSL/TLS
4. Set resource limits
5. Use external volumes for persistence

```yaml
# Example production settings
services:
  sfppy-studio:
    restart: always
    deploy:
      resources:
        limits:
          memory: 2G
```
