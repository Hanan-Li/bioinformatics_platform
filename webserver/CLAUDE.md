# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture Overview

This is a bioinformatics platform with a full-stack web application architecture:

### Backend (Flask API)
- **Location**: `/api/` directory
- **Framework**: Flask with Blueprint-based modular structure
- **Database**: SQLite with custom schema in `api/data/schema.sql`
- **Authentication**: JWT-based authentication with custom password hashing (SHA-512 + salt)
- **External Integration**: Seqera Platform API for workflow execution

**Key Modules**:
- `api/__init__.py`: Flask app factory with configuration
- `api/auth.py`: User registration, login, JWT token management  
- `api/workflow.py`: Workflow creation, instantiation, and Seqera integration
- `api/db.py`: SQLite database connection and initialization
- `api/common.py`: JWT token validation decorator
- `api/pipeline_defs/`: Static pipeline definitions (JSON format)

### Frontend (Next.js)
- **Location**: `/bioinformatics_frontend/` directory  
- **Framework**: Next.js 15 with TypeScript
- **Styling**: Tailwind CSS
- **Fonts**: Geist font family

### Database Schema
The application uses SQLite with these main tables:
- `user`: User accounts with SHA-512 hashed passwords
- `project`: User projects 
- `workflow`: Workflow definitions
- `instantiation`: Workflow execution instances linked to Seqera

### External Services
- **Seqera Platform**: Cloud-based workflow execution platform
- **API Integration**: Uses Bearer token authentication to launch workflows

## Development Commands

### Frontend Development
```bash
cd bioinformatics_frontend
npm run dev      # Start development server (localhost:3000)
npm run build    # Build for production
npm run start    # Start production server
npm run lint     # Run ESLint
```

### Backend Development
```bash
# From the root webserver directory
export FLASK_APP=api
flask run        # Start Flask development server

# Database initialization
flask init-db    # Initialize SQLite database with schema
```

## Configuration

### Flask Configuration (`api/__init__.py`)
- `SECRET_KEY`: JWT signing key (currently hardcoded for development)
- `DATABASE`: SQLite database path  
- `SEQERA_TOKEN`: Authentication token for Seqera API
- `PIPELINE_DEFINITIONS_DIR`: Directory containing pipeline JSON definitions

### Environment Setup
The application expects these environment variables or Flask config:
- `SEQERA_TOKEN`: Required for workflow execution via Seqera Platform

## API Endpoints

### Authentication (`/api/auth`)
- `POST /api/auth/register`: User registration
- `POST /api/auth/login`: User login (returns JWT token)

### Workflows (`/api/workflow`)  
- `GET /api/workflow/get`: List available pipeline definitions
- `POST /api/workflow/instantiate`: Launch workflow on Seqera Platform
- `GET /api/workflow/create`: Workflow creation endpoint (protected)

All workflow endpoints require JWT token in `x-access-token` header.

## Development Notes

- The codebase uses custom password hashing instead of werkzeug's built-in functions
- JWT tokens have 30-minute expiration
- Pipeline definitions are stored as static JSON files in `api/pipeline_defs/`
- Frontend and backend run as separate services during development
- The application integrates with Nextflow workflows via Seqera Platform