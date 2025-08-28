#!/bin/bash

# Protein Interaction Analysis Server - Deployment Script
# This script helps deploy the application to various platforms

set -e

echo "🚀 Protein Interaction Analysis Server - Deployment Helper"
echo "=========================================================="

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check if production requirements file exists
check_prod_requirements() {
    if [ ! -f "requirements-prod.txt" ]; then
        print_warning "requirements-prod.txt not found. Creating optimized production requirements..."
        cp requirements.txt requirements-prod.txt
        # Remove development dependencies from production file
        sed -i '/pytest/d' requirements-prod.txt
        sed -i '/black/d' requirements-prod.txt
        sed -i '/flake8/d' requirements-prod.txt
        sed -i '/psutil/d' requirements-prod.txt
        print_status "Created requirements-prod.txt for deployment"
    fi
}

# Check if required files exist
check_requirements() {
    print_status "Checking deployment requirements..."

    required_files=("server.py" "requirements-prod.txt" ".streamlit/config.toml")
    for file in "${required_files[@]}"; do
        if [ ! -f "$file" ]; then
            print_error "Required file missing: $file"
            exit 1
        fi
    done

    print_status "All required files found ✓"
}

# Test local run
test_local() {
    print_status "Testing local application..."

    # Check if streamlit is installed
    if ! command -v streamlit &> /dev/null; then
        print_error "Streamlit is not installed. Please install it first."
        exit 1
    fi

    # Try to run the app briefly to check for syntax errors
    timeout 10s streamlit run server.py --server.headless true --server.port 8502 > /dev/null 2>&1 || true

    print_status "Local test completed ✓"
}

# Deploy to Streamlit Cloud
deploy_streamlit_cloud() {
    print_status "Preparing for Streamlit Cloud deployment..."

    echo ""
    echo "📋 Streamlit Cloud Deployment Steps:"
    echo "====================================="
    echo "1. Make sure your GitHub repository is public"
    echo "2. Go to https://share.streamlit.io"
    echo "3. Connect your GitHub account"
    echo "4. Select this repository"
    echo "5. Set main file path to: server.py"
    echo "6. Click 'Deploy'"
    echo ""
    echo "Your app will be available at: https://your-app-name.streamlit.app"
    echo ""
    echo "💡 The app uses requirements-prod.txt for faster deployment"
}

# Deploy to Heroku
deploy_heroku() {
    print_status "Preparing for Heroku deployment..."

    # Check if Heroku CLI is installed
    if ! command -v heroku &> /dev/null; then
        print_error "Heroku CLI is not installed."
        echo "Please install it from: https://devcenter.heroku.com/articles/heroku-cli"
        exit 1
    fi

    echo ""
    echo "📋 Heroku Deployment Steps:"
    echo "==========================="
    echo "1. Create a Heroku account at https://heroku.com"
    echo "2. Login to Heroku CLI: heroku login"
    echo "3. Create a new app: heroku create your-protein-app"
    echo "4. Set the stack to container:"
    echo "   heroku stack:set container -a your-protein-app"
    echo "5. Push to Heroku:"
    echo "   git push heroku main"
    echo ""
    echo "Your app will be available at: https://your-protein-app.herokuapp.com"
}

# Deploy using Docker
deploy_docker() {
    print_status "Preparing Docker deployment..."

    # Check if Docker is installed
    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed."
        echo "Please install Docker from: https://docs.docker.com/get-docker/"
        exit 1
    fi

    echo ""
    echo "📋 Docker Deployment Steps:"
    echo "==========================="
    echo "1. Build the Docker image:"
    echo "   docker build -t protein-interaction-app ."
    echo ""
    echo "2. Run locally for testing:"
    echo "   docker run -p 8501:8501 protein-interaction-app"
    echo ""
    echo "3. Deploy to your preferred platform:"
    echo "   - AWS: Use ECS or Elastic Beanstalk"
    echo "   - Google Cloud: Use Cloud Run"
    echo "   - Azure: Use Container Instances"
    echo "   - DigitalOcean: Use App Platform"
    echo ""
    echo "Example for local testing:"
    echo "docker run -p 8501:8501 -e STREAMLIT_SERVER_ADDRESS=0.0.0.0 protein-interaction-app"
}

# Main menu
show_menu() {
    echo ""
    echo "Choose deployment option:"
    echo "1) Streamlit Cloud (Free & Easy)"
    echo "2) Heroku (Professional)"
    echo "3) Docker (Full Control)"
    echo "4) Test Local Application"
    echo "5) Check Requirements"
    echo "6) Exit"
    echo ""
    read -p "Enter your choice (1-6): " choice

    case $choice in
        1)
            deploy_streamlit_cloud
            ;;
        2)
            deploy_heroku
            ;;
        3)
            deploy_docker
            ;;
        4)
            test_local
            ;;
        5)
            check_requirements
            ;;
        6)
            print_status "Goodbye! 👋"
            exit 0
            ;;
        *)
            print_error "Invalid choice. Please select 1-6."
            show_menu
            ;;
    esac
}

# Run checks first
check_prod_requirements
check_requirements

# Show menu
show_menu
