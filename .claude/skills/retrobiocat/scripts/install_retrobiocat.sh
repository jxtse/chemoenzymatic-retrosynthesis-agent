#!/bin/bash
# RetroBioCat 2.0 Installation Script (using uv)

set -e

echo "=== RetroBioCat 2.0 Installation Script ==="
echo ""

# Detect OS
OS="$(uname -s)"
ARCH="$(uname -m)"

echo "Detected OS: $OS"
echo "Architecture: $ARCH"
echo ""

# Check if uv is installed
if ! command -v uv &> /dev/null; then
    echo "Error: uv not found. Please install uv first:"
    echo "  curl -LsSf https://astral.sh/uv/install.sh | sh"
    exit 1
fi

echo "Using uv package manager"
echo ""

# Check for ARM Mac
if [[ "$OS" == "Darwin" ]] && [[ "$ARCH" == "arm64" ]]; then
    echo ""
    echo "Detected ARM Mac (Apple Silicon)"
    echo "Installing HDF5 dependencies..."

    # Check if homebrew is installed
    if ! command -v brew &> /dev/null; then
        echo "Error: Homebrew not found. Please install Homebrew first:"
        echo "  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
        exit 1
    fi

    brew install hdf5 c-blosc

    echo "Setting environment variables..."
    export HDF5_DIR=/opt/homebrew/opt/hdf5
    export BLOSC_DIR=/opt/homebrew/opt/c-blosc

    echo "Installing pytables separately..."
    uv pip install tables
fi

# Install RetroBioCat 2.0 into current uv environment
echo ""
echo "Installing RetroBioCat 2.0..."
uv pip install git+https://github.com/willfinnigan/RetroBioCat-2.git

# Verify installation
echo ""
echo "Verifying installation..."
uv run python -c "import rbc2; print('RetroBioCat 2.0 imported successfully!')" || {
    echo "Error: Installation verification failed"
    exit 1
}

# Test basic functionality
echo ""
echo "Testing basic functionality..."
uv run python << EOF
from rbc2 import RetroBioCatExpander
print("Initializing RetroBioCat expander (will download data on first run)...")
expander = RetroBioCatExpander()
print("Testing with simple molecule...")
reactions = expander.get_reactions('CCO')
print(f"Success! Found {len(reactions)} reactions for ethanol")
EOF

echo ""
echo "=== Installation Complete! ==="
echo ""
echo "Quick start (using uv):"
echo "  uv run python -c \""
echo "    from rbc2 import MCTS, get_expanders"
echo "    expanders = get_expanders(['retrobiocat'])"
echo "    mcts = MCTS('CCCC=O', expanders)"
echo "    mcts.run()"
echo "  \""
echo ""
echo "See skill documentation for more examples:"
echo "  .claude/skills/retrobiocat/SKILL.md"
