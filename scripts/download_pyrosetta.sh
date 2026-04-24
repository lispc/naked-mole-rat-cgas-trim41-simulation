#!/bin/bash
# Download PyRosetta with resume support using curl -C -
# This script can be re-run to resume interrupted downloads

URL="https://conda.rosettacommons.org/osx-arm64/pyrosetta-2025.06+release.e5e4b278be-py312_0.tar.bz2"
OUTDIR="/Users/zhangzhuo/repos/personal/naked-mole-rat-cgas-trim41-simulation/tmp_downloads"
OUTFILE="$OUTDIR/pyrosetta-2025.06+release.e5e4b278be-py312_0.tar.bz2"

echo "=== PyRosetta Download with Resume Support ==="
echo "Target: $OUTFILE"
echo ""

# Check existing file size
if [ -f "$OUTFILE" ]; then
    EXISTING_SIZE=$(stat -f%z "$OUTFILE" 2>/dev/null || stat -c%s "$OUTFILE" 2>/dev/null)
    echo "Existing file size: $EXISTING_SIZE bytes"
else
    echo "No existing file found, starting fresh download"
fi

# Download with resume
echo "Starting/resuming download..."
curl -C - -L -o "$OUTFILE" "$URL" --progress-bar

EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✅ Download completed successfully"
    echo "File: $OUTFILE"
    ls -lh "$OUTFILE"
    
    # Verify SHA256
    echo ""
    echo "Verifying SHA256..."
    EXPECTED_SHA="50f2cbdde1876df1c8c3271bc14657a0794e20032c17adfd482e3e4b9b4980ec"
    ACTUAL_SHA=$(shasum -a 256 "$OUTFILE" | awk '{print $1}')
    if [ "$EXPECTED_SHA" = "$ACTUAL_SHA" ]; then
        echo "✅ SHA256 verified"
    else
        echo "⚠️ SHA256 mismatch!"
        echo "Expected: $EXPECTED_SHA"
        echo "Actual:   $ACTUAL_SHA"
        exit 1
    fi
elif [ $EXIT_CODE -eq 33 ]; then
    # curl exit code 33 means "HTTP server doesn't seem to support byte ranges"
    echo ""
    echo "⚠️ Server doesn't support resume, restarting download..."
    curl -L -o "$OUTFILE" "$URL" --progress-bar
else
    echo ""
    echo "❌ Download interrupted (curl exit code: $EXIT_CODE)"
    echo "Run this script again to resume downloading"
    exit $EXIT_CODE
fi
