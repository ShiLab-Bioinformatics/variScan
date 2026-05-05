export CGO_ENABLED=0

# Variables
BINARY_DIR=bin
MODULE_NAME=github.com/username/project-name
GO_FILES=cmd/find-best-align.go cmd/match-two-ends.go
GOOPT=-ldflags="-s -w"

# Default target: builds everything
all: build

# Create the binary directory and build both tools
build:
	@mkdir -p $(BINARY_DIR)
	go build $(GOOPT) -o $(BINARY_DIR)/find-best-align cmd/find-best-align.go
	go build $(GOOPT) -o $(BINARY_DIR)/match-two-ends cmd/match-two-ends.go
	go build $(GOOPT) -o $(BINARY_DIR)/build-xlsx cmd/build-xlsx.go

# Clean up binaries
clean:
	rm -rf $(BINARY_DIR)

# This ensures your dependencies are synced before building
tidy:
	go mod tidy
