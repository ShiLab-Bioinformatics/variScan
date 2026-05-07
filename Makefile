export CGO_ENABLED=0

# Variables
VARIVER=0.1.0
BINARY_DIR=bin
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

release:
	rm -rf ~/temp/VariantAlign-$(VARIVER) ~/temp/VariantAlign-$(VARIVER)-binary-Linux-x86_64.tar.gz 
	mkdir ~/temp/VariantAlign-$(VARIVER)
	cp -r bin ~/temp/VariantAlign-$(VARIVER)
	cp README.md LICENSE run_VariantAlign.sh  ~/temp/VariantAlign-$(VARIVER)
