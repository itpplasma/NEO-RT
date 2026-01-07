CONFIG ?= Release
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all configure reconfigure build test install clean ctest pytest deps-debian deps doc golden

all: build

deps:
	@. /etc/os-release; \
	case $$ID in \
		debian|ubuntu) $(MAKE) deps-debian ;; \
		*) echo "Unsupported OS: $$ID". Please install dependencies manually. && exit 1 ;; \
	esac

deps-debian:
	sudo apt-get update
	sudo apt-get install -y build-essential gfortran cmake ninja-build python3-numpy
	sudo apt-get install -y libblas-dev liblapack-dev libsuitesparse-dev libnetcdf-dev libnetcdff-dev

$(BUILD_NINJA):
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

configure: $(BUILD_NINJA)

reconfigure:
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG)

test: ctest pytest

ctest: build
	ctest --test-dir $(BUILD_DIR) --output-on-failure

pytest: build
	cd test && python -m pytest -v

golden:
	python3 test/golden_record/ensure_golden.py

doc:
	@if ! command -v ford >/dev/null; then echo "FORD executable not found. Install it with 'pip install ford'."; exit 1; fi
	mkdir -p build/doc
	ford doc.md

clean:
	rm -rf $(BUILD_DIR)

include examples/Makefile
