CONFIG ?= Release
USE_THICK_ORBITS ?= OFF
BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all configure reconfigure build test install clean ctest pytest deps-debian deps

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
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON -DUSE_THICK_ORBITS=$(USE_THICK_ORBITS)

configure: $(BUILD_NINJA)

reconfigure:
	cmake -S . -B$(BUILD_DIR) -GNinja -DCMAKE_BUILD_TYPE=$(CONFIG) -DCMAKE_COLOR_DIAGNOSTICS=ON -DUSE_THICK_ORBITS=$(USE_THICK_ORBITS)

build: configure
	cmake --build $(BUILD_DIR) --config $(CONFIG)

test: ctest pytest

ctest: build
	cd $(BUILD_DIR) && ctest

pytest: build
	cd test/ripple_plateau && python3 test_ripple_plateau.py

doc: configure
	cmake --build --preset default --target doc

clean:
	rm -rf $(BUILD_DIR)
