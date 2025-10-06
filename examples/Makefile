# Example helper targets for running bundled NEO-RT cases.
# Add new example metadata by following the EXAMPLE_<name>_* pattern below.

EXAMPLES_ROOT := $(patsubst %/,%,$(dir $(lastword $(MAKEFILE_LIST))))
NEO_RT_BIN := $(abspath $(BUILD_DIR)/neo_rt.x)

# Metadata for available examples
EXAMPLE_base_DIR := $(EXAMPLES_ROOT)/base
EXAMPLE_base_RUN := driftorbit

EXAMPLE_TARGETS := base

# Helpers for looking up example metadata by name.
example_var = EXAMPLE_$(1)_$(2)
example_value = $($(call example_var,$1,$2))

.PHONY: examples-list examples-%

examples-list:
	@printf "Available examples: %s\n" "$(EXAMPLE_TARGETS)"

examples-%: build | $(NEO_RT_BIN)
	$(call run_example,$*)

# Internal helper: run a single example by name.
define run_example
	@ex_dir="$(call example_value,$1,DIR)"; \
	ex_run="$(call example_value,$1,RUN)"; \
	if [ -z "$$ex_dir" ]; then \
		echo "Unknown example '$1'."; \
		exit 1; \
	fi; \
	echo "Running example '$1' with run name '$$ex_run'"; \
	cd "$$ex_dir" && $(NEO_RT_BIN) "$$ex_run"
endef
