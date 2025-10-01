# Convenience functions
define nl


endef

define isoneof
$(if $(filter $(1),$(2)),true,false)
endef

# Convenience variables to track target
kernel := $(shell uname -s)
ifeq ($(kernel), Darwin)
	using_mac = true
else
	# Assuming no other support
	using_linux = true
endif

# Validate external options from env vars

dmd_allowed_opts := dmd ldc2
ifeq ($(call isoneof, $(DMD), $(dmd_allowed_opts)), false)
    $(warning DMD=$(DMD) is not available as an option for setting the D compiler.)
    $(warning Available options are:${nl}    $(dmd_allowed_opts))
    $(error Exiting unsuccessfully.${nl}Please fix DMD selection on command line or as environment variable setting.)
else ifeq ($(DMD), dmd)
	using_dmd := true
else ifeq ($(DMD), ldc2)
	using_ldc := true
endif

# Variables that can be set by the user
imports ?= 
import_paths ?=
d_previews ?=
d_versions ?=

string_imports ?=

llvm_features ?=
llvm_features += \
	$(if $(filter 1,$(WITH_THREAD_SANITIZER)),sanitize=thread) \
	$(if $(filter 1,$(WITH_ADDRESS_SANITIZER)),sanitize=address) \
	$(if $(filter profile,$(FLAVOUR)),profile-generate) \
	$(if $(filter fast,$(FLAVOUR)),fast-math)

# Default values
# -w:   enable warnings as errors
# -g:   add symbolic debug info
# NOTE: Ensure that we use *lazy* (recursive) variables here.
compiler_flags = -w -g

ifeq ($(call isoneof,$(FLAVOUR),fast profile), true)
	compiler_flags += -release -boundscheck=off
    # 2025-01-15: -O2 messes with the piston-in-tube example on MacOS so we back off to -O1
	compiler_flags += $(if $(using_dmd), -O, $(if $(using_mac), -O1, -O2))
	compiler_flags += $(if $(using_dmd), -inline, -enable-inlining)
    # Full link-time optimization does not play nicely on macOS
	llvm_features += $(if $(using_linux), lto=full)
endif

ifeq ($(FLAVOUR), profile)
	compiler_flags += $(if $(using_dmd),-profile)
endif

ifeq ($(FLAVOUR), debug)
	d_versions += enable_fp_exceptions
	compiler_flags += $(if $(using_dmd),-debug,-d-debug)
endif

# Default values
linker_libs = dl $(if $(using_mac),d_classic)
linker_paths = $(LIBLUAPATH)

linker_flags = $(if $(using_linux),-E)
linker_flags += $(addprefix -L,$(linker_paths))
linker_flags += $(addprefix -l,$(linker_libs))

compiler_flags += $(addprefix -I=,$(import_paths))
compiler_flags += $(addprefix -i=,$(imports))
compiler_flags += $(addprefix -J=,$(string_imports))
compiler_flags += $(addprefix -L,$(linker_flags))
compiler_flags += $(addprefix -preview=,$(d_previews))
compiler_flags += $(addprefix $(if $(using_dmd),-version=,-d-version=),$(d_versions))

ifeq ($(DMD), ldc2)
	compiler_flags += $(addprefix -f,$(llvm_features))
endif
