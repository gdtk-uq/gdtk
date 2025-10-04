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

INSTALL= install -D -p
INSTALL_EXEC= $(INSTALL) -m 0755
INSTALL_DATA= $(INSTALL) -m 0644

# ----------------------------------------------------------------------------
#                   Validate external options from env vars
# ----------------------------------------------------------------------------
# Currently this is DMD and FLAVOUR

DMD ?= dmd
allowed_DMD_opts := dmd ldc2
ifeq ($(call isoneof, $(DMD), $(allowed_DMD_opts)), false)
    $(warning DMD=$(DMD) is not available as an option for setting the D compiler.)
    $(warning Available options are:${nl}    $(allowed_DMD_opts))
    $(error Exiting unsuccessfully.${nl}Please fix DMD selection on command line or as environment variable setting.)
else ifeq ($(DMD), dmd)
	using_dmd := true
else ifeq ($(DMD), ldc2)
	using_ldc := true
endif

FLAVOUR ?= debug
allowed_FLAVOUR_opts := debug fast profile
ifeq ($(call isoneof, $(FLAVOUR), $(allowed_FLAVOUR_opts)), false)
    $(warning FLAVOUR=$(FLAVOUR) is not available as an option for building.)
    $(warning Available options are: ${nl}    $(allowed_FLAVOUR_opts))
    $(error Exiting unsuccessfully.${nl}Please fix FLAVOUR selection.)
endif

# ----------------------------------------------------------------------------
#                    Variables that can be set by the user
# ----------------------------------------------------------------------------

target_type ?= executable
allowed_target_types := executable staticLibrary dynamicLibrary
ifeq ($(call isoneof, $(target_type), $(allowed_target_types)), false)
    $(warning target_type=$(target_type) is not available as an option for building.)
    $(warning Available options are: ${nl}    $(allowed_target_types))
    $(error Exiting unsuccessfully.${nl}Please fix target_type selection.)
endif

imports ?= 
import_paths ?=
d_previews ?=
d_versions ?=

string_imports ?=

dmd_features ?=
dmd_flags ?=

ldc_features ?=
ldc_flags ?=

# ----------------------------------------------------------------------------
#                     Construction of remaining variables
# ----------------------------------------------------------------------------

dmd_features += \
	$(if $(filter dynamicLibrary,$(target_type)),PIC)

ldc_features += \
	$(if $(filter 1,$(WITH_THREAD_SANITIZER)),sanitize=thread) \
	$(if $(filter 1,$(WITH_ADDRESS_SANITIZER)),sanitize=address) \
	$(if $(filter profile,$(FLAVOUR)),profile-generate) \
	$(if $(filter fast,$(FLAVOUR)),fast-math)

ldc_flags += \
	$(if $(filter dynamicLibrary,$(target_type)),-relocation-model=pic)

# Default values
# -w:   enable warnings as errors
# -g:   add symbolic debug info
# NOTE: Ensure that we use *lazy* (recursive) variables here.
compiler_flags ?= -w -g \
	$(if $(filter dynamicLibrary,$(target_type)),-shared -op) \
	$(if $(filter staticLibrary,$(target_type)),-c -op) \

ifeq ($(call isoneof,$(FLAVOUR),fast profile), true)
	compiler_flags += -release -boundscheck=off
    # 2025-01-15: -O2 messes with the piston-in-tube example on MacOS so we back off to -O1
	compiler_flags += $(if $(using_dmd), -O, $(if $(using_mac), -O1, -O2))
	compiler_flags += $(if $(using_dmd), -inline, -enable-inlining)
    # Full link-time optimization does not play nicely on macOS
	ldc_features += $(if $(using_linux), lto=full)
endif

ifeq ($(FLAVOUR), profile)
	compiler_flags += $(if $(using_dmd),-profile)
endif

ifeq ($(FLAVOUR), debug)
	d_versions += enable_fp_exceptions
	ldc_flags += -link-defaultlib-debug
	compiler_flags += $(if $(using_dmd),-debug,-d-debug)
endif

# Default values
linker_libs = dl $(if $(using_mac),d_classic)
linker_paths ?=

linker_flags = $(if $(using_linux),-E)
linker_flags += $(addprefix -L,$(linker_paths))
linker_flags += $(addprefix -l,$(linker_libs))

compiler_flags += $(addprefix -I=,$(import_paths))
compiler_flags += $(addprefix -i=,$(imports))
compiler_flags += $(addprefix -J=,$(string_imports))
compiler_flags += $(addprefix -L,$(linker_flags))
compiler_flags += $(addprefix -preview=,$(d_previews))

ifeq ($(DMD), dmd)
	compiler_flags += $(addprefix -version=,$(d_versions))
	compiler_flags += $(dmd_flags) $(addprefix -f,$(dmd_features))
else ifeq ($(DMD), ldc2)
	compiler_flags += $(addprefix -d-version=,$(d_versions))
	compiler_flags += $(ldc_flags) $(addprefix -f,$(ldc_features))
endif

# ========================== LINE-BY-LINE EXPLANATION =========================
# This is more convoluted than I'd like, and could definitely be done cleaner
# with another language, but I want to avoid adding even more files.
#
# 1. Skip the first line and read from the second line onwards
# 2. Ensure all lines end in ' \', by removing (if it exists) and appending
#    a new ' \'. Yes, it's sed. Yes, I hate the syntax. So much escaping...
#
# 3. Filter out lines referring to stdlib (for readability)
# 4. Sort lines (for readability)
# 5. Prepend with '<target>: .EXTRA_PREREQS = \'
#    This allows us to check the dependencies, without including them
#    in the automatic variables, which clutters our output.
#
# 6. Write out to a temporary file
# 7. Overwrite the generated deps file with our clean one!
# ========================== LINE-BY-LINE EXPLANATION =========================

define cleanup-deps
@echo "Cleaning up generated dependency file"
@tail --lines=+2 $@.deps \
| sed 's| \\$$||;s|$$| \\|' \
| grep --invert-match '$(STDLIB_DIR)' \
| sort \
| { echo '$@: .EXTRA_PREREQS = \'; cat; } \
> $@.deps.tmp
@mv $@.deps.tmp $@.deps
endef

# Remove the stdlib from our deps, for readability
STDLIB_DIR = $(realpath $(dir $(shell which $(DMD))))

