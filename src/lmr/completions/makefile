user_shell := $(shell echo $$SHELL)
shell_name := $(notdir $(user_shell))
install_file := install -Dm644

guess:
	$(MAKE) --no-print-directory $(shell_name)

bash: lmr.bash
	$(install_file) $< ~/.local/share/bash-completion/completions/lmr

fish: lmr.fish
	$(install_file) $< ~/.config/fish/completions
	@echo "Please restart your shell for the changes to take effect"
	@echo "Or alternatively run 'source lmr.fish'"

# If using oh-my-zsh, this will direct to the relevant directory
# Else, we install to the default
ZSH ?= ~/.zsh
zsh: lmr.zsh
	$(install_file) $< $(ZSH)/completions/_lmr
	@echo "Please restart your shell for the changes to take effect"
	@echo "Or alternatively run 'compinit'"

remove:
	# Removing all options to be sure
	- rm -f ~/.local/share/bash-completion/completions/lmr
	- rm -f ~/.config/fish/completions/lmr.fish
	- rm -f $(ZSH)/completions/_lmr
	@echo "Please restart your shell for the changes to take effect"
ifeq ($(shell_name),bash)
	@echo "Or alternatively run 'complete -r lmr'"
else ifeq ($(shell_name),fish)
	@echo "Or alternatively run 'complete -c lmr -e'"
endif

.PHONY: guess bash fish zsh remove
