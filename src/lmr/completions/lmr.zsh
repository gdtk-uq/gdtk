#compdef lmr

_lmr() {
  local curcontext="$curcontext" state line
  typeset -A opt_args

  _arguments -C \
    '1: :->subcommand' \
    '*::arg:->args'

  case $state in
    subcommand)
      # Check if lmr command is available
      if (( $+commands[lmr] )); then
        local -a subcommands
        local desc_output
        
        # Get subcommands with descriptions
        desc_output=$(lmr completions -d 2>/dev/null)
        
        # Parse tab-separated subcommands and descriptions
        while IFS=$'\t' read -r cmd desc; do
          [[ -n "$cmd" ]] && subcommands+=("$cmd:$desc")
        done <<< "$desc_output"
        
        # Offer completions
        (( ${#subcommands[@]} > 0 )) && _describe 'lmr subcommand' subcommands
      fi
      ;;
    args)
      # For subcommands, complete files (default behavior)
      _files
      ;;
  esac
}

_lmr "$@"
