function __fish_lmr_completions -a args
    # Don't fail if lmr doesn't exist
    if not command -q lmr
        return
    end
    lmr completions $args 2>/dev/null
end

complete -c lmr -n "command -q lmr" -f
complete -c lmr -n "not __fish_seen_subcommand_from (__fish_lmr_completions)" -ka "(__fish_lmr_completions --describe)"
complete -c lmr -n "__fish_seen_subcommand_from (__fish_lmr_completions)" -F

complete -c lmr-debug -w lmr
