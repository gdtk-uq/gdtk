_lmr_completions() {
    local cur
    cur="${COMP_WORDS[COMP_CWORD]}"
    
    if [ "${COMP_CWORD}" -eq 1 ]; then
        local opts=$(lmr completions 2>/dev/null)
        COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
    else
        COMPREPLY=()
    fi
    
    return 0
}

# If no completions, default to completing files
complete -F _lmr_completions -o default lmr
