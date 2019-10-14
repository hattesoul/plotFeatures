# ~/.bashrc: executed by bash(1) for non-login shells.
# see /usr/share/doc/bash/examples/startup-files (in the package bash-doc)
# for examples

# If not running interactively, don't do anything
case $- in
    *i*) ;;
      *) return;;
esac

# don't put duplicate lines or lines starting with space in the history.
# See bash(1) for more options
HISTCONTROL=ignoreboth

# append to the history file, don't overwrite it
shopt -s histappend

# for setting history length see HISTSIZE and HISTFILESIZE in bash(1)
HISTSIZE=1000
HISTFILESIZE=2000

# check the window size after each command and, if necessary,
# update the values of LINES and COLUMNS.
shopt -s checkwinsize

# If set, the pattern "**" used in a pathname expansion context will
# match all files and zero or more directories and subdirectories.
#shopt -s globstar

# make less more friendly for non-text input files, see lesspipe(1)
[ -x /usr/bin/lesspipe ] && eval "$(SHELL=/bin/sh lesspipe)"

# set variable identifying the chroot you work in (used in the prompt below)
if [ -z "${debian_chroot:-}" ] && [ -r /etc/debian_chroot ]; then
    debian_chroot=$(cat /etc/debian_chroot)
fi

# set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

# uncomment for a colored prompt, if the terminal has the capability; turned
# off by default to not distract the user: the focus in a terminal window
# should be on the output of commands, not on the prompt
#force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
	# We have color support; assume it's compliant with Ecma-48
	# (ISO/IEC-6429). (Lack of such support is extremely rare, and such
	# a case would tend to support setf rather than setaf.)
	color_prompt=yes
    else
	color_prompt=
    fi
fi

# store colors
# RESET="\[$(tput sgr0)\]"
RESET="\[\e[0m\]"
BOLD="\[\e[1m\]"
MAGENTA="\[\033[0;35m\]"
YELLOW="\[\033[01;33m\]"
BLUE="\[\033[00;34m\]"
LIGHT_GRAY="\[\033[0;37m\]"
CYAN="\[\033[0;36m\]"
GREEN="\[\033[00;32m\]"
RED="\[\033[0;31m\]"
VIOLET="\[\033[01;35m\]"
MYYELLOW="\[\033[38;5;11m\]"
MYCYAN="\[\033[38;5;14m\]"

FGDEFAULT="\[\e[38;5;39m\]"
FGBLACK="\[\e[38;5;0m\]"
FGGRAY="\[\e[38;5;242m\]"
FGMAGENTA="\[\e[38;5;104m\]"
FGBLUE="\[\e[38;5;75m\]"
FGGREEN="\[\e[38;5;77m\]"
FGYELLOW="\[\e[38;5;185m\]"
FGRED="\[\e[38;5;197m\]"

BGDEFAULT="\[\e[48;5;49m\]"
BGBLACK="\[\e[48;5;40m\]"
BGMAGENTA="\[\e[48;5;104m\]"
BGBLUE="\[\e[48;5;75m\]"
BGGREEN="\[\e[48;5;77m\]"
BGYELLOW="\[\e[48;5;185m\]"
BGRED="\[\e[48;5;197m\]"

function color_my_prompt {
  # local __user_and_host="$GREEN\u@\h"
  # local __cur_location="$BLUE\W"           # capital 'W': current directory, small 'w': full file path
  local __git_branch_color="$FGGREEN"
  # local __prompt_tail="$VIOLET$"
  # local __user_input_color="$GREEN"
  # local __git_branch='$(__git_ps1)';

  # colour branch name depending on state
  if [[ "$(__git_ps1)" =~ "*" ]]; then     # if repository is dirty
      __git_branch_color="$FGRED"
  elif [[ "$(__git_ps1)" =~ "$" ]]; then   # if there is something stashed
      __git_branch_color="$FGYELLOW"
  elif [[ "$(__git_ps1)" =~ "%" ]]; then   # if there are only untracked files
      __git_branch_color="$FGGRAY"
  elif [[ "$(__git_ps1)" =~ "+" ]]; then   # if there are staged files
      __git_branch_color="$FGMAGENTA"
  fi

  # Build the PS1 (Prompt String)
  # PS1="$__user_and_host $__cur_location$__git_branch_color$__git_branch $__prompt_tail$__user_input_color "
  git_branch() {
      MYBRANCH=$(git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/\1/')
      if test ${#MYBRANCH} -gt 0; then
        MYBRANCH="[$__git_branch_color$MYBRANCH$RESET]"
        # MYBRANCH=$MYNEWBRANCH
      fi
  }

  if [ "$color_prompt" = yes ]; then
      # PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
      # PS1="\d \t \[$(tput bold)\]\[$(tput sgr0)\]\[\033[38;5;11m\]\u\[$(tput sgr0)\]\[$(tput sgr0)\]\[\033[38;5;15m\]@\[$(tput sgr0)\]\[\033[38;5;11m\]\h\[$(tput sgr0)\]\[\033[38;5;15m\]:\[$(tput sgr0)\]\[\033[38;5;14m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\]\\$\[$(tput sgr0)\]"
      # export PS1="\[\033[38;5;0m\]\[\033[48;5;97m\]\d\[$(tput sgr0)\]\[\033[38;5;15m\]\[\033[48;5;-1m\]>\[$(tput sgr0)\]\[\033[38;5;0m\]\[\033[48;5;32m\]\t\[$(tput sgr0)\]\[\033[38;5;15m\]\[\033[48;5;-1m\]>\[$(tput sgr0)\]\[\033[38;5;0m\]\[\033[48;5;41m\]\u@\h\[$(tput sgr0)\]\[\033[38;5;15m\]\[\033[48;5;-1m\]>\[$(tput sgr0)\]\[\033[38;5;0m\]\[\033[48;5;221m\]\w\[$(tput sgr0)\]\[\033[38;5;15m\]\[\033[48;5;-1m\]>\\$ \$?\[$(tput sgr0)\]"
      git_branch
      # PS1="📅$RESET\d 🕓\t \[$(tput bold)\]$MYYELLOW\u$RESET@$MYYELLOW\h$RESET:$MYCYAN\w$RESET$MYBRANCH$RESET\$ "
      # PS1="$RESET$BGMAGENTA$FGBLACK▌📆\d$BGBLUE$FGMAGENTA▌$FGBLACK🕓\t$BGGREEN$FGBLUE▌$FGBLACK🙂$BOLD\u$RESET$BGGREEN$FGBLACK@\h$RESET$FGGREEN▌\n$BGYELLOW$FGBLACK▌📂\w$RESET$FGYELLOW▌$RESET$MYBRANCH$RESET\$ "
      PS1="$RESET$BGMAGENTA$FGBLACK▌📆\d$BGBLUE$FGMAGENTA▌$FGBLACK🕓\t$BGGREEN$FGBLUE▌$FGBLACK🤓$BOLD\u$RESET$BGYELLOW$FGGREEN▌$FGBLACK💻\h$BGRED$FGYELLOW▌$FGBLACK📂\w$RESET$FGRED▌$RESET$MYBRANCH\$ "
  else
      PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
  fi
}

# configure PROMPT_COMMAND which is executed each time before PS1
export PROMPT_COMMAND=color_my_prompt

#unset color_prompt force_color_prompt

# if .git-prompt.sh exists, set options and execute it
if [ -f ~/.git-prompt.sh ]; then
  GIT_PS1_SHOWDIRTYSTATE=true
  GIT_PS1_SHOWSTASHSTATE=true
  GIT_PS1_SHOWUNTRACKEDFILES=true
  GIT_PS1_SHOWUPSTREAM="auto"
  GIT_PS1_HIDE_IF_PWD_IGNORED=true
  GIT_PS1_SHOWCOLORHINTS=true
  . ~/.git-prompt.sh
fi

# If this is an xterm set the title to user@host:dir
case "$TERM" in
xterm*|rxvt*)
    PS1="\[\e]0;${debian_chroot:+($debian_chroot)}\u@\h: \w\a\]$PS1"
    ;;
*)
    ;;
esac

# enable color support of ls and also add handy aliases
if [ -x /usr/bin/dircolors ]; then
    test -r ~/.dircolors && eval "$(dircolors -b ~/.dircolors)" || eval "$(dircolors -b)"
    alias ls='ls --color=auto'
    #alias dir='dir --color=auto'
    #alias vdir='vdir --color=auto'

    alias grep='grep --color=auto'
    alias fgrep='fgrep --color=auto'
    alias egrep='egrep --color=auto'
fi

# colored GCC warnings and errors
# export GCC_COLORS='error=01;31:warning=01;35:note=01;36:caret=01;32:locus=01:quote=01'

# some more ls aliases
alias ll='ls -alF'
alias la='ls -A'
alias l='ls -CF'

# Add an "alert" alias for long running commands.  Use like so:
#   sleep 10; alert
alias alert='notify-send --urgency=low -i "$([ $? = 0 ] && echo terminal || echo error)" "$(history|tail -n1|sed -e '\''s/^\s*[0-9]\+\s*//;s/[;&|]\s*alert$//'\'')"'

# Alias definitions.
# You may want to put all your additions into a separate file like
# ~/.bash_aliases, instead of adding them here directly.
# See /usr/share/doc/bash-doc/examples in the bash-doc package.

if [ -f ~/.bash_aliases ]; then
    . ~/.bash_aliases
fi

# enable programmable completion features (you don't need to enable
# this, if it's already enabled in /etc/bash.bashrc and /etc/profile
# sources /etc/bash.bashrc).
if ! shopt -oq posix; then
  if [ -f /usr/share/bash-completion/bash_completion ]; then
    . /usr/share/bash-completion/bash_completion
  elif [ -f /etc/bash_completion ]; then
    . /etc/bash_completion
  fi
fi

# add cellranger path
export PATH="$PATH:/home/hattesoul/cellranger/cellranger-3.1.0"

# start xbindkeys
xbindkeys &
