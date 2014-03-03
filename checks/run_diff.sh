# This file contains a function used in the various example/test
# scripts. It relies on bashisms to get the substrings when running
# the echo commands.

run_diff()
{
    # This function is run after each check. It needs two arguments:
    # $1: first file to compare
    # $2: second file to compare
    # $3: start message to print on the output line before OK or FAILED
    # $4-: any arguments to pass to the diff command
    file1=$1
    file2=$2
    name=""
    if [  ${#} -ge 3 ]; then
        name=$3
        shift 3
        args=$@
    fi

    blanks="                                                                      "

    if diff $args "$file1" "$file2"; then
        echo -e "${name}${blanks:${#name}} OK"
    else
        echo -e "${name}${blanks:${#name}} FAILED"
        exit 1
    fi
}
