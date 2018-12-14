# USAGE:
#   awk -v repl=replacement.txt -f script.awk sequences.fa

BEGIN {
    if (repl == "") {
        printf("repl is unset.\n") >"/dev/stderr"
        err = 1
        exit
    }
    err = 0
}

function print_nicely(string) {
    # Prints the string and wraps it at 70 characters.
    m = 1
    while (m <= length(string)) {
        print substr(string, m, 70)
        m += 70
    }
}

function split_and_print(string) {
    # Splits the string on N+ and prints each bit preceeded by
    # the next sequence header read from the repl file.

    n = split(string, a, "N+")
    for (i = 1; i <= n; ++i) {
        if ((getline <repl) != 1) {
            printf("Can not read sequence name from %s\n", repl) >"/dev/stderr"
            err = 1
            exit
        }
        printf(">%s\n", $2 "|pilon")
        print_nicely(a[i])
    }
}

/^>/ {
    # We found a sequence header. Print the current sequence
    # (if we have any) and continue with the next line.
    if (sequence != "") {
        split_and_print(sequence)
        sequence = ""
    }
    next
}

{
    # Accumulate sequence.
    sequence = sequence $0
}

END {
    # Unless we are here due to an error, and if we still have sequence
    # data that is not outputted (the last sequence in the input),
    # output the sequence.
    if (!err && sequence != "")
        split_and_print(sequence)

    exit err
}

