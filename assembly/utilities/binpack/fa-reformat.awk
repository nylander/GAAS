/^>/ && NR > 1 && seq != "" { printf("%s\n", seq); seq = "" }
/^>/    {
    printf("%s\n", $0) >"/dev/stderr"
    print
    next
}

{
    seq = seq $0

    while (length(seq) >= maxlen) {
        printf("%s\n", substr(seq, 1, maxlen))
        seq = substr(seq, maxlen + 1)
    }
}

END {
    if (seq != "")
        printf("%s\n", seq)
}
