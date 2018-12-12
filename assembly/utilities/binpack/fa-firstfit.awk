{
    seq_size[$1] = $2
    seq[nseq++] = $1
    totsize += $2
}

END {
    maxsize = int(totsize/maxcontigs + padding*nseq/maxcontigs)
    maxbin = 1000

    while (maxbin >= maxcontigs) {
        delete sizes
        delete names
        delete bins
        maxbin=0

        for (i = 0; i < nseq; ++i) {
            bin = 0
            while (sizes[bin] + seq_size[seq[i]] + (sizes[bin] == 0 ? 0 : padding) > maxsize) {
                if (sizes[bin] == 0)
                    break
                ++bin
            }

            sizes[bin] += seq_size[seq[i]] + (sizes[bin] == 0 ? 0 : padding)
            names[bin] = (names[bin] == "" ? seq[i] : sprintf("packed_contig-%03d", bin + 1))
            bins[seq[i]] = bin

            maxbin = (bin > maxbin ? bin : maxbin)
        }

        printf("max contig size = %d --> %d contigs\n", maxsize, maxbin + 1) >"/dev/stderr"
        if (maxbin >= maxcontigs)
            maxsize += int(maxsize/100)
    }
}

END {
    for (i in bins)
        printf("%s\t%s\t%d\n", names[bins[i]], i, seq_size[i])

    printf("Produced %d contigs from %d sequences\n", maxbin + 1, NR) >"/dev/stderr"
}
