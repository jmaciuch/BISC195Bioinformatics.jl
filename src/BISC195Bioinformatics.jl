module Assignment07

export normalizeDNA

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

export composition

function composition(sequence)
    sequence = normalizeDNA(sequence) # make uppercase string, check invalid bases
    base_count = Dict('A' => 0,
                      'T' => 0,
                      'G' => 0,
                      'C' => 0,
                      'N' => 0)

    for base in sequence
        if base == 'A'
            base_count['A'] += 1
        elseif base == 'C'
            base_count['C'] += 1
        elseif base == 'G'
            base_count['G'] += 1
        elseif base == 'T'
            base_count['T'] += 1
        elseif base == 'N' 
            base_count['N'] += 1
        end
    end
    return base_count
end

export gc_content

function gc_content(sequence)
    ## Convert sequence to uppercase string
     sequence = normalizeDNA(sequence)
    
    ## throw an error if the string contains anything other than ACGT
    if any(c-> !in(c, ['A','C','G','T','N']), sequence)
        throw(ArgumentError("Sequence must only contain ACGTN"))
    end
    
    ## Determine length of string
    seqlength = length(sequence)
    
    ## Pull number of G's and C's from `composition()` and add values together
    gc_count = composition(sequence)['G'] + composition(sequence)['C']
    
    ## return GC content as number of G's and C's out of sequence length
    return gc_count / seqlength
end

export complement 

function complement(sequence::AbstractString)
    sequence = normalizeDNA(sequence)

    complements = Dict("A" => 'T',
                       "T" => 'A',
                       "G" => 'C',
                       "C" => 'G',
                       "N" => 'N')
    
    complement_seq = ""
    for base in sequence 
        base = uppercase(string(base))
    
        !(base in keys(complements)) && error("Invalid base $base")

        complement_seq = complement_seq * complements[base]
    end
    return complement_seq
end

export reverse_complement

function reverse_complement(sequence)
    sequence = normalizeDNA(sequence)
    comp_seq = complement(sequence)
    rev_comp_seq = reverse(comp_seq)

    return rev_comp_seq
end

export parse_fasta 

function parse_fasta(path)
    headers = []
    sequences = []
    current_seq = ""

    for line in eachline(path)                
        if startswith(line, '>')                     ##When line is a header
            line = chop(line, head = 1, tail = 0)    ##Chop '>' character
            push!(headers, line)                     ##Push line to headers

            if current_seq != ""                     ##Push whatever is being held in current_seq unless it is first line
                push!(sequences, current_seq)
                current_seq = ""                                   
            end
        else                                         ##If line is not header, concatenate line to current_seq
            line = normalizeDNA(line)
            current_seq = current_seq * line
        end
    end
    push!(sequences, current_seq)                    ##Push final current_seq to sequences
    return tuple(headers, sequences)
end
# Your code here.
# Don't forget to export your functions!


end # module Assignment07
