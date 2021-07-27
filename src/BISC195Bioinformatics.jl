module BISC195Bioinformatics

export normalizeDNA

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    #for base in seq
        # note: `N` indicates an unknown base
        #occursin(base, "AGCTNR") || error("invalid base $base")
    #end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

export composition
"""
    composition(::AbstractString)

Counts number of each base type in a sequence.
('N' for unknown bases).
Returns a Tuple.
"""
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
"""
    gc_content(::AbstractString)

Calculates GC content for given sequence.
Returns a Float. 
"""
function gc_content(sequence)
    ## Convert sequence to uppercase string
     sequence = normalizeDNA(sequence)
    
    ## Determine length of string
    seqlength = length(sequence)
    
    ## Pull number of G's and C's from `composition()` and add values together
    gc_count = composition(sequence)['G'] + composition(sequence)['C']
    
    ## return GC content as number of G's and C's out of sequence length
    return gc_count / seqlength
end

export complement 
"""
    complement(::AbstractString)

Returns complement sequence of given sequence. 
"""
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
"""
    reverse_complement(::AbstractString)

Returns reverse complement of given sequence.
"""
function reverse_complement(sequence)
    sequence = normalizeDNA(sequence)
    comp_seq = complement(sequence)
    rev_comp_seq = reverse(comp_seq)

    return rev_comp_seq
end

export parse_fasta 
"""
    parse_fasta(path)

Given path to a fasta file, stores headers and sequences in separate vectors. 
Returns a tuple. 
"""
function parse_fasta(path)
    headers = []
    sequences = []
    current_seq = ""

    for line in eachline(path)                
        if startswith(line, '>')
            line = chop(line, head = 1, tail = 0)
            push!(headers, line)

            if current_seq != ""               
                push!(sequences, current_seq)
                current_seq = ""                                   
            end
        else
            line = normalizeDNA(line)
            
            current_seq = current_seq * line
        end
    end
    push!(sequences, current_seq)             
    return tuple(headers, sequences)
end

export mean_and_std_fasta
"""
    mean_and_std_fasta(path)

Given path to a fasta file, returns mean and standard deviation of sequence length and gc_content.

Returns tuple (mean length, std length, mean gc content, std gc content)
"""

using Statistics 

function mean_and_std_fasta(path)

    sequences = parse_fasta(path)[2]

    seq_lengths = []
    gc_contents = []

    for sequence in sequences
        push!(seq_lengths, length(sequence))
        push!(gc_contents, gc_content(sequence))
    end

    mean_length = mean(seq_lengths)
    std_length = std(seq_lengths)
    mean_gc_content = mean(gc_contents)
    std_gc_content = std(gc_contents)
    
    tuple(mean_length, std_length, mean_gc_content, std_gc_content)
end

export sequence_lengths
"""
    sequence_lengths(path)
Given path to a fasta file, return a vector with length of each sequence.

"""
function sequence_lengths(path)
    seq_lengths = []
    sequences = parse_fasta(path)[2]

    for sequence in sequences 
        push!(seq_lengths, length(sequence))
    end

    return seq_lengths
end

end # module BISC195Bioinformatics
