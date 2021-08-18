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
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
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
            line = uppercase(line)
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

export morethan25k
"""
    morethan25k(headers, sequences)

Given a tuple of headers and sequences (output from parse_fasta function), returns tuple containing only sequences, headers, and sequence lengths for sequences with more than 25k bases.
"""
function morethan25k(headers, sequences)
    seq_lengths = [length(sequence) for sequence in sequences] ## pushes sequence lengths for each sequence to new array
    
    keep = findall(x -> x > 25000, seq_lengths) ##keep sequence if sequence length more than 25000 bases
    
    headers_25k = headers[keep]
    sequences_25k = sequences[keep]

    return tuple(headers_25k, sequences_25k, seq_lengths)
end

export uniquekmers
"""
    uniquekmers(sequence, k)

Given a sequence and an integer k, returns all unique sequences of length k. 
Does not recognize kmers containing invalid bases. 

Example: 
julia> seq = "ATGCGATXGTAC";

julia> uniquekmers(seq, 4)
Set{Any} with 5 elements:
  "GTAC"
  "TGCG"
  "ATGC"
  "CGAT"
  "GCGA"
"""
function  uniquekmers(sequence, k)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence") ##Throws error if k is larger than sequence length

    kmers = Set() ## initialize set
     
    stopindex = length(sequence) - k + 1

    for i in 1:stopindex
        kmer = sequence[i:(i+k-1)]
        kmer = uppercase(kmer) 
        
        keep = true 

        for base in kmer       ## Changes flag to false if invalid base encountered
            if !occursin(base, "AGCT") 
                keep = false
                break          ## Stop evalulating kmer once the first invalid base is encountered
            end
        end

        keep && push!(kmers, kmer) ## If keep still true, push kmer to set
    end
    return kmers
end

export splice_fasta
"""
    splice_fasta(path, variants, k)

Given a path to a fasta file and array of pre-defined variant arrays, outputs array of tuples containing k number of sequences 
from each variant in variants and the variant name. 

Example: 
julia> variant_dict = Dict("Alpha" => "Alpha",
                           "B.1.1.7" => "Alpha",
                           "Beta" => "Beta",
                           "B.1.351" => "Beta",
                           "B.1.351.2" => "Beta",
                           "B.1.351.3" => "Beta",
                           ...)

julia> splice_fasta("data/Analysis1_test.fasta", variant_dict, 2)
8-element Vector{Any}:
("Beta", "TTTGCGTTTTTAAAGCGCCCCGATAAGCTAGATCGATCGCGTAGCGCTCAGCTAGCTTAG")
⋮
("Alpha", "CCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACA")

julia> splice_fasta("data/Analysis1_test.fasta", variant_dict, 100)
ERROR: Dataset contains less than 100 entries for Beta variant

"""
function splice_fasta(path, variant_dict, k)
    headers, sequences = parse_fasta(path)

    pangolins = []
    
    for header in headers                     
        split_header = split(header, "|")
        push!(pangolins, split_header[2])                                          
    end

    var_seq = []                                 ## Array contains tuples of variant name and sequences

    i = 1                                        ## If panglolin in dict, pull variant name and pair with sequence
    for pangolin in pangolins
        haskey(variant_dict, pangolin) && push!(var_seq, tuple(variant_dict[pangolin], sequences[i]))
        i += 1
    end  

    variant_set = Set(values(variant_dict))      ## Get a set of all variant names without repeats

    indices = []

    for variant in variant_set                   ## find indices of k # of tuples for each variant
        current_var = Set()
        
        k_count = 0                              ## Count amount of sequences available for each variant to use in error eval
        for tup in var_seq                       
            if tup[1] == variant
                k_count += 1
            end
        end

        (k_count < k) && error("Dataset contains less than $k entries for $variant variant")

        while length(current_var) < k
            push!(current_var, rand(findall(tup -> tup[1] == variant, var_seq)))
        end
        
        union!(indices, current_var)
    end

    return var_seq[indices]
end

end # module BISC195Bioinformatics