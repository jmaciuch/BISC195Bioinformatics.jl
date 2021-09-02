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
end ##normalizeDNA

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
end ## composition

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
end ##gc_content

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
end ##complement

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
end ## reverse_complement

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
end ## parse_fasta

export mean_and_std_fasta

using Statistics

"""
    mean_and_std_fasta(path)

Given path to a fasta file, returns mean and standard deviation of sequence length and gc_content.

Returns tuple (mean length, std length, mean gc content, std gc content)
"""
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
end ## mean_and_std_fasta

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
end ## sequence_lengths

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
end ## morethan25k

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
end ## uniquekmers

export slice_fasta_var
"""
    slice_fasta_var(path, variants, k)

Given a path to a fasta file and array of pre-defined variant arrays, outputs array of tuples containing k number of sequences 
from each variant in variants and the variant name. 

Example: 
julia> variant_dict = Dict("Alpha" => "Alpha",
                           "B.1.1.7" => "Alpha",
                           "Beta" => "Beta",
                           "B.1.351" => "Beta",
                           etc...)

julia> test_slice = slice_fasta_var("data/Analysis1_test.fasta", variant_dict, 2)
8-element Vector{Any}:
("Beta", "TTTGCGTTTTTAAAGCGCCCCGATAAGCTAGATCGATCGCGTAGCGCTCAGCTAGCTTAG")
⋮
("Alpha", "CCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACA")

julia> slice_fasta_var("data/Analysis1_test.fasta", variant_dict, 100)
ERROR: Dataset contains less than 100 entries for Beta variant

"""
function slice_fasta_var(path, variant_dict, k)
    headers, sequences = parse_fasta(path)

    pangolins = []
    
    for header in headers                     
        split_header = split(header, "|")
        push!(pangolins, split_header[2])                                          
    end

    var_seq = []                                 ## Array contains tuples of variant name and sequences
    variant_set = Set()                          ## Set containing all variant names present in data set w/o repeats

    i = 1                                        
    for pangolin in pangolins                    
        if haskey(variant_dict, pangolin)                                ## If panglolin in variant_dict
            push!(var_seq, tuple(variant_dict[pangolin], sequences[i]))  ## Tuple & push variant name and sequence
            push!(variant_set, variant_dict[pangolin])       ## Push variant name into set to keep track of variants in dataset
            i += 1
        end
    end  

    indices = []

    for variant in variant_set                   ## finding indices for k # of tuples for each variant
        current_var = Set()                      ## avoiding repeat indices for current variant being evaluated
        
        k_count = 0                              ## Count amount of sequences available for each variant to use in error eval
        for tup in var_seq                       
            if tup[1] == variant
                k_count += 1
            end
        end

        (k_count < k) && error("Dataset contains less than $k entries for $variant variant")

        while length(current_var) < k            ## Pick random number out of all indices for variant
            push!(current_var, rand(findall(tup -> tup[1] == variant, var_seq)))
        end
        
        union!(indices, current_var)             ## Add current index picks to indices array
    end

    return var_seq[indices]
end ## slice_fasta_var

export uniquekmer_mean_and_std
"""
    uniquekmer_mean_and_std(category_sequence, k)

Given an array of tuples of the form ("Category", "Sequence") and a kmer length (k), returns array of tuples of the form 
("Category", mean # unique kmers, std dev). 

Example:
julia> test_slice = slice_fasta_var("data/Analysis1_test.fasta", variant_dict, 2);

julia> uniquekmer_mean_and_std(test_slice, 3)
4-element Vector{Any}:
("Beta", 31.0, 4.242640687119285)
("Gamma", 28.5, 0.7071067811865476)
("Delta", 35.5, 9.192388155425117)
("Alpha", 30.5, 13.435028842544403)
"""
function uniquekmer_mean_and_std(category_sequence, k)
    categories = []                                        ## array containing category names (i.e. "Alpha", "Beta", etc.)
    uniquekmer_counts = []

    current_category = []
    for tup in category_sequence

        if !(tup[1] in categories)                          ## If function encounters new category name
            push!(categories, tup[1])                       ## Push current category to list of all categories

            if length(current_category) != 0                ## If there's something stored in current_category
                push!(uniquekmer_counts, current_category)                ## Push current category to uniquekmer_counts
                current_category = []                                     ## Re-initialize empty set
                push!(current_category, length(uniquekmers(tup[2], k)))   ## Push first kmer count of new category to current_category
            else                                            
                push!(current_category, length(uniquekmers(tup[2], k)))   ## If this is the first iter, just push kmer count to current
            end
        else    
            push!(current_category, length(uniquekmers(tup[2], k))) ## If function is in the middle of a category, proceed pushing count 
        end
    end
    push!(uniquekmer_counts, current_category)              ## Push last category's kmer counts
    
    category_mean_and_std = []

    for i in eachindex(categories)              
        mean_kmer = mean(uniquekmer_counts[i])            
        std_dev_kmer = std(uniquekmer_counts[i])

        push!(category_mean_and_std, tuple(categories[i], mean_kmer, std_dev_kmer))
    end
    return category_mean_and_std
end ## uniquekmer_mean_and_std

export format_date

using Dates

"""
    format_date(date)

Given a date in the format "YYYY-MM-DD", returns date in the formatting for Date package, i.e. Date(YYYY, MM, DD).

Example:
julia> format_date("2021/07/25")
2021-07-25

julia> format_date("20-04-21")
ERROR: Date is not in YYYY-MM-DD format
"""
function format_date(date)
    length(date) == 10 || error("Date is not in YYYY-MM-DD format")  ## Checking that input is properly formatted

    for i in eachindex(date)
        if i != 5 && i != 8             
            occursin(date[i], "0123456789") || error("Date is not in YYYY-MM-DD format")
        end

        if date[i] == '/'
            date = replace(date, "/" => "-")
        end
    end
    
    date_split = split(date, "-")
    
    YYYY = parse(Int64, date_split[1]) ## Converting substrings into numbers for Date function
    MM = parse(Int64, date_split[2])
    DD = parse(Int64, date_split[3])
    
    return Date(YYYY, MM, DD)
end ##format_date

export date_diff
"""
    date_diff(date1, date2)

Returns the length of time (in days) between two dates of the form "YYYY-MM-DD".

Example:
julia> date1 = "2020-01-10";

julia> date2 = "2020-01-02";

julia> test = date_diff(date1, date2)
8

julia> typeof(test)
Int64
"""
function date_diff(date1, date2)
    typeof(date1) == Date || (date1 = format_date(date1))          ##Formatting dates for Dates package operation
    typeof(date2) == Date || (date2 = format_date(date2))

    diff = abs(date2 - date1)
    return Dates.value(diff)                                       ## Convert Date type into integer 
end ## date_diff

export slice_fasta_date
"""
    slice_fasta_date(path, k)

Given a path to a fasta file and a sample size k, provides an array of tuples of the form (collection date, sequence).

Example:
julia> path = "data/Analysis2_test.fasta"

julia> test = slice_fasta_date(path, 3);

julia> test[1]
3-element Vector{Any}:
 2021-04-25
 2021-04-23
 2021-03-03

 julia> test[2]
 3-element Vector{Any}:
 "NGVKGFNCYFPLQSY..."
"""
function slice_fasta_date(path, k)
    headers, sequences = parse_fasta(path)

    length(headers) >= k || error("Data set contains less than $k entries")

    dates = []
    
    for header in headers
        split_header = split(header, "|")
        date = split_header[3]

        if length(date) == 7                 ## If collection date only gives month and year i.e. 2020-07
            date = date * "-01"              ## Add default day (1st of the month)
            push!(dates, date)
        else
            push!(dates, date)
        end                                       
    end

    dates_formatted = []

    for i in eachindex(dates)                ## Format dates
        date = format_date(dates[i])
        push!(dates_formatted, date)
    end
        
    indices = []

    while length(indices) < k               ## Pick k (non-repeating) index numbers out of dates array
        index = rand(1:length(dates_formatted), 1)
        !(index in indices) && union!(indices, index)
    end

    return tuple(dates_formatted[indices], sequences[indices])
end ## slice_fasta_date

export time_vs_align_score

using BioAlignments
"""
    time_vs_align_score(path, k)

Given a path to a fasta file and a sample size (k), produces array of data points with length of time between collection dates and
protein alignment score. 
Produces one data point for every match pair within fasta file.
    
Score model only accommodates fasta files containing protein sequences. 

Example:
julia> time_vs_align_score(path, 3)
(Any[2.0, 51.0, 53.0], Any[-32.0, -30.0, -29.0])
"""
function time_vs_align_score(path, k)
    dates, sequences = slice_fasta_date(path, k)

    align_mat = zeros(length(sequences), length(sequences))       ## Initializing matrices for alignment and date comparison
    time_mat = zeros(length(dates), length(dates))

    scoremodel = AffineGapScoreModel(BLOSUM62, gap_open=-10, gap_extend=-1) ## Score model for alignment algorithm

    times = []
    scores = []

    for i in 1:length(dates)         
        for j in 1:length(dates)
            if i < j                ## for upper half of each matrix triangle, calculate and push time and score
                align_mat[i, j] = score(pairalign(GlobalAlignment(), sequences[i], sequences[j], scoremodel))
                time_mat[i, j] = date_diff(dates[i], dates[j])
                push!(times, time_mat[i, j])
                push!(scores, align_mat[i, j])
            end
        end
    end
    return tuple(times, scores)
end ## time_vs_align_score

end # module BISC195Bioinformatics