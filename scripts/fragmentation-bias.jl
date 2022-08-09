#!/usr/bin/env julia

using ArgParse
using BioAlignments
using BioSequences
using CodecZlib
using FASTX
using Printf
using Statistics

arg_parser = ArgParseSettings()

@add_arg_table arg_parser begin
    "--bam"
        help = "BAM alignment file."
    "--fastq"
        default = ""
        help = "FASTQ.GZ file."
    "--ref"
        help = "Reference sequence in FASTA format."
    "--rna"
        help = "Use 'U' instead of 'T'."
        action = :store_true
    "--left-bp"
        arg_type = Int
        default = 25
        help = "Number base-pairs to the right of DNA cut."
    "--right-bp"
        arg_type = Int
        default = 25
        help = "Number base-pairs to the right of DNA cut."
    "--ignore-R1"
        help = "Ignore R1 reads."
        action = :store_true
    "--ignore-R2"
        help = "Ignore R2 reads."
        action = :store_true
    "--output-format"
        default = "csv"
        help = "Output format: CSV, transfac (for weblogo)."
end

parsed_args = parse_args(ARGS, arg_parser)

function increment_nucl_count!(nucl_count, current_seq, left_bp, right_bp, total_nucl_counter)
    for pos in [1:length(current_seq);]
        if pos > left_bp + right_bp
            break
        end

        if string(current_seq[pos]) == "A"
            nucl_count["A"][pos] += 1
        elseif string(current_seq[pos]) == "C"
            nucl_count["C"][pos] += 1
        elseif string(current_seq[pos]) == "G"
            nucl_count["G"][pos] += 1
        elseif string(current_seq[pos]) == "T"
            nucl_count["T"][pos] += 1
        end

        # Increments total counter.
        if @isdefined total_nucl_counter
           total_nucl_counter[pos] += 1
        end
    end
end

function main()
    # Argument parsing.
    options = Dict()
    for (key,val) in parsed_args
        options[key] = val
    end

    bam = options["bam"]
    fastq_gz = options["fastq"]
    ref = options["ref"]
    left_bp = options["left-bp"]
    right_bp = options["right-bp"]
    use_rna = options["rna"]
    ignore_r1 = options["ignore-R1"]
    ignore_r2 = options["ignore-R2"]
    output_format = options["output-format"]

    nucl_count = Dict()
    nucl_count["A"] = zeros(left_bp + right_bp)
    nucl_count["C"] = zeros(left_bp + right_bp)
    nucl_count["G"] = zeros(left_bp + right_bp)
    nucl_count["T"] = zeros(left_bp + right_bp)

    total_nucl_count  = zeros(left_bp + right_bp)

    if fastq_gz != ""
        # Reading fastq.gz.
        reader = BioSequences.FASTQ.Reader(GzipDecompressorStream(open(fastq_gz)))
        for record in reader
            increment_nucl_count!(nucl_count, uppercase(string(sequence(record))), 0, right_bp, total_nucl_count)
        end
        close(reader)
    else
        # Read reference sequences.
        references = BioSequences.FASTA.Reader(Base.open(ref, "r"), index=index="$ref.fai")

        references_seq = Dict()
        references_length = Dict()

        # Reading bam.
        bam_reader = open(BAM.Reader, bam)

        ref_nucl_count = Dict()
        ref_nucl_count["A"] = 0
        ref_nucl_count["C"] = 0
        ref_nucl_count["G"] = 0
        ref_nucl_count["T"] = 0

        total_ref_nucl_count  = 0

        for reference in references
            for nucl in sort(collect(keys(nucl_count)))
                # Calculates nucleotide frequency in reference sequences.
                current_count = length(collect(eachmatch(Regex(nucl),
                                                         string(BioSequences.FASTA.sequence(reference)))))
                ref_nucl_count[nucl] += current_count
                total_ref_nucl_count += current_count
            end
        end

        ref_nucl_freq = Dict()
        for nucl in sort(collect(keys(nucl_count)))
            ref_nucl_freq[nucl] = ref_nucl_count[nucl] / total_ref_nucl_count
        end

        bam_record = BAM.Record()
        i = 1

        while !eof(bam_reader)
            read!(bam_reader, bam_record)

            if BAM.ismapped(bam_record)
                # Determines how to modify nucleotide frequency from record position,
                # direction and length.
                left_pos = BAM.position(bam_record)
                right_pos = BAM.rightposition(bam_record)

                # Flag positions:
                # 1  - read paired
                # 2  - read mapped in proper pair
                # 3  - read unmapped
                # 4  - mate unmapped
                # 5  - read reverse strand
                # 6  - mate reverse strand
                # 7  - first in pair
                # 8  - second in pair
                # 9  - not primary alignment
                # 10 - read fails platform/vendor quality checks
                # 11 - read is PCR or optical duplicate
                # 12 - supplementary alignment
                flags_binary = reverse(string(BAM.flag(bam_record), base=2, pad=12))

                # Retrieve sequence from the reference.
                ref_name = BAM.refname(bam_record)

                # Parse data that is related to references.
                if ! haskey(references_seq, ref_name)
                   references_seq[ref_name] = sequence(references[ref_name])
                   references_length[ref_name] = length(string(sequence(references[ref_name])))
                end

                # Positions near chromosome starts and ends are ignored.
                if right_pos - right_bp + 1 < 0 ||
                   left_pos - left_bp < 0 ||
                   right_pos + left_bp > references_length[ref_name] ||
                   left_pos+right_bp-1 > references_length[ref_name]
                    continue
                end

                # If reverse.
                if flags_binary[5] == '1'
                    seq = uppercase(string(reverse_complement(DNASequence(string(sequence(references[ref_name])[right_pos-right_bp+1:right_pos+left_bp])))))
                else
                    seq = uppercase(string(sequence(references[ref_name])[left_pos-left_bp:left_pos+right_bp-1]))
                end

                increment_nucl_count!(nucl_count, seq, left_bp, right_bp, total_nucl_count)
            end
        end

        close(bam_reader)
    end

    if output_format == "transfac"
        if use_rna
            println("P0\tA\tC\tG\tU")
        else
            println("P0\tA\tC\tG\tT")
        end
        for pos in [1:length(nucl_count["A"]);]
            if pos - left_bp > 0
               corrected_pos = pos-left_bp
            else
               corrected_pos = pos-left_bp - 1
            end
            print(string(corrected_pos) * "\t" *
                  string(nucl_count["A"][pos]) * "\t" *
                  string(nucl_count["C"][pos]) * "\t" *
                  string(nucl_count["G"][pos]) * "\t" *
                  string(nucl_count["T"][pos]) * "\n")
        end
    else
        println("nucl,pos,freq,avg_freq")
        for nucl in sort(collect(keys(nucl_count)))
            if use_rna && nucl == "T"
                nucl_visual = "U"
            else
                nucl_visual = nucl
            end

            # Calculates nucleotide frequency in reads.
            for pos in [1:length(nucl_count[nucl]);]
                if pos - left_bp > 0
                   corrected_pos = pos-left_bp
                else
                   corrected_pos = pos-left_bp - 1
                end

                if fastq_gz != ""
                    print(nucl_visual * "," * string(corrected_pos) * "," *
                          string(round(nucl_count[nucl][pos] / total_nucl_count[pos],
                                       digits=4)) * "\n")
                else
                    print(nucl_visual * "," * string(corrected_pos) * "," *
                          string(round(nucl_count[nucl][pos] / total_nucl_count[pos],
                                       digits=4)) *
                          "," * string(round(ref_nucl_freq[nucl], digits=4)) * "\n")
                end
            end
        end
    end
end

main()
