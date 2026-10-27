###############################################################################
# PER-SEQUENCE CHARGE O/E ANALYSIS
#
# Goal:
#   Position-specific observed/expected enrichment for charge-changing
#   mutations.
#
# UNIT OF ANALYSIS:
#   ONE SEQUENCE = ONE OBSERVATION.
#
# IMPORTANT:
#   There is NO clone collapsing.
#   Every row in myDb is retained independently.
#
# REGION:
#   IMGT_V = 312 nt = 104 codons.
#
# PRIMARY EVENT CLASSES:
#   positive_gain : neutral -> H/K/R
#   positive_loss : H/K/R -> neutral
#   negative_gain : neutral -> D/E
#   negative_loss : D/E -> neutral
#   reversal      : H/K/R <-> D/E
#
# PRIMARY EVENTS:
#   Exactly ONE nucleotide difference within a codon.
#
# MULTI-NUCLEOTIDE EVENTS:
#   Retained as QC only.
#   Excluded from primary O/E.
#
# EXPECTED MODEL:
#   SHazaM HH_S5F targeting model.
#
# OBSERVED / EXPECTED UNIVERSE:
#   A sequence/codon contributes only when BOTH the observed sequence and
#   germline contain a complete A/C/G/T codon at that IMGT coordinate.
#
# SHORT GERMLINES:
#   germline_alignment may contain fewer than 312 characters.
#
#   Positions beyond the actual germline alignment are treated as MISSING.
#   They are NOT padded with N.
#
#   Missing positions contribute to neither observed nor expected.
#
#   An actual N in germline_alignment remains N and is also non-informative
#   for primary observed events / codon-level expected event paths.
#
# REQUIRED myDb COLUMNS:
#   - sequence_alignment
#   - germline_alignment
###############################################################################


###############################################################################
# 0. PACKAGES
###############################################################################

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(alakazam)
    library(shazam)
    library(Biostrings)
})


###############################################################################
# 1. INPUT
###############################################################################

# Uncomment if myDb is not already loaded:

myDb <- read_tsv(
    "sequences_lambda.tsv",
    show_col_types = FALSE
)

if (!exists("myDb")) {
    stop("myDb is not loaded. Load your Change-O database first.")
}

sequence_col <- "sequence_alignment"
germline_col <- "germline_alignment"

required_columns <- c(
    sequence_col,
    germline_col
)

missing_columns <- setdiff(
    required_columns,
    names(myDb)
)

if (length(missing_columns) > 0) {
    stop(
        "Missing required columns: ",
        paste(missing_columns, collapse = ", ")
    )
}


###############################################################################
# 2. ANALYSIS PARAMETERS
###############################################################################

n_pos <- IMGT_V@seqLength
n_codons <- n_pos / 3

stopifnot(
    n_pos == 312,
    n_codons == 104
)

positive_aas <- c("H", "K", "R")
negative_aas <- c("D", "E")

neutral_aas <- c(
    "A", "C", "F", "G", "I", "L", "M",
    "N", "P", "Q", "S", "T", "V", "W", "Y"
)

event_names <- c(
    "positive_gain",
    "positive_loss",
    "negative_gain",
    "negative_loss",
    "reversal"
)

nucleotide_names <- c(
    "A", "C", "G", "T"
)


###############################################################################
# 3. INPUT-DATASET QC
###############################################################################

cat("\n")
cat("============================================\n")
cat("INPUT DATASET QC\n")
cat("============================================\n")

n_input <- nrow(myDb)

cat(
    "Input sequences:",
    n_input,
    "\n"
)

if (n_input == 0) {
    stop("myDb contains zero sequences.")
}

cat(
    "Unique nucleotide sequences:",
    length(
        unique(
            as.character(
                myDb[[sequence_col]]
            )
        )
    ),
    "\n"
)

cat(
    "Unique germline sequences:",
    length(
        unique(
            as.character(
                myDb[[germline_col]]
            )
        )
    ),
    "\n"
)


###############################################################################
# 4. NORMALIZE ALIGNMENTS TO THE 312-NT ANALYSIS COORDINATE SYSTEM
#
# IMPORTANT:
#
# We retain two separate concepts:
#
#   alignment sequence:
#       the actual characters supplied in the database
#
#   coordinate presence:
#       whether that position actually exists in the supplied alignment
#
# A short sequence/germline is NOT padded with N.
#
# Missing suffix positions are represented internally by "." ONLY after
# recording the original length. These "." positions are explicitly treated
# as missing and cannot contribute to observed or expected calculations.
#
# Existing N remains N.
###############################################################################

normalize_alignment <- function(
    x,
    target_length
) {

    x <- as.character(x)

    if (anyNA(x)) {
        stop(
            "NA values found in alignment column."
        )
    }

    x <- toupper(x)

    # Treat "-" as an alignment gap.
    x <- chartr(
        "-",
        ".",
        x
    )

    # Keep only the first target_length coordinates.
    x <- substr(
        x,
        1,
        target_length
    )

    # Pad SHORT alignments with ".".
    #
    # These positions are explicitly missing.
    x <- vapply(
        x,
        function(seq) {

            len <- nchar(seq)

            if (len < target_length) {

                paste0(
                    seq,
                    paste0(
                        rep(
                            ".",
                            target_length - len
                        ),
                        collapse = ""
                    )
                )

            } else {

                seq
            }
        },
        character(1)
    )

    x
}


sequences <- normalize_alignment(
    myDb[[sequence_col]],
    n_pos
)

germlines <- normalize_alignment(
    myDb[[germline_col]],
    n_pos
)


###############################################################################
# 5. ORIGINAL ALIGNMENT LENGTHS
###############################################################################

get_original_length <- function(x) {

    x <- as.character(x)

    x[is.na(x)] <- ""

    nchar(
        chartr(
            "-",
            ".",
            x
        )
    )
}


sequence_original_length <- get_original_length(
    myDb[[sequence_col]]
)

germline_original_length <- get_original_length(
    myDb[[germline_col]]
)


cat("\n")
cat("============================================\n")
cat("ALIGNMENT LENGTH QC\n")
cat("============================================\n")

cat(
    "Sequences shorter than",
    n_pos,
    "nt:",
    sum(
        sequence_original_length < n_pos
    ),
    "\n"
)

cat(
    "Germlines shorter than",
    n_pos,
    "nt:",
    sum(
        germline_original_length < n_pos
    ),
    "\n"
)

cat(
    "Sequences exactly",
    n_pos,
    "nt:",
    sum(
        sequence_original_length == n_pos
    ),
    "\n"
)

cat(
    "Germlines exactly",
    n_pos,
    "nt:",
    sum(
        germline_original_length == n_pos
    ),
    "\n"
)


###############################################################################
# 6. ALPHABET QC
#
# Allowed:
#   A C G T N .
#
# "." represents an alignment gap OR an internally padded missing suffix.
#
# We keep track of coordinate presence separately, so these cases cannot
# accidentally contribute to observed/expected calculations.
###############################################################################

bad_sequence <- unique(
    unlist(
        regmatches(
            sequences,
            gregexpr(
                "[^ACGTN.]",
                sequences
            )
        )
    )
)

bad_germline <- unique(
    unlist(
        regmatches(
            germlines,
            gregexpr(
                "[^ACGTN.]",
                germlines
            )
        )
    )
)

if (length(bad_sequence) > 0) {

    stop(
        "Unexpected characters in ",
        sequence_col,
        ": ",
        paste(
            bad_sequence,
            collapse = ", "
        )
    )
}

if (length(bad_germline) > 0) {

    stop(
        "Unexpected characters in ",
        germline_col,
        ": ",
        paste(
            bad_germline,
            collapse = ", "
        )
    )
}


cat(
    "Sequence rows containing N:",
    sum(
        grepl(
            "N",
            sequences,
            fixed = TRUE
        )
    ),
    "\n"
)

cat(
    "Germline rows containing N:",
    sum(
        grepl(
            "N",
            germlines,
            fixed = TRUE
        )
    ),
    "\n"
)

cat(
    "Sequence rows containing gaps:",
    sum(
        grepl(
            "\\.",
            sequences
        )
    ),
    "\n"
)

cat(
    "Germline rows containing gaps:",
    sum(
        grepl(
            "\\.",
            germlines
        )
    ),
    "\n"
)


###############################################################################
# 7. COORDINATE-PRESENCE MASKS
#
# A position is "present" when it lies inside the original alignment.
#
# A "." occurring INSIDE the alignment is still a gap and therefore not a
# valid nucleotide.
#
# For a short alignment, all coordinates after the original length are also
# absent.
###############################################################################

sequence_present <- matrix(
    FALSE,
    nrow = n_input,
    ncol = n_pos
)

germline_present <- matrix(
    FALSE,
    nrow = n_input,
    ncol = n_pos
)

for (i in seq_len(n_input)) {

    if (sequence_original_length[i] > 0) {

        sequence_present[
            i,
            seq_len(
                min(
                    sequence_original_length[i],
                    n_pos
                )
            )
        ] <- TRUE
    }

    if (germline_original_length[i] > 0) {

        germline_present[
            i,
            seq_len(
                min(
                    germline_original_length[i],
                    n_pos
                )
            )
        ] <- TRUE
    }
}


###############################################################################
# 8. NUCLEOTIDE-VALIDITY MASKS
###############################################################################

sequence_valid_nt <- matrix(
    FALSE,
    nrow = n_input,
    ncol = n_pos
)

germline_valid_nt <- matrix(
    FALSE,
    nrow = n_input,
    ncol = n_pos
)

for (i in seq_len(n_input)) {

    sequence_chars <- strsplit(
        sequences[i],
        "",
        fixed = TRUE
    )[[1]]

    germline_chars <- strsplit(
        germlines[i],
        "",
        fixed = TRUE
    )[[1]]

    sequence_valid_nt[
        i,
        ] <- (
            sequence_present[i, ] &
            sequence_chars %in% nucleotide_names
        )

    germline_valid_nt[
        i,
        ] <- (
            germline_present[i, ] &
            germline_chars %in% nucleotide_names
        )
}


###############################################################################
# 9. CHARGE CLASSIFICATION
###############################################################################

classify_charge <- function(aa) {

    out <- rep(
        NA_character_,
        length(aa)
    )

    out[
        aa %in% positive_aas
    ] <- "positive"

    out[
        aa %in% negative_aas
    ] <- "negative"

    out[
        aa %in% neutral_aas
    ] <- "neutral"

    out
}


event_matches <- function(
    from_class,
    to_class,
    event
) {

    if (event == "positive_gain") {

        return(
            from_class == "neutral" &
            to_class == "positive"
        )
    }

    if (event == "positive_loss") {

        return(
            from_class == "positive" &
            to_class == "neutral"
        )
    }

    if (event == "negative_gain") {

        return(
            from_class == "neutral" &
            to_class == "negative"
        )
    }

    if (event == "negative_loss") {

        return(
            from_class == "negative" &
            to_class == "neutral"
        )
    }

    if (event == "reversal") {

        return(
            (
                from_class == "positive" &
                to_class == "negative"
            ) |
            (
                from_class == "negative" &
                to_class == "positive"
            )
        )
    }

    stop(
        "Unknown event: ",
        event
    )
}


###############################################################################
# 10. CODON-LEVEL INFORMATIVE MASK
#
# A sequence/codon is informative ONLY when:
#
#   1. all three sequence positions exist,
#   2. all three germline positions exist,
#   3. all six nucleotides are A/C/G/T.
#
# This SAME mask is used for:
#   - observed events
#   - expected targeting
#
# Thus observed and expected have the same sequence-level universe.
###############################################################################

cat("\n")
cat("============================================\n")
cat("CALCULATING CODON-LEVEL COVERAGE\n")
cat("============================================\n")

valid_codon <- matrix(
    FALSE,
    nrow = n_input,
    ncol = n_codons
)

for (j in seq_len(n_codons)) {

    idx <- (
        3 * j - 2
    ):(3 * j)

    valid_codon[, j] <- (
        rowSums(
            sequence_valid_nt[
                ,
                idx,
                drop = FALSE
            ]
        ) == 3 &
        rowSums(
            germline_valid_nt[
                ,
                idx,
                drop = FALSE
            ]
        ) == 3
    )
}

coverage_codon <- colSums(
    valid_codon
)

cat(
    "Total informative sequence-codon observations:",
    sum(
        valid_codon
    ),
    "\n"
)

cat(
    "Codons with zero informative observations:",
    sum(
        coverage_codon == 0
    ),
    "\n"
)


###############################################################################
# 11. FAST HH_S5F TARGETING
#
# This is evaluated on an individual germline sequence.
#
# IMPORTANT:
#   Missing suffix positions remain missing.
#
#   The targeting calculation is performed on the actual existing germline
#   nucleotide sequence, not on an artificial 312-nt N-padded sequence.
#
#   Targeting is then mapped back to the original 312-position coordinate
#   system.
###############################################################################
###############################################################################
# 11. FAST HH_S5F TARGETING
#
# The analysis coordinate system is fixed at 312 IMGT nucleotide positions.
#
# IMPORTANT:
#
#   * germline_alignment is the source sequence.
#   * Existing alignment gaps "." are NOT removed and concatenated across.
#   * Positions outside the actual germline length are missing.
#   * A 5-mer is evaluated only when all five required positions exist and
#     are usable for the SHazaM targeting lookup.
#
# This preserves the original IMGT coordinate system while preventing a
# short germline from being artificially padded with N or from creating
# artificial 5-mers across gaps.
###############################################################################

s5f <- HKL_S5F@targeting

stopifnot(
    nrow(s5f) == 5,
    ncol(s5f) == 3125
)

nucleotide_names <- c(
    "A",
    "C",
    "G",
    "T"
)


###############################################################################
# 11A. PREPARE ONE GERMLINE FOR TARGETING
###############################################################################

prepare_germline_for_targeting <- function(
    germline_seq,
    n_pos = 312
) {

    germline_seq <- toupper(
        as.character(germline_seq)
    )

    germline_seq <- chartr(
        "-",
        ".",
        germline_seq
    )

    if (nchar(germline_seq) > n_pos) {

        germline_seq <- substr(
            germline_seq,
            1,
            n_pos
        )
    }

    if (nchar(germline_seq) < n_pos) {

        germline_seq <- paste0(
            germline_seq,
            paste0(
                rep(
                    ".",
                    n_pos - nchar(germline_seq)
                ),
                collapse = ""
            )
        )
    }

    chars <- strsplit(
        germline_seq,
        "",
        fixed = TRUE
    )[[1]]

    if (length(chars) != n_pos) {

        stop(
            "Germline preparation produced ",
            length(chars),
            " positions instead of ",
            n_pos,
            "."
        )
    }

    # SHazaM ambiguity conversion.
    #
    # Existing ambiguous IUPAC nucleotides become N.
    #
    # "." remains "." and therefore represents an unavailable aligned
    # position for purposes of the targeting calculation.
    chars <- gsub(
        "[MRWSYKVHDB]",
        "N",
        chars
    )

    chars
}


###############################################################################
# 11B. FAST HH_S5F TARGETING FOR ONE GERMLINE
#
# Output:
#
#   rows    = possible mutant nucleotide A/C/G/T
#   columns = IMGT nucleotide positions 1..312
#
# NA means that the position does not have a usable HH_S5F targeting value.
###############################################################################

fast_targeting <- function(
    germline_seq,
    n_pos = 312
) {

    chars <- prepare_germline_for_targeting(
        germline_seq,
        n_pos = n_pos
    )

    targeting <- matrix(
        NA_real_,
        nrow = 4,
        ncol = n_pos,
        dimnames = list(
            nucleotide_names,
            seq_len(n_pos)
        )
    )

    ###########################################################################
    # A 5-mer centered at position p requires:
    #
    #   p-2, p-1, p, p+1, p+2
    #
    # to exist in the alignment.
    #
    # Therefore positions 1, 2, 311, and 312 cannot receive an ordinary
    # complete 5-mer.
    ###########################################################################

    if (n_pos < 5) {
        return(targeting)
    }

    for (p in 3:(n_pos - 2)) {

        five_positions <- (
            p - 2
        ):(p + 2)

        five_mer <- paste0(
            chars[
                five_positions
            ],
            collapse = ""
        )

        #######################################################################
        # A complete 5-mer must have no alignment gaps.
        #
        # N is allowed as an ambiguity character because the SHazaM targeting
        # matrix contains N-containing 5-mers.
        #######################################################################

        if (
            grepl(
                "\\.",
                five_mer
            )
        ) {
            next
        }

        if (
            !grepl(
                "^[ACGTN]{5}$",
                five_mer
            )
        ) {
            next
        }

        #######################################################################
        # The center nucleotide is the germline nucleotide being mutated.
        #
        # For an exact A/C/G/T center, use the four possible alternative
        # nucleotides only.
        #
        # For an N center, HH_S5F's N row can be used if present.
        #######################################################################

        center <- chars[p]

        if (
            center %in%
            nucleotide_names
        ) {

            lookup <- s5f[
                nucleotide_names,
                five_mer,
                drop = FALSE
            ]

            targeting[
                ,
                p
            ] <- lookup[, 1]

            # Remove the self-substitution.
            self_row <- match(
                center,
                nucleotide_names
            )

            targeting[
                self_row,
                p
            ] <- NA_real_

        } else if (
            center == "N"
        ) {

            ###################################################################
            # For an ambiguous germline center, the exact germline nucleotide
            # is unknown. We therefore retain the four A/C/G/T rows rather
            # than inventing a self-substitution exclusion.
            ###################################################################

            if (
                "N" %in%
                rownames(s5f)
            ) {

                lookup <- s5f[
                    nucleotide_names,
                    five_mer,
                    drop = FALSE
                ]

                targeting[
                    ,
                    p
                ] <- lookup[, 1]
            }
        }
    }

    targeting
}


###############################################################################
# 12. EVENT PATHS FOR ONE GERMLINE
#
# Output:
#
#   rows    = possible mutant nucleotide A/C/G/T
#   columns = IMGT nucleotide positions 1..312
#
# TRUE means:
#
#   "Changing the germline nucleotide at this position to this new
#    nucleotide produces this charge event."
#
# Only complete A/C/G/T germline codons are used to construct event paths.
###############################################################################

make_event_paths <- function(
    germline_seq,
    n_pos = 312
) {

    chars <- prepare_germline_for_targeting(
        germline_seq,
        n_pos = n_pos
    )

    paths <- lapply(
        event_names,
        function(event) {

            matrix(
                FALSE,
                nrow = 4,
                ncol = n_pos,
                dimnames = list(
                    nucleotide_names,
                    seq_len(n_pos)
                )
            )
        }
    )

    names(paths) <- event_names


    ###########################################################################
    # Evaluate every codon independently.
    ###########################################################################

    for (j in seq_len(n_codons)) {

        idx <- (
            3 * j - 2
        ):(3 * j)

        codon <- paste0(
            chars[idx],
            collapse = ""
        )

        #######################################################################
        # An event path requires an exact germline codon.
        #
        # N, gap, or missing positions therefore produce no event path.
        #######################################################################

        if (
            !grepl(
                "^[ACGT]{3}$",
                codon
            )
        ) {
            next
        }

        germline_aa <- unname(
            Biostrings::GENETIC_CODE[
                codon
            ]
        )

        from_class <- classify_charge(
            germline_aa
        )

        if (is.na(from_class)) {
            next
        }

        germline_nt <- strsplit(
            codon,
            "",
            fixed = TRUE
        )[[1]]


        #######################################################################
        # Test all three nucleotide positions in the codon.
        #######################################################################

        for (p in seq_len(3)) {

            position <- idx[p]

            original_nt <- germline_nt[p]

            for (new_nt in nucleotide_names) {

                # A nucleotide cannot mutate to itself.
                if (
                    new_nt ==
                    original_nt
                ) {
                    next
                }

                mutant_nt <- germline_nt

                mutant_nt[p] <- new_nt

                mutant_codon <- paste0(
                    mutant_nt,
                    collapse = ""
                )

                mutant_aa <- unname(
                    Biostrings::GENETIC_CODE[
                        mutant_codon
                    ]
                )

                to_class <- classify_charge(
                    mutant_aa
                )

                if (is.na(to_class)) {
                    next
                }

                for (event in event_names) {

                    paths[[event]][
                        new_nt,
                        position
                    ] <- isTRUE(
                        event_matches(
                            from_class,
                            to_class,
                            event
                        )
                    )
                }
            }
        }
    }

    paths
}


###############################################################################
# 13. EXPECTED TARGETING FOR ONE GERMLINE
#
# For every nucleotide position:
#
#   expected_event_targeting =
#       sum(
#           HH_S5F targeting for possible mutant nucleotide
#           x
#           whether that mutation belongs to the event
#       )
#
# Positions without a valid 5-mer contribute zero.
###############################################################################

calculate_expected_targeting <- function(
    germline_seq
) {

    targeting <- fast_targeting(
        germline_seq,
        n_pos = n_pos
    )

    paths <- make_event_paths(
        germline_seq,
        n_pos = n_pos
    )

    expected <- matrix(
        0,
        nrow = length(event_names),
        ncol = n_pos,
        dimnames = list(
            event_names,
            seq_len(n_pos)
        )
    )

    for (event in event_names) {

        event_path <- paths[[event]]

        #######################################################################
        # NA targeting values cannot contribute.
        #######################################################################

        product <- targeting *
            event_path

        product[
            !is.finite(product)
        ] <- NA_real_

        expected[event, ] <- colSums(
            product,
            na.rm = TRUE
        )
    }

    expected
}


###############################################################################
# 14. VALIDATE FAST TARGETING AGAINST SHAZAM
#
# Use a FULL-LENGTH germline for the primary validation because the goal here
# is to establish that the fast lookup reproduces SHazaM's HH_S5F targeting
# calculation.
#
# The sequence is passed as inputSeq as well, matching the original workflow.
###############################################################################

cat("\n")
cat("============================================\n")
cat("VALIDATING FAST HH_S5F TARGETING\n")
cat("============================================\n")

validation_i <- which(
    germline_original_length >= n_pos
)[1]

if (is.na(validation_i)) {

    stop(
        "No germline of at least ",
        n_pos,
        " nt is available for HH_S5F validation."
    )
}

validation_germline <- toupper(
    chartr(
        "-",
        ".",
        substr(
            as.character(
                myDb[[germline_col]][
                    validation_i
                ]
            ),
            1,
            n_pos
        )
    )
)

validation_sequence <- toupper(
    chartr(
        "-",
        ".",
        substr(
            as.character(
                myDb[[sequence_col]][
                    validation_i
                ]
            ),
            1,
            n_pos
        )
    )
)


###############################################################################
# SHazaM reference calculation
###############################################################################

slow_targeting <- shazam:::calculateTargeting(
    germlineSeq = validation_germline,
    inputSeq = validation_sequence,
    targetingModel = HKL_S5F,
    regionDefinition = IMGT_V
)


###############################################################################
# Our fast calculation
###############################################################################

fast_targeting_result <- fast_targeting(
    validation_germline,
    n_pos = n_pos
)


###############################################################################
# Compare only the A/C/G/T rows.
#
# Both matrices are 4 x 312 for these rows.
###############################################################################

comparison <- (
    slow_targeting[
        nucleotide_names,
        ,
        drop = FALSE
    ] -
    fast_targeting_result
)

finite_difference <- comparison[
    is.finite(comparison)
]

if (length(finite_difference) == 0) {

    stop(
        "No finite targeting values were available for validation."
    )
}

max_targeting_difference <- max(
    abs(
        finite_difference
    )
)

cat(
    "Validation sequence:",
    validation_i,
    "\n"
)

cat(
    "Validation germline length:",
    germline_original_length[
        validation_i
    ],
    "\n"
)

cat(
    "Maximum HH_S5F targeting difference:",
    max_targeting_difference,
    "\n"
)


if (
    !is.finite(
        max_targeting_difference
    ) ||
    max_targeting_difference > 1e-8
) {

    stop(
        "Fast HH_S5F targeting does not reproduce SHazaM targeting. ",
        "Maximum difference = ",
        max_targeting_difference
    )
}

cat(
    "Fast HH_S5F targeting validation: PASSED\n"
)


###############################################################################
# 15. VALIDATE EVENT PATHS
###############################################################################

cat("\n")
cat("============================================\n")
cat("VALIDATING EVENT PATHS\n")
cat("============================================\n")

validation_paths <- make_event_paths(
    validation_germline,
    n_pos = n_pos
)

path_counts <- vapply(
    validation_paths,
    function(x) {

        sum(
            x,
            na.rm = TRUE
        )
    },
    numeric(1)
)

print(
    path_counts
)

if (
    any(
        !is.finite(
            path_counts
        )
    )
) {

    stop(
        "Invalid values found in event-path matrices."
    )
}

cat(
    "Event-path validation: PASSED\n"
)


###############################################################################
# 16. UNIQUE-GERMLINE EXPECTED-TARGETING CACHE
#
# IMPORTANT:
#
# This does NOT collapse observations.
#
# If 100 rows have the same germline, the cached expected targeting is reused
# 100 times.
#
# The sequence-specific informative-codon mask is applied later, so a short
# germline/codon is excluded for the appropriate individual observations.
###############################################################################

cat("\n")
cat("============================================\n")
cat("BUILDING GERMLINE EXPECTED-TARGETING CACHE\n")
cat("============================================\n")

unique_germlines <- unique(
    germlines
)

n_unique_germlines <- length(
    unique_germlines
)

germline_index <- match(
    germlines,
    unique_germlines
)

cat(
    "Input sequences:",
    n_input,
    "\n"
)

cat(
    "Unique germlines:",
    n_unique_germlines,
    "\n"
)

cat(
    "Every input sequence remains an independent observation.",
    "\n"
)


###############################################################################
# Cache:
#
#   germline x event x nucleotide_position
###############################################################################

expected_targeting_cache <- array(
    0,
    dim = c(
        n_unique_germlines,
        length(event_names),
        n_pos
    ),
    dimnames = list(
        germline = seq_len(
            n_unique_germlines
        ),
        event = event_names,
        nucleotide_position = seq_len(
            n_pos
        )
    )
)


for (u in seq_len(n_unique_germlines)) {

    expected_targeting_cache[
        u,
        ,
    ] <- calculate_expected_targeting(
        unique_germlines[u]
    )

    if (
        u %% 500 == 0 ||
        u == n_unique_germlines
    ) {

        cat(
            "\rGermlines processed:",
            u,
            "/",
            n_unique_germlines,
            "   "
        )
    }
}

cat("\n")


###############################################################################
# 16A. EXPECTED TARGETING ACROSS THE SEQUENCE-LEVEL UNIVERSE
#
# For each sequence i and codon j:
#
#   if valid_codon[i,j] == TRUE:
#
#       add the cached germline targeting for sequence i.
#
#   otherwise:
#
#       add nothing.
#
# This is the critical step that ensures expected and observed calculations
# use exactly the same sequence/codon universe.
###############################################################################

cat("\n")
cat("============================================\n")
cat("CALCULATING SEQUENCE-LEVEL EXPECTED TARGETING\n")
cat("============================================\n")

expected_event_nt <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_pos,
    dimnames = list(
        event_names,
        seq_len(n_pos)
    )
)


for (j in seq_len(n_codons)) {

    informative_sequences <- which(
        valid_codon[, j]
    )

    if (
        length(
            informative_sequences
        ) == 0
    ) {
        next
    }

    for (i in informative_sequences) {

        u <- germline_index[i]

        expected_event_nt <-
            expected_event_nt +
            expected_targeting_cache[
                u,
                ,
                ,
                drop = TRUE
            ]
    }

    if (
        j %% 10 == 0 ||
        j == n_codons
    ) {

        cat(
            "\rCodons processed:",
            j,
            "/",
            n_codons,
            "   "
        )
    }
}

cat(
    "\nExpected targeting calculation complete.\n"
)


###############################################################################
# 16B. COLLAPSE EXPECTED TARGETING TO CODON LEVEL
###############################################################################

expected_event_codon <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

for (j in seq_len(n_codons)) {

    idx <- (
        3 * j - 2
    ):(3 * j)

    expected_event_codon[
        ,
        j
    ] <- rowSums(
        expected_event_nt[
            ,
            idx,
            drop = FALSE
        ]
    )
}


###############################################################################
# 16C. EXPECTED-TARGETING QC
###############################################################################

expected_QC <- data.frame(
    aa_position = seq_len(n_codons),
    coverage = coverage_codon
)

for (event in event_names) {

  expected_QC[[paste0(event, "_expected_targeting")]] <-
    expected_event_codon[event, ]

  expected_QC[[paste0(event, "_expected_zero")]] <-
    expected_event_codon[event, ] <= 0
}

write_tsv(
    expected_QC,
    "charge_OE_expected_targeting_QC.tsv"
)

cat(
    "Expected-targeting cache and sequence-level expected values: PASSED\n"
)




###############################################################################
# 17. EXPECTED TARGETING AT CODON LEVEL
#
# Sum the three nucleotide positions belonging to each codon.
#
# IMPORTANT:
#   Expected values are SUMMED, not averaged.
###############################################################################

expected_event_codon <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

for (j in seq_len(n_codons)) {

    idx <- (
        3 * j - 2
    ):(3 * j)

    expected_event_codon[
        ,
        j
    ] <- rowSums(
        expected_event_nt[
            ,
            idx,
            drop = FALSE
        ]
    )
}


###############################################################################
# 18. OBSERVED EVENT COUNTS
###############################################################################

cat("\n")
cat("============================================\n")
cat("CALCULATING OBSERVED EVENTS\n")
cat("============================================\n")

observed_event_nt <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_pos,
    dimnames = list(
        event_names,
        seq_len(n_pos)
    )
)

observed_event_codon <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

multi_event_codon <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

single_event_totals <- setNames(
    numeric(
        length(event_names)
    ),
    event_names
)

multi_event_totals <- setNames(
    numeric(
        length(event_names)
    ),
    event_names
)


###############################################################################
# 19. OBSERVED EVENTS, CODON BY CODON
#
# Every sequence is evaluated independently.
###############################################################################

for (j in seq_len(n_codons)) {

    idx <- (
        3 * j - 2
    ):(3 * j)

    informative <- valid_codon[, j]

    if (!any(informative)) {
        next
    }

    g_codons <- substring(
        germlines[informative],
        idx[1],
        idx[3]
    )

    s_codons <- substring(
        sequences[informative],
        idx[1],
        idx[3]
    )

    g_split <- strsplit(
        g_codons,
        "",
        fixed = TRUE
    )

    s_split <- strsplit(
        s_codons,
        "",
        fixed = TRUE
    )

    g1 <- vapply(
        g_split,
        `[`,
        character(1),
        1
    )

    g2 <- vapply(
        g_split,
        `[`,
        character(1),
        2
    )

    g3 <- vapply(
        g_split,
        `[`,
        character(1),
        3
    )

    s1 <- vapply(
        s_split,
        `[`,
        character(1),
        1
    )

    s2 <- vapply(
        s_split,
        `[`,
        character(1),
        2
    )

    s3 <- vapply(
        s_split,
        `[`,
        character(1),
        3
    )

    difference_1 <- g1 != s1
    difference_2 <- g2 != s2
    difference_3 <- g3 != s3

    n_nt_changes <- (
        difference_1 +
        difference_2 +
        difference_3
    )


    ###########################################################################
    # SINGLE-NUCLEOTIDE EVENTS
    ###########################################################################

    single_nt <- (
        n_nt_changes == 1
    )

    if (any(single_nt)) {

        g_single <- g_codons[
            single_nt
        ]

        s_single <- s_codons[
            single_nt
        ]

        germline_aa <- unname(
            Biostrings::GENETIC_CODE[
                g_single
            ]
        )

        observed_aa <- unname(
            Biostrings::GENETIC_CODE[
                s_single
            ]
        )

        from_class <- classify_charge(
            germline_aa
        )

        to_class <- classify_charge(
            observed_aa
        )

        single_indices <- which(
            single_nt
        )

        for (event in event_names) {

            hit <- event_matches(
                from_class,
                to_class,
                event
            )

            hit[is.na(hit)] <- FALSE

            if (!any(hit)) {
                next
            }

            event_indices <- single_indices[
                hit
            ]

            ###################################################################
            # Count the event at the exact nucleotide position.
            ###################################################################

            if (
                any(
                    difference_1[
                        event_indices
                    ]
                )
            ) {

                observed_event_nt[
                    event,
                    idx[1]
                ] <-
                    observed_event_nt[
                        event,
                        idx[1]
                    ] +
                    sum(
                        difference_1[
                            event_indices
                        ]
                    )
            }

            if (
                any(
                    difference_2[
                        event_indices
                    ]
                )
            ) {

                observed_event_nt[
                    event,
                    idx[2]
                ] <-
                    observed_event_nt[
                        event,
                        idx[2]
                    ] +
                    sum(
                        difference_2[
                            event_indices
                        ]
                    )
            }

            if (
                any(
                    difference_3[
                        event_indices
                    ]
                )
            ) {

                observed_event_nt[
                    event,
                    idx[3]
                ] <-
                    observed_event_nt[
                        event,
                        idx[3]
                    ] +
                    sum(
                        difference_3[
                            event_indices
                        ]
                    )
            }

            n_events <- sum(
                hit
            )

            observed_event_codon[
                event,
                j
            ] <-
                observed_event_codon[
                    event,
                    j
                ] +
                n_events

            single_event_totals[event] <-
                single_event_totals[event] +
                n_events
        }
    }


    ###########################################################################
    # MULTI-NUCLEOTIDE EVENTS
    #
    # QC ONLY.
    #
    # These are classified by their resulting amino-acid change, but are NOT
    # included in the primary observed O/E calculation.
    ###########################################################################

    multi_nt <- (
        n_nt_changes > 1
    )

    if (any(multi_nt)) {

        g_multi <- g_codons[
            multi_nt
        ]

        s_multi <- s_codons[
            multi_nt
        ]

        germline_aa_multi <- unname(
            Biostrings::GENETIC_CODE[
                g_multi
            ]
        )

        observed_aa_multi <- unname(
            Biostrings::GENETIC_CODE[
                s_multi
            ]
        )

        from_multi <- classify_charge(
            germline_aa_multi
        )

        to_multi <- classify_charge(
            observed_aa_multi
        )

        for (event in event_names) {

            hit_multi <- event_matches(
                from_multi,
                to_multi,
                event
            )

            hit_multi[is.na(hit_multi)] <- FALSE

            n_multi_events <- sum(
                hit_multi
            )

            if (n_multi_events == 0) {
                next
            }

            multi_event_codon[
                event,
                j
            ] <-
                multi_event_codon[
                    event,
                    j
                ] +
                n_multi_events

            multi_event_totals[event] <-
                multi_event_totals[event] +
                n_multi_events
        }
    }

    if (
        j %% 10 == 0 ||
        j == n_codons
    ) {

        cat(
            "\rCodons processed:",
            j,
            "/",
            n_codons,
            "   "
        )
    }
}

cat("\nObserved event calculation complete.\n")


###############################################################################
# 20. INTERNAL OBSERVED-COUNT CHECK
###############################################################################
#
# Every primary observed event is exactly one nucleotide event.
#
# Therefore:
#
#   sum(observed_event_nt across 312 positions)
#     ==
#   sum(observed_event_codon across 104 codons)
#
###############################################################################

observed_event_codon_from_nt <- matrix(
    0,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

for (j in seq_len(n_codons)) {

    idx <- (
        3 * j - 2
    ):(3 * j)

    observed_event_codon_from_nt[
        ,
        j
    ] <- rowSums(
        observed_event_nt[
            ,
            idx,
            drop = FALSE
        ]
    )
}

if (
    !isTRUE(
        all.equal(
            observed_event_codon,
            observed_event_codon_from_nt,
            tolerance = 0
        )
    )
) {

    stop(
        "Internal observed-event aggregation mismatch."
    )
}


###############################################################################
# 21. EVENT QC
###############################################################################

cat("\n")
cat("============================================\n")
cat("EVENT QC\n")
cat("============================================\n")

event_qc <- data.frame(
    event = event_names,
    single_nt = as.numeric(
        single_event_totals
    ),
    multi_nt = as.numeric(
        multi_event_totals
    )
)

event_qc$total_events <- (
    event_qc$single_nt +
    event_qc$multi_nt
)

event_qc$multi_fraction <- ifelse(
    event_qc$total_events > 0,
    event_qc$multi_nt /
        event_qc$total_events,
    NA_real_
)

print(
    event_qc
)

write_tsv(
    event_qc,
    "charge_OE_event_QC.tsv"
)


###############################################################################
# 22. OBSERVED AND EXPECTED TOTALS
###############################################################################

cat("\n")
cat("============================================\n")
cat("OBSERVED / EXPECTED TOTALS\n")
cat("============================================\n")

for (event in event_names) {

    cat(
        event,
        "\n",
        "  observed:",
        sum(
            observed_event_codon[
                event,
                ],
            na.rm = TRUE
        ),
        "\n",
        "  expected targeting:",
        sum(
            expected_event_codon[
                event,
                ],
            na.rm = TRUE
        ),
        "\n"
    )
}


###############################################################################
# 23. POSITIONAL OBSERVED / EXPECTED FREQUENCIES
#
# O/E is a positional enrichment:
#
#       O_j / sum(O)
# O/E = -------------
#       E_j / sum(E)
#
# This is NOT:
#
#       observed_count / expected_count
#
# because the goal is enrichment of the POSITIONAL DISTRIBUTION of events.
###############################################################################

observed_frequency_event <- matrix(
    NA_real_,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

expected_frequency_event <- matrix(
    NA_real_,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)

OE_event <- matrix(
    NA_real_,
    nrow = length(event_names),
    ncol = n_codons,
    dimnames = list(
        event_names,
        seq_len(n_codons)
    )
)


for (event in event_names) {

    observed <- observed_event_codon[
        event,
        ]

    expected <- expected_event_codon[
        event,
        ]

    total_observed <- sum(
        observed,
        na.rm = TRUE
    )

    total_expected <- sum(
        expected,
        na.rm = TRUE
    )

    if (total_observed > 0) {

        observed_frequency_event[
            event,
            ] <- observed /
            total_observed
    }

    if (total_expected > 0) {

        expected_frequency_event[
            event,
            ] <- expected /
            total_expected
    }

    valid <- (
        is.finite(
            observed_frequency_event[
                event,
                ]
        ) &
        is.finite(
            expected_frequency_event[
                event,
                ]
        ) &
        expected_frequency_event[
            event,
            ] > 0
    )

    OE_event[
        event,
        valid
    ] <-
        observed_frequency_event[
            event,
            valid
        ] /
        expected_frequency_event[
            event,
            valid
        ]
}


###############################################################################
# 24. REGION LABELS
###############################################################################

region_nt <- as.character(
    IMGT_V@boundaries
)

if (
    length(region_nt) !=
    n_pos
) {

    region_nt <- rep(
        NA_character_,
        n_pos
    )
}

region_codon <- vapply(
    seq_len(n_codons),
    function(j) {

        idx <- (
            3 * j - 2
        ):(3 * j)

        values <- unique(
            region_nt[idx]
        )

        values <- values[
            !is.na(values) &
            values != ""
        ]

        if (length(values) == 0) {
            return(NA_character_)
        }

        paste(
            values,
            collapse = "/"
        )
    },
    character(1)
)


###############################################################################
# 25. COVERAGE / MISSINGNESS QC
###############################################################################

codon_missing_sequence <- (
    n_input -
    rowSums(
        matrix(
            valid_codon,
            nrow = n_input,
            ncol = n_codons
        )
    )
)

# More useful per-codon summary:
codon_qc <- data.frame(
    aa_position = seq_len(n_codons),
    region = region_codon,
    informative_sequences = coverage_codon,
    noninformative_sequences =
        n_input - coverage_codon,
    informative_fraction =
        coverage_codon / n_input
)

write_tsv(
    codon_qc,
    "charge_OE_codon_coverage_QC.tsv"
)


###############################################################################
# 26. FINAL 104-CODON RESULT TABLE
#
# Naming convention:
#
#   <event>_observed
#       observed number of primary single-NT events
#
#   <event>_expected_targeting
#       summed HH_S5F expected targeting for the same sequence/codon universe
#
#   <event>_observed_frequency
#       observed positional frequency
#
#   <event>_expected_frequency
#       expected positional frequency
#
#   <event>_OE
#       observed_frequency / expected_frequency
#
#   <event>_multi_nt
#       number of multi-NT QC events
#
#   <event>_expected_zero
#       TRUE when expected targeting at this codon is zero
###############################################################################
result <- data.frame(
  aa_position = seq_len(n_codons),
  region = region_codon,
  coverage = coverage_codon
)

for (event in event_names) {

  result[[paste0(event, "_observed")]] <-
    observed_event_codon[event, ]

  result[[paste0(event, "_expected_targeting")]] <-
    expected_event_codon[event, ]

  result[[paste0(event, "_observed_frequency")]] <-
    observed_frequency_event[event, ]

  result[[paste0(event, "_expected_frequency")]] <-
    expected_frequency_event[event, ]

  result[[paste0(event, "_OE")]] <-
    OE_event[event, ]

  result[[paste0(event, "_multi_nt")]] <-
    multi_event_codon[event, ]

  result[[paste0(event, "_expected_zero")]] <-
    expected_event_codon[event, ] == 0
}

###############################################################################
# 27. O/E QC FILTER
#
# Purpose:
#   Flag positions where the O/E ratio is unstable because either:
#
#   1. expected targeting is extremely low, or
#   2. observed frequency is extremely low.
#
# The thresholds are calculated separately for each event.
#
# IMPORTANT:
#   This does NOT modify the primary O/E calculation.
#
#   It only creates QC flags.
###############################################################################

calculate_OE_QC <- function(
    observed,
    expected,
    OE
) {

    total_observed <- sum(
        observed,
        na.rm = TRUE
    )

    observed_frequency <- if (
        total_observed > 0
    ) {

        observed /
            total_observed

    } else {

        rep(
            NA_real_,
            length(observed)
        )
    }

    ###########################################################################
    # Expected-targeting threshold
    ###########################################################################

    expected_values <- expected[
        is.finite(expected)
    ]

    if (length(expected_values) > 0) {

        expected_Q1 <- quantile(
            expected_values,
            probs = 0.25,
            names = FALSE
        )

        expected_Q3 <- quantile(
            expected_values,
            probs = 0.75,
            names = FALSE
        )

        expected_IQR <- (
            expected_Q3 -
            expected_Q1
        )

        expected_lower_fence <- (
            expected_Q1 -
            1.5 * expected_IQR
        )

        # Expected targeting cannot be negative.
        expected_threshold <- max(
            0,
            expected_lower_fence
        )

    } else {

        expected_threshold <- NA_real_
    }


    ###########################################################################
    # Observed-frequency threshold
    ###########################################################################

    observed_values <- observed_frequency[
        is.finite(
            observed_frequency
        )
    ]

    if (length(observed_values) > 0) {

        observed_Q1 <- quantile(
            observed_values,
            probs = 0.25,
            names = FALSE
        )

        observed_Q3 <- quantile(
            observed_values,
            probs = 0.75,
            names = FALSE
        )

        observed_IQR <- (
            observed_Q3 -
            observed_Q1
        )

        observed_lower_fence <- (
            observed_Q1 -
            1.5 * observed_IQR
        )

        observed_threshold <- max(
            0,
            observed_lower_fence
        )

    } else {

        observed_threshold <- NA_real_
    }


    ###########################################################################
    # Flags
    ###########################################################################

    low_expected <- (
        is.finite(expected) &
        expected <= expected_threshold
    )

    low_observed <- (
        is.finite(observed_frequency) &
        observed_frequency <= observed_threshold
    )

    expected_zero <- (
        !is.finite(expected) |
        expected <= 0
    )

    OE_unstable <- (
        expected_zero |
        low_expected |
        low_observed
    )

    list(
        observed_frequency = observed_frequency,
        expected_threshold = expected_threshold,
        observed_threshold = observed_threshold,
        low_expected = low_expected,
        low_observed = low_observed,
        expected_zero = expected_zero,
        OE_unstable = OE_unstable
    )
}


###############################################################################
# 28. ADD O/E QC COLUMNS
###############################################################################
for (event in event_names) {

  qc <- calculate_OE_QC(
    observed = observed_event_codon[event, ],
    expected = expected_event_codon[event, ],
    OE = OE_event[event, ]
  )

  result[[paste0(event, "_low_expected")]] <-
    qc$low_expected

  result[[paste0(event, "_low_observed")]] <-
    qc$low_observed

  result[[paste0(event, "_OE_unstable")]] <-
    qc$OE_unstable
}

###############################################################################
# 29. PRIMARY O/E OUTPUT
#
# This is the clean 104-position table.
###############################################################################

write_tsv(
    result,
    "charge_OE_per_sequence_104_codons.tsv"
)


###############################################################################
# 30. OPTIONAL FILTERED O/E TABLE
#
# Here we retain the O/E value only when:
#
#   - expected targeting > 0
#   - expected targeting is not flagged low
#   - observed frequency is not flagged low
#
# The raw result table above remains untouched.
###############################################################################

filtered_result <- result

for (event in event_names) {

    OE_column <- paste0(
        event,
        "_OE"
    )

    unstable_column <- paste0(
        event,
        "_OE_unstable"
    )

    filtered_result[
        ,
        OE_column
    ][
        filtered_result[
            ,
            unstable_column
        ]
    ] <- NA_real_
}


write_tsv(
    filtered_result,
    "charge_OE_per_sequence_104_codons_QC_filtered.tsv"
)


###############################################################################
# 31. EVENT SUMMARY
###############################################################################

event_summary <- data.frame(
    event = event_names,
    observed_primary_events = vapply(
        event_names,
        function(event) {

            sum(
                observed_event_codon[
                    event,
                    ],
                na.rm = TRUE
            )
        },
        numeric(1)
    ),
    expected_targeting = vapply(
        event_names,
        function(event) {

            sum(
                expected_event_codon[
                    event,
                    ],
                na.rm = TRUE
            )
        },
        numeric(1)
    ),
    multi_nt_events = vapply(
        event_names,
        function(event) {

            sum(
                multi_event_codon[
                    event,
                    ],
                na.rm = TRUE
            )
        },
        numeric(1)
    )
)

event_summary$observed_fraction <- ifelse(
    event_summary$observed_primary_events > 0,
    event_summary$observed_primary_events /
        sum(
            event_summary$observed_primary_events
        ),
    NA_real_
)

event_summary$expected_fraction <- ifelse(
    event_summary$expected_targeting > 0,
    event_summary$expected_targeting /
        sum(
            event_summary$expected_targeting
        ),
    NA_real_
)

write_tsv(
    event_summary,
    "charge_OE_event_summary.tsv"
)


###############################################################################
# 32. FINAL QC SUMMARY
###############################################################################

cat("\n")
cat("============================================\n")
cat("FINAL QC SUMMARY\n")
cat("============================================\n")

cat(
    "Input sequences:",
    n_input,
    "\n"
)

cat(
    "Unique germlines:",
    n_unique_germlines,
    "\n"
)

cat(
    "312-nt analysis positions:",
    n_pos,
    "\n"
)

cat(
    "104 codons:",
    n_codons,
    "\n"
)

cat(
    "Informative sequence-codon observations:",
    sum(
        valid_codon
    ),
    "\n"
)

cat(
    "Non-informative sequence-codon observations:",
    sum(
        !valid_codon
    ),
    "\n"
)

cat("\n")

print(
    event_summary
)

cat("\n")
cat("Output files:\n")
cat(
    "  charge_OE_event_QC.tsv\n"
)
cat(
    "  charge_OE_codon_coverage_QC.tsv\n"
)
cat(
    "  charge_OE_per_sequence_104_codons.tsv\n"
)
cat(
    "  charge_OE_per_sequence_104_codons_QC_filtered.tsv\n"
)
cat(
    "  charge_OE_event_summary.tsv\n"
)

cat("\n")
cat("Analysis complete.\n")
