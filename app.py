
import streamlit as st
import tempfile
import subprocess
import os
import pandas as pd

st.set_page_config(
    page_title="Bacterial Identification Pipeline",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Bacterial Identification Pipeline")
st.subheader("16S rRNA Gene Sequence Analysis")

st.write(
    "Upload a bacterial 16S rRNA FASTA sequence to analyze "
    "its nucleotide composition and identify the closest "
    "matching organism."
)

st.divider()


# ==============================
# FASTA READER
# ==============================

def read_fasta(filepath):

    with open(filepath, "r") as file:
        data = file.read()

    lines = data.strip().splitlines()

    if not lines:
        return None, None

    header = lines[0]

    if not header.startswith(">"):
        return None, None

    sequence = "".join(lines[1:]).upper()

    return header, sequence


# ==============================
# DNA ANALYSIS
# ==============================

def analyze_dna(sequence, sequence_id="Unknown"):

    sequence = sequence.upper().replace(" ", "").replace("\n", "")

    valid_bases = set("ATGCN")
    invalid_bases = set(sequence) - valid_bases

    if not sequence:
        return {"Error": "DNA sequence is empty."}

    if invalid_bases:
        return {
            "Error": "Invalid DNA characters found.",
            "Invalid characters": str(invalid_bases)
        }

    length = len(sequence)

    A = sequence.count("A")
    T = sequence.count("T")
    G = sequence.count("G")
    C = sequence.count("C")
    N = sequence.count("N")

    known_bases = A + T + G + C

    if known_bases == 0:
        return {"Error": "No valid A/T/G/C bases found."}

    gc_percentage = ((G + C) / known_bases) * 100

    return {
        "Sequence ID": sequence_id,
        "Sequence Length (bp)": length,
        "A Count": A,
        "T Count": T,
        "G Count": G,
        "C Count": C,
        "N Count": N,
        "GC Percentage": round(gc_percentage, 2)
    }


# ==============================
# BLAST FUNCTION
# ==============================

def run_blast(fasta_file):

    db_path = "/content/blastdb/16S_ribosomal_RNA"
    output_file = "/tmp/blast_results.tsv"

    if not os.path.exists(db_path + ".nhr"):
        return None, "BLAST database not found."

    command = [
        "blastn",
        "-query", fasta_file,
        "-db", db_path,
        "-task", "megablast",
        "-max_target_seqs", "10",
        "-outfmt",
        "6 qseqid sacc pident length qcovs evalue bitscore stitle",
        "-out", output_file
    ]

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=120
        )

    except subprocess.TimeoutExpired:

        return None, "BLAST search timed out after 120 seconds."

    if result.returncode != 0:

        return None, result.stderr

    if not os.path.exists(output_file):

        return None, "BLAST output file was not created."

    if os.path.getsize(output_file) == 0:

        return None, "No BLAST matches were found."

    columns = [
        "Query_ID",
        "Accession",
        "Percent_Identity",
        "Alignment_Length",
        "Query_Coverage",
        "E_Value",
        "Bit_Score",
        "Description"
    ]

    blast_df = pd.read_csv(
        output_file,
        sep="\t",
        header=None,
        names=columns
    )

    return blast_df, None


# ==============================
# ORGANISM IDENTIFICATION
# ==============================

def identify_organism(blast_df):

    if blast_df is None or blast_df.empty:
        return None

    best_hit = blast_df.iloc[0]

    description = str(best_hit["Description"])

    words = description.split()

    predicted_organism = "Unknown"

    if len(words) >= 2:

        if (
            words[0][0].isupper()
            and words[1][0].islower()
        ):
            predicted_organism = (
                words[0] + " " + words[1]
            )

    return {
        "Predicted Organism": predicted_organism,
        "Percent Identity": float(
            best_hit["Percent_Identity"]
        ),
        "Query Coverage": float(
            best_hit["Query_Coverage"]
        ),
        "E-value": best_hit["E_Value"],
        "Accession": best_hit["Accession"]
    }


# ==============================
# USER INTERFACE
# ==============================

st.header("📁 Upload FASTA File")

uploaded_file = st.file_uploader(
    "Choose a FASTA file",
    type=["fasta", "fa", "fna"]
)


if uploaded_file is not None:

    st.success("✅ FASTA file uploaded successfully.")

    with tempfile.NamedTemporaryFile(
        delete=False,
        suffix=".fasta"
    ) as temp_file:

        temp_file.write(
            uploaded_file.getvalue()
        )

        fasta_path = temp_file.name


    header, sequence = read_fasta(fasta_path)


    if header is None or sequence is None:

        st.error("❌ Invalid or empty FASTA file.")


    else:

        results = analyze_dna(
            sequence,
            header
        )


        if "Error" in results:

            st.error(
                "❌ DNA sequence validation failed."
            )

            for key, value in results.items():

                st.write(
                    f"**{key}:** {value}"
                )


        else:

            # ==============================
            # SEQUENCE INFORMATION
            # ==============================

            st.header("🧬 Sequence Information")

            col1, col2 = st.columns(2)

            with col1:

                st.metric(
                    "Sequence Length",
                    f"{results['Sequence Length (bp)']} bp"
                )

            with col2:

                st.metric(
                    "GC Percentage",
                    f"{results['GC Percentage']}%"
                )


            # ==============================
            # NUCLEOTIDE COMPOSITION
            # ==============================

            st.header("🧪 Nucleotide Composition")

            col1, col2, col3, col4, col5 = st.columns(5)

            with col1:
                st.metric("A", results["A Count"])

            with col2:
                st.metric("T", results["T Count"])

            with col3:
                st.metric("G", results["G Count"])

            with col4:
                st.metric("C", results["C Count"])

            with col5:
                st.metric("N", results["N Count"])


            # ==============================
            # COMPOSITION PLOT
            # ==============================

            st.header("📊 Nucleotide Composition Plot")

            chart_data = {
                "A": results["A Count"],
                "T": results["T Count"],
                "G": results["G Count"],
                "C": results["C Count"],
                "N": results["N Count"]
            }

            st.bar_chart(chart_data)


            # ==============================
            # BLAST ANALYSIS
            # ==============================

            st.header("🔬 BLAST Analysis")

            st.info(
                "Searching the NCBI 16S ribosomal RNA "
                "reference database..."
            )

            blast_results, blast_error = run_blast(
                fasta_path
            )


            if blast_error:

                st.error(
                    "❌ BLAST search failed."
                )

                st.code(blast_error)


            else:

                st.success(
                    "✅ BLAST search completed successfully."
                )


                # ==============================
                # IDENTIFICATION
                # ==============================

                identification = identify_organism(
                    blast_results
                )


                if identification:

                    st.header(
                        "🧬 Predicted Organism"
                    )

                    st.success(
                        f"**{identification['Predicted Organism']}**"
                    )


                    col1, col2, col3 = st.columns(3)

                    with col1:

                        st.metric(
                            "Identity",
                            f"{identification['Percent Identity']:.2f}%"
                        )

                    with col2:

                        st.metric(
                            "Query Coverage",
                            f"{identification['Query Coverage']:.2f}%"
                        )

                    with col3:

                        st.metric(
                            "E-value",
                            identification["E-value"]
                        )


                    st.write(
                        "**Best-hit accession:**",
                        identification["Accession"]
                    )


                # ==============================
                # TOP BLAST MATCHES
                # ==============================

                st.header(
                    "🔎 Top 10 BLAST Matches"
                )

                display_df = blast_results[
                    [
                        "Accession",
                        "Percent_Identity",
                        "Query_Coverage",
                        "E_Value",
                        "Alignment_Length",
                        "Description"
                    ]
                ].copy()

                display_df.columns = [
                    "Accession",
                    "Identity (%)",
                    "Coverage (%)",
                    "E-value",
                    "Alignment Length (bp)",
                    "Description"
                ]

                st.dataframe(
                    display_df,
                    use_container_width=True,
                    hide_index=True
                )


            # ==============================
            # COMPLETION
            # ==============================

            st.success(
                "✅ Pipeline analysis completed."
            )
