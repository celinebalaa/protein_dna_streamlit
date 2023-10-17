"""Website app built using the streamlit library.

It translates a given DNA sequence into a protein
and color codes the amino acids based on their properties.
"""

import streamlit as st

frames = ["Frame 1", "Frame 2", "Frame 3"]
frames_int = {"Frame 1": 0, "Frame 2": 1, "Frame 3": 2}
# Frame will be changed to the according position of the first nucleotide: int.

genetic_code = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        # _ is the stop codon
}

allowed_bases = ["A", "C", "G", "T"]
# Only the codon containing these four nucleotides will be translated.

hydrophobic = ["G", "A", "V", "L", "I", "P", "M"]
polar = ["S", "T", "N", "Q", "Y", "F", "C", "W"]
basic = ["K", "R", "H"]
acid = ["D", "E"]
# Four lists dividing the amino acids based on their properties.


def read_fasta(dna_sequence_fasta, sequence_type_five):
    """Take a fasta sequence inputted by the user and its type and format it.

    Parameters
    ----------
    dna_sequence_fasta : str
        Fasta sequence inputted by the user for translation.
    sequence_type : str
        Type of the sequence inputted: RNA or DNA.

    Returns
    -------
    dict
        Dictionary containing the sequence accession number as keys
        and the corresponding sequences as values.
    """
    dna_sequences = {}
    dna_sequence_fields = dna_sequence_fasta.split("\n")
    dna_id = ""
    for line in dna_sequence_fields:
        if line.startswith(">"):
            dna_id = line[1:].split()[0]
            dna_sequences[dna_id] = ""
            # The key will be the sequence's identifier.
        else:
            if sequence_type_five == "RNA":
                line = line.replace("U", "T")
            line = line.upper().strip().replace(" ", "")
            dna_sequences[dna_id] += line
    return dna_sequences


def read_fasta_three(dna_sequence_fasta, sequence_type_three):
    """Take a fasta sequence and its type inputted and change its direction.

    Parameters
    ----------
    dna_sequence_fasta : string
        Fasta sequence inputted by the user for translation.
    sequence_type : string
        Type of the sequence inputted: RNA or DNA.

    Returns
    -------
    dictionnary
        Dictionary containing the sequence accession number as keys
        and the corresponding sequences as values.
    """
    dna_sequences_three = {}
    dna_sequence_fields = dna_sequence_fasta.split("\n")
    dna_id = ""
    for line in dna_sequence_fields:
        if line.startswith(">"):
            dna_id = line[1:].split()[0]
            dna_sequences_three[dna_id] = ""
        else:
            if sequence_type_three == "RNA":
                line = line.replace("U", "T")
            line = line.upper().strip().replace(" ", "")
            line = line[::-1]
            # Reverse the sequence hence translation will be from 3' to 5'.
            dna_sequences_three[dna_id] += line
    return dna_sequences_three


def dna_protein(dna_seq, reading_frame):
    """Translate a dna sequence into a protein.

    Three reading frame can be chosen (1,2,3)
    If one is chosen the first nucleotide will be at position 0.

    Parameters
    ----------
    dna_seq : str
        DNA sequence in a fasta format.
    reading_frame : int
        Reading frame chosen by the user.

    Returns
    -------
    string
        Protein sequence, one letter, fasta format.
    """
    protein_sequence = ""
    for i in range(reading_frame, len(dna_seq)-reading_frame, 3):
        # Sequence will be read from the nucleotide specified to its end.
        # Three nucleotides aka the present codon will be skiped each time.
        codon = dna_seq[i:i+3]
        # A codon is formed by three successive nucleotides.
        not_base = False
        if len(codon) != 3:
            break
        # The codon length has to be equal to three.
        for base in codon:
            # Check if the codon only has A,C,T or G.
            if base not in allowed_bases:
                not_base = True
        if not_base:
            continue
        protein_sequence += genetic_code.get(codon) + ''
        # Gets the amino acid corresponding to the given codon.
    return protein_sequence


def translate_sequence(dna_dict):
    """Translate multiple fasta dna sequences into protein sequences.

    Parameters
    ----------
    dna_dict : dict
        Dictionnary containing dna sequences and their identifiers.

    Returns
    -------
    dict
        Dictionnary containing the accession number as keys and
        the corresponding protein sequences as values.
    """
    protein_dict = {}
    for sequence_id_prot in dna_dict.keys():
        # Each sequence in the dictionnary will be translated into a protein.
        protein_dict[sequence_id_prot] = dna_protein(
            dna_dict[sequence_id_prot], frame_output)
        # The translated sequence will be stored along with its identifier.
    return protein_dict


def aminoacids_properties(aminoacid):
    """Return a color depending on the amino acid given.

    The color code is a follow: yellow for hydrophobic,
    green for non-charged polar, blue for basic and red for acid.

    Parameters
    ----------
    aminoacid : str
        Amino acid from a protein sequence.

    Returns
    -------
    string
        Amino acid name colored depending on its properties.
    """
    color_aa = "black"
    # For stop codon it will be a _ in black
    if aminoacid in hydrophobic:
        color_aa = "#ffba00"
    elif aminoacid in polar:
        color_aa = "#279f27"
    elif aminoacid in basic:
        color_aa = "red"
    elif aminoacid in acid:
        color_aa = "#005bff"
    # Format returned is in css and html.
    return f"<font color='{color_aa}'>{aminoacid}"


st.title("Translation app using streamlit")
st.header("DNA to Protein Translation")
st.markdown("Input a DNA sequence in a **FASTA** or ** Multi-FASTA** format.")

sequence_type = st.radio("Choose the type of the sequence", ["DNA", "RNA"])
# Sequence can be read 5'3' (default) or 3'5'.
five_orientation = st.checkbox("Forward (5'3')", True)
three_orientation = st.checkbox("Reverse (3'5')", False)

frame_input = st.selectbox("Select a reading frame.", frames)
frame_output = frames_int.get(frame_input)

dna_sequence_input = st.text_area("DNA sequence: ", key="DNA Sequence")

st.markdown("You can upload a Multi-FASTA file,"
            "for an example click the button below!")
if st.button("Show multifasta example"):
    # st.code was used in order to have the copy button.
    st.code(""">NC_001802.1:336-4642 Human immunodeficiency virus 1
ATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAG
GGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAA
TCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACA
GGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGA

>NC_045512.2 Severe acute respiratory syndrome coronavirus 2
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA
CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC
TAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTG
TTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTC""")

if dna_sequence_input != "":
    # Will prevent the program from running if no sequence is inputted.

    st.success("The protein sequence generated is in one letter FASTA format.")
    st.header("Protein Sequence")

    if five_orientation:
        st.markdown("### Orientation 5'3' ")
        dna_sequence_output_five = read_fasta(
            dna_sequence_input, sequence_type)
        # Formated sequences are returned and ready for translation.
        protein_sequence_five = translate_sequence(dna_sequence_output_five)
        # Sequence are translated and stored in a dictionary.

        for sequence_id, sequence_adn in protein_sequence_five.items():
            # Protein sequences are read one by one.
            amino_acid_str = ""
            for amino_acid in sequence_adn:
                # Check the properties of each amino acid and color code it.
                amino_acid_str += aminoacids_properties(amino_acid)
            st.write(f"\>{sequence_id} \n <p>{amino_acid_str}</p>",
                     unsafe_allow_html=True)
            # Sequence's identifier and the colored amino acids are displayed.

    if three_orientation:
        # Same process but with a reversed orientation.
        st.markdown("### Orientation 3'5' ")
        dna_sequence_output_three = read_fasta_three(
            dna_sequence_input, sequence_type)
        protein_sequence_three = translate_sequence(dna_sequence_output_three)

        for sequence_id, sequence_adn in protein_sequence_three.items():
            amino_acid_str = ""
            for amino_acid in sequence_adn:
                amino_acid_str += aminoacids_properties(amino_acid)
            st.write(f"\>{sequence_id} \n <p>{amino_acid_str}</p>",
                     unsafe_allow_html=True)
