# VCF-Variant-Annotation-and-Analysis
Analyzes VCF variants using GFF and FASTA to classify mutation impact as synonymous, non-synonymous, or non-coding. Developed as coursework for MSc Bioinformatics 2024/25, University of Glasgow.


### **Overview**

This Python script analyzes genomic variants from a VCF file by integrating annotation data from GFF and FASTA files. It identifies whether each variant lies within a coding region and evaluates its impact at the protein level by translating affected codons. The tool classifies mutations as synonymous, non-synonymous, or non-coding, and outputs both a detailed TSV file and a summary bar plot of mutation types.


### **Features**

* Parses and filters variants from a VCF file based on quality.
* Identifies coding regions using gene annotations from a GFF file.
* Extracts reference sequences from a FASTA file.
* Translates codons to determine amino acid changes.
* Classifies each variant as:

  * Synonymous
  * Non-synonymous
  * Non-coding
* Outputs a table of annotated variants.
* Generates a bar plot showing proportions of mutation types.

---

### **Requirements**

* Python 3.6 or higher

Install dependencies using:

```bash
pip install PyVCF gffutils biopython pandas matplotlib seaborn
```

---

### **Usage**

```bash
python script.py \
  --vcf path/to/input.vcf \
  --gff path/to/annotations.gff \
  --fasta path/to/genome.fasta \
  --log path/to/log.txt \
  --output path/to/output.tsv
```

**Arguments:**

| Argument   | Description                             |
| ---------- | --------------------------------------- |
| `--vcf`    | Input VCF file with variant calls       |
| `--gff`    | GFF file containing gene annotations    |
| `--fasta`  | FASTA file of the reference genome      |
| `--log`    | Log file for runtime messages/errors    |
| `--output` | Output TSV file with annotated variants |

---

### **Output**

1. **Annotated TSV file**
   Each row includes:

   * Chromosome
   * Position
   * Reference and Alternate Allele
   * Mutation Type
   * Transcript ID (if applicable)
   * Protein Location (1-based)
   * Reference and Alternate Amino Acid

2. **Mutation Type Bar Plot**
   Saves a plot as `mut_prop_barplot.png` showing the proportion of each mutation type.

**Example TSV Output:**

```
Chrom  Pos     Ref  Alt  Type            Transcript    Protein Location  Ref AA  Alt AA
chr1   10523   C    T    Synonymous      AT1G12345.1    117               L       NA
chr2   34210   A    G    Non-synonymous  AT2G67890.2    98                D       G
chr3   789000  T    C    Non-coding      NA             NA                NA      NA
```

---

### **Notes**

* Designed for genomic data where VCF, GFF, and FASTA chromosome IDs match.
* Assumes no alternative splicing or indel handling.
* Protein location is approximate based on linear CDS position.
* The GFF database is created on first run and reused for performance.
