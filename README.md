# Lacertid genomes — analysis scripts

Scripts accompanying the manuscript:

> **Four lacertid lizard genome assemblies reveal a conserved micro/macrochromosome architecture**

## Repository structure

| Folder | Contents |
|---|---|
| `01_macrosynteny/` | Pairwise macrosynteny inference (JCVI/MCScan), block parsing, and minspan sensitivity analysis |
| `02_hic/` | Hi-C interchromosomal contact extraction, PCA, radar plots, and PC1 vs chromosome size |
| `03_gc_gene_density/` | Per-chromosome GC content and gene density computation and visualisation |
| `04_orthofinder_agora/` | OrthoFinder orthogroup inference and AGORA ancestral genome reconstruction |

## Dependencies

- **Macrosynteny**: [JCVI/MCScan](https://github.com/tanghaibao/jcvi), [DIAMOND](https://github.com/bbuchfink/diamond)
- **Hi-C**: [cooler](https://github.com/open2c/cooler), R (`data.table`, `ggplot2`, `fmsb`)
- **GC/gene density**: Python (`pandas`, `matplotlib`, `scipy`)
- **OrthoFinder/AGORA**: [OrthoFinder v3](https://github.com/davidemms/OrthoFinder), [AGORA](https://github.com/sib-swiss/AGORA)

## Notes

Input files (genome FASTA, GFF3, `.cool` Hi-C files) are available from the European Nucleotide Archive under project accession PRJEB[XXXX].
Path variables at the top of each script should be set to match your local directory layout.
