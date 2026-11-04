<div align="center">

<img src="images/header02.jpg" width="100%">

# Integrated Data Analysis Uncovers New COVID-19 Related Genes and Potential Drug Re-Purposing Candidates

**Alexandros Xenos\*, Noël Malod-Dognin\*, Carme Zambrana and Nataša Pržulj**

<sub>\* These authors contributed equally to this work</sub>

Corresponding author: Prof. N. Pržulj, natasa.przulj@mbzuai.ac.ae

</div>

This repository contains the Python scripts and data for the paper "Integrated data analysis uncovers new COVID-19 related genes and potential drug re-purposing candidates".

In this work, to better understand the COVID-19 molecular basis and design therapeutic strategies, we build upon the recently proposed concept of an integrated cell, iCell, fusing three omics, tissue-specific human molecular interaction networks. We apply this methodology to construct infected and control iCells using gene expression data from patient samples and three cell lines. We find large differences between patient-based and cell line-based iCells, suggesting that cell lines are ill-suited to study this disease. We compare patient-based infected and control iCells and uncover genes whose functioning and network wiring patterns are altered by the disease. We validate in the literature that 18 out of the top 20 of the most rewired genes are indeed COVID-19 related. Since only three of these genes are targets of approved drugs, we apply another data fusion step to predict drugs for re-purposing. Our most interesting prediction is Artenimol, an antimalarial agent targeting ZFP62, one of the genes that we newly identified as COVID-19 related. This drug is a derivative of Artemisinin drugs that are already under clinical investigation for their potential role in the treatment of COVID-19. Our results demonstrate further applicability of the iCell framework in integrative comparative studies of human diseases.

## Downloads

- [Data](https://github.com/przulj-lab/covid-icell/releases/download/v1.0/Data.zip)

Download and unzip the data in the repository folder.

## Usage

### Generating the iCells

Run `Main_create_iCells.py` with an integer between 0 and 7 that selects the tissue:

```bash
python Main_create_iCells.py <tissue>
```

| Parameter | Tissue |
|---|---|
| `0` | A549_Control |
| `1` | A549_Infected |
| `2` | NHBE_Control |
| `3` | NHBE_Infected |
| `4` | CALU_Control |
| `5` | CALU_Infected |
| `6` | Patient_Control |
| `7` | Patient_Infected |

### Drug re-purposing

To perform the data fusion step that predicts drugs for re-purposing, run:

```bash
python Main_drug_repurposing.py
```

It takes as input the G matrix factor of the patients' data and fuses it with the known drug-target interactions.

## Citation

A. Xenos, N. Malod-Dognin, C. Zambrana and N. Pržulj. Integrated data analysis uncovers new COVID-19 related genes and potential drug re-purposing candidates. *International Journal of Molecular Sciences* 24(2), 1431 (2023). [PMID: 36674947](https://pubmed.ncbi.nlm.nih.gov/36674947/)
