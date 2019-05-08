Day 11: Structural Bioinformatics (Pt1)
================

Revisit the PDB
---------------

**Question 1: Download a CSV file from the PDB site (accessible from “Analyze” -&gt; “PDB Statistics” &gt; “by Experimental Method and Molecular Type”. Move this CSV file into your RStudio project and determine the percentage of structures solved by X-Ray and Electron Microscopy. Also can you determine what proportion of structures are protein? Aim to have a rendered GitHub document with working code that yields your answers**
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Download a CSV file from: <https://www.rcsb.org/stats/summary>

``` r
# percent of xray structures 
db <- read.csv("Data Export Summary.csv", row.names = 1)
head(db)
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 127027          2012               6549     8 135596
    ## NMR                    11064          1279                259     8  12610
    ## Electron Microscopy     2296            31                803     0   3130
    ## Other                    256             4                  6    13    279
    ## Multi Method             131             5                  2     1    139

``` r
sum(db$Total)
```

    ## [1] 151754

``` r
rowSums(db)
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##              271192               25220                6260 
    ##               Other        Multi Method 
    ##                 558                 278

``` r
(db$Total/sum(db$Total)) * 100 
```

    ## [1] 89.35250471  8.30950090  2.06254860  0.18385018  0.09159561

``` r
# proportion of protein 
(sum(db$Proteins)/sum(db$Total)) *100
```

    ## [1] 92.76461

Percent Xray : 89.35%
---------------------

Percent Electron Microscopy: 2.06%
----------------------------------

Proportion of Protein : 92.76461%
---------------------------------

### **Q2: Type HIV in the PDB website search box on the home page and determine how many HIV-1 protease structures are in the current PDB?**

1,157 protease structures
-------------------------

``` r
library(bio3d)
library(bio3d.view)
library(rgl)
pdb <- read.pdb(file = "1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
pdb$seqres
```

    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "PRO" "GLN" "ILE" "THR" "LEU" "TRP" "GLN" "ARG" "PRO" "LEU" "VAL" "THR" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "ILE" "LYS" "ILE" "GLY" "GLY" "GLN" "LEU" "LYS" "GLU" "ALA" "LEU" "LEU" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "ASP" "THR" "GLY" "ALA" "ASP" "ASP" "THR" "VAL" "LEU" "GLU" "GLU" "MET" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "SER" "LEU" "PRO" "GLY" "ARG" "TRP" "LYS" "PRO" "LYS" "MET" "ILE" "GLY" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "GLY" "ILE" "GLY" "GLY" "PHE" "ILE" "LYS" "VAL" "ARG" "GLN" "TYR" "ASP" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "GLN" "ILE" "LEU" "ILE" "GLU" "ILE" "CYS" "GLY" "HIS" "LYS" "ALA" "ILE" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "GLY" "THR" "VAL" "LEU" "VAL" "GLY" "PRO" "THR" "PRO" "VAL" "ASN" "ILE" 
    ##     A     A     A     A     A     A     A     A     A     A     A     A 
    ## "ILE" "GLY" "ARG" "ASN" "LEU" "LEU" "THR" "GLN" "ILE" "GLY" "CYS" "THR" 
    ##     A     A     A     B     B     B     B     B     B     B     B     B 
    ## "LEU" "ASN" "PHE" "PRO" "GLN" "ILE" "THR" "LEU" "TRP" "GLN" "ARG" "PRO" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "LEU" "VAL" "THR" "ILE" "LYS" "ILE" "GLY" "GLY" "GLN" "LEU" "LYS" "GLU" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "ALA" "LEU" "LEU" "ASP" "THR" "GLY" "ALA" "ASP" "ASP" "THR" "VAL" "LEU" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "GLU" "GLU" "MET" "SER" "LEU" "PRO" "GLY" "ARG" "TRP" "LYS" "PRO" "LYS" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "MET" "ILE" "GLY" "GLY" "ILE" "GLY" "GLY" "PHE" "ILE" "LYS" "VAL" "ARG" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "GLN" "TYR" "ASP" "GLN" "ILE" "LEU" "ILE" "GLU" "ILE" "CYS" "GLY" "HIS" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "LYS" "ALA" "ILE" "GLY" "THR" "VAL" "LEU" "VAL" "GLY" "PRO" "THR" "PRO" 
    ##     B     B     B     B     B     B     B     B     B     B     B     B 
    ## "VAL" "ASN" "ILE" "ILE" "GLY" "ARG" "ASN" "LEU" "LEU" "THR" "GLN" "ILE" 
    ##     B     B     B     B     B     B 
    ## "GLY" "CYS" "THR" "LEU" "ASN" "PHE"

``` r
ca.ins <- atom.select(pdb, "calpha")

atom.select(pdb, resno = 10)
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, resno = 10)
    ## 
    ##    Atom Indices#: 16  ($atom)
    ##    XYZ  Indices#: 48  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
pdb$atom[1:2, c("eleno", "elety", "x","y","z")]
```

    ##   eleno elety      x      y     z
    ## 1     1     N 29.361 39.686 5.862
    ## 2     2    CA 30.307 38.663 5.319

Select protein only and write out a PDB file titled "1HSE\_Protein.pdb" Select ligand/drug and write out a PDB file titled "1HSE\_Ligand.pdb"

``` r
protein <- atom.select(pdb, "protein", value = T)
protein
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
ligand <- atom.select(pdb, "ligand", value = T)
ligand
```

    ## 
    ##  Call:  trim.pdb(pdb = pdb, sele)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
write.pdb(protein, file = "1HSE_Protein.pdb")
write.pdb(ligand, file="1HSE_Ligand.pdb")
```

``` r
view(pdb, "overview", col="sse")
```

    ## Computing connectivity from coordinates...

``` r
view(ligand)
```

    ## Computing connectivity from coordinates...

Section 6. Working with multiple PDBs
-------------------------------------

``` r
ids <- c("1TND_B","1AGR_A","1TAG_A","1GG2_A","1KJY_A","4G5Q_A")
files <- get.pdb(ids, split = TRUE)
```

    ## Warning in get.pdb(ids, split = TRUE): ./1TND.pdb exists. Skipping download

    ## Warning in get.pdb(ids, split = TRUE): ./1AGR.pdb exists. Skipping download

    ## Warning in get.pdb(ids, split = TRUE): ./1TAG.pdb exists. Skipping download

    ## Warning in get.pdb(ids, split = TRUE): ./1GG2.pdb exists. Skipping download

    ## Warning in get.pdb(ids, split = TRUE): ./1KJY.pdb exists. Skipping download

    ## Warning in get.pdb(ids, split = TRUE): ./4G5Q.pdb exists. Skipping download

    ## 
      |                                                                       
      |                                                                 |   0%
      |                                                                       
      |===========                                                      |  17%
      |                                                                       
      |======================                                           |  33%
      |                                                                       
      |================================                                 |  50%
      |                                                                       
      |===========================================                      |  67%
      |                                                                       
      |======================================================           |  83%
      |                                                                       
      |=================================================================| 100%

``` r
pdbs <- pdbaln(files, fit = T)
```

    ## Reading PDB files:
    ## ./split_chain/1TND_B.pdb
    ## ./split_chain/1AGR_A.pdb
    ## ./split_chain/1TAG_A.pdb
    ## ./split_chain/1GG2_A.pdb
    ## ./split_chain/1KJY_A.pdb
    ## ./split_chain/4G5Q_A.pdb
    ## .....   PDB has ALT records, taking A only, rm.alt=TRUE
    ## .
    ## 
    ## Extracting sequences
    ## 
    ## pdb/seq: 1   name: ./split_chain/1TND_B.pdb 
    ## pdb/seq: 2   name: ./split_chain/1AGR_A.pdb 
    ## pdb/seq: 3   name: ./split_chain/1TAG_A.pdb 
    ## pdb/seq: 4   name: ./split_chain/1GG2_A.pdb 
    ## pdb/seq: 5   name: ./split_chain/1KJY_A.pdb 
    ## pdb/seq: 6   name: ./split_chain/4G5Q_A.pdb 
    ##    PDB has ALT records, taking A only, rm.alt=TRUE

``` r
pc.xray <- pca(pdbs)
plot(pc.xray)
```

![](Day11_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
pc1 <- mktrj(pc.xray, pc=1, file="pc_1.pdb")
library(bio3d.view)
# Structural displacements captured by PC1
view(pc1)
```

    ## Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()

``` r
# The rglwidget() function from the rgl
# package will show output in your Rmd
# notebook and rendered html_output
# documents
library(rgl)
#rglwidget(pc1)
```
