# Fetal Immune Atlas figure colors

Colors for the assay markers and tissue strips in Figure 1b.

| Element | Color |
|---|---|
| scRNA-seq and sc TCR/BCR-seq | `#DA3B35` |
| Spectral flow cytometry | `#ECA243` |
| Plasma proteomics | `#6E68A8` |
| Stimulation proteomics | `#6E68A8` |
| Stimulation bulk RNA-seq | `#43AC71` |
| Organ and cell-count strips | `#5BB9BB` |
| Frames and text | `#202020` |

## Figure 3c differential incoming signaling

The same pathway uses the same color in every cell group. Colors are read from the `M03c` rows in `palette_registry.tsv`.

| Pathway | Color |
|---|---|
| MHC-I | `#E06E52` |
| CLEC | `#5C8D76` |
| CD99 | `#8EC9E5` |
| ICAM | `#B98C39` |
| SELPLG | `#775D86` |
| LCK | `#A85B6A` |
| ITGB2 | `#99C356` |

## Figure 3d CellChat mechanism bubble plot

Early and late stages use the same communication-probability color scale. Colors are read from the `M03d` rows in `palette_registry.tsv`; points are fully opaque.

| Element | Color |
|---|---|
| Communication probability very low | `#2166AC` |
| Communication probability low | `#518BBD` |
| Communication probability midpoint | `#F0E68C` |
| Communication probability high | `#D25E55` |
| Communication probability very high | `#B2182B` |
| Frames and text | `#202020` |
| Vertical guide lines | `#C9C9C9` |

## Figure 5i IGH D-J heatmap

Colors for IGH D-J frequency and correlation direction.

| Element | Color |
|---|---|
| Heatmap low | `#246BAE` |
| Heatmap midpoint | `#FFFFFF` |
| Heatmap warm | `#EDAE11` |
| Heatmap high | `#C71000` |
| Positive correlation | `#8F1336` |
| Negative correlation | `#3D3B25` |

## Figure 2c local-correlation heatmaps

The Figure 2c diverging palette is centered at zero. Dense correlation tiles are rasterized while labels and legends remain vector elements in PDF output.

| Element | Color |
|---|---|
| Negative local correlation | `#343A9A` |
| Zero local correlation | `#18A6B4` |
| Positive local correlation | `#F4D943` |
| Frames and text | `#202020` |

## Figure 2d pathway enrichment

The Figure 2d sequential palette encodes increasing `-log10(BH P)` from light blue to deep violet.

| Element | Color |
|---|---|
| Lower significance | `#8BB8D3` |
| Intermediate significance | `#6959A6` |
| Higher significance | `#292467` |
| Frames and text | `#202020` |

The M1/M7 module gene-change companion uses the Figure 2e early and late stage endpoints, with a neutral gray for genes retained in both stages.

| Element | Color |
|---|---|
| Early-only module genes | `#D8A33E` |
| Shared module genes | `#8A8A8A` |
| Late-only module genes | `#236F68` |

The Early-to-Late module correspondence matrix uses the same sequential family as Figure 2d. The light endpoint represents zero overlap, and the fixed scale spans 0% to 100% of each Early module.

| Element | Color |
|---|---|
| Zero overlap | `#F2F3F3` |
| Lower overlap | `#8BB8D3` |
| Intermediate overlap | `#6959A6` |
| Higher overlap | `#292467` |
| M4-STEM context outline | `#686868` |
| Frames and text | `#202020` |

## Figure 2e temporal transcript programs

The Figure 2e heatmap uses a cool-neutral-warm scale. Age-bin and transcript-program annotations use the fixed colors below.

| Element | Color |
|---|---|
| Heatmap low | `#1498B5` |
| Heatmap midpoint | `#F7F7F2` |
| Heatmap high | `#E85D3F` |
| 16-22 pcw | `#D8A33E` |
| 22-28 pcw | `#D86B3A` |
| 28-34 pcw | `#5A8C73` |
| 34-40 pcw | `#236F68` |
| Program 1 | `#167D8D` |
| Program 2 | `#C678A6` |
| Program 3 | `#D55C8C` |
| Program 4 | `#149D77` |
| Program 5 | `#A77568` |

## Supplementary Figure 11d clone cell number

Cell-type colors for the top-clone cell-count plot.

| Cell type | Color |
|---|---|
| DP(P) T | `#E06E3A` |
| DN(Q) T | `#C73934` |
| Treg | `#008EA0` |
| Naïve CD8 T | `#894491` |
| Naïve CD4 T | `#5A9599` |
| Th17-like innate T | `#DE624E` |
| Cycling Treg | `#88C9D6` |
| Sample strip | `#DAD9D9` |
| Frames and text | `#181818` |

## Figure 7b plasma protein age associations

Protein colors follow manuscript Figure 7b and are read by protein name from the `M07b` rows in `palette_registry.tsv`. Points and regression lines use the same color; confidence bands are gray.

| Protein | Color |
|---|---|
| CD8A | `#49A055` |
| CCL11 | `#1BAA96` |
| IL18 | `#F2BE7D` |
| CD6 | `#9BD2E3` |
| CD5 | `#D93A45` |
| TRAIL | `#EBA042` |
| OSM | `#333366` |
| FGF19 | `#4A2963` |
| IL22RA1 | `#9A3333` |
| FGF5 | `#393C8F` |
| Flt3L | `#CACC4F` |
| CD40 | `#009966` |
| IL8 | `#989A49` |
| OPG | `#CD393A` |
| MMP10 | `#CB9965` |
| MCP1 | `#CD393A` |
| CCL3 | `#2F337B` |

## Graphical abstract

The README graphical abstract is stored in `docs/figures/graphical_abstract.png` and `.pdf`, relative to the repository root.

| Element | Color |
|---|---|
| Earlier development | `#D8A33E` |
| Later maturation | `#236F68` |
| Protein features | `#6E68A8` |
| Receptor details | `#246BAE` |
| Outlines and text | `#202020` |
| Background | `#FFFFFF` |

## Figure 8 multiomics integration

View colors are stored in the `M08_views` rows of `palette_registry.tsv` and read by the Figure 8b script. The HSC/MPP-to-BCR order follows the manuscript. Factor 3 expression uses `#334555`, white and `#8E328A`; GO enrichment uses `#058786`, `#7EB5B4` and `#224767`.

### Figure 8b assay categories

The category strip reads the `M08_modality` rows of `palette_registry.tsv`. Adjacent categories are drawn as continuous color blocks.

| Assay | Color |
|---|---|
| Single-cell RNA-seq | `#E41A1C` |
| TCR | `#11A579` |
| BCR | `#F6CF71` |
| Olink | `#7570B3` |
| Stimulation bulk RNA-seq | `#4DAF4A` |
| Spectral flow cytometry | `#FF7F00` |
